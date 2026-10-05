"""Surface energy-balance flux terms, transcribed from sflux.F90.

Fortran source: sflux.F90 (9 442 bytes, 221 lines, 1 subroutine). See
docs_transcricao/doc-transcricao-sflux.docx for the full transcription
rationale, the call-site argument mapping, and the pending validation work
referenced throughout this module as "the transcription doc".

What it computes
------------------
Given the radiative inputs, the turbulent-transfer coefficients, and the
current thermodynamic state of one surface node, this assembles every term
of the surface energy balance -- downwelling/upwelling solar, downwelling/
upwelling longwave, precipitation heat, sensible heat, latent heat, ground
conduction, and a lumped "water/vapour flow + phase change" term -- their
sum ``net``, and the two Jacobian sums ``dnet1``/``dnet2`` that
soil_tmp.F90 consumes in its implicit temperature solve. It is called from
surfenergy.F90 (twice) and new_profile.F90 (three times).

Sign convention (Fortran lines 75-77): ``> 0`` -> towards the surface,
``< 0`` -> away from it; positive z is upward relative to sea level.

``vflag`` -- which surface, and which temperature is T1
-------------------------------------------------------
``vflag == 1``  -> low vegetation: ``vf = +1``, and the routine reads
``ftemp`` (the foliage temperature) from state for ``tt1``; at both
observed ``vflag == 1`` call sites the caller also passes
``temp1 = ftemp`` (doc point 9), so the longwave term
``vf*cc1*(temp1**4 - ftemp**4)`` is *identically zero* on those calls
while its linearisation is not.

``vflag == 0``  -> bare soil / snow: ``vf = -1``, ``tt1 = temp1`` (the
argument). The longwave-derivative branch at Fortran line 165 keys on
``vflag == 0`` specifically. The Fortran header comment records the
argument convention as: ``vflag == 1`` -> T1 = ftemp, T2 = toptemp;
``vflag == 0`` -> T1 = toptemp, T2 = ftemp.

Not pure -- ``state: FasstState``
---------------------------------
Like the impure functions in fasst/functions.py, ``sflux`` reads
module-global state alongside its declared arguments: the scalars
``sigfl, ftemp, epf, emis, albf, albedo_fasst, hfol_tot, hsaccum, ptemp,
deltat_fasst, mflag`` and the node-indexed arrays ``flowu, flowl, ice,
wvc, fv1, soil_moist, phead`` (all indexed by ``sn``). It therefore takes
``state`` as an explicit first parameter. Promoting these to plain
value arguments is phase-2 work (doc PENDENCIAS), not done here.

Returning a dataclass, not 12 out-arguments -- doc point 1
----------------------------------------------------------
The Fortran declares 12 ``intent(out)`` arguments
(``radswd, radswu, radlwd, radlwu, phw, shw, lhw, chw, chw1, net, dnet1,
dnet2``) and additionally computes ``dchw`` -- routed through the same
``anint(x*1d20)*1d-20`` rounding block (Fortran lines 203, 207) -- and
``dchw1`` -- initialised and never assigned -- both of which are then
discarded. ``dchw`` is the conduction-term Jacobian ``d(chw)/dT1``; it is
*not* added into ``dnet1`` (Fortran line 215), because the Fortran comment
says that contribution is supplied by soil_tmp.F90 itself. This is the
same pattern sp_humid.F90 shows (13 computed-and-discarded quantities, 6
still rounded). Per the doc's recommendation ("Recomendo o dataclass --
custo nulo, ganho de diagnóstico") :func:`sflux` returns a single
:class:`SurfaceFluxes` carrying:

  * the 12 Fortran ``intent(out)`` values (unchanged arithmetic);
  * ``dchw`` and ``dchw1`` -- the computed-and-discarded pair above;
  * ``cheati, cheatv, cheatw, fheat, fchw`` -- the five intermediates that
    sum into ``chw1``. These are genuine Fortran locals but the Fortran
    never routes them out; exposing them is a TRANSCRIPTION-TIME DIAGNOSTIC
    ADDITION, not a preserved discard. They are ``0.0`` whenever
    ``vflag == 1`` or ``hsaccum > eps``.

Dead Fortran locals NOT reproduced: ``W, gamma, dgdT, t1`` (declared,
zeroed, never read -- ``gamma``/``dgdT`` name a surface-tension term that
was never implemented, doc point 2).

Known, UNCORRECTED points (all doc-verified; do not silently "fix")
------------------------------------------------------------------------
* doc point 3 -- ``2500775.6`` (Fortran line 123) is a hand-written latent
  heat of vaporisation [J/kg] used in the vapour phase-change term. There
  is no named constant for it. It is inconsistent with the module's own
  pair by ``LHSUB - LHFUS = 2.838e6 - 3.335e5 = 2.5045e6`` -- a 3.7e3 J/kg
  (1.5e-3 relative) gap. Kept literal.
* doc point 4 -- the vapour phase-change term evaluates ``spheats`` and
  ``dense`` with phase ``d1i = 1`` (WATER), not ``0`` (WATER_VAPOR); the
  water content change is priced with liquid-water thermodynamic
  properties. ``dense``/``spheats`` are declared as bare ``real(dp)``
  externals in the Fortran (no interface), so nothing checks this at the
  call site. Transcribed exactly; flagged for the author.
* doc point 5 -- three rounding idioms in one routine: ``anint(x*1d20)*
  1d-20`` on almost every line, ``anint(x*1d18)*1d-18`` on ``dnet1``/
  ``dnet2`` (Fortran 216, 218), and ``aint(x*1d18)*1d-18`` -- truncation
  toward zero, not rounding -- on ``net`` (Fortran 213). Plus
  ``aint(|x - mflag|*1d5)*1d-5`` as the missing-data test for ``swu`` and
  ``lwu`` (Fortran 143, 155). All four preserved verbatim via
  :func:`_anint` / :func:`_aint`; see fasst/functions.py for why the idiom
  is never removed during transcription.
* doc point 8 -- ``radlwu`` (Fortran 160-163) mixes ``tt1*tt1*tt1*tt1``
  (repeated multiply -> strict left-to-right rounding) with
  ``temp1**4d0`` / ``ftemp**4d0`` (a *real* exponent -> ``pow``) in one
  expression. Both are preserved as written. The single-line Fortran guard
  ``if(temp1 > eps.and.ftemp > eps)`` covers ONLY the ``radlwu``
  assignment, not the ``dradlwu1``/``dradlwu2`` lines below it -- so when
  ``temp1 <= eps`` or ``ftemp <= eps`` the flux stays ``0`` while its
  derivatives are still computed.
* doc point 9 -- veg-branch Jacobian vs. a zero flux term. At both
  observed ``vflag == 1`` call sites the caller passes ``temp1 == ftemp``,
  so ``vf*cc1*(temp1**4 - ftemp**4) == 0`` contributes nothing to
  ``radlwu`` -- yet ``dradlwu1`` carries ``-vf*cc1*4*tt1**3`` (Fortran
  169, sign opposite the ground branch's line 166) and ``dradlwu2``
  carries ``+vf*cc1*4*temp2**3`` (Fortran 170). The ground branch (line
  166-167) is self-consistent with its flux term; at those call sites the
  veg branch's ``cc1`` linearisation has no matching non-zero flux.
  Transcribed exactly; open question for the author.
* doc point 10 -- unguarded divisions, no ``eps`` check in the Fortran:
  ``f3 = h/(deltat_fasst*s)`` (line 113) and ``dlhw``/``dlhwfg``'s
  ``.../d1`` (lines 193, 195). Fortran turns a zero denominator into
  +/-inf or nan and keeps going. Routed through :func:`_fdiv` here (plain
  IEEE-754, no exception), mirroring sp_humid.py. The conduction
  divisions ``kth/h`` (lines 201-207) sit inside ``if(h /= 0d0)`` -- an
  EXACT non-zero test, not an ``eps`` guard -- so they stay plain ``/``.
* doc point 7 -- ``sgl`` vs ``sgl1``: the solar terms use
  ``sgl1 = max(0, 1 - sigfl)``; the longwave and conduction terms use the
  unclamped ``sgl = 1 - sigfl``. They differ only if ``sigfl > 1``. The
  Fortran carries commented-out ``sgl`` variants of the three solar lines
  (140, 144, 147) and a ``*5d-1)`` fragment on the clamp (line 93),
  evidence the split was introduced by hand.
* doc point "s per caller" -- ``s`` (Fortran ``s``) is passed as the
  literal ``1`` from surfenergy.F90 (``d2i = 1``) but as ``step`` from
  new_profile.F90, where ``step = timstep*3600/deltat_fasst`` is the
  sub-step count of the current met interval. So ``f3 = h/(deltat_fasst*s)``
  spreads the phase-change energy over one sub-step from one caller and
  over the whole met interval from the other. Preserved; see doc.
* doc point "inert spheats arg" -- ``cheati``/``cheatv`` call
  ``spheats(tt1, ...)`` while ``cheatw`` calls ``spheats(temp1, ...)``.
  The block only runs when ``vflag != 1``, where ``tt1 == temp1``, so this
  is inert today -- but a latent trap if the block ever runs with
  ``vflag == 1``.
"""

import math
from dataclasses import dataclass

import numpy as np

from .constants import EPS, GRAV, KVEG, LHFUS, SIGMA, TREF
from .functions import Phase, dense, spheats
from .state import FasstState

__all__ = ["SurfaceFluxes", "sflux"]

# Hand-written latent heat of vaporisation [J/kg], Fortran sflux.F90 line 123.
# No named constant exists for it in fasst_global.F90. Kept as a module
# literal (not folded into a call) so it stays greppable -- see doc point 3.
_LHVAP_LITERAL = 2500775.6


def _anint(x: float, p: int) -> float:
    """Fortran ``ANINT(x * 10**p) * 10**-p`` -- round half **away from zero**.

    NOT ``round`` / ``np.rint`` (both round half to even). Mirrors
    ``fasst/functions.py::_anint``; duplicated so each transcribed file
    stays a self-contained record of its own source's idioms. NaN/inf pass
    through unchanged, as ``ANINT`` does.

    Never the identity: ``x * 10**p`` and the closing ``* 10**-p`` are each
    an inexact float64 multiply. For ``p = 20`` and ``|x| >~ 4.5e-5`` the
    ``anint`` step is a no-op and only ~1 ULP of double-rounding remains;
    for smaller ``|x|`` it quantises with an absolute ~1e-20 step. See
    fasst/functions.py's module docstring for the full analysis.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0 ** p)
    return math.copysign(math.floor(abs(xs) + 0.5), xs) * (10.0 ** -p)


def _aint(x: float, p: int) -> float:
    """Fortran ``AINT(x * 10**p) * 10**-p`` -- truncate **toward zero**.

    ``AINT`` differs from ``ANINT`` only in dropping the ``+ 0.5``: it
    chops the fractional part rather than rounding it. Fortran uses it on
    ``net`` (``aint(net*1d18)*1d-18``, sflux.F90 line 213) and on the
    missing-data comparisons for ``swu``/``lwu`` (``aint(|.|*1d5)*1d-5``,
    lines 143, 155). NaN/inf pass through unchanged.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0 ** p)
    return math.copysign(math.floor(abs(xs)), xs) * (10.0 ** -p)


def _round20(x: float) -> float:
    """Fortran ``anint(x*1d20)*1d-20`` -- the routine's default idiom."""
    return _anint(x, 20)


def _round18(x: float) -> float:
    """Fortran ``anint(x*1d18)*1d-18`` -- used only on ``dnet1``/``dnet2``."""
    return _anint(x, 18)


def _trunc18(x: float) -> float:
    """Fortran ``aint(x*1d18)*1d-18`` -- used only on ``net`` (truncation)."""
    return _aint(x, 18)


def _fdiv(x: float, y: float) -> float:
    """Plain IEEE-754 float division: ``x / 0.0`` yields +/-inf or nan
    instead of raising, matching what the Fortran arithmetic silently does
    at the call sites that carry no ``eps`` guard (``f3``, ``dlhw``,
    ``dlhwfg`` -- see module docstring, doc point 10). Mirrors
    ``fasst/sp_humid.py::_fdiv``.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


@dataclass
class SurfaceFluxes:
    """Return value of :func:`sflux`.

    All 12 Fortran ``intent(out)`` values, plus the computed-and-discarded
    ``dchw``/``dchw1`` pair, plus the five ``chw1`` intermediates (the
    latter are a transcription-time diagnostic addition -- see module
    docstring, "Returning a dataclass"). Every field is in W/m^2 except the
    derivative fields (W/m^2 per K) as noted.

    Sign convention: ``> 0`` towards the surface, ``< 0`` away from it.
    """

    # -- Fortran intent(out), same order --------------------------------------
    radswd: float   # incoming / downwelling solar (>= 0)
    radswu: float   # upwelling / reflected solar (<= 0)
    radlwd: float   # incoming longwave / IR (>= 0)
    radlwu: float   # upwelling / emitted longwave / IR (<= 0)
    phw: float      # precipitation heat flux (+/-)
    shw: float      # sensible heat flux (+/-)
    lhw: float      # latent heat flux (+ condensation / - evaporation)
    chw: float      # ground conductive heat flux (+/-)
    chw1: float     # lumped water/vapour flow + phase-change heat flux (+/-)
    net: float      # sum of all flux terms (+/-) -- AINT-truncated at 1e18
    dnet1: float    # d(net)/dT1  [W/m^2/K] -- w.r.t. layer of interest (soil_tmp.F90)
    dnet2: float    # d(net)/dT2  [W/m^2/K] -- w.r.t. the other layer (soil_tmp.F90)

    # -- computed by the Fortran, routed through its rounding block, then
    # discarded (doc point 1) ----------------------------------------------
    dchw: float     # d(chw)/dT1 [W/m^2/K] -- NOT added into dnet1 (see docstring)
    dchw1: float    # d(chw1)/dT1 -- Fortran initialises it and never assigns it (always 0.0)

    # -- transcription-time diagnostic addition: the five terms summed into
    # chw1 (Fortran chw1 = -cheati + cheatv + cheatw + fheat + fchw). All
    # five are 0.0 when vflag == 1. When vflag == 0: fheat/fchw are always
    # computed (fchw needs hfol_tot > eps); cheati/cheatv/cheatw need
    # additionally hsaccum <= eps and their own content guard. -------------
    cheati: float   # ice/water phase-change heat (+/-)
    cheatv: float   # vapour/water phase-change heat (+/-), includes state.fv1[sn]
    cheatw: float   # water-flow heat (+/-)
    fheat: float    # heat carried by water+vapour flow through the node (+/-)
    fchw: float     # conductive heat from foliage to the ground surface (+/-)


def sflux(
    state: FasstState,
    vflag: int,
    s: int,
    sn: int,
    swu: float,
    swd: float,
    lwu: float,
    lwd: float,
    sfac: float,
    temp1: float,
    cc1: float,
    taf: float,
    ph: float,
    sh: float,
    lh: float,
    sphumidaf: float,
    sph: float,
    kth: float,
    temp2: float,
    h: float,
    r1: float,
    r2: float,
    iceo: float,
    wvco: float,
    smo: float,
    dqdt1: float,
    dqdt2: float,
    d1: float,
) -> SurfaceFluxes:
    """Assemble the surface energy-balance flux terms.

    Fortran: sflux.F90, subroutine ``sflux``. No subroutines called; uses
    the functions ``dense`` and ``spheats``.

    Args are in Fortran declaration order (the 12 ``intent(out)`` values
    are returned in :class:`SurfaceFluxes` instead of being args):

        state: the model state; see module docstring for the fields read.
        vflag: 1 = low vegetation, 0 = bare soil / snow. See module
            docstring, "``vflag``".
        s: sub-step divisor for the phase-change term -- ``1`` from
            surfenergy.F90, ``step`` from new_profile.F90 (doc point
            "s per caller").
        sn: node index into ``state.flowu / flowl / ice / wvc / fv1 /
            soil_moist / phead``.
        swu, swd: upwelling / downwelling solar flux [W/m^2]. ``swu`` equal
            to ``state.mflag`` (within an ``aint(|.|*1e5)*1e-5`` tolerance)
            means "missing" -> the reflected term falls back to
            ``swd * albedo``.
        lwu, lwd: upwelling / downwelling longwave flux [W/m^2]. ``lwu ==
            state.mflag`` (same tolerance) means "missing" -> the emitted
            term is modelled from ``tt1**4`` instead.
        sfac: fraction of shortwave penetrating the surface (bare/snow
            only; the veg branch forces ``1.0``).
        temp1: T1 [K]. With ``vflag == 0`` this is ``tt1`` directly; with
            ``vflag == 1`` the caller passes ``state.ftemp`` here and
            ``tt1`` is read from ``state.ftemp`` too.
        cc1: canopy longwave-exchange coefficient [W/m^2/K^4].
        taf: air/foliage reference temperature for the sensible term [K].
        ph: precipitation heat-transfer coefficient [W/m^2/K]; also the
            raw ``d(phw)/dT1``.
        sh, lh: sensible / latent transfer coefficients [W/m^2 per unit
            driving difference].
        sphumidaf: reference specific humidity for the latent term.
        sph: surface saturation specific humidity.
        kth: thermal conductivity across the conduction path [W/m/K].
        temp2: T2 [K] -- the "other" layer (see module docstring).
        h: conduction path length [m]; also the numerator of ``f3``.
            ``h == 0.0`` (exact) disables the conduction term.
        r1, r2: relative-humidity factors for the latent term and its two
            derivatives. Swapped between the ground call
            (``rhsurf, rpp``) and the veg call (``rpp, rhsurf``) -- this is
            the T1/T2 convention, not a defect.
        iceo, wvco, smo: previous-step ice / water-vapour / soil-moisture
            content at the node, for the phase-change deltas. From
            new_profile.F90 these are genuine previous-step values
            (``ioo/woo/smoo``); from surfenergy.F90 they are the CURRENT
            ``ice(nnodes)/wvc(nnodes)/soil_moist(nnodes)`` -- see doc
            PENDENCIAS on auditing whether ``sn == nnodes`` there makes the
            deltas identically zero.
        dqdt1, dqdt2: ``d(specific humidity)/dT`` for T1 / T2. Also swapped
            between the ground and veg calls.
        d1: appears only in the latent-derivative terms, both as a
            subtrahend and as a divisor (Fortran lines 193, 195). Its
            meaning is not documented anywhere in the Fortran source -- see
            doc PENDENCIAS.

    Returns:
        A :class:`SurfaceFluxes`.
    """
    # -- branch-independent setup ------------------------------------------
    # Fortran zeroes ~50 locals here; only the ones actually read below are
    # kept. vf: +1 for vegetation, -1 for ground. f1/f2 weight sigfl in the
    # sensible / latent derivatives and swap between the two modes.
    fchw = 0.0    # conductive heat from foliage to ground surface [W/m^2]
    fheat = 0.0   # heat carried by water+vapour flow [W/m^2]
    cheati = 0.0  # ice/water phase-change heat [W/m^2]
    cheatv = 0.0  # vapour/water phase-change heat [W/m^2]
    cheatw = 0.0  # water-flow heat [W/m^2]
    chw1 = 0.0
    radlwu = 0.0    # Fortran line 63 -- stays 0.0 if the lwu-missing branch's
                    # `temp1 > eps .and. ftemp > eps` guard is false (doc point 8)
    dradlwu1 = 0.0  # d(radlwu)/dT1
    dradlwu2 = 0.0  # d(radlwu)/dT2
    dchw = 0.0      # d(chw)/dT1 -- computed, never returned into dnet1
    dchw1 = 0.0     # d(chw1)/dT1 -- Fortran never assigns this; stays 0.0

    if vflag == 1:  # low vegetation
        vf = 1
        alb = state.albf                     # fallback surface reflectance -- IS read below
                                             # when swu == mflag (new_profile.F90 can pass
                                             # the missing sentinel as swu)
        sgl = state.sigfl                    # low-veg density
        sgl1 = state.sigfl
        sfc = 1.0                             # veg branch forces full penetration
        tt1 = state.ftemp                    # T1 [K]
        ems = state.epf                      # surface emissivity
        f1 = 6e-1
        f2 = 1e-1
    else:  # bare soil / snow
        vf = -1
        alb = state.albedo_fasst
        sgl = 1.0 - state.sigfl
        sgl1 = max(0.0, 1.0 - state.sigfl)   # Fortran carries a commented `*5d-1)` here
        sfc = sfac
        tt1 = temp1
        ems = state.emis
        f1 = 1e-1
        f2 = 6e-1

        # foliage conduction
        if state.hfol_tot > EPS:
            fchw = sgl * KVEG * (state.ftemp - temp1) / state.hfol_tot
            fchw = _round20(fchw)

        # water + vapour flow conduction
        fheat = (state.flowu[sn] - state.flowl[sn]) * abs(temp1 - TREF)
        fheat = _round20(fheat)

        # phase change; warming/cooling of pore-space materials
        if state.hsaccum <= EPS:
            f3 = _fdiv(h, state.deltat_fasst * s)  # unguarded in Fortran -- see doc point 10

            # ice -- Fortran d1i = 2 (ICE); spheats at tt1, dense at temp1
            if state.ice[sn] + iceo > EPS:
                cheati = (
                    (LHFUS + spheats(tt1, Phase.ICE) * abs(temp1 - TREF))
                    * dense(temp1, 0.0, Phase.ICE)
                    * (state.ice[sn] - iceo)
                    * f3
                )
            # vapour -- Fortran d1i = 1 (WATER, not WATER_VAPOR); see doc point 4
            if state.wvc[sn] + wvco > EPS:
                cheatv = (
                    (_LHVAP_LITERAL + spheats(tt1, Phase.WATER) * abs(temp1 - TREF))
                    * dense(temp1, 0.0, Phase.WATER)
                    * ((state.wvc[sn] - wvco) * f3)
                    + state.fv1[sn]
                )
            # water -- Fortran d1i = 1 (WATER); spheats at temp1 (not tt1) here
            if state.soil_moist[sn] + smo > EPS:
                cheatw = (
                    dense(temp1, 0.0, Phase.WATER)
                    * (
                        spheats(temp1, Phase.WATER) * abs(temp1 - TREF)
                        + GRAV * state.phead[sn]
                    )
                    * (state.soil_moist[sn] - smo)
                    * f3
                )

        chw1 = -cheati + cheatv + cheatw + fheat + fchw
        chw1 = _round20(chw1)

    # -- solar / short-wave ----------------------------------------------------
    radswd = sgl1 * swd * sfc                                    # Fortran: commented alt uses sgl
    radswd = _round20(radswd)
    if _aint(abs(swu - state.mflag), 5) > EPS:                   # swu is a real measurement
        radswu = -sgl1 * swu * sfc
    else:                                                        # swu missing -> reflect swd
        radswu = -sgl1 * swd * alb * sfc
    radswu = _round20(radswu)

    # -- IR / long-wave ------------------------------------------------------
    radlwd = sgl * lwd
    radlwd = _round20(radlwd)
    if _aint(abs(lwu - state.mflag), 5) > EPS:                   # lwu is a real measurement
        radlwu = -sgl * lwu
        dradlwu1 = 0.0
        dradlwu2 = 0.0
    else:
        # Fortran: this assignment ONLY is guarded by temp1/ftemp > eps;
        # the dradlwu* lines below are NOT (doc point 8). tt1**4 is a
        # repeated multiply; temp1/ftemp use a real exponent (`** 4d0`).
        if temp1 > EPS and state.ftemp > EPS:
            radlwu = (
                -sgl * (1.0 - ems) * lwd
                - sgl * ems * SIGMA * (tt1 * tt1 * tt1 * tt1)
                + vf * cc1 * (temp1 ** 4.0 - state.ftemp ** 4.0)
            )

        if vflag == 0:  # ground
            dradlwu1 = (-sgl * ems * SIGMA + vf * cc1) * (4.0 * tt1 * tt1 * tt1)
            dradlwu2 = -vf * cc1 * 4.0 * state.ftemp * state.ftemp * state.ftemp
        else:  # vegetation -- see doc point 9 (zero flux term, non-zero Jacobian)
            dradlwu1 = (-sgl * ems * SIGMA - vf * cc1) * (4.0 * tt1 * tt1 * tt1)
            dradlwu2 = vf * cc1 * 4.0 * temp2 * temp2 * temp2

    radlwu = _round20(radlwu)
    dradlwu1 = _round20(dradlwu1)
    dradlwu2 = _round20(dradlwu2)

    # -- precipitation ------------------------------------------------------
    phw = ph * (state.ptemp - tt1)
    phw = _round20(phw)
    dphw = _round20(ph)

    # -- sensible ---------------------------------------------------------------
    shw = sh * (taf - tt1)
    shw = _round20(shw)
    dshw = sh * (f1 * state.sigfl - 1.0)                        # w.r.t. T1
    dshw = _round20(dshw)
    dshwfg = sh * f2 * state.sigfl                              # w.r.t. T2
    dshwfg = _round20(dshwfg)

    # -- latent ---------------------------------------------------------------
    lhw = lh * (sphumidaf - r1 * sph)
    lhw = _round20(lhw)
    dlhw = _fdiv(lh * (f1 * state.sigfl - d1) * r1 * dqdt1, d1)  # w.r.t. T1 -- unguarded /d1
    dlhw = _round20(dlhw)
    dlhwfg = _fdiv(lh * f2 * state.sigfl * r2 * dqdt2, d1)       # w.r.t. T2 -- unguarded /d1
    dlhwfg = _round20(dlhwfg)

    # -- ground conduction -- `h /= 0d0` is an EXACT test, not an eps guard ---
    chw = 0.0
    if h != 0.0:
        if state.hsaccum > EPS:
            chw = sgl * kth * (temp2 - tt1) / h
            chw = _round20(chw)
            dchw = -_round20(sgl * (kth / h))                   # minus OUTSIDE the round
        else:
            chw = kth * (temp2 - tt1) / h
            chw = _round20(chw)
            dchw = -_round20(kth / h)                           # minus OUTSIDE the round

    # -- net and Jacobian sums ----------------------------------------------
    net = radswd + radswu + radlwd + radlwu + shw + lhw + chw + chw1 + phw
    net = _trunc18(net)                                         # AINT, not ANINT

    dnet1 = dradlwu1 + dshw + dlhw - dphw                       # w.r.t. layer of interest
    dnet1 = _round18(dnet1)
    dnet2 = dradlwu2 + dshwfg + dlhwfg                          # w.r.t. non-layer of interest
    dnet2 = _round18(dnet2)

    return SurfaceFluxes(
        radswd=radswd,
        radswu=radswu,
        radlwd=radlwd,
        radlwu=radlwu,
        phw=phw,
        shw=shw,
        lhw=lhw,
        chw=chw,
        chw1=chw1,
        net=net,
        dnet1=dnet1,
        dnet2=dnet2,
        dchw=dchw,
        dchw1=dchw1,
        # float(...) only strips the numpy scalar wrapper the `state.<array>[sn]`
        # reads introduce -- float64 -> float is the identity on the value.
        # These five are NOT rounded (the Fortran does not round them).
        cheati=float(cheati),
        cheatv=float(cheatv),
        cheatw=float(cheatw),
        fheat=float(fheat),
        fchw=float(fchw),
    )
