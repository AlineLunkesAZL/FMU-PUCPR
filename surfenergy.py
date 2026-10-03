"""Surface sensible/latent/precip heat coefficients and surface energy
balance, transcribed from surfenergy.F90.

Fortran source: surfenergy.F90 (26 158 bytes, 706 lines). The file holds
three program units:

  * ``subroutine surfenergy`` (lines 1-406) -- transcribed here as
    :func:`surfenergy`.
  * ``real(dp) function drag`` (lines 413-515) -- Louis (1979) / Jacobson
    (2005) / Mascart (1995) drag coefficient. **NOT transcribed** -- it is
    dead code (its only two call sites, surfenergy.F90:252 and :269, are
    commented out; ``newdrag`` replaced it). See "Dead code: ``drag``"
    below and the transcription doc.
  * ``real(dp) function newdrag`` (lines 522-705) -- Louis (1979) / Li
    (2010) / Beljaars & Holtslag (1991) iterative drag coefficient.
    Transcribed here as :func:`newdrag`.

See docs_transcricao/doc-transcricao-surfenergy.docx for the full
transcription rationale, the call-site argument mapping, and the pending
validation work referenced throughout this module as "the transcription
doc".

What :func:`surfenergy` does
-----------------------------
For one surface node it (1) calls :func:`~fasst.sp_humid.sp_humid` twice
-- once for the *ground* saturation quantities, once for the *air*
quantities; (2) sets the surface relative humidity ``rhsurf`` from a
branch tree keyed on melt/snow/water/ice/pond state; (3) sets the foliage
reference temperature ``ftemp`` and, on the vegetated path, calls
``lowveg_met`` (see below); (4) computes the drag coefficients via
:func:`newdrag`; (5) builds the latent, sensible and precipitation heat
transfer coefficients; (6) computes the shortwave-penetration fraction
``sfac``; and (7) calls :func:`~fasst.sflux.sflux` once (ground) or twice
(ground + foliage) to assemble the two "surface" energy sums ``isurfg``
and ``isurff``.

The Hughes et al. (1993) reference in the Fortran header
(Int. J. Remote Sensing 14(7):1383-1412, pp.1389-1390) is the source of
the sensible/latent formulation.

``lowveg_met`` -- an untranscribed dependency (doc PENDENCIAS)
-------------------------------------------------------------
``module_lowveg.F90``'s ``lowveg_met`` (root uptake, precip interception,
transpiration for low vegetation) is NOT transcribed. :func:`surfenergy`
takes it as an injected callable ``lowveg_met`` (keyword-only). It is only
reached on the ``sigfl > eps and icase != 4`` path; on the bare-soil /
snow path (``sigfl <= eps or icase == 4``) the argument is unused and may
be left ``None``. Expected contract, matching the Fortran call at
surfenergy.F90:211-213::

    lowveg_met(state, ii, iter_, pt1, pt2, mixrgr, wetness, taf, dmet,
               wetbulba) -> LowVegMetResult

where ``pt1 = int(met[iw, PT])``, ``pt2 = int(met[iw, PT2])``,
``mixrgr = mixrgrs`` (ground saturation mixing ratio), ``wetness = rhsurf``.
It is ``intent(inout)`` on ``dmet`` (it MAY mutate the array in place) and
on ``wetbulba`` (it OVERWRITES it -- so ``state.ptemp`` at the end of
:func:`surfenergy` is :func:`sp_humid`'s wet-bulb value only on the
bare-soil path; on the vegetated path it is whatever ``lowveg_met``
returned). It also sets ``state.uaf`` (module_lowveg.F90:108) and
``state.chnf`` is read by :func:`newdrag`'s caller-side blend.

Not pure -- reads AND writes ``state``
--------------------------------------
Reads: ``single_multi_flag, freq_id, nnodes, iw, met, icase, sigfl, epf,
emis, hm, hi, hsaccum, newsd, dsnow, hfol_tot, lail, iheightn, elev,
rough, z0l, vegl_type, idens, node_type, ice, soil_moist, stt, nsoilp,
phead, wvc, sdens, tmelt, ptemp`` and the constants below.
**Writes** (as the Fortran does, via module globals):
``state.ftemp, state.ptemp, state.uaf, state.chnf, state.meltfl,
state.error_code, state.error_type``. A caller that reads only the
returned :class:`SurfaceEnergy` would miss every one of those -- they are
enumerated in :func:`surfenergy`'s docstring.

Known, UNCORRECTED points (all doc-verified; do not silently "fix")
------------------------------------------------------------------------
* doc point 1 -- dead code ``drag`` (fan-out 0, both call sites commented).
  NOT transcribed, exactly as ``fasst/functions.py`` omits ``dddt`` /
  ``USCS_sandsiltclay``. Behavioural note for anyone tempted to resurrect
  it: ``drag`` blends ``(1-sigfl)*Cpng0 + sigfl*chnf`` *inside* the
  function (line 511), whereas ``newdrag`` returns a bare coefficient and
  the caller does that blend (lines 255, 272) -- naively swapping ``drag``
  back in would double-apply the vegetation blend.
* doc point 2 -- ``2500775.6`` reappears here (line 229) as the intercept
  of ``lh = 2500775.6 - 2369.729*(Tair - Tref)`` [J/kg], the WES latent
  heat of vaporisation fit. This is the SAME hand-written constant flagged
  in sflux.py's doc point 2 (there it is ``_LHVAP_LITERAL``); the codebase
  now has two independent copies. Both kept literal (``_LHEVAP_INTERCEPT``
  / ``_LHEVAP_SLOPE`` here).
* doc point 3 -- ``stempt < 35.9`` (line 132): ``35.9`` is a magic Kelvin
  sanity floor (a real surface temperature can never be that cold) with no
  named constant. When tripped, the Fortran does ``write(10, ...)`` --
  hardcoded unit 10, guarded by ``single_multi_flag == 0`` -- substitutes
  ``stempt = dmet(4)`` (air temperature for surface temperature,
  silently), and sets ``error_code = 1, error_type = 3``. Transcribed as a
  ``warnings.warn`` plus the same state writes; ``_STEMPT_FLOOR = 35.9``.
* doc point 4 -- negative ``sfac`` for open water without ice. Line 366:
  ``sfac = 1d0 - dexp(0.4d0)`` evaluates to ~-0.492 (note the POSITIVE
  exponent, unlike every sibling line which uses ``dexp(-t2)``). A
  negative shortwave-penetration fraction flows into
  :func:`~fasst.sflux.sflux` as ``sfac`` and makes ``radswd`` negative.
  Very likely a sign typo for ``dexp(-0.4d0)``. Transcribed exactly.
* doc point 5 -- competing ``beta`` (surface-wetness) formulations. Active
  (line 173, Lee & Pielke 1992): ``(1 - cos(pi*sm/por))**2``. Commented
  (Maneta & Silverman 2013): ``(1 - cos((pi*sm/por)**2))`` -- different
  parenthesisation, different function. Same class as module_snow's
  point 11.
* doc point 6 -- five rounding idioms: ``anint(x*1d20)*1d-20`` (most),
  ``anint(x*1d10)*1d-10`` on ``lheatg`` (line 279), ``aint(x*1d10)*1d-10``
  -- truncation -- on ``evapcm`` (line 282), ``anint(x*1d18)*1d-18`` on
  ``isurff`` (line 399) while ``isurfg`` at line 385 uses ``1d20`` for the
  structurally identical sum, and bare ``aint(x)`` to truncate the precip
  type code (lines 310, 314). All preserved via :func:`_anint` /
  :func:`_aint`.
* doc point 7 -- ``sp_humid`` positional-out mapping. The Fortran passes 10
  positional out-arguments into reused locals whose names do NOT match the
  callee's parameter names. Verified against sp_humid.F90's argument list
  order ``mixr, dmrdt, vpress, wetbulb, rhov, rhoda, vpsat, thvc, dthvdt,
  dthvdh`` (NOT the declaration grouping, which lists ``vpsat`` before
  ``rhov, rhoda``). :func:`sp_humid` returns a
  :class:`~fasst.sp_humid.MoistAirState`, so this transcription maps BY
  NAME and the positional hazard does not apply -- but the doc records
  which local caught which field.
* doc point 8 -- ``sflux``'s ``chw``, ``chw1`` and ``net`` outputs are
  DISCARDED here (received into scratch locals ``t1``/``t3``/``net`` and
  never read). ``isurfg``/``isurff`` are rebuilt from only 7 of the 9
  ``sflux`` flux components (``radswd + radswu + radlwd + radlwu + shw +
  lhw + phw``), i.e. the surface balance *excluding* ground conduction and
  the lumped phase-change term -- those are handled downstream in
  soil_tmp.F90.
* doc point 9 -- ``d1i`` at line 281 (``dense(stempt, 0d0, d1i)``). ``d1i``
  is a reused local phase selector (the fasst_functions "bare integer
  phase" hazard). Traced on both paths: it is deterministically ``0``
  (WATER_VAPOR) at line 281 -- on the ``iheightn <= hm`` path nothing
  reassigns it after the line-50 init; on the other path line 270 sets it
  to ``0`` unconditionally. Safe here, but fragile.
* doc point 10 -- ``sheats``/``spheatr``/``pdens1`` are read at lines
  321-325 while possibly still zero-initialised: if ``int(met[iw, PT])``
  is not in {2, 3, 4} and ``dmet[7] <= 0``, neither precip-type branch
  runs and ``pheatg`` collapses to ``0`` via
  ``spheatr*pdens1*dmet[6] + sheats*pdens*dmet[7]`` (all-zero). Benign
  (no precip => no precip heat, and ``dmet[6/7]`` are ~0 too) but it is
  the sp_humid-point-2 defect class -- a silent zero, not an assert.
* doc point 11 -- ``newdrag``'s convergence loop. The result is assigned
  inside a ``do while(counter <= 5)`` fixed-point iteration on ``ustar``;
  the returned value is the last iteration's. The ``else`` branch (calm +
  neutral, lines 691-693) does ``ustar1 = ustar`` while ``ustar`` is still
  ``0`` from init, then the exit test compares ``0`` to ``0`` -- so that
  branch always exits after one pass. Inert (removing the line changes
  nothing: the exit test still fires), but the statement order is
  suspicious. There is NO non-convergence error flag (unlike
  module_snow's ``bottom``): a run that never converges in 6 passes just
  returns the last value.
* doc point 12 -- dead locals NOT reproduced: ``mixrgr`` (no ``s`` --
  declared, zeroed, never used; distinct from the live ``mixrgrs``), and
  the ``drag``-function locals (whole function omitted). Also note the
  name ``beta`` carries two unrelated meanings in the Fortran: the Lee &
  Pielke wetness factor in ``surfenergy`` (lines 122/173) and
  ``log(ht*f2)`` inside ``newdrag`` (line 534). Both are function-local
  here so there is no collision, but a grep for ``beta`` finds two
  quantities.
* doc point 13 -- the header comment "uses the function:
  dense,spheats,thconds,soilhumid,drag" is wrong on two counts:
  ``soilhumid`` appears only in commented-out lines (170-171) and ``drag``
  is never called. Actually used: ``dense``, ``spheats``, ``thconds``,
  ``newdrag``, plus the subroutines ``sp_humid``, ``lowveg_met``,
  ``sflux``.
* doc point 14 -- ``dmet`` local vs ``DerivedMet``. In the caller
  (new_profile.F90:482-497) ``dmet`` columns 1-5 and 8-13 are copied
  straight from ``dmet1`` (so :class:`~fasst.indices.DerivedMet` names
  apply), but columns 6 and 7 are the *slope-projected* precipitation
  ``delvar1(6/7)*cos(sloper)`` rounded at ``1d15`` -- NOT the raw
  ``dmet1`` columns. A reader mapping ``dmet`` onto ``DerivedMet``
  wholesale gets the two precip columns subtly wrong.
* doc point 15 -- ``chnf`` is the ``fasst_global`` module variable, not a
  local of ``surfenergy``. It is reset to ``0`` only on the bare-soil path
  (line 193). On the vegetated path it is (re)assigned by :func:`newdrag`
  only when ``sigfl > eps and icase in (2, 3) and iheightn > hm``; for
  ``icase == 1`` (still vegetated) with ``iheightn > hm`` the drag blends
  at lines 255/272 read whatever ``state.chnf`` was left at by a PRIOR
  call. Transcribed literally as ``state.chnf`` throughout.
* doc point 16 -- unguarded divisions. ``_fdiv`` (plain IEEE-754, no
  exception) is used for ``dsqrt(dsnow)`` in the ``sfac`` block and
  ``hfol_tot/hm`` on the ``icase == 4`` path, where a zero denominator is
  plausible in normal operation and the Fortran recovers via the
  ``if(t2 > 5d1)`` clamp / a large ``ftemp``. ``vpressa/vpsatgr`` and
  ``soil_moist/porosity`` (``rhsurf`` branch) and every division inside
  :func:`newdrag` are left as plain ``/``: their zero cases are corrupt
  input (zero air pressure, unset soil parameters) outside the validation
  grid. The ``sfac`` snow branch also collapses the Fortran's four
  identical ``if(t2 > 5d1) t2 = 5d1`` / ``sfac = 1 - dexp(-0.8)*dexp(-t2)``
  tails into one -- value-identical.
"""

import math
import warnings
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from .constants import EPS, GRAV, IDENS, KVEG, LHSUB, PI, SDENSW, SIGMA, TREF, VK, VPSAT0
from .functions import Phase, dense, spheats, thconds
from .indices import DerivedMet as D
from .indices import MetCol
from .sflux import SurfaceFluxes, sflux
from .sp_humid import sp_humid
from .state import FasstState

__all__ = ["SurfaceEnergy", "LowVegMetResult", "surfenergy", "newdrag"]

# -- hand-written literals, kept greppable (see module docstring) -------------
_LHEVAP_INTERCEPT = 2500775.6   # J/kg -- WES lhevap fit intercept (line 229); == sflux._LHVAP_LITERAL
_LHEVAP_SLOPE = 2369.729        # J/(kg*K) -- WES lhevap fit slope (line 229)
_STEMPT_FLOOR = 35.9            # K -- magic surface-temperature sanity floor (line 132)
_SNOW_SPHEAT = 2.09e3           # J/(kg*K) -- snow specific heat, hand-written (line 317)

_CMAX = 5                       # newdrag: maximum fixed-point iterations (line 537)
_PRANDTL = 9.5e-1               # newdrag: Prandtl number R (line 538)


def _anint(x: float, p: int) -> float:
    """Fortran ``ANINT(x * 10**p) * 10**-p`` -- round half **away from zero**.

    Mirrors ``fasst/functions.py::_anint`` / ``fasst/sflux.py::_anint``;
    duplicated so each transcribed file stays a self-contained record of
    its own source's idioms. NaN/inf pass through unchanged. Never the
    identity -- the two multiplies are inexact float64 ops.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0 ** p)
    return math.copysign(math.floor(abs(xs) + 0.5), xs) * (10.0 ** -p)


def _aint(x: float, p: int) -> float:
    """Fortran ``AINT(x * 10**p) * 10**-p`` -- truncate **toward zero**.

    ``p = 0`` reproduces a bare ``aint(x)`` (drop the fractional part),
    used to coerce the float precip-type code to an integer at
    surfenergy.F90:310, 314.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0 ** p)
    return math.copysign(math.floor(abs(xs)), xs) * (10.0 ** -p)


def _fdiv(x: float, y: float) -> float:
    """Plain IEEE-754 float division: ``x / 0.0`` yields +/-inf or nan
    instead of raising, matching the Fortran arithmetic at the unguarded
    divisions this file uses where a zero denominator is plausible in
    normal operation -- ``dsqrt(dsnow)`` in the ``sfac`` block (``dsnow``
    is 0 before the snow model has run a metamorphism step) and
    ``hfol_tot/hm`` on the ``icase == 4`` path. Fortran turns those into
    ``+inf`` and the immediately-following ``if(t2 > 5d1) t2 = 5d1`` clamp
    recovers; Python's ``/`` would raise first. Mirrors
    ``fasst/sp_humid.py::_fdiv``. Other unguarded divisions in this file
    (``vpressa/vpsatgr``, ``sm/porosity``, and :func:`newdrag`'s internals)
    are left as plain ``/`` -- their zero cases are corrupt input outside
    the validation grid; see the transcription doc, doc point 16.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


def _round20(x: float) -> float:
    """Fortran ``anint(x*1d20)*1d-20`` -- this file's default idiom."""
    return _anint(x, 20)


def _round18(x: float) -> float:
    """Fortran ``anint(x*1d18)*1d-18`` -- used only on ``isurff`` (line 399)."""
    return _anint(x, 18)


def _round10(x: float) -> float:
    """Fortran ``anint(x*1d10)*1d-10`` -- used only on ``lheatg`` (line 279)."""
    return _anint(x, 10)


def _trunc10(x: float) -> float:
    """Fortran ``aint(x*1d10)*1d-10`` -- used only on ``evapcm`` (line 282)."""
    return _aint(x, 10)


@dataclass
class LowVegMetResult:
    """Expected return shape of the injected ``lowveg_met`` callable.

    Mirrors ``module_lowveg.F90``'s ``lowveg_met`` ``intent(out)`` list
    (``rhoaf, cf, pheatf, rpp, mixra, mixrf, dqdtf``) plus the updated
    ``wetbulba`` (``intent(inout)`` in the Fortran). The callable may also
    mutate the ``dmet`` array it is handed. This class is provided only so
    a caller (or a test stub) has a name to build; :func:`surfenergy` does
    not import ``module_lowveg``.
    """

    rhoaf: float
    cf: float
    pheatf: float
    rpp: float
    mixra: float
    mixrf: float
    dqdtf: float
    wetbulba: float


@dataclass
class SurfaceEnergy:
    """Return value of :func:`surfenergy`.

    The 25 Fortran ``intent(out)`` arguments, plus ``pdens`` (declared
    ``intent(inout)`` but unconditionally zeroed at surfenergy.F90:84 and
    only written afterward, so effectively an output), plus a small tier of
    transcription-time diagnostic intermediates.

    Remember (module docstring, "Not pure"): :func:`surfenergy` ALSO writes
    ``state.ftemp / ptemp / uaf / chnf / meltfl / error_code / error_type``
    -- those are not fields here.
    """

    # -- Fortran intent(out), subroutine-argument order --------------------
    evapcm: float     # surface evaporation/sublimation rate [m/s] -- AINT-truncated at 1e10
    rhsurf: float     # surface relative humidity [0..1]
    pheatf: float     # foliage precipitation heat coefficient [W/m^2/K] (from lowveg_met)
    pheatg: float     # ground precipitation heat coefficient [W/m^2/K]
    mixrf: float      # foliage mixing ratio (from lowveg_met)
    mixraf: float     # blended air/foliage/ground reference mixing ratio
    mixrgrs: float    # ground saturation mixing ratio (sp_humid call 1)
    mixra: float      # air mixing ratio (sp_humid call 2, or lowveg_met)
    lhtf: float       # foliage latent heat transfer coefficient [W/m^2]
    lheatg: float     # ground latent heat transfer coefficient [W/m^2] -- rounded at 1e10
    sheatf: float     # foliage sensible heat transfer coefficient [W/m^2/K]
    sheatg: float     # ground sensible heat transfer coefficient [W/m^2/K]
    dqdtf: float      # d(specific humidity)/dT, foliage (from lowveg_met; 0 on bare path)
    dqdtg: float      # d(specific humidity)/dT, ground (sp_humid call 1)
    sfac: float       # shortwave-penetration fraction [unitless] -- can be < 0, see doc point 4
    d1: float         # latent-term denominator 1 - sigfl*(0.6*(1-rpp) + 0.1*(1-rhsurf))
    rpp: float        # foliage precip fraction (from lowveg_met; 0 on bare path)
    cc1: float        # canopy longwave-exchange coefficient [W/m^2/K^4]
    isurfg: float     # ground surface energy sum (7 sflux terms) [W/m^2] -- rounded at 1e20
    isurff: float     # foliage surface energy sum (7 sflux terms) [W/m^2] -- rounded at 1e18; 0 unless icase in {2,3}
    disurfg: float    # d(isurfg)/dT1  (sflux dnet1, ground call)
    disurff: float    # d(isurff)/dT1  (sflux dnet1, foliage call); 0 unless icase in {2,3}
    disurffg: float   # d(isurff)/dT2  (sflux dnet2, foliage call); 0 unless icase in {2,3}
    disurfgf: float   # d(isurfg)/dT2  (sflux dnet2, ground call)
    lh: float         # latent heat of evaporation or sublimation [J/kg]

    # -- intent(inout), effectively an output ------------------------------
    pdens: float      # precipitation (new-snow) density [kg/m^3]; 0 unless snow

    # -- transcription-time diagnostic additions (genuine Fortran locals
    # feeding the outputs, but never routed out by the Fortran) -----------
    ldragcoeff: float  # latent-flux drag coefficient (feeds lheatg, evapcm)
    sdragcoeff: float  # sensible-flux drag coefficient (feeds sheatg)
    rhoag: float       # ground moist-air density [kg/m^3] = rhov + rhoda (sp_humid call 1)
    rhoaf: float       # foliage air density [kg/m^3] (from lowveg_met; 0 on bare path)
    wetbulba: float    # air wet-bulb temperature [K]; becomes state.ptemp. From sp_humid on
                       # the bare path, from lowveg_met on the vegetated path (doc point 7)


# ===========================================================================
def newdrag(state: FasstState, vflag: int, taf: float, stempt: float, dps: float) -> float:
    """Iterative neutral-to-stable drag coefficient.

    Fortran: surfenergy.F90 lines 522-705, ``function newdrag``. Method:
    Louis (1979) for the first ``ustar``, Li (2010) for ``zeta = z/L``,
    Hogstrom (1995) (unstable) or Beljaars & Holtslag (1991) (stable) for
    the stability functions, Li (2012) for the final coefficient, iterated
    up to ``_CMAX`` times on ``ustar``.

    Args:
        state: reads ``iheightn, hm, uaf, rough, z0l`` (+ constants
            ``VK, GRAV, PI``).
        vflag: 0 -> ground roughness (``rough`` if ``hm <= eps`` else
            ``7.775e-3``); anything else -> foliage roughness ``z0l``.
        taf, stempt: reference and surface temperature [K] for the bulk
            Richardson number.
        dps: molecular diffusion coefficient (``dv`` for latent, ``dh`` for
            sensible) used to set the scalar roughness ``z0p``.

    Returns:
        The drag coefficient (the last fixed-point iterate). See module
        docstring, doc point 11, on the loop's quirks.
    """
    ht = state.iheightn - state.hm
    ug = state.uaf
    if ug < 5e-2:
        ug = 0.0

    if vflag == 0:
        z0g = state.rough if state.hm <= EPS else 7.775e-3
    else:
        z0g = state.z0l

    # bulk Richardson number: < 0 unstable, 0 neutral, > 0 stable
    if abs(taf - stempt) <= EPS or ug <= EPS:
        rib = 0.0
    else:
        rib = 2.0 * GRAV * ht * (taf - stempt) / (ug * ug * (taf + stempt))

    f1 = 1.0 / z0g
    cdng0 = VK * VK / (math.log(ht * f1) * math.log(ht * f1))

    if rib < 0.0:  # unstable
        c = 7.4 * cdng0 * 9.4 * math.sqrt(ht * f1)
        gammam = 1.0 - 9.4 * rib / (1.0 + c * math.sqrt(abs(rib)))
    elif abs(rib) < EPS:  # neutral
        gammam = 1.0
    else:  # rib > eps, stable
        gammam = 1.0 / ((1.0 + 4.7 * rib) * (1.0 + 4.7 * rib))

    ustar1 = math.sqrt(cdng0 * gammam) * ug

    ustar = 0.0
    result = 0.0
    counter = 0
    while counter <= _CMAX:
        if ustar1 <= EPS:
            z0p = z0g
        else:
            z0p = dps / (4e-1 * ustar1)  # FAM2dEd
            if z0p >= z0g:
                z0p = 0.9 * z0g
        f2 = 1.0 / z0p

        alpha = math.log(ht * f1)
        beta = math.log(ht * f2)

        psim1 = psim2 = psih1 = psih2 = 0.0
        if rib <= 0.0:  # unstable and neutral
            if rib < -5.0:  # unstable, Li eqn. breaks down -> interpolate
                rib2 = -2.0
                zeta2 = 4.5e-2 * alpha * rib2 * rib2 + rib2 * (
                    (3e-3 * beta + 5.9e-3) * alpha * alpha
                    + (-8.28e-2 * beta + 8.845e-1) * alpha
                    + (1.739e-1 * beta * beta - 9.213e-1 * beta - 1.057e-1)
                )
                rib5 = -5.0
                zeta5 = 4.5e-2 * alpha * rib5 * rib5 + rib5 * (
                    (3e-3 * beta + 5.9e-3) * alpha * alpha
                    + (-8.28e-2 * beta + 8.845e-1) * alpha
                    + (1.739e-1 * beta * beta - 9.213e-1 * beta - 1.057e-1)
                )
                zeta = ((zeta5 - zeta2) / (rib5 - rib2)) * (rib + rib2) + zeta2
            elif rib < 0.0 and rib >= -5.0:  # unstable, Li eqn. works
                zeta = 4.5e-2 * alpha * rib * rib + rib * (
                    (3e-3 * beta + 5.9e-3) * alpha * alpha
                    + (-8.28e-2 * beta + 8.845e-1) * alpha
                    + (1.739e-1 * beta * beta - 9.213e-1 * beta - 1.057e-1)
                )
            else:  # neutral (rib == 0)
                zeta = 0.0

            # Hogstrom (1995)
            zeta0 = z0p * zeta / ht
            x1 = max(0.0, (1.0 - 1.9e1 * zeta)) ** 2.5e-1
            psim1 = math.log((1.0 + x1 * x1) * (1.0 + x1) * (1.0 + x1) / 8.0) - 2.0 * math.atan(x1) + PI * 5e-1
            x2 = max(0.0, (1.0 - 1.9e1 * zeta0)) ** 2.5e-1
            psim2 = math.log((1.0 + x2 * x2) * (1.0 + x2) * (1.0 + x2) / 8.0) - 2.0 * math.atan(x2) + PI * 5e-1

            y1 = max(0.0, (1.0 - 1.16e1 * zeta)) ** 5e-1
            psih1 = 2.0 * math.log((1.0 + y1) * 5e-1)
            y2 = max(0.0, (1.0 - 1.16e1 * zeta0)) ** 5e-1
            psih2 = math.log((1.0 + y2) * 5e-1)

        elif rib > 0.0:  # stable
            if rib <= 2e-1:  # weakly stable
                zeta = rib * rib * (
                    (5.738e-1 * beta - 4.399e-1) * alpha + (-4.901e0 * beta + 5.25e1)
                ) + rib * (
                    (-5.39e-2 * beta + 1.54e0) * alpha + (-6.69e-1 * beta - 3.282e0)
                )
                zeta = max(0.0, zeta)
            else:  # strongly stable
                zeta = (7.529e-1 * alpha + 1.494e1) * rib + 1.569e-1 * alpha - 3.091e-1 * beta - 1.303e0

            # Beljaars & Holtslag (1991)
            zeta0 = z0p * zeta / ht
            a = 1.0
            b = 6.667e-1
            c = 5.0
            d = 3.5e-1
            psim1 = -a * zeta - b * (zeta - c / d) * math.exp(-d * zeta) - b * c / d
            psim2 = -a * zeta0 - b * (zeta0 - c / d) * math.exp(-d * zeta0) - b * c / d
            psih1 = -(1.0 + 2.0 * a * zeta / 3.0) ** 1.5 - b * (zeta - c / d) * math.exp(-d * zeta) - b * c / d + 1.0
            psih2 = -(1.0 + 2.0 * a * zeta0 / 3.0) ** 1.5 - b * (zeta0 - c / d) * math.exp(-d * zeta0) - b * c / d + 1.0

        # Li (2012): ustar and the coefficient
        result = 0.0
        if ug > 1e-3 or abs(rib) >= EPS:
            ustar = max(0.0, ug * VK / (math.log(ht * f2) - psim1 + psim2))
            result = max(0.0, (ustar / ug) * (VK / _PRANDTL) / (math.log(ht * f2) - psih1 + psih2))
        else:
            ustar1 = ustar  # doc point 11: ustar is still 0.0 here on the first pass
            result = (VK / _PRANDTL) / (math.log(ht * f2) - psih1 + psih2)
        result = _round20(result)

        if abs(ustar - ustar1) <= 1e-3:
            counter = _CMAX + 1
        else:
            counter = counter + 1
            ustar1 = ustar

    return result


# ===========================================================================
def surfenergy(
    state: FasstState,
    ii: int,
    iter_: int,
    sn: int,
    kave: float,
    sthick: float,
    pdens: float,
    dmet,
    *,
    lowveg_met: Optional[Callable[..., LowVegMetResult]] = None,
) -> SurfaceEnergy:
    """Surface energy-balance coefficients and sums for one node.

    Fortran: surfenergy.F90 lines 1-406, ``subroutine surfenergy``. Calls
    ``sp_humid`` (x2), ``lowveg_met`` (vegetated path only), ``sflux``
    (x1 or x2); uses ``dense``, ``spheats``, ``thconds``, ``newdrag``.

    Args:
        state: model state -- see module docstring for the fields read and
            **written**.
        ii, iter_: pass-through indices for ``lowveg_met`` (Fortran ``ii``,
            ``iter``). Not otherwise used.
        sn: surface node index. ``stempt = state.stt[sn]``; ``sflux`` is
            handed ``state.stt[sn - 1]`` as its ``temp2``.
        kave: average thermal conductivity across the ground conduction
            path [W/m/K] -- passed to the ground ``sflux`` call as ``kth``.
        sthick: ground conduction path length [m] -- ``sflux`` ``h``.
        pdens: incoming value ignored (the Fortran zeroes it immediately);
            returned in :class:`SurfaceEnergy`.
        dmet: 13-element derived-meteorology array, 1-based (phantom slot
            0), :class:`~fasst.indices.DerivedMet` layout EXCEPT columns
            6/7 (slope-projected precip -- doc point 14). May be mutated in
            place by ``lowveg_met``.
        lowveg_met: injected callable, required only on the
            ``sigfl > eps and icase != 4`` path. See module docstring for
            the contract. ``None`` raises :class:`NotImplementedError` if
            that path is reached.

    Returns:
        A :class:`SurfaceEnergy`. Also mutates ``state`` (docstring).
    """
    icase = state.icase
    nnodes = state.nnodes
    iw = state.iw

    # Fortran zeroes these near the top (smr line 121, taf line 101); every
    # branch below reassigns them, but the pre-init keeps a later edit from
    # turning a missed branch into an UnboundLocalError where Fortran would
    # silently use 0d0 (cf. sflux.py's radlwu).
    smr = 0.0
    taf = 0.0

    # -- surface temperature, with the 35.9 K sanity floor (doc point 3) --
    stempt = state.stt[sn]
    if stempt < _STEMPT_FLOOR:
        if state.single_multi_flag == 0:
            day = hour = -1
            if state.met is not None:
                day = int(state.met[iw, MetCol.DOY])
                hour = int(state.met[iw, MetCol.HR])
            warnings.warn(
                f"freq_id {state.freq_id} toptemp = {stempt:.3f} airtemp = {dmet[D.TAIR_K]:.3f} "
                f"Top temp too low in surfenergy; day {day} hour {hour}",
                RuntimeWarning,
                stacklevel=2,
            )
        stempt = dmet[D.TAIR_K]
        state.error_code = 1
        state.error_type = 3

    # -- sp_humid call 1: GROUND saturation quantities (rh2 = 1) ----------
    # Fortran positional out -> local: mixrgrs<-mixr, dqdtg<-dmrdt,
    # (t1<-vpress, t2<-wetbulb discarded), rhov<-rhov, rhoda<-rhoda,
    # vpsatgr<-vpsat, (t3/t4/t5<-thvc/dthvdt/dthvdh discarded). doc point 7.
    f1i = 1 if state.hm > EPS else 0
    g = sp_humid(f1i, dmet[D.AIR_PRESS], stempt, 1.0, state.phead[nnodes])
    mixrgrs = g.mixr
    dqdtg = g.dmrdt
    vpsatgr = g.vpsat
    rhoag = g.rhov + g.rhoda  # kg/m^3

    # -- sp_humid call 2: AIR quantities (rh2 = dmet[5]% / 100) ----------
    # RH_PCT -> fraction: this is the unit trap flagged in fasst/indices.py;
    # the Fortran gets it right here (`c2 = dmet(5)*1d-2`).
    a = sp_humid(f1i, dmet[D.AIR_PRESS], dmet[D.TAIR_K], dmet[D.RH_PCT] * 1e-2, 0.0)
    mixra = a.mixr
    vpressa = a.vpress
    wetbulba = a.wetbulb
    vpsata = a.vpsat

    # -- surface relative humidity rhsurf, and smr -----------------------
    node_top = state.node_type[nnodes]
    if state.hm > EPS or node_top in ("SN", "WA"):
        rhsurf = 1.0
        smr = 1.0
    elif state.hpond > EPS or state.met[iw, MetCol.PREC] + state.met[iw, MetCol.PREC2] > EPS:
        rhsurf = 1.0
        smr = 1.0
    elif state.stt[nnodes] <= TREF or state.ice[nnodes] > EPS:
        rhsurf = 1.0
        smr = 1.0
    else:
        # Lee & Pielke (1992). Commented Maneta & Silverman variant differs
        # in parenthesisation -- doc point 5.
        beta = min(
            1.0,
            2.5e-1 * (1.0 - math.cos(PI * state.soil_moist[nnodes] / state.nsoilp[nnodes, 9])) ** 2.0,
        )
        rhsurf = beta + (1.0 - beta) * vpressa / vpsatgr
        smr = state.soil_moist[nnodes] / state.nsoilp[nnodes, 9]
    if state.vegl_type == 8:
        rhsurf = smr

    # -- foliage reference temperature ftemp, taf, and the vegetated path -
    # NOTE: `chnf` in the Fortran is the module global (fasst_global), never
    # a local -- so it is `state.chnf` here, and it is only RESET to 0 on
    # the bare-soil path. On the vegetated path with icase == 1 (or with
    # `iheightn <= hm`) it keeps its prior value, which then feeds the drag
    # blends at lines 255/272. See doc point 15.
    epf1 = 0.0
    rpp = 0.0
    sheatf = 0.0
    pheatf = 0.0
    lhtf = 0.0
    cc1 = 0.0
    mixrf = 0.0
    rhoaf = 0.0
    dqdtf = 0.0
    cf = 0.0

    if state.sigfl <= EPS or icase == 4:
        state.ftemp = dmet[D.TAIR_K]  # air temperature
        if icase == 4:
            state.ftemp = state.stt[nnodes] + _fdiv(state.hfol_tot, state.hm) * (
                state.stt[nnodes + 1] - state.stt[nnodes]
            )
            state.ftemp = _round20(state.ftemp)
        taf = dmet[D.TAIR_K]
        state.uaf = dmet[D.WIND]
        state.chnf = 0.0
        d1 = 1.0
        mixraf = mixra
    else:
        if icase == 2:
            state.ftemp = state.stt[nnodes + 1]
        elif icase == 3:
            state.ftemp = state.stt[nnodes + 2]
        epf1 = state.epf
        taf = (1.0 - 7e-1 * state.sigfl) * dmet[D.TAIR_K] + state.sigfl * (
            6e-1 * state.ftemp + 1e-1 * stempt
        )  # foliage temp at the ground surface [K]

        if lowveg_met is None:
            raise NotImplementedError(
                "surfenergy's vegetated path (sigfl > eps and icase != 4) needs the "
                "injected `lowveg_met` callable; module_lowveg.F90 is not transcribed. "
                "See fasst/surfenergy.py module docstring."
            )
        lv = lowveg_met(
            state, ii, iter_,
            int(state.met[iw, MetCol.PT]), int(state.met[iw, MetCol.PT2]),
            mixrgrs, rhsurf, taf, dmet, wetbulba,
        )
        rhoaf, cf, pheatf, rpp = lv.rhoaf, lv.cf, lv.pheatf, lv.rpp
        mixra, mixrf, dqdtf, wetbulba = lv.mixra, lv.mixrf, lv.dqdtf, lv.wetbulba

        eps1 = epf1 + state.emis - epf1 * state.emis
        if eps1 > EPS:
            cc1 = epf1 * state.emis * state.sigfl * SIGMA / eps1

        d1 = 1.0 - state.sigfl * (6e-1 * (1.0 - rpp) + 1e-1 * (1.0 - rhsurf))
        mixraf = (
            (1.0 - 7e-1 * state.sigfl) * mixra
            + 6e-1 * state.sigfl * rpp * mixrf
            + 1e-1 * state.sigfl * rhsurf * mixrgrs
        )
        mixraf = mixraf / d1

    d1 = _round20(d1)
    cc1 = _round20(cc1)
    mixraf = _round20(mixraf)

    # -- latent heat of evaporation / sublimation ------------------------
    state.meltfl = "m"
    lh = _LHEVAP_INTERCEPT - _LHEVAP_SLOPE * (dmet[D.TAIR_K] - TREF)  # J/kg
    if vpressa <= VPSAT0 and rhsurf * vpsatgr <= VPSAT0:
        if state.uaf < 0.5 and state.tmelt[iw] < TREF:
            state.meltfl = "s"
            lh = LHSUB

    # -- drag coefficients ----------------------------------------------------
    dv = 0.0
    if state.iheightn <= state.hm or state.iheightn <= EPS:
        sdragcoeff = 2e-3 + 6e-3 * (state.elev * 2e-4)  # unitless
        if state.hsaccum > EPS or state.hi > EPS:
            sdragcoeff = 1.35e-3  # Rachel
        ldragcoeff = sdragcoeff
    else:
        # dv etc. from Jacobson, "Fundamentals of Atmospheric Modeling" 2nd ed.
        if dmet[D.TAIR_K] > EPS:
            dv = 2.11e-5 * ((dmet[D.TAIR_K] / TREF) ** 1.94) * (1013.25 / dmet[D.AIR_PRESS])

        if state.sigfl > EPS and icase in (2, 3):
            state.chnf = newdrag(state, 1, dmet[D.TAIR_K], state.ftemp, dv)
        ldragcoeff = newdrag(state, 0, dmet[D.TAIR_K], stempt, dv)
        ldragcoeff = max(0.0, (1.0 - state.sigfl) * ldragcoeff + state.sigfl * state.chnf)

        kv = thconds(dmet[D.TAIR_K], Phase.WATER_VAPOR)  # W/m/K, water vapour
        kd = thconds(dmet[D.TAIR_K], Phase.DRY_AIR)      # W/m/K, dry air
        ka = kd * (1.0 - (1.17 - 1.02 * kv / kd) * dmet[D.RH_PCT] * 1e-2)  # moist-air k
        dh = ka / (rhoag * ((1.0 + 0.87 * mixra) * spheats(dmet[D.TAIR_K], Phase.DRY_AIR)))

        if state.sigfl > EPS and icase in (2, 3):
            state.chnf = newdrag(state, 1, dmet[D.TAIR_K], state.ftemp, dh)
        sdragcoeff = newdrag(state, 0, dmet[D.TAIR_K], stempt, dh)
        sdragcoeff = max(0.0, (1.0 - state.sigfl) * sdragcoeff + state.sigfl * state.chnf)

    # -- latent heat loss/gain, surface evaporation/sublimation -----------
    lheatg = ldragcoeff * lh * state.uaf * rhoag
    lheatg = _round10(lheatg)
    # d1i == 0 (WATER_VAPOR) deterministically here -- doc point 9
    evapcm = lheatg * (mixraf - rhsurf * mixrgrs) / (dense(stempt, 0.0, Phase.WATER_VAPOR) * lh)  # m/s
    evapcm = _trunc10(evapcm)

    lhtf = state.lail * cf * lh * state.uaf * rhoaf  # W/m^2
    lhtf = _round20(lhtf)

    # -- sensible heat loss/gain ----------------------------------------------
    spheat = spheats(taf, Phase.DRY_AIR)

    shmax = 2.0 if state.hm > EPS else 0.0
    sheatg = max(shmax, sdragcoeff * spheat * state.uaf * rhoag)  # W/m^2/K
    sheatg = _round20(sheatg)

    if abs(state.lail) <= EPS:
        sheatf = 0.0
    else:
        shmax = 2.0 if state.hm > EPS else 0.0
        sheatf = max(shmax, 1.1 * state.lail * cf * spheat * state.uaf * rhoaf)  # foliage W/m^2/K
    sheatf = _round20(sheatf)

    # -- heat loss due to precipitation (Jordan, CRREL Rep. 91-16) --------
    pdens = 0.0
    pdens1 = 0.0
    spheatr = 0.0
    sheats = 0.0
    pt_code = _aint(state.met[iw, MetCol.PT], 0)
    if pt_code in (2.0, 4.0):  # rain, freezing rain
        pdens1 = dense(wetbulba, 0.0, Phase.WATER)  # kg/m^3
        spheatr = spheats(wetbulba, Phase.WATER)    # J/kg/K
    elif pt_code == 3.0 or dmet[D.SNOW_M] > 0.0:  # snow
        pdens = dense(wetbulba, dmet[D.WIND], Phase.SNOW)
        sheats = _SNOW_SPHEAT  # J/kg/K
    state.ptemp = wetbulba

    if pdens > EPS and dmet[D.SNOW_M] <= EPS:
        pheatg = sheats * pdens * dmet[D.PRECIP_M]
    else:
        # doc point 10: all three of spheatr/pdens1/sheats can still be 0 here
        pheatg = spheatr * pdens1 * dmet[D.PRECIP_M] + sheats * pdens * dmet[D.SNOW_M]

    pheatg = (1.0 - state.sigfl) * pheatg
    pheatg = _round20(pheatg)

    # -- shortwave-penetration fraction sfac -----------------------------
    if state.hm <= EPS:  # no snow or ice
        sfac = 1.0
    elif state.hsaccum > EPS or state.newsd > EPS:  # snow
        if state.newsd <= EPS:
            if iw >= 2:
                t2 = _fdiv(state.hsaccum * state.sdens[iw - 1] * 3.795e-3, math.sqrt(state.dsnow))
            else:
                t2 = state.hsaccum * SDENSW * 3.795e-3 / math.sqrt(1e-3)
        else:
            if iw >= 2:
                t2 = _fdiv(state.newsd * state.sdens[iw - 1] * 3.795e-3, math.sqrt(state.dsnow))
            else:
                t2 = state.hsaccum * SDENSW * 3.795e-3 / math.sqrt(1e-3)
        if t2 > 5e1:
            t2 = 5e1
        sfac = 1.0 - math.exp(-0.8) * math.exp(-t2)
    elif state.hi > EPS and state.hsaccum <= EPS and state.newsd <= EPS:  # ice, no snow
        t2 = 0.8 + state.hi * IDENS * 3.795e-3 / math.sqrt(1e-3)
        if t2 > 5e1:
            t2 = 5e1
        sfac = 1.0 - math.exp(-t2)
    elif node_top == "WA":
        if state.ice[nnodes] > EPS:
            t2 = 0.8 + state.met[iw, MetCol.HI] * IDENS * 3.795e-3 / math.sqrt(1e-3)
            if t2 > 5e1:
                t2 = 5e1
            sfac = 1.0 - math.exp(-t2)
        else:
            sfac = 1.0 - math.exp(0.4)  # doc point 4: POSITIVE exponent -> sfac ~ -0.492
    else:
        sfac = 0.0  # Fortran: sfac keeps its 0d0 initializer on this fall-through
    sfac = _round20(sfac)

    # -- solar + incoming IR terms via sflux ----------------------------
    # sflux out -> local: radsd<-radswd, radsu<-radswu, radld<-radlwd,
    # radlu<-radlwu, pht1<-phw, sh<-shw, lhw<-lhw, (chw/chw1/net discarded),
    # disurfg<-dnet1, disurfgf<-dnet2. doc point 8.
    fg = sflux(
        state, 0, 1, sn,
        dmet[D.UPSOL], dmet[D.TSOL_NET], dmet[D.IR_UP], dmet[D.IR_DOWN], sfac,
        stempt, cc1, taf, pheatg, sheatg, lheatg, mixraf, mixrgrs,
        kave, state.stt[sn - 1], sthick, rhsurf, rpp, state.ice[nnodes],
        state.wvc[nnodes], state.soil_moist[nnodes], dqdtg, dqdtf, d1,
    )
    isurfg = fg.radswd + fg.radswu + fg.radlwd + fg.radlwu + fg.shw + fg.lhw + fg.phw
    isurfg = _round20(isurfg)
    disurfg = fg.dnet1
    disurfgf = fg.dnet2

    if icase in (2, 3):
        ff = sflux(
            state, 1, 1, sn,
            dmet[D.UPSOL], dmet[D.TSOL_NET], dmet[D.IR_UP], dmet[D.IR_DOWN], 1.0,
            state.ftemp, cc1, taf, pheatf, sheatf, lhtf, mixraf, mixrf,
            KVEG, stempt, state.hfol_tot, rpp, rhsurf, state.ice[nnodes],
            state.wvc[nnodes], state.soil_moist[nnodes], dqdtf, dqdtg, d1,
        )
        isurff = ff.radswd + ff.radswu + ff.radlwd + ff.radlwu + ff.shw + ff.lhw + ff.phw
        isurff = _round18(isurff)
        disurff = ff.dnet1
        disurffg = ff.dnet2
    else:
        isurff = 0.0
        disurff = 0.0
        disurffg = 0.0

    # float(...) strips the numpy-scalar wrapper that `dmet[...]` /
    # `state.<array>[...]` reads introduce (float64 -> float is the identity
    # on the value); the _round* helpers already return plain float.
    return SurfaceEnergy(
        evapcm=float(evapcm),
        rhsurf=float(rhsurf),
        pheatf=float(pheatf),
        pheatg=float(pheatg),
        mixrf=float(mixrf),
        mixraf=float(mixraf),
        mixrgrs=float(mixrgrs),
        mixra=float(mixra),
        lhtf=float(lhtf),
        lheatg=float(lheatg),
        sheatf=float(sheatf),
        sheatg=float(sheatg),
        dqdtf=float(dqdtf),
        dqdtg=float(dqdtg),
        sfac=float(sfac),
        d1=float(d1),
        rpp=float(rpp),
        cc1=float(cc1),
        isurfg=float(isurfg),
        isurff=float(isurff),
        disurfg=float(disurfg),
        disurff=float(disurff),
        disurffg=float(disurffg),
        disurfgf=float(disurfgf),
        lh=float(lh),
        pdens=float(pdens),
        ldragcoeff=float(ldragcoeff),
        sdragcoeff=float(sdragcoeff),
        rhoag=float(rhoag),
        rhoaf=float(rhoaf),
        wetbulba=float(wetbulba),
    )
