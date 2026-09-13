"""Kinematic melt-water routing through a homogeneous snow pack.

Fortran source: module_snow.F90 (88 670 bytes, 2 268 lines, 12 subroutines:
1 public entry point + 11 helpers inside the module ``contains``). See
docs_transcricao/2026-09-01-doc-transcricao-module_snow.docx for the full
transcription rationale, the numeric verification of every suspected
defect, and the pending validation work referenced throughout this module
as "the doc".

What it does
------------
Each time step, :func:`snow_model` (1) updates snow grain geometry
(dendricity, sphericity, grain diameter), density and effective porosity
(:func:`inivarivals`); (2) computes top melt (radiation/turbulence), bottom
melt (soil heat flux), wind-drift loss and precipitation gain
(:func:`getwave`); (3) packages the step's liquid water as a new "wave"
propagating vertically by a power-law permeability/saturation relation
(kinematic gravity flow, Colbeck 1972 / Jordan 1991); and (4) decides
whether that wave collides with earlier waves, partly refreezes, or leaves
the base of the pack as runoff (:func:`shoveout`, :func:`satvolume`).
Algorithm attributed in the Fortran header to Mary Albert / G. Krajeski
(BASIC) and G. Koenig (Fortran, 1996); grain metamorphism and wind drift
cite Jordan (SNTHERM, 1991), Jordan & Andreas (1999) and Crocus (Vionnet
et al., 2012).

Structure of this transcription
-------------------------------
The Fortran keeps every ``save::`` variable of ``snow_model`` (and two more
in helpers, see below) as module-scoped persistent state, and the eleven
helpers receive most of that state as ``intent(inout)`` arguments -- i.e.
the arguments are *aliases* of the saved variables. That literal semantics
is reproduced here as a :class:`SnowState` dataclass whose fields the
helpers **mutate in place**; only genuine ``intent(in)``/``intent(out)``
locals are passed/returned by value (``bottom`` -> ``bot``, ``depth`` ->
``d``, ``predsndp`` -> ``(sumf, d1)``, ``getwave`` -> the five melt/gain
terms, ``shoveout`` -> ``v``). All twelve functions keep their Fortran
names and 1:1 correspondence so the doc's per-line validation stays
possible. Module-global reads (``met``, ``dmet1``, ``stt``, ``toptemp``,
``ptemp``, ``sloper``, ...) go through an explicit :class:`~fasst.state.FasstState`
argument, exactly as in fasst/functions.py and fasst/sp_humid.py.

External calls: ``dense`` and ``spheats`` from fasst_functions.F90 (the
Fortran declares them as bare ``real(dp)`` locals -- the F77 idiom already
documented for sp_humid.F90 -- and here they are simply imported).

Persistence census -- 29 persistent values, only 17 survive a checkpoint
-----------------------------------------------------------------------
``snow_model`` has 27 ``save::`` names (26 scalars + ``wave(100,5)``).
Two more live in helpers and are *not* in ``snow_model``'s zero-out block:
``dcold`` (:func:`inivarivals`, live state -- read at Fortran line 577,
written 586) and ``newdense`` (:func:`predsndp`, inert today: only read
inside ``case(1)`` after a fresh assignment). All 29 are :class:`SnowState`
fields.

When ``infer_test == 1`` (restart from a checkpoint, fasst_driver.F90),
``snow_model`` restores 15 scalars + 2 integers from the global
``sn_stat(15)`` / ``sn_istat(2)`` arrays. The mapping (this closes the
fasst_global.F90 doc's open "catalog sn_stat/sn_istat" item):

    sn_stat : 1 sdold  2 sd  3 rho  4 maxwet  5 maxwetold  6 vd  7 vw
              8 tottim  9 temphist  10 phie  11 swi  12 swe  13 rain2
              14 dend  15 sph
    sn_istat: 1 tims   2 smwaves

The other 12 persistent values (``n, xf, a, at, waterheld, kappa, wave,
d, k, timestp`` + ``dcold, newdense`` -- the "in-transit" state: the water
waves and the current permeability) are **not** restored; they are
re-zeroed in the same ``if(first_snow == 0)`` block. This is a deliberate
partial restore, replicated verbatim (doc point 1). Consequence: a snow
pack resumed from a checkpoint reappears with correct depth / density /
temperature history but an *empty* wave array, ``kappa = 0`` and no free
water -- physically inconsistent with an evolved pack. A
``SnowState`` seeded from ``sn_stat``/``sn_istat`` is therefore NOT
equivalent to an uninterrupted run; see PENDENCIAS in the doc (open
decision: widen the checkpoint format to cover the missing 12).

Numeric-idiom helpers
---------------------
* :func:`_anint` / :func:`_aint` -- Fortran ``ANINT`` (round half away from
  zero -- NOT numpy's round-half-to-even) and ``AINT`` (truncate toward
  zero), each with an explicit power-of-ten scale. NaN/inf pass through
  unchanged, matching Fortran.
* :func:`_fdiv` / :func:`_fpow` / :func:`_fsqrt` -- plain IEEE-754
  semantics (silent ``inf``/``nan``), extending the ``_fdiv`` precedent
  from fasst/sp_humid.py. Needed because Python ``x / 0.0`` raises,
  ``(-0.3) ** 0.7`` returns a *complex* number (poisoning downstream
  comparisons), and ``math.sqrt(-1)`` raises -- where gfortran yields
  ``nan``/``inf`` and keeps going. Every ``**`` and ``dsqrt`` in the wave
  physics (``bottom``, ``depth``, ``collide``, ``shoveout``, ``satvolume``)
  goes through them, because Newton overshoot routinely drives ``tau``
  negative into a real exponent.
* ``f1`` sign trap: :func:`bottom` sets ``f1 = 1/(1-n)`` (negative);
  :func:`depth` and :func:`collide` set ``f1 = 1/(n-1)`` (positive). Same
  name, same ``n``, inverted. Not factored into a shared helper.

I/O: this module performs NO file I/O. Fortran ``snow_model`` writes a
diagnostic line to unit 55 when ``sprint == 1`` (lines 236-246); ``dates``
and ``wet`` are computed only for that write, and line 240 is already
commented out upstream. The ``sprint`` parameter is kept for signature
parity and the write is skipped.

``dmet1`` columns are used by raw index with an inferred-meaning comment
(3 = wind speed, 4 = air temperature [K], 6 = rain SWE [cm-equiv], 7 = snow
SWE) -- cataloguing ``dmet1`` as an enum is left open pending
fasst_main.F90 (doc PENDENCIAS).
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import (
    EPS,
    GRAV,
    HT_MINM,
    ITHCOND,
    LHFUS,
    LHSUB,
    PI,
    SDENSD,
    SDENSW,
    SPFLAG,
    TREF,
)
from .functions import Phase, dense, spheats
from .indices import MetCol
from .state import FasstState

__all__ = ["SnowState", "snow_model", "PERMEABILITY_EXPONENT"]

# doc point 12: only value `n` is ever given. Documented calibration range 2.7-3.3.
PERMEABILITY_EXPONENT = 3.3

# doc point 5: von Karman constant LOCAL to getwave (Fortran line 653).
# Diverges from fasst.constants.VK (0.40, the value the rest of the project
# uses); (0.35/0.40)**2 = 0.7656, so Cdng0 is 23.4 % lower here.
VK_LOCAL = 0.35

# Fortran ip_* met-column pointers (fasst_global.F90) == fasst.indices.MetCol.
IP_YEAR = MetCol.YEAR
IP_DOY = MetCol.DOY
IP_HR = MetCol.HR
IP_MIN = MetCol.MIN
IP_AP = MetCol.AP
IP_WS = MetCol.WS
IP_PT = MetCol.PT
IP_PT2 = MetCol.PT2
IP_SD = MetCol.SD
IP_TSOIL = MetCol.TSOIL
IP_HI = MetCol.HI


# Newton-iteration census for `bottom`, split by branch (*1 = no-frozen-water
# loop, *2 = frozen-water loop). Not part of the model: read + reset by
# tests/ref/check_snow_ref.py.
_BOTTOM_CENSUS = {"iters1": 0, "calls1": 0, "iters2": 0, "calls2": 0}


def _reset_bottom_census() -> None:
    for k in _BOTTOM_CENSUS:
        _BOTTOM_CENSUS[k] = 0


# ---------------------------------------------------------------------------
# numeric-idiom helpers
# ---------------------------------------------------------------------------
def _anint(x: float, p: int = 0) -> float:
    """Fortran ``ANINT(x*10**p)*10**-p`` -- round half *away from zero*.

    NOT ``np.rint`` (which rounds half to even); at the 1e10/1e15 scales
    used in this file the difference is reachable. NaN/inf pass through
    unchanged, as ``ANINT`` does.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def _aint(x: float, p: int = 0) -> float:
    """Fortran ``AINT(x*10**p)*10**-p`` -- truncate toward zero.

    NaN/inf pass through unchanged, as ``AINT`` does.
    """
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _fdiv(x: float, y: float) -> float:
    """Plain IEEE-754 float division: ``x / 0.0`` -> +-inf / nan, no raise.

    See fasst/sp_humid.py ``_fdiv``; duplicated so each transcribed file is
    a self-contained record of its own source's idioms.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


def _fpow(x: float, y: float) -> float:
    """Plain IEEE-754 real power: negative base with non-integer exponent
    -> ``nan`` (matching gfortran ``x**y``), ``0.0**negative`` -> ``inf``.

    Python's ``**`` returns a *complex* number for ``(-0.3) ** 0.7``, which
    silently poisons the ``>=`` comparisons in the Newton loops; this
    keeps the arithmetic real, with ``nan`` ending the iteration exactly
    as a Fortran ``nan`` would (``nan >= x`` is false).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float_power(np.float64(x), np.float64(y)))


def _fsqrt(x: float) -> float:
    """Plain IEEE-754 sqrt: ``sqrt(-1)`` -> ``nan`` instead of raising,
    matching gfortran ``dsqrt`` of a negative argument.
    """
    with np.errstate(invalid="ignore"):
        return float(np.sqrt(np.float64(x)))


# ---------------------------------------------------------------------------
# persistent state
# ---------------------------------------------------------------------------
@dataclass
class SnowState:
    """One snow pack's worth of Fortran ``save::`` state (module_snow.F90).

    Instantiate one per grid point. The eleven helper functions mutate this
    in place -- that is the literal meaning of the Fortran ``intent(inout)``
    arguments, which alias the saved module variables.

    Fields are grouped by whether ``snow_model`` restores them from the
    global ``sn_stat`` / ``sn_istat`` checkpoint arrays when
    ``state.infer_test == 1``. Seeding this object from a checkpoint is NOT
    equivalent to continuing an uninterrupted run: the second group (the
    water waves, permeability and free water) is always re-zeroed, never
    restored. See the module docstring, "Persistence census".
    """

    # -- restored from sn_stat / sn_istat on checkpoint (doc point 1) --------
    sdold: float = 0.0        # sn_stat(1)  previous-step snow depth (cm)
    sd: float = 0.0           # sn_stat(2)  snow depth (cm)
    rho: float = 0.0          # sn_stat(3)  ice density in the pack (g/cm^3)
    maxwet: float = 0.0       # sn_stat(4)  depth of wet snow (cm)
    maxwetold: float = 0.0    # sn_stat(5)  previous-step wet-snow depth (cm)
    vd: float = 0.0           # sn_stat(6)  mean crystal volume, dry part (cm^3)
    vw: float = 0.0           # sn_stat(7)  mean crystal volume, wet part (cm^3)
    tottim: float = 0.0       # sn_stat(8)  total model time (s)
    temphist: float = 0.0     # sn_stat(9)  temperature history for refreeze (C)
    phie: float = 0.0         # sn_stat(10) effective porosity (unitless)
    swi: float = 0.0          # sn_stat(11) irreducible water saturation (unitless)
    swe: float = 0.0          # sn_stat(12) snow water equivalent (cm)
    rain2: float = 0.0        # sn_stat(13) accumulated non-runoff rain (cm)
    dend: float = 0.0         # sn_stat(14) grain dendricity (unitless)
    sph: float = 0.0          # sn_stat(15) grain sphericity (unitless)
    tims: int = 0             # sn_istat(1) number of steps in refreezing
    smwaves: int = 0          # sn_istat(2) number of active melt waves

    # -- NOT restored on checkpoint: the "in-transit" state -----------------
    n: float = 0.0            # permeability/saturation exponent (always 3.3)
    xf: float = 0.0           # amount of the wave that freezes this step (cm)
    a: float = 0.0            # current melt depth (cm)
    at: float = 0.0           # net surface balance, freeze/melt units (cm)
    waterheld: float = 0.0    # water held in the pack (cm)
    kappa: float = 0.0        # numerical factor in the flow relation
    d: float = 0.0            # mean crystal diameter (cm)
    k: float = 0.0            # intrinsic permeability (cm^2)
    timestp: int = 0          # model time step (s)
    # wave(100,5) in Fortran; rows 1..100, cols 1..5 used. Row 0 / col 0 are
    # phantom (1-based indexing, see fasst.state docstring). Row 101 is spare
    # headroom for removewave's `wave(ic+1,:)` shift (doc point 8).
    wave: np.ndarray = field(default_factory=lambda: np.zeros((102, 6), dtype=np.float64))

    # -- helper-local save:: (never in snow_model's zero-out; also not restored)
    dcold: float = 0.0        # inivarivals: previous crystal diameter (cm), live state
    newdense: float = 0.0     # predsndp: last new-snow density (g/cm^3), inert today


# ---------------------------------------------------------------------------
# 1/12  snow_model
# ---------------------------------------------------------------------------
def snow_model(
    state: FasstState,
    snow: SnowState,
    sprint: int,
    water_flag: int,
    first_snow: int,
    sdensi: float,
    sdf: float,
):
    """Run one time step of the snow melt-water model.

    Fortran: module_snow.F90 subroutine ``snow_model`` (lines 5-326). Called
    from fasst_main.F90 (STEP 5, twice: ice-water and soil cases) and
    open_water.F90 (open-water surface).

    Args:
        state: model global state (``met``, ``dmet1``, ``stt``, ``sn_stat``, ...).
        snow: persistent pack state (mutated in place).
        sprint: 1 -> Fortran writes a unit-55 diagnostic line; the write is
            omitted here (see module docstring), the flag is kept for parity.
        water_flag: 0 -> land (soil heat flux from the profile); else -> ice.
        first_snow: 0 -> (re)initialize the pack this call; returned as 1.
        sdensi: initial snow density (kg/m^3), used only when ``first_snow == 0``.
        sdf: snow depth (m) in / out.

    Returns:
        ``(first_snow, sdf, phie1)`` -- the three Fortran ``intent(inout)`` /
        ``intent(out)`` arguments. ``phie1`` is the effective porosity that
        new_profile.F90 consumes across time steps (doc point 9).
    """
    iw = state.iw
    met = state.met
    nnodes = state.nnodes

    phie1 = 0.0
    state.vsmelt = 0.0

    if first_snow == 0:
        snow.n = 3.3
        snow.phie = 6.5e-1
        if abs(state.hsaccum - SPFLAG) > EPS and abs(state.iswe - SPFLAG) > EPS:
            snow.sd = state.hsaccum * 1e2
            snow.swe = state.iswe * 1e2
            if state.hsaccum > EPS and state.iswe > EPS:
                snow.rho = state.iswe / state.hsaccum
            else:
                snow.rho = sdensi * 1e-3
            if snow.swe <= EPS and snow.sd > EPS:
                snow.swe = snow.sd / 3.0
        elif abs(state.hsaccum - SPFLAG) > EPS and abs(state.iswe - SPFLAG) <= EPS:
            snow.sd = state.hsaccum * 1e2
            snow.rho = sdensi * 1e-3
            snow.swe = snow.sd / 3.0
        elif abs(state.hsaccum - SPFLAG) <= EPS and abs(state.iswe - SPFLAG) > EPS:
            snow.swe = state.iswe * 1e2
            snow.rho = sdensi * 1e-3
            snow.sd = state.iswe * 3.0
        # No else: if both are the missing-data flag, sd/swe/rho keep their
        # previous (persisted) values -- Fortran has no case(default) here.
        snow.swi = 2.88e-2
        snow.xf = 0.0
        snow.a = 0.0
        snow.at = 0.0
        snow.waterheld = 0.0
        snow.maxwet = 0.0
        snow.tottim = 0.0
        snow.maxwetold = 0.0
        snow.vd = 0.0
        snow.vw = 0.0
        snow.sdold = sdf * 1e2
        snow.kappa = 0.0
        for i in range(1, 51):  # Fortran zeros wave(1:50,:) only, though it is wave(100,5)
            for j in range(1, 6):
                snow.wave[i, j] = 0.0
        snow.rain2 = 0.0
        snow.smwaves = 0
        snow.timestp = int(state.timstep * 3.6e2)  # seconds; int() truncates toward zero

        snow.temphist = met[iw, IP_TSOIL] - TREF  # C
        snow.d = 0.0
        snow.k = 0.0
        snow.tims = 0
        snow.dend = 0.0
        snow.sph = 0.0

        if state.infer_test == 1 and sdf > EPS:
            snow.sdold = state.sn_stat[1]
            snow.sd = state.sn_stat[2]
            snow.rho = state.sn_stat[3]
            snow.maxwet = state.sn_stat[4]
            snow.maxwetold = state.sn_stat[5]
            snow.vd = state.sn_stat[6]
            snow.vw = state.sn_stat[7]
            snow.tottim = state.sn_stat[8]
            snow.temphist = state.sn_stat[9]
            snow.phie = state.sn_stat[10]
            snow.swi = state.sn_stat[11]
            snow.swe = state.sn_stat[12]
            snow.rain2 = state.sn_stat[13]
            snow.dend = state.sn_stat[14]
            snow.sph = state.sn_stat[15]
            snow.tims = int(state.sn_istat[1])
            snow.smwaves = int(state.sn_istat[2])

    if abs(snow.sd) <= EPS:
        snow.smwaves = 0
        snow.tottim = 0.0

    timestp = snow.timestp  # persists between calls when first_snow != 0

    nsnow = 0
    if _aint(met[iw, IP_PT2]) == 3 or _aint(met[iw, IP_PT]) == 3:
        nsnow = 1
    snow.tottim = snow.tottim + float(timestp)  # total model run time (s)

    inivarivals(state, snow, timestp, nsnow)

    # --- main physics ---
    if water_flag == 0:
        qb_f1 = _fdiv(
            (state.grthcond[nnodes] + state.grthcond[nnodes - 1]) * 5e-1,
            (state.nz[nnodes] - state.nz[nnodes - 1]) * 5e-1,
        )  # W/m^2 K -- unguarded in Fortran (nz is monotone by construction); parity
    else:
        qb_f1 = 0.0
        if met[iw, IP_HI] > EPS:
            qb_f1 = ITHCOND / met[iw, IP_HI]
        if met[iw, IP_HI] > EPS:
            qb_f1 = ITHCOND * state.km[1] / (state.km[1] * met[iw, IP_HI] + ITHCOND * snow.sd)

    abot, meta, add, dwind, atop = getwave(state, snow, water_flag, timestp, qb_f1)

    if abot > EPS:
        abot1 = abot * snow.rho  # cm swe
    else:
        abot1 = 0.0

    v = shoveout(state, snow, timestp, snow.kappa, snow.n, snow.maxwet, snow.d)

    # sprint == 1: Fortran writes met/depth/melt terms to unit 55 here.
    # Omitted (see module docstring); no state depends on it.

    if snow.sd < HT_MINM * 1e2:
        abot1 = abot1 + snow.sd * snow.rho
        snow.sd = 0.0
        snow.rho = 0.0

    # doc point 10: snow_model rounds at 1e10; every other subroutine at 1e15.
    snow.n = _anint(snow.n, 10)
    snow.phie = _anint(snow.phie, 10)
    snow.sd = _anint(snow.sd, 10)
    snow.swe = _anint(snow.swe, 10)
    snow.swi = _anint(snow.swi, 10)
    snow.xf = _anint(snow.xf, 10)
    snow.a = _anint(snow.a, 10)
    snow.at = _anint(snow.at, 10)
    snow.waterheld = _anint(snow.waterheld, 10)
    snow.maxwet = _anint(snow.maxwet, 10)
    snow.rho = _anint(snow.rho, 10)
    snow.maxwetold = _anint(snow.maxwetold, 10)
    snow.vd = _anint(snow.vd, 10)
    snow.vw = _anint(snow.vw, 10)
    snow.kappa = _anint(snow.kappa, 10)
    snow.rain2 = _anint(snow.rain2, 10)
    snow.d = _anint(snow.d, 10)
    snow.k = _anint(snow.k, 10)
    snow.tottim = _anint(snow.tottim, 10)
    snow.temphist = _anint(snow.temphist, 10)
    snow.dend = _anint(snow.dend, 10)
    snow.sph = _anint(snow.sph, 10)
    abot1 = _anint(abot1, 10)

    for i in range(1, 51):  # rows 1..50 only (see the wave zero-out above)
        for j in range(1, 6):
            snow.wave[i, j] = _anint(snow.wave[i, j], 10)

    sdf = snow.sd * 1e-2  # m
    sdf = _anint(sdf, 10)
    state.atopf = atop * 1e-2  # m
    state.atopf = _anint(state.atopf, 10)
    state.sdens[iw] = snow.rho * 1e3  # kg/m^3
    first_snow = 1
    state.dsnow = snow.d * 1e-2  # m
    phie1 = snow.phie  # unitless

    # melting needs energy -> cools the surface
    state.refreeze = 0.0  # m
    if snow.sdold > EPS:
        state.refreeze = 1e-2 * (atop + abot) * snow.rho
    elif snow.sdold <= EPS and snow.sd > EPS:
        state.refreeze = 1e-2 * (atop + abot) * snow.rho
    state.refreeze = _aint(state.refreeze, 15)  # AINT (truncate), not ANINT -- doc point 10

    state.vsmelt = (v + abot1) * 1e-2  # m
    state.vsmelt = _anint(state.vsmelt, 10)

    state.slushy[iw] = 0
    if snow.sd > EPS and snow.maxwet / snow.sd > 0.15:
        state.slushy[iw] = 1

    state.sn_stat[1] = snow.sdold
    state.sn_stat[2] = snow.sd
    state.sn_stat[3] = snow.rho
    state.sn_stat[4] = snow.maxwet
    state.sn_stat[5] = snow.maxwetold
    state.sn_stat[6] = snow.vd
    state.sn_stat[7] = snow.vw
    state.sn_stat[8] = snow.tottim
    state.sn_stat[9] = snow.temphist
    state.sn_stat[10] = snow.phie
    state.sn_stat[11] = snow.swi
    state.sn_stat[12] = snow.swe
    state.sn_stat[13] = snow.rain2
    state.sn_stat[14] = snow.dend
    state.sn_stat[15] = snow.sph
    state.sn_istat[1] = snow.tims
    state.sn_istat[2] = snow.smwaves

    return first_snow, sdf, phie1


# ---------------------------------------------------------------------------
# 2/12  inivarivals -- grain geometry, density, permeability
# ---------------------------------------------------------------------------
def inivarivals(state: FasstState, snow: SnowState, timestp: int, nsnow: int) -> None:
    """Initialize the snow grain / density / permeability parameters.

    Fortran: module_snow.F90 lines 329-589. Grain growth: Brun (1989)
    Annals of Glaciology 13; grain shape / size: Crocus (Vionnet et al.
    2012). ``dcold`` is a persistent (``save::``) local -> ``snow.dcold``.
    """
    iw = state.iw
    met = state.met
    nnodes = state.nnodes
    pi = PI  # Fortran: parameter pi = 3.141592654d0 (== fasst.constants.PI)

    l = 0.0
    ke = 0.0
    vav = 0.0
    f1 = 0.0
    f2 = 0.0
    dc = 0.0
    dend_new = 0.0
    sph_new = 0.0
    pwet = 0.0
    tgrad = 0.0
    add = 0.0

    dtime = float(timestp) / (6e1 * 2.4e1)  # time step in days

    if snow.sd > EPS:
        tgrad = abs(state.toptemp - state.stt[nnodes]) / (snow.sd * 1e-2)  # K/m
    tave = 5e-1 * (state.toptemp + state.stt[nnodes])  # K
    if snow.sd > EPS:
        pwet = 1e2 * snow.maxwet * snow.swi / snow.sd  # unitless

    if abs(snow.rho) <= EPS:
        snow.rho = 5e-1 * (SDENSD + SDENSW) * 1e-3  # middle of the allowed range (g/cm^3)
    vo = 5e-5

    # Shift crystals between the wet and dry portions based on the change in
    # snow depth and wet depth since the previous step (no real growth here).
    if snow.maxwet > EPS:
        f1 = 1.0 / snow.maxwet
    if snow.sd - snow.maxwet > EPS:
        f2 = 1.0 / (snow.sd - snow.maxwet)

    if snow.maxwet > snow.maxwetold:  # deeper wet snow
        snow.vw = snow.vw * (snow.maxwetold * f1) + snow.vd * ((snow.maxwet - snow.maxwetold) * f1)
        if snow.sd > snow.sdold:  # deeper pack this step
            if snow.sd > snow.maxwet:  # pack not entirely wet
                snow.vd = snow.vd * ((snow.sdold - snow.maxwet) * f2) + vo * (
                    (snow.sd - snow.sdold) * f2
                )
    elif snow.maxwet < snow.maxwetold:  # shallower wet snow
        if snow.sd > snow.sdold:
            if snow.sd > snow.maxwet:
                snow.vd = (
                    vo * ((snow.sd - snow.sdold) * f2)
                    + snow.vd * ((snow.sdold - snow.maxwetold) * f2)
                    + snow.vw * ((snow.maxwetold - snow.maxwet) * f2)
                )
        else:
            if snow.sd > snow.maxwet:
                snow.vd = snow.vd * ((snow.sd - snow.maxwetold) * f2) + snow.vw * (
                    (snow.maxwetold - snow.maxwet) * f2
                )
    else:
        if snow.sd > snow.sdold:
            if snow.sd > snow.maxwet:
                snow.vd = vo * ((snow.sd - snow.sdold) * f2) + snow.vd * (
                    (snow.sdold - snow.maxwet) * f2
                )

    if snow.maxwet >= snow.sd:
        snow.vd = vo  # maxwet should not exceed sd
    if snow.maxwet <= EPS:
        snow.vw = 0.0
    snow.sdold = snow.sd
    snow.maxwetold = snow.maxwet

    v1 = 1.28e-11
    snow.vd = snow.vd + v1 * float(timestp)  # cm^3
    if pwet > EPS:
        l = pwet
        v1 = 1.28e-11 + 4.22e-13 * (l**3.0)  # cm^3/s
        snow.vw = snow.vw + v1 * float(timestp)  # cm^3

    if snow.sd > EPS:
        f1 = 1.0 / snow.sd
        vav = snow.vw * (snow.maxwet * f1) + snow.vd * ((snow.sd - snow.maxwet) * f1)
        if vav < 0.0:
            vav = vo
    else:
        vav = vo

    if vav > EPS:
        snow.d = 2.0 * (((3.0 / (4.0 * pi)) * vav) ** (1.0 / 3.0))

    if 7.8 * snow.rho > 5e1 or snow.d <= 0.0:
        snow.k = 7.7e-2
    else:
        snow.k = 7.7e-2 * (snow.d**2.0) * math.exp(-7.8 * snow.rho)

    if abs(snow.n) > 0.0:
        f1 = 1.0 / snow.n
    if snow.k > EPS and abs(f1) > EPS:
        ke = _fdiv(snow.k**f1, snow.phie)
    if abs(f1) > EPS:
        snow.kappa = snow.n * ke * (5.47e4**f1)  # constant used in the model

    # grain shape -- Crocus (Vionnet et al. 2012, GMD 5, pp.773-791)
    if nsnow == 1:
        dend_new = min(max(1.29 - 1.7e-1 * met[iw, IP_WS], 2e-1), 1.0)
        sph_new = min(max(8e-2 * met[iw, IP_WS] + 3.8e-1, 5e-1), 9e-1)
        add = (state.dmet1[iw, 7] + state.dmet1[iw, 6]) * 1e2 * math.cos(state.sloper)

    f1 = math.exp(_fdiv(-6e3, tave))  # tave ~ 2*273 K in practice; _fdiv guards tave == 0
    dend_old = snow.dend
    sph_old = snow.sph
    if pwet > EPS:  # wet snow
        if dend_old > EPS:
            if snow.sph < 1.0:  # dendritic
                snow.dend = dend_old - dtime * pwet * pwet * pwet / 1.6e1
                snow.sph = sph_old
            else:
                snow.dend = dend_old
                snow.sph = sph_old + dtime * pwet * pwet * pwet / 1.6e1
        else:  # non-dendritic
            if snow.sph < 1.0:
                snow.sph = sph_old + dtime * pwet * pwet * pwet / 1.6e1
            else:
                snow.sph = sph_old
    else:  # dry snow
        if dend_old > EPS:
            if tgrad <= 5.0:
                snow.dend = dend_old - 2e8 * f1 * dtime
                snow.sph = sph_old + 1e9 * f1 * dtime
            else:
                snow.dend = dend_old - 2e8 * f1 * (tgrad**4e-1) * dtime
                snow.sph = sph_old + 2e8 * f1 * (tgrad**4e-1) * dtime
        else:
            snow.dend = 0.0
            if tgrad <= 5.0:
                snow.sph = sph_old + 1e9 * f1 * dtime
            else:
                if sph_old > EPS:
                    snow.sph = sph_old + 2e8 * f1 * (tgrad**4e-1) * dtime
                else:
                    snow.sph = 0.0
    snow.dend = min(max(snow.dend, 0.0), 1.0)
    snow.sph = min(max(snow.sph, 0.0), 1.0)

    if snow.sd + add > EPS:
        # finding #9: guard is sd+add, denominator is sd*rho+add; and the sph
        # numerator drops the *rho the dend numerator keeps. Literal; NaN can
        # propagate through the clamps when sd > 0, rho == 0, add == 0.
        snow.dend = _fdiv(
            snow.dend * snow.sd * snow.rho + dend_new * add, snow.sd * snow.rho + add
        )
        snow.sph = _fdiv(snow.sph * snow.sd + sph_new * add, snow.sd * snow.rho + add)
    else:
        snow.dend = 0.0
        snow.sph = 0.0

    snow.dend = min(max(snow.dend, 0.0), 1.0)
    snow.dend = _anint(snow.dend, 15)
    snow.sph = min(max(snow.sph, 0.0), 1.0)
    snow.sph = _anint(snow.sph, 15)

    # grain size change -- Crocus
    tempt = state.toptemp - TREF

    if tempt < -4e1:
        f = 0.0
    elif tempt >= -4e1 and tempt > -2.2e1:
        f = 1.1e-2 * (tempt + 4e1)
    elif tempt >= -2.2e1 and tempt > -6e0:
        f = 2e-1 + 5e-2 * (tempt + 2.2e1)
    else:
        f = 1.0 - 5e-2 * tempt

    if snow.rho <= 1.5e-1:
        h = 1.0
    elif snow.rho > 1.5e-1 and snow.rho <= 4e-1:
        h = 1.0 - 4e-3 * (snow.rho * 1e3 - 1.5e2)
    else:
        h = 0.0

    if tgrad < 1.5e1:
        g = 0.0
    elif tgrad >= 1.5e1 and tgrad < 2.5e1:
        g = 1e-2 * (tgrad - 1.5e1)
    elif tgrad >= 2.5e1 and tgrad < 4e1:
        g = 1e-1 + 3.7e-2 * (tgrad - 2.5e1)
    elif tgrad >= 4e1 and tgrad < 5e1:
        g = 6.5e-1 + 2e-2 * (tgrad - 4e1)
    elif tgrad >= 5e1 and tgrad < 7e1:
        g = 8.5e-1 + 7.5e-3 * (tgrad - 5e1)
    else:
        g = 1.0

    if snow.sd <= EPS:
        dc = 0.0
    else:
        if snow.dcold <= EPS:
            snow.dcold = snow.d
        if snow.sph < 1.0 and snow.dend <= EPS:
            dc = snow.dcold + f * h * g * 1.0417e-9 * float(timestp) * 1e2  # cm
            snow.d = dc
        else:
            dc = snow.d  # cm
    dc = _anint(dc, 15)
    snow.dcold = dc


# ---------------------------------------------------------------------------
# 3/12  getwave -- step physics: precip, melt (top/bottom), wind drift, new wave
# ---------------------------------------------------------------------------
def getwave(
    state: FasstState,
    snow: SnowState,
    water_flag: int,
    timestp: int,
    qb_f1: float,
):
    """Energy / mass budget of the step's melt wave and its interaction with
    existing waves.

    Fortran: module_snow.F90 lines 592-879. Wind transport after Jordan &
    Andreas (1999) JGR 104(C4) pp.7785-7806, parameters for snow on sea ice.

    Returns ``(abot, meta, add, dwind, atop)`` -- all cm:
        abot  bottom melt depth
        meta  metamorphic/settling depth loss (negative of the predsndp d1)
        add   new-snow depth gain
        dwind wind-drift depth loss
        atop  top melt (or sublimation) depth
    """
    iw = state.iw
    met = state.met
    dmet1 = state.dmet1
    nnodes = state.nnodes

    atop = 0.0
    abot = 0.0
    dwind = 0.0
    meta = 0.0
    add = 0.0
    atopw = 0.0
    precamtr = 0.0
    precamts = 0.0
    ustarc = 0.0
    b4 = 0.0
    G = 0.0
    sdo = snow.sd

    if _aint(met[iw, IP_PT2]) == 3 or _aint(met[iw, IP_PT]) == 3:
        precamts = (dmet1[iw, 7] + dmet1[iw, 6]) * 1e2 * math.cos(state.sloper)  # cm swe
    elif _aint(met[iw, IP_PT]) == 2:
        precamtr = dmet1[iw, 6] * 1e2 * math.cos(state.sloper)  # cm swe

    # add snow or rain
    if precamtr > EPS:
        _sumf, d1 = predsndp(state, snow, 2, timestp, precamtr, 0.0)
        rain = _anint(d1, 15)  # noqa: F841 -- computed, never read again (matches Fortran)
    elif precamts > EPS:
        _sumf, d1 = predsndp(state, snow, 1, timestp, precamts, 0.0)
        add = _anint(d1, 15)

    # top melt depth
    f3 = _fdiv(1.0, LHFUS * snow.rho * 1e3)  # rho can be 0 -> Fortran 1/0 = inf; parity
    tave = (state.toptemp + met[iw, IP_TSOIL]) * 5e-1
    cold = spheats(tave, Phase.ICE) * 1e3 * (1e-2 * snow.sd * snow.rho) * (TREF - tave)  # J/m^2
    cold = 1e2 * cold * f3  # noqa: F841 -- computed, never read again (matches Fortran)

    if snow.sd > 0.0 and state.melt[iw] > 0.0:
        f1 = 1.0
        if state.meltfl.strip() == "s":
            atop = _fdiv(1e2 * float(timestp) * state.melt[iw], LHSUB * snow.rho * 1e3)  # cm sublimation
            atop = min(snow.sd, atop, f1 * float(timestp) / 3.6e2)
            atopw = 0.0
        else:
            state.meltfl = "m"
            atop = 1e2 * float(timestp) * state.melt[iw] * f3  # cm top melt
            atop = min(snow.sd, atop, f1 * float(timestp) / 3.6e2)
            atopw = atop * snow.rho  # cm swe
        atop = _anint(atop, 15)
        snow.sd = max(0.0, snow.sd - atop)

    snow.a = atopw + precamtr
    snow.at = atopw

    # bottom melt depth
    if snow.sd > 0.0:
        if water_flag == 0 and met[iw, IP_TSOIL] > TREF:
            # active form: internal gradient between the last two soil nodes.
            # Fortran also carries a commented qb_f1*(met(iw,ip_tsoil)-Tref) form (doc point 11).
            qbot = max(0.0, qb_f1 * (state.stt[nnodes - 1] - state.stt[nnodes]))  # W/m^2
        else:
            qbot = 0.0

        f1 = 1.0
        abot = min(snow.sd, 1e2 * float(timestp) * qbot * f3, f1 * float(timestp) / 3.6e2)  # cm
        abot = _anint(abot, 15)  # cm
        snow.sd = max(0.0, snow.sd - abot)

    # wind transport
    if snow.sd > EPS:
        n1 = 12.91 - 24.52 * snow.phie + 11.88 * snow.phie * snow.phie
        ustart = 0.1 + 4.0 * _fsqrt((1.0 - snow.phie) * n1 * 0.2)  # m/s

        ht = state.iheightn - (state.hsaccum + state.hi + state.newsd)
        while ht <= 2.0:
            ht = ht + 5e-1
        z0g = 7.775e-3  # (0.05 - 1.5 mm)

        # Richardson number
        if abs(dmet1[iw, 4] - state.toptemp) <= EPS or dmet1[iw, 3] <= EPS:
            Rib = 0.0
        else:
            Rib = (
                2.0
                * GRAV
                * ht
                * (dmet1[iw, 4] - state.toptemp)
                / (dmet1[iw, 3] * dmet1[iw, 3] * (dmet1[iw, 4] + state.toptemp))
            )

        # z0h, z0q; Louis (1979) for ustar
        f2 = ht / z0g
        Cdng0 = VK_LOCAL * VK_LOCAL / (math.log(f2) * math.log(f2))  # doc point 5: local 0.35

        if Rib < 0.0:  # unstable
            c = 7.4 * Cdng0 * 9.4 * _fsqrt(f2)
            Gammam = 1.0 - 9.4 * Rib / (1.0 + c * _fsqrt(abs(Rib)))
        elif abs(Rib) < EPS:  # neutral
            Gammam = 1.0
        else:  # stable
            Gammam = 1.0 / ((1.0 + 4.7 * Rib) * (1.0 + 4.7 * Rib))

        ustars = _fsqrt(Cdng0 * Gammam) * dmet1[iw, 3]

        rhoa = 3.48e-3 * (met[iw, IP_AP] * 1e2 / dmet1[iw, 4])  # kg/m^3
        taut = ustart * ustart * rhoa  # kg/m*s^2
        tau = ustars * ustars * rhoa  # kg/m*s^2
        G = (5.0 / GRAV) * (ustars - ustart) * (tau - taut)  # kg/m*s

        if snow.rho * 1e3 < 4e2 or abs(snow.rho * 1e3 - 4e2) > 5e1:
            b4 = 1.0
        else:
            b4 = math.exp(-4.6e-2 * (snow.rho * 1e3 - 4e2))  # s/m

        # from Crocus
        f1 = 1.25 - 4.2e-3 * (max(SDENSD, snow.rho * 1e3) - SDENSD)
        if snow.dend > EPS:
            Mo = 3.4e-1 * (7.5e-1 * snow.dend - 5e-1 * snow.sph) + 6.6e-1 * f1
        else:
            Mo = 3.4e-1 * (-5.83e-1 * (snow.d * 1e1) - 8.33e-1 * snow.sph + 8.33e-1) + 6.6e-1 * f1
        Mo = max(Mo, 0.0)
        # doc point 6: active version uses +2.868 (commented: -2.868). With +,
        # driftp > eps for ANY wind/Mo, so `if driftp > eps` never blocks.
        driftp = 2.868 * math.exp(-8.5e-2 * met[iw, IP_WS]) + 1.0 + Mo
        if driftp > EPS:
            ustarc = -math.log((1.0 + Mo) / 2.868) / 8.5e-2

        t1 = 4e-2 * dmet1[iw, 4] + 8.84e-2 * (snow.rho * 1e3) * GRAV * 1e-2
        if (met[iw, IP_WS] >= ustarc and driftp > EPS) and t1 < 5e1:
            # active factor 2.66e-2 (commented in Fortran: 2.66e-3) -- doc point 11
            dwind = snow.sd * 2.66e-2 * b4 * G * math.exp(-t1)  # m/s
        else:
            dwind = 0.0
        dwind = 1e2 * dwind * float(timestp)  # cm
        dwind = _anint(dwind, 15)

        snow.sd = max(0.0, snow.sd - dwind)
    snow.sd = _anint(snow.sd, 15)

    # reported snow depth, else predict metamorphism / settling
    if _aint(abs(met[iw, IP_SD] - state.mflag) * 1e5) * 1e-5 > EPS:  # reported depth
        snow.sd = met[iw, IP_SD] * 1e2
        if abs(snow.sd) <= EPS and sdo > EPS:
            snow.sd = sdo
        snow.swe = snow.sd * snow.rho
    else:  # depth missing -> predict
        _sumf, d1 = predsndp(state, snow, 0, timestp, precamts, 0.0)
        meta = -_anint(d1, 15)

    # freezing / refreezing
    satvolume(
        state, snow, timestp, snow.rho, snow.tottim, snow.at, snow.kappa, snow.n, snow.sdold
    )

    # irreducible liquid water
    predictswi(state, snow, 1, timestp, 0.0)

    if snow.a > 1e-4:  # water in the pack
        snow.smwaves = snow.smwaves + 1
        assert snow.smwaves < 100, (
            "snow melt-wave array overflow: Fortran relies on -fcheck=all here (doc point 8)"
        )
        snow.wave[snow.smwaves, 1] = -snow.tottim
        snow.wave[snow.smwaves, 2] = snow.sd
        snow.wave[snow.smwaves, 3] = snow.a

        if snow.sd > 0.0:  # snow on the ground
            vol = snow.a
            wval = snow.smwaves
            bot = bottom(state, snow, wval, timestp, vol, snow.kappa, snow.n, snow.sd)
            snow.wave[snow.smwaves, 4] = bot  # time for the wave to reach the base
        else:
            snow.wave[snow.smwaves, 4] = 0.0

        snow.wave[snow.smwaves, 5] = -1.0  # -1 flags "infinite time to reach the base"
        if snow.smwaves > 1:
            checkwaves(state, snow, timestp, snow.kappa, snow.n, snow.tottim, snow.sd)
    elif snow.a > 0.0 and snow.a <= 1e-4 and snow.smwaves >= 2:
        snow.wave[snow.smwaves - 1, 3] = snow.wave[snow.smwaves - 1, 3] + snow.a

    return abot, meta, add, dwind, atop


# ---------------------------------------------------------------------------
# 4/12  checkwaves -- does the new wave overtake earlier waves?
# ---------------------------------------------------------------------------
def checkwaves(
    state: FasstState,
    snow: SnowState,
    timestp: int,
    kappa: float,
    n: float,
    tottim: float,
    sd: float,
) -> None:
    """Combine successive waves so the numbering stays top-to-bottom.

    Fortran: module_snow.F90 lines 882-931. Header comment claims "no
    subroutines called" -- WRONG: it calls :func:`collide` (line 926). This
    is a real defect in the Fortran header, flagged in the doc.
    """
    d4 = 1
    if snow.smwaves > 1:
        while (
            snow.wave[snow.smwaves - 1, 4] - snow.wave[snow.smwaves - 1, 1]
            > snow.wave[snow.smwaves, 4] - snow.wave[snow.smwaves, 1]
        ):
            w1 = snow.smwaves - 1
            w2 = snow.smwaves
            collide(state, snow, w1, w2, d4, timestp, kappa, n, tottim, sd)
            if snow.smwaves <= 1:
                break


# ---------------------------------------------------------------------------
# 5/12  shoveout -- water leaving the base of the pack this step
# ---------------------------------------------------------------------------
def shoveout(
    state: FasstState,
    snow: SnowState,
    timestp: int,
    kappa: float,
    n: float,
    maxwet: float,
    d: float,
) -> float:
    """Integrate the volume of water flowing out of the pack this step.

    Fortran: module_snow.F90 lines 934-1099. Returns ``v`` (cm), the outflow
    volume. Mutates ``snow`` in place (waves removed/merged, ``tottim``
    temporarily advanced then restored, ``rho``/``swe``/... via ``predsndp``).
    """
    wave = snow.wave
    wd = 0
    smwavesout = 0
    timend = 0
    flipper = 0
    tb = 0.0
    t1 = 0.0
    h1 = 0.0
    h2 = 0.0
    p1 = 0.0
    p2 = 0.0
    p3 = 0.0
    p4 = 0.0
    p5 = 0.0
    v = 0.0

    if (snow.smwaves > 0) and (snow.sd > 0.0):  # waves and snow present
        tb = snow.tottim  # start time

        smwaves_entry = snow.smwaves  # Fortran DO trip count is fixed at entry
        for ctr in range(1, smwaves_entry + 1):
            if wave[ctr, 4] - wave[ctr, 1] < snow.tottim + float(timestp):
                smwavesout = smwavesout + 1

        while snow.tottim < (tb + float(timestp)) and smwavesout > wd:
            wd = wd + 1  # wave being processed

            if snow.tottim >= (wave[1, 4] - wave[1, 1]):
                t1 = snow.tottim
            else:
                t1 = wave[1, 4] - wave[1, 1]

            if wd == smwavesout:  # last wave out
                if abs(wave[1, 5] + 1.0) <= EPS:  # lasts forever
                    snow.tottim = tb + float(timestp)
                else:
                    if (wave[1, 5] - wave[1, 1]) > (tb + float(timestp)):
                        snow.tottim = tb + float(timestp)
                    else:
                        snow.tottim = wave[1, 5] - wave[1, 1]
                        timend = -1
                flipper = 0
            else:  # not the last wave
                if abs(wave[1, 5] + 1.0) <= EPS:  # no freezing
                    snow.tottim = wave[2, 4] - wave[2, 1]  # time the 2nd wave exits
                else:  # freezing
                    if wave[1, 5] - wave[1, 1] > wave[2, 4] - wave[2, 1]:
                        snow.tottim = wave[2, 4] - wave[2, 1]
                    else:
                        snow.tottim = wave[1, 5] - wave[1, 1]
                        timend = -1
                flipper = -1  # waves will collide

            if (abs(kappa) > EPS and wave[1, 2] > 0.0) and abs(n - 1.0) > EPS:
                p1 = _fpow(wave[1, 2] / kappa, n / (n - 1.0))
            p2 = 1.0 - n
            p3 = snow.tottim + wave[1, 1]
            if p3 <= EPS:
                p3 = 1.0  # stops the program bombing
            p4 = t1 + wave[1, 1]
            if p4 <= EPS:
                p4 = 1.0  # ditto
            if abs(n - 1.0) > EPS:
                p5 = 1.0 / (1.0 - n)
            if p3 > EPS and abs(p5) > EPS:
                h1 = p1 * p2 * _fpow(p3, p5)
            if p4 > EPS and abs(p5) > EPS:
                h2 = p1 * p2 * _fpow(p4, p5)  # h2 = -wave(#,3)
            v = v + h1 - h2  # water flowing out of the pack

            d2 = 1
            if timend != 0:
                removewave(snow, d2)
            if flipper != 0:
                if timend == 0:
                    d3 = 2
                    d4 = 2
                    collide(state, snow, d2, d3, d4, timestp, kappa, n, snow.tottim, snow.sd)
            flipper = 0
            timend = 0
        snow.tottim = tb
    elif snow.smwaves != 0:  # waves but no snow
        smwaves_entry = snow.smwaves
        for ctr in range(1, smwaves_entry + 1):
            if wave[ctr, 4] - wave[ctr, 1] < snow.tottim:
                if (abs(kappa) > EPS and wave[ctr, 2] > 0.0) and abs(n - 1.0) > EPS:
                    p1 = _fpow(wave[ctr, 2] / kappa, n / (n - 1.0))
                p2 = 1.0 - n
                p3 = snow.tottim + wave[ctr, 1]
                if p3 <= EPS:
                    p3 = 1.0  # stops the program bombing
                if abs(n - 1.0) > EPS:
                    p5 = 1.0 / (1.0 - n)
                if p3 > EPS and abs(p5) > EPS:
                    h1 = p1 * p2 * _fpow(p3, p5)
                v = v - h1
            else:
                v = v + wave[ctr, 3]
            removewave(snow, ctr)

    t1 = 0.0
    if v > 0.0:
        predsndp(state, snow, 3, timestp, t1, v)

    return v


# ---------------------------------------------------------------------------
# 6/12  predsndp -- snow depth / SWE: settling, new snow, rain runoff
# ---------------------------------------------------------------------------
def predsndp(
    state: FasstState,
    snow: SnowState,
    smode: int,
    timestp: int,
    newsnow: float,
    meltswe: float,
):
    """Predict snow depth and snow water equivalent.

    Fortran: module_snow.F90 lines 1102-1381. ``newsnow`` is a SWE-equivalent
    depth. Uses ``dense`` (Jordan et al. 1999). ``newdense`` is a persistent
    (``save::``) local -> ``snow.newdense``.

    smode: 0 = metamorphic densification + overburden settling;
           1 = add new snow;  2 = add rain water;  3 = rain runoff / SWE loss.

    Returns ``(sumf, d1)`` -- ``sumf`` (new-snow density, only set for
    smode 1) is read by no caller; ``d1`` is the depth change (cm).
    """
    iw = state.iw
    dmet1 = state.dmet1
    nnodes = state.nnodes

    sumf = 0.0
    d1 = 0.0
    meta = 0.0
    over = 0.0
    overc = 0.0
    eda = 0.0
    avgsnpress = 0.0
    f1 = 0.0
    rave = 0.0  # doc point 3: only assigned in case(0); stays 0.0 in case(1)

    # local recompute of the slope angle (doc point 7: duplicates state.sloper,
    # minus the anint(...*1d20) rounding the two other call sites apply).
    sloper = max(0.0, min(1.57, state.slope_fasst * PI / 1.8e2))

    if smode == 0:
        # --- metamorphic densification + overburden settling ---
        if snow.sd > EPS:  # snow present
            f1 = 1.0 / snow.sd
            rave = snow.rho
            critdense = 1e-1
            c1 = 2.778e-6  # 1/s
            c3 = min(1.0, max(0.0, math.exp(-4.6e1 * (rave - critdense))))
            c4 = max(1.0, min(2.0, 1.0 + snow.maxwet * f1))
            ttemp = (state.toptemp + state.stt[nnodes]) * 5e-1  # K

            if abs(4e-2 * (TREF - ttemp)) < 5e1:
                meta = min(0.0, -c1 * c3 * c4 * math.exp(-4e-2 * (TREF - ttemp)))  # 1/s
            else:
                meta = 0.0
            meta = _anint(meta, 15)

            # sntherm overburden (result `over` -- the one actually used)
            c5 = 1e-1  # 1/C
            c6 = 2.3e1  # cm^3/g
            eda0 = 9e5  # kg*s/m^2
            if abs(c5 * (TREF - ttemp) + c6 * rave) < 5e1:
                avgsnpress = 1e1 * snow.swe * (2.0 / 3.0) * math.cos(sloper)  # kg/m^2
                eda = eda0 * math.exp(c5 * (TREF - ttemp) + c6 * rave)  # kg/m*s
                over = min(0.0, -avgsnpress / eda)  # 1/s
            else:
                over = 0.0
            over = _anint(over, 15)

            # crocus overburden (result `overc` computed and DISCARDED -- doc point 11)
            c5 = 1e-1
            c6 = 2.3e1
            eda0 = 7.62237e6  # kg/m*s
            if abs(c5 * (TREF - ttemp) + c6 * rave) < 5e1:
                avgsnpress = 9.81e1 * snow.swe * (2.0 / 3.0) * math.cos(sloper)
                f1c = 1.0 / (1.0 + 6e1 * snow.maxwet * snow.swi * f1)
                f2c = min(4.0, math.exp(1e1 * min(4e-1, snow.d * 1e1 - 2e-1)))
                edac = eda0 * f1c * f2c * rave / 2.5e-1
                eda = edac * math.exp(c5 * (TREF - ttemp) + c6 * rave)
                overc = min(0.0, -avgsnpress / eda)
            else:
                overc = 0.0
            overc = _anint(overc, 15)  # noqa: F841 -- discarded (doc point 11)

            delta = snow.sd * (over + meta) * float(timestp)  # cm
            d1 = delta  # cm

            snow.sd = snow.sd + delta  # cm

            splus = snow.swi * snow.swe  # water to satisfy the swi requirement
            snow.rho = _fdiv(snow.swe - splus, snow.sd)  # sd may be <= 0 after delta; parity
            if snow.rho > SDENSW * 1e-3:
                snow.rho = SDENSW * 1e-3
                snow.swi = 0.05
                snow.swe = snow.rho * snow.sd / (1.0 - snow.swi)
            snow.phie = (1.0 - snow.rho) * (1.0 - snow.swi)

    elif smode == 1:
        # --- add new snow to the pack ---
        snow.newdense = 0.0
        # Jordan et al. (1999) JGR 104(C4) p.7785-7806
        snow.newdense = dense(state.ptemp, dmet1[iw, 3], Phase.SNOW) * 1e-3  # g/cm^3
        # Fortran keeps four alternative new-snow-density schemes commented
        # here (Crocus; Anderson 1976; Hedstrom & Pomeroy 1998; Eta/NWP) -- doc point 11.
        f1 = _fdiv(1.0, snow.newdense)  # dense() clamps to [0.05, 0.55] -> safe; parity
        sumf = snow.newdense
        critdense = 1e-1
        c1 = 2.778e-6  # 1/s
        c3 = min(1.0, max(0.0, math.exp(-4.6e1 * (sumf - critdense))))
        c4 = 1.0

        if abs(4e-2 * (TREF - state.ptemp)) < 5e1:
            meta = min(0.0, -c1 * c3 * c4 * math.exp(-4e-2 * (TREF - state.ptemp)))  # 1/s
        else:
            meta = 0.0
        meta = _anint(meta, 15)

        # sntherm
        c5 = 1e-1  # 1/C
        c6 = 2.3e1  # cm^3/g
        eda0 = 9e5  # kg*s/m^2
        if abs(c5 * (TREF - state.ptemp) + c6 * sumf) < 5e1:
            avgsnpress = 1e1 * newsnow * (2.0 / 3.0) * math.cos(sloper)  # kg/m^2
            # doc point 3: guard above uses c6*sumf, this line uses c6*rave --
            # and rave is 0.0 in this branch, so the density term vanishes.
            eda = eda0 * math.exp(c5 * (TREF - state.ptemp) + c6 * rave)  # kg/m*s
            over = min(0.0, -avgsnpress / eda)  # 1/s
        else:
            over = 0.0
        over = _anint(over, 15)

        delta = newsnow * (over + meta) * float(timestp) * f1  # cm
        snow.sd = snow.sd + newsnow * f1 + delta
        snow.rho = _fdiv(snow.rho * snow.sd + newsnow, snow.sd + newsnow * f1)  # parity
        snow.phie = (1.0 - snow.rho) * (1.0 - snow.swi)

        d1 = newsnow * f1
        snow.swe = snow.swe + newsnow

    elif smode == 2:
        # --- add liquid water from rainfall ---
        snow.rain2 = newsnow  # Fortran has `rain2 + newsnow` commented out

    elif smode == 3:
        # --- rain runoff / SWE loss ---
        if snow.rain2 > meltswe:  # runoff
            snow.rain2 = snow.rain2 - meltswe
        else:  # rain held in the pack -> SWE grows
            if snow.swe > EPS:
                snow.swe = snow.swe + snow.rain2 - meltswe
                snow.rain2 = 0.0
        if snow.swe < 0.0:
            d1 = -snow.sd
            snow.swe = 0.0
            snow.sd = 0.0
            snow.rho = 0.0

    else:
        state.error_code = 1
        state.error_type = 4

    if snow.sd < 0.0:
        snow.sd = 0.0
        snow.swe = 0.0
        snow.rho = 0.0

    return sumf, d1


# ---------------------------------------------------------------------------
# 7/12  satvolume -- reach irreducible saturation; refreezing
# ---------------------------------------------------------------------------
def satvolume(
    state: FasstState,
    snow: SnowState,
    timestp: int,
    rho: float,
    tottim: float,
    at: float,
    kappa: float,
    n: float,
    sdold: float,
) -> None:
    """Decide whether the pack has reached irreducible water saturation, how
    much of the current wave is used getting there, and any refreezing.

    Fortran: module_snow.F90 lines 1384-1598. Mutates ``snow`` in place;
    sets ``snow.a`` to the remaining melt volume on exit.
    """
    wave = snow.wave
    p1 = 0.0
    p2 = 0.0
    p3 = 0.0
    p5 = 0.0
    volgone = 0.0
    heatneed = 0.0
    heatgot = 0.0
    upneed = 0.0
    xfold = 0.0
    exdelta = 0.0

    if rho > EPS:
        kfs = 0.023 + (7.75e-5 * rho * 1e3 + 1.105e-6 * ((rho * 1e3) ** 2.0)) * (2.29 - 0.023)  # W/m*K
        kfs = kfs * 1e-2
    else:
        kfs = 0.0045
    l = 333.05
    xfold = snow.xf
    vol = snow.a
    t1 = state.dmet1[state.iw, 4] - TREF

    if vol <= 0.0:
        if t1 < 0.0:  # no melt, freezing conditions
            snow.temphist = (snow.temphist * snow.tims - t1) / (snow.tims + 1.0)
            snow.tims = snow.tims + 1

            # sqrt(l*swi) is the latent heat to freeze snow at swi level
            snow.xf = _fsqrt((2.0 * kfs * snow.temphist * timestp * snow.tims) / (l * snow.swi))
            if snow.xf > xfold:
                exdelta = snow.xf - xfold
                predictswi(state, snow, 3, timestp, exdelta)
            xfold = snow.xf

            if snow.xf >= snow.sd:  # freezing to the base -> drop all waves that would exit
                smwaves_entry = snow.smwaves
                for _ctr in range(1, smwaves_entry + 1):
                    removewave(snow, _ctr)
            else:  # freezing partway -> drop waves accordingly
                smwaves_entry = snow.smwaves
                for ctr in range(smwaves_entry, 0, -1):
                    parm1 = max(0.0, tottim + wave[ctr, 1])
                    ctr1 = ctr
                    depthval = depth(state, snow, ctr1, parm1, kappa, n, snow.sd)  # depth in the pack

                    if (snow.xf - depthval - snow.sd + wave[ctr, 2]) > 0.0:
                        removewave(snow, ctr)
                    else:
                        # time to freeze the wave; if that beats time-to-base, drop it
                        if (
                            wave[ctr, 2] + snow.xf - snow.sd > 0.0
                            and abs(kappa) > EPS
                            and abs(n - 1.0) > EPS
                        ):
                            p1 = _fpow((wave[ctr, 2] + snow.xf - snow.sd) / kappa, n / (n - 1.0))
                            p2 = n - 1.0
                            p3 = tottim + wave[ctr, 1]
                            if p3 <= EPS:
                                p3 = 1.0
                            p5 = 1.0 / (1.0 - n)
                            if abs(p5) > EPS:
                                volgone = p1 * p2 * _fpow(p3, p5)
                            if volgone > EPS:
                                wave[ctr, 5] = _fpow((n - 1.0) / volgone, n - 1.0) * _fpow(
                                    wave[ctr, 2] / kappa, n
                                )
                            if wave[ctr, 5] < wave[ctr, 4]:
                                removewave(snow, ctr)
                        break

        vol = 0.0
    elif vol > 0.0:
        if snow.maxwet - snow.sd > EPS and snow.sd > 0.0:  # wet pack, no frozen snow
            snow.maxwet = snow.sd
            if (snow.tims != 0) or (snow.temphist > EPS):
                snow.temphist = 0.0
                snow.tims = 0
                snow.xf = 0.0
                xfold = 0.0
        else:
            # heat to raise the frozen pack to melt level
            if (snow.tims != 0) or (snow.temphist > EPS):
                heatneed = 1.0405 * rho * snow.xf * snow.temphist
                heatgot = vol * 3330464.0
                if heatgot > heatneed:  # brings snow to melt temp
                    vol = (1.0 - (heatneed / heatgot)) * vol
                    snow.xf = 0.0
                    snow.tims = 0
                    snow.temphist = 0.0
                else:
                    vol = 0.0
                    snow.temphist = 0.9611 * _fdiv(heatneed - heatgot, rho * snow.xf)  # parity

            if snow.maxwet < snow.sd * snow.swi:  # pack not fully at swi level
                upneed = (snow.sd - snow.maxwet) * snow.swi
            else:
                upneed = 0.0

            if vol >= upneed:  # how much of the pack can reach swi level
                predictswi(state, snow, 2, timestp, upneed)
                vol = vol - upneed  # remainder leaves the pack
                if (snow.tims != 0) or (snow.temphist > EPS):
                    snow.temphist = 0.0
                    snow.tims = 0
                    snow.xf = 0.0
                    xfold = 0.0
            else:
                predictswi(state, snow, 2, timestp, vol)
                vol = 0.0
                if (snow.tims != 0) or (snow.temphist > EPS):
                    snow.temphist = (snow.temphist * snow.tims - t1) / (snow.tims + 1.0)
                    snow.tims = snow.tims + 1
                    if snow.temphist < EPS:
                        snow.temphist = 0.0
                        snow.tims = 0
                        snow.xf = 0.0
                        xfold = 0.0

    snow.a = vol


# ---------------------------------------------------------------------------
# 8/12  predictswi -- irreducible water saturation and wet-snow depth
# ---------------------------------------------------------------------------
def predictswi(
    state: FasstState,
    snow: SnowState,
    wmode: int,
    timestp: int,
    exdelta: float,
) -> None:
    """Predict irreducible water saturation ``swi`` and wet-snow depth
    ``maxwet``.

    Fortran: module_snow.F90 lines 1601-1734. ``exdelta`` (cm) is the water
    that freezes this step (wmode 2/3).

    wmode: 1 = predict swi, adjust maxwet;  2 = add water;  3 = refreeze.

    doc point 4: of the seven ``*const`` coefficients only ``tempconst``,
    ``rainconst``, ``runoffconst``, ``refreezeconst`` are ever assigned;
    ``heatbalconst``, ``wetconst`` and ``swiconst`` stay 0, so the ``at``,
    ``wetuse`` and ``swi`` feedback terms of ``swidelta`` are structurally
    zero. Preserved verbatim.
    """
    iw = state.iw

    rainuse = 0.0
    freezeuse = 0.0
    wetuse = 0.0
    swiold = 0.0
    f1 = 0.0
    # doc point 4: NEVER assigned anywhere in this routine.
    heatbalconst = 0.0
    wetconst = 0.0
    swiconst = 0.0

    if wmode == 1:  # predict swi and adjust maxwet
        if snow.sd > EPS:
            if _aint(state.met[iw, IP_PT]) == 2:
                rainuse = state.dmet1[iw, 6] * 1e2

            if snow.xf > EPS and snow.sd > EPS:
                freezeuse = min(1.0, max(0.0, snow.xf / snow.sd))  # frac. of snow with frozen water

            if snow.maxwet > EPS and snow.sd > EPS:
                wetuse = min(1.0, max(0.0, snow.maxwet / snow.sd))  # frac. at irreducible sat.

            f1 = 1.0 / 3.6e2
            tempconst = -5e-6 * f1
            rainconst = -4e-2 * f1
            runoffconst = -1.5e-3 * f1
            refreezeconst = -1.5e-4 * f1

            swidelta = (
                tempconst * (state.ptemp - TREF)
                + rainconst * rainuse
                + refreezeconst * freezeuse
                + runoffconst * snow.a
                + heatbalconst * snow.at  # == 0 (doc point 4)
                + wetconst * wetuse  # == 0 (doc point 4)
                + swiconst * snow.swi  # == 0 (doc point 4)
            )
            swidelta = swidelta * timestp
            swiold = snow.swi
            snow.swi = max(5e-2, min(1e-1, snow.swi + swidelta))

            if (snow.sd > snow.sdold) and (snow.tottim > (1.5 * float(timestp))):
                snow.swi = max(
                    5e-2, min(1e-1, (snow.swi * snow.sdold + 5e-2 * (snow.sd - snow.sdold)) / snow.sd)
                )

            snow.maxwet = max(0.0, min(snow.waterheld / snow.swi, snow.sd))

            # current melt depth
            if (rainuse > 0.0 or snow.maxwet > snow.sd - snow.swi) or (
                snow.a < snow.at or snow.swi < swiold
            ):
                if (snow.maxwet - snow.sd) * snow.swi > snow.at - snow.a:
                    snow.a = snow.at
                else:
                    snow.a = max(0.0, snow.a + snow.maxwet - (snow.sd - snow.swi))

            if snow.maxwet - (snow.sd - snow.swi) > EPS:  # wet-snow depth & water held
                snow.maxwet = max(0.0, snow.sd - snow.swi)
                snow.waterheld = snow.maxwet * snow.swi
        else:
            snow.swi = 5e-2
            snow.waterheld = 0.0
            snow.maxwet = 0.0

    elif wmode == 2:  # add water to the pack
        snow.waterheld = snow.waterheld + exdelta
        snow.maxwet = max(0.0, min(snow.waterheld / snow.swi, snow.sd))
        if snow.maxwet - (snow.sd - snow.swi) > EPS:
            snow.maxwet = max(0.0, snow.sd - snow.swi)

    elif wmode == 3:  # refreezing -> remove water
        omw = snow.maxwet
        snow.maxwet = max(0.0, snow.maxwet - exdelta)  # remove frozen water
        if (snow.maxwet + snow.xf) - snow.sd > EPS:  # physical consistency
            snow.maxwet = snow.sd - snow.xf
        snow.waterheld = max(0.0, snow.waterheld - (omw - snow.maxwet) * snow.swi)

    else:
        state.error_code = 1
        state.error_type = 4


# ---------------------------------------------------------------------------
# 9/12  bottom -- time for a wave to reach the base (Newton iteration)
# ---------------------------------------------------------------------------
def bottom(
    state: FasstState,
    snow: SnowState,
    wavenumber: int,
    timestp: int,
    volume: float,
    kappa: float,
    n: float,
    sd: float,
) -> float:
    """Time (s) for a wave to reach the base of the pack, by Newton's method
    when more than one wave is present.

    Fortran: module_snow.F90 lines 1737-1984. Returns ``bot``.

    ``f1 = 1/(1-n)`` here -- OPPOSITE sign to :func:`depth` / :func:`collide`
    which use ``1/(n-1)``. Not a typo; not factored out.

    Non-convergence sets ``state.error_code = 1`` / ``error_type = 4``. The
    first Newton branch checks ``iter >= 100`` (fires); the second checks
    ``iter > 100`` (Fortran line 1953) -- the loop exits at ``iter == 100``
    exactly, so that check is UNREACHABLE. Kept verbatim (finding, not in doc).
    """
    wave = snow.wave
    asign = 0
    ticker = 0
    it = 0
    wnt = wavenumber
    tau = 0.0
    d = 0.0
    d1 = 0.0
    d2 = 0.0
    d3 = 0.0
    d4 = 0.0
    p1 = 0.0
    p2 = 0.0
    p5 = 0.0
    diff = 0.0
    bwn = 0.0
    bwn1 = 0.0
    dta = 0.0
    bot = 0.0
    f1 = 0.0
    f2 = 0.0
    if abs(1.0 - n) > EPS:
        f1 = 1.0 / (1.0 - n)  # NOTE: 1/(1-n), negative -- opposite to depth()/collide()
    if n > 0.0:
        f2 = 1.0 / n

    if wavenumber == 1:
        # single wave: closed form
        if (abs(n - 1.0) > EPS and volume > EPS) and (wave[wnt, 2] > EPS and kappa > EPS):
            bot = _fpow((n - 1.0) / volume, n - 1.0) * _fpow(wave[wnt, 2] / kappa, n)
    else:
        # more than one wave; the previous wave may have lost water to freezing
        asign = -1

        if wave[wnt, 2] > 0.0:  # snow present
            if abs(wave[wnt - 1, 5] + 1.0) <= EPS:  # no frozen water
                tau = timestp * 5e-1  # half a time step
                dta = wave[wnt - 1, 1] - wave[wnt, 1]

                it = 0
                d = depth(state, snow, wnt, tau, kappa, n, sd)
                while abs(d - wave[wnt, 2]) >= 5e-3 * wave[wnt, 2] and it < 100:
                    d1 = d
                    d2 = 0.0
                    d3 = 0.0
                    if (tau > EPS and abs(f1) > EPS) and abs(tau + dta) > EPS:
                        d2 = (f2 - 1.0) * _fpow(_fpow(tau, f1) - _fpow(tau + dta, f1), asign)
                        d3 = f1 * (_fpow(tau, n * f1) - _fpow(tau + dta, n * f1))
                    d4 = d1 * d2 * d3
                    if abs(d4) > EPS:
                        tau = tau - (d - wave[wnt, 2]) / d4
                    if tau < 0.0 and abs(ticker - 1) > 0:  # tau cannot be negative -- retry
                        tau = timestp * _fpow(4.0, ticker - 1.0)
                        ticker = ticker - 1
                    d = depth(state, snow, wnt, tau, kappa, n, sd)
                    it = it + 1
                bot = tau
                _BOTTOM_CENSUS["iters1"] += it   # instrumentation
                _BOTTOM_CENSUS["calls1"] += 1

                if it >= 100:
                    state.error_code = 1
                    state.error_type = 4

                if (volume > EPS and abs(n - 1.0) > EPS) and (
                    wave[wnt, 2] > EPS and abs(kappa) > EPS
                ):
                    if bot > _fpow((n - 1.0) / volume, n - 1.0) * _fpow(wave[wnt, 2] / kappa, n):
                        state.error_code = 1
                        state.error_type = 4

            elif abs(wave[wnt - 1, 5] + 1.0) > EPS:  # frozen water present
                if (volume > EPS and abs(n - 1.0) > EPS) and (
                    wave[wnt, 2] > EPS and abs(kappa) > EPS
                ):
                    bot = _fpow((n - 1.0) / volume, n - 1.0) * _fpow(wave[wnt, 2] / kappa, n)
                if bot - wave[wnt, 1] < wave[wnt - 1, 5] - wave[wnt - 1, 1]:
                    tau = timestp * 5e-1
                    dta = wave[wnt - 1, 1] - wave[wnt, 1]
                    if (wave[wnt - 1, 2] > EPS and abs(kappa) > EPS) and abs(f1) > EPS:
                        p1 = _fpow(wave[wnt - 1, 2] / kappa, n * f1)
                    p2 = n - 1.0
                    p5 = f1

                    # water missing in the previous wave due to freezing
                    if abs(wave[wnt - 1, 5]) > EPS and abs(p5) > EPS:
                        diff = p1 * p2 * _fpow(wave[wnt - 1, 5], p5)

                    it = 0
                    if (
                        (abs(1.0 - f2) > EPS and f2 > EPS)
                        and (abs(tau) > EPS and abs(tau + dta) > EPS)
                    ) and (volume > EPS and abs(diff) > EPS):
                        bwn = kappa * _fpow(volume * (-f1), 1.0 - f2) * _fpow(tau, f2)
                        bwn1 = kappa * _fpow(diff * (-f1), 1.0 - f2) * _fpow(tau + dta, f2)
                        if bwn > bwn1:  # overtakes the previous wave
                            if wave[wnt, 3] > diff:
                                d1 = kappa * _fpow((wave[wnt, 3] - diff) * (-f1), 1.0 - f2) * _fpow(
                                    tau, f2
                                )
                                d2 = _fpow(1.0 - _fpow(tau / (tau + dta), -f1), f2 - 1.0)
                                d = d1 * d2
                            else:
                                d = bwn
                        else:
                            d = bwn

                    while abs(d - wave[wnt, 2]) >= 5e-3 * wave[wnt, 2] and it < 100:
                        d1 = d
                        d2 = 0.0
                        d3 = 0.0
                        if (abs(tau) > EPS and abs(tau + dta) > EPS) and abs(f1) > EPS:
                            d2 = (f2 - 1.0) * _fpow(_fpow(tau, f1) - _fpow(tau + dta, f1), asign)
                            d3 = f1 * (_fpow(tau, n * f1) - _fpow(tau + dta, n * f1))
                        d4 = d1 * d2 * d3
                        if abs(d4) > EPS:
                            tau = tau - (d - wave[wnt, 2]) / d4
                        if tau < 0.0 and abs(ticker - 1) > 0:
                            tau = timestp * _fpow(2.0, float(ticker - 1))
                            ticker = ticker - 1

                        if (
                            (abs(1.0 - f2) > EPS and f2 > EPS)
                            and (abs(tau) > EPS and abs(tau + dta) > EPS)
                        ) and (volume > EPS and abs(diff) > EPS):
                            bwn = kappa * _fpow(volume * (-f1), 1.0 - f2) * _fpow(tau, f2)
                            bwn1 = kappa * _fpow(diff * (-f1), 1.0 - f2) * _fpow(tau + dta, f2)
                            if bwn > bwn1:  # overtakes the previous wave
                                if wave[wnt, 3] > diff:
                                    d1 = kappa * _fpow(
                                        (wave[wnt, 3] - diff) * (-f1), 1.0 - f2
                                    ) * _fpow(tau, f2)
                                    d2 = _fpow(1.0 - _fpow(tau / (tau + dta), -f1), f2 - 1.0)
                                    d = d1 * d2
                                else:
                                    d = bwn
                            else:
                                d = bwn

                        it = it + 1
                    bot = tau
                    _BOTTOM_CENSUS["iters2"] += it   # instrumentation
                    _BOTTOM_CENSUS["calls2"] += 1

                    if it > 100:  # UNREACHABLE: loop exits at it == 100 (finding, not in doc)
                        state.error_code = 1
                        state.error_type = 4

                    if (volume > EPS and abs(n - 1.0) > EPS) and (
                        wave[wnt, 2] > EPS and abs(kappa) > EPS
                    ):
                        if bot > _fpow((n - 1.0) / volume, n - 1.0) * _fpow(
                            wave[wnt, 2] / kappa, n
                        ):
                            state.error_code = 1
                            state.error_type = 4
        else:
            bot = 0.0

    return bot


# ---------------------------------------------------------------------------
# 10/12  collide -- merge two waves that have met
# ---------------------------------------------------------------------------
def collide(
    state: FasstState,
    snow: SnowState,
    w1: int,
    w2: int,
    mode: int,
    timestp: int,
    kappa: float,
    n: float,
    tottim: float,
    sd: float,
) -> None:
    """Combine wave ``w2`` into the earlier wave ``w1``.

    Fortran: module_snow.F90 lines 1987-2100. mode 1 = they meet inside the
    pack; mode 2 = they meet at the base.

    ``f1 = 1/(n-1)`` here (positive) -- opposite sign to :func:`bottom`.
    """
    wave = snow.wave
    p1 = 0.0
    p2 = 0.0
    p5 = 0.0
    diff = 0.0
    vol = 0.0
    bot = 0.0
    x1 = 0.0
    vc1 = 0.0
    f1 = 0.0
    if abs(n - 1.0) > EPS:
        f1 = 1.0 / (n - 1.0)  # 1/(n-1), positive -- opposite to bottom()

    delta = wave[w1, 1] - wave[w2, 1]  # noqa: F841 -- Fortran computes this and never reads it; kept for parity

    if mode == 1:  # collide inside the pack
        wave[w1, 1] = wave[w2, 1]
        wave[w1, 2] = wave[w2, 2]

        if abs(wave[w1, 5] + 1.0) <= EPS:  # no freezing
            wave[w1, 3] = wave[w1, 3] + wave[w2, 3]  # combine water
        else:
            if (wave[w1, 2] > EPS and abs(kappa) > EPS) and abs(f1) > EPS:
                p1 = _fpow(wave[w1, 2] / kappa, n * f1)
            p2 = n - 1.0
            p5 = -f1
            if wave[w1, 5] > EPS and abs(p5) > EPS:
                diff = p1 * p2 * _fpow(wave[w1, 5], p5)
            else:
                diff = 0.0
            wave[w1, 3] = max(0.0, wave[w1, 3] + wave[w2, 3] - diff)  # minus frozen water
        wave[w1, 5] = wave[w2, 5]
        removewave(snow, w2)

        if wave[w1, 2] > 0.0 and wave[w1, 3] > EPS:
            vol = wave[w1, 3]
            wval = w1
            bot = bottom(state, snow, wval, timestp, vol, kappa, n, sd)
            wave[w1, 4] = bot  # time for the wave to reach the base
        else:
            wave[w1, 4] = 0.0

    elif mode == 2:  # collide at the base
        if (wave[w1, 2] > EPS and abs(kappa) > EPS) and abs(f1) > EPS:
            x1 = _fpow(wave[w1, 2] / kappa, n * f1)
            vc1 = x1 * (n - 1.0) * _fpow(tottim + wave[w1, 1], -f1)
            if wave[w1, 5] != -1.0:  # frozen water present
                p1 = _fpow(wave[w1, 2] / kappa, n * f1)
                p2 = n - 1.0
                p5 = -f1
                if abs(p5) > EPS:
                    diff = p1 * p2 * _fpow(wave[w1, 5], p5)
                vc1 = max(0.0, vc1 - diff)  # remove frozen water

        wave[w1, 1] = wave[w2, 1]
        wave[w1, 2] = wave[w2, 2]
        wave[w1, 3] = vc1 + wave[w2, 3]
        wave[w1, 4] = tottim + wave[w1, 1]
        wave[w1, 5] = wave[w2, 5]
        removewave(snow, w2)


# ---------------------------------------------------------------------------
# 11/12  removewave -- drop a wave from wave(100,5) and shift the rest up
# ---------------------------------------------------------------------------
def removewave(snow: SnowState, wn: int) -> None:
    """Remove wave ``wn`` (exited or merged) and shift the following waves.

    Fortran: module_snow.F90 lines 2103-2142.
    """
    wave = snow.wave
    for ic in range(wn, snow.smwaves + 1):
        for icc in range(1, 6):
            wave[ic, icc] = wave[ic + 1, icc]
    snow.smwaves = snow.smwaves - 1


# ---------------------------------------------------------------------------
# 12/12  depth -- depth of a wave at a given instant
# ---------------------------------------------------------------------------
def depth(
    state: FasstState,
    snow: SnowState,
    wn: int,
    t: float,
    kappa: float,
    n: float,
    sd: float,
) -> float:
    """Depth (cm) of wave ``wn`` at time ``t`` (s).

    Fortran: module_snow.F90 lines 2145-2266. Single wave, multiple waves,
    and partial freezing of the previous wave are handled separately.

    ``f1 = 1/(n-1)`` here (positive) -- opposite sign to :func:`bottom`.
    ``d < c*0.99`` sets ``state.error_code = 1`` / ``error_type = 4``.
    """
    wave = snow.wave
    aa = 0
    dta = 0.0
    d1 = 0.0
    d2 = 0.0
    c = 0.0
    p1 = 0.0
    p2 = 0.0
    p5 = 0.0
    diff = 0.0
    d = 0.0
    f1 = 0.0
    f2 = 0.0
    if abs(n - 1.0) > EPS:
        f1 = 1.0 / (n - 1.0)  # 1/(n-1), positive -- opposite to bottom()
    if n > EPS:
        f2 = 1.0 / n

    if wn == 1:  # single wave in the pack
        if (abs(t) > EPS and abs(1.0 - f2) > EPS) and abs(wave[1, 3] * f1) > EPS:
            d = kappa * _fpow(wave[1, 3] * f1, 1.0 - f2) * _fpow(t, f2)
    else:  # more than one wave
        aa = -1

        if abs(wave[wn - 1, 5] + 1.0) <= EPS:  # no freezing
            dta = wave[wn - 1, 1] - wave[wn, 1]

            if (abs(t) > EPS and abs(1.0 - f2) > EPS) and abs(wave[wn, 3] * f1) > EPS:
                d1 = kappa * _fpow(wave[wn, 3] * f1, 1.0 - f2) * _fpow(t, f2)
                if dta <= 0.0:
                    d = sd
                else:
                    d2 = _fpow(1.0 - _fpow(t / (t + dta), f1), f2 - 1.0)
                    d = d1 * d2
                c = kappa * _fpow(wave[wn, 3] * f1, 1.0 - f2) * _fpow(t, f2)

            # c is the depth a single wave could have; with multiple waves d must exceed it
            if d < (c * 0.99):
                state.error_code = 1
                state.error_type = 4
        else:  # freezing associated with the previous wave
            if (abs(kappa) > EPS and wave[wn - 1, 2] > EPS) and abs(f1) > EPS:
                p1 = _fpow(wave[wn - 1, 2] / kappa, n * f1)
            p2 = n - 1.0
            p5 = f1
            if abs(p5) > EPS:
                diff = p1 * p2 * _fpow(wave[wn - 1, 5], p5)
            dta = wave[wn - 1, 1] - wave[wn, 1]

            if (wave[wn, 3] > EPS and abs(diff) > EPS) and (
                abs(n - 1.0) > EPS and abs(_fpow(wave[wn, 3] / diff, n - 1.0) - 1.0) > EPS
            ):
                if (
                    t < _fpow(_fpow(wave[wn, 3] / diff, n - 1.0) - 1.0, aa) * dta
                    or diff > wave[wn, 3]
                ):
                    if (abs(t) > EPS and abs(1.0 - f2) > EPS) and abs(wave[wn, 3] * f1) > EPS:
                        d = kappa * _fpow(wave[wn, 3] * f1, 1.0 - f2) * _fpow(t, f2)
                else:
                    if (abs(t) > EPS and abs(1.0 - f2) > EPS) and abs(wave[wn, 3] * f1) > EPS:
                        d1 = kappa * _fpow((wave[wn, 3] - diff) * f1, 1.0 - f2) * _fpow(t, f2)
                        if abs(dta) <= EPS:
                            d = sd
                        else:
                            d2 = _fpow(1.0 - _fpow(t / (t + dta), f1), f2 - 1.0)
                            d = d1 * d2

    return d
