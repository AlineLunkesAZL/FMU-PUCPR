"""Constitutive-relation function library, transcribed from fasst_functions.F90.

Fortran source: fasst_functions.F90 (21 025 bytes, 680 lines, 11 functions).
See docs_transcricao/2026-08-24-doc-transcricao-fasst_functions.docx for the
full transcription rationale, numeric verification, and pending validation
work referenced throughout this module as "the transcription doc".

Scope -- 9 of 11 Fortran functions transcribed, 2 excluded as dead code
-------------------------------------------------------------------------
``dddt`` (derivative of ``dense`` w.r.t. temperature) and
``USCS_sandsiltclay`` (soil classification by the Casagrande plasticity
chart) are called by none of the 30 Fortran source files (fan-out 0) and
are deliberately NOT transcribed here (doc point 1, PENDENCIAS). Both carry
known defects that would need fixing before any use: ``dddt`` disagrees
with ``dense``'s saturation clamps and silently returns 0 for phase 4
(snow, doc point 2); ``USCS_sandsiltclay`` compares a plasticity *sum*
against a line defined for a plasticity *difference*, and has an
unreachable branch (``sand <= gravel`` instead of ``sand >= gravel``, doc
point 11). If either is ever needed, transcribe from the Fortran source
directly and fix both defects first -- do not resurrect this omission
silently.

The `phase` argument -- Fortran `f1`
-------------------------------------
``dense``, ``spheats``, and ``thconds`` select a material phase via an
integer passed as a bare local variable in the Fortran call sites (doc
point "Convenção do argumento f1"). Promoted here to the :class:`Phase`
IntEnum for readability; the underlying integer values (0-4) are
unchanged, so this is a value-preserving rename, not a behavior change.

Kelvin vs. Celsius within the same `select case` -- read before editing
-------------------------------------------------------------------------
Several functions use the *absolute* temperature (`temp1`, Kelvin) in one
phase branch and the *shifted* temperature (`t1 = temp1 - 273.15`, Celsius)
in another, within the same function (doc point 7). This is preserved
exactly; each branch below is commented with which one it uses. Do not
"normalize" these to a single unit without re-deriving the fit constants
from source (Jordan 1991 CRREL SR 91-16 p.17 for the Kelvin-linear ice
specific heat in `spheats`, in particular).

Polynomials: repeated multiplication, coefficient first
------------------------------------------------------------
The Fortran writes every polynomial term as `coef*temp1*temp1*temp1*...`
(coefficient first, then N-1 multiplies), NOT `coef*temp1**N`. gfortran at
`-O0 -ffp-contract=off` evaluates that strictly left to right, a different
rounding sequence from `pow(temp1, N)`. The transcription below keeps the
Fortran token order verbatim (`c * t * t * t * t`); an earlier draft used
`c * t**N`, which diverged 3-9 ULP from the reference driver
(tests/ref/gen_functions_ref.F90, step-4 validation 2026-09-03). Only
*integer*-power terms were converted; real exponents (`temp1**(-1.00336)`,
`(278.15 - temp1)**(-1.15)`, `t3**1.7`) stay as `**` -- Fortran uses `pow`
for those too.

Global-state parameters -- `state: FasstState`
-------------------------------------------------
``head``, ``soilhumid``, ``vap_press``, and ``maxinfiltrate`` are not pure
in the Fortran source: alongside their declared arguments they read
module-global arrays (`nsoilp`, `pheadmin`, `ntype`, `ice`, `stt`, `nnodes`,
`hm`) indexed by the node number `i` (doc "Dependência oculta de estado
global"). Since this codebase already replaced Fortran's module-COMMON
state with an explicit :class:`~fasst.state.FasstState` instance (see
fasst/state.py module docstring), the literal transcription of "read the
global" is "read `state.<field>`" -- so each of these four functions takes
`state` as an explicit first parameter. This is phase 1 (preserve the
original call shape, minus the invisible global reads); promoting them to
fully pure functions with every field passed by value is deferred to phase
2 (doc PENDENCIAS) and NOT done here.

The `anint(x*1d20)*1d-20` idiom
----------------------------------
Appears at the end of nearly every function in the Fortran source. It is
**never** guaranteed to be the identity: `x*1e20` and the closing `*1e-20`
are each an inexact float64 multiply, so the round trip perturbs the value
by up to ~1-2 ULP at *any* magnitude (measured: `_round20(962.0)` ==
961.9999999999999). What varies with magnitude is only the `anint` step:

  * `|x| >~ 2^52/1e20 = 4.5e-5`:  `x*1e20 >= 2^52`, which a float64 already
    holds with no fractional bits, so `anint` is a no-op on the scaled
    value -- the residual ~1 ULP is pure double-rounding.
  * `|x| <~ 4.5e-5`:  `x*1e20` keeps a fractional part and `anint` really
    rounds it, so the result is quantized with an ABSOLUTE step of ~1e-20
    (half-step error bound ~5e-21). This is the regime the small (~1e-5)
    sp_humid outputs live in.
  * `|x| <~ 5e-21`:  `anint(x*1e20)` rounds to 0 -- the idiom flushes the
    value to exactly 0.0.

(doc point 3.) The idiom never improves the result and sometimes worsens
it -- but it is preserved verbatim as :func:`_round20` because bit-exact
comparison against the Fortran reference output is the only reliable
transcription-error detector, and removing the idiom changes the last bit
of small outputs. Do not remove; see doc point 3 and PENDENCIAS ("Teste de
regressão do idioma anint") for the isolated, measured removal this is
deferred to.

One function -- :func:`met_date` -- uses `anint(x*1d10)*1d-10` in the
Fortran source (fasst_functions.F90 line 330), i.e. it quantizes to 1e-10,
not 1e-20 like every other function. For a decimal year near 2000 that
quantization is *active* (2000 * 1e10 = 2e13, which keeps a fractional part
below 2^52), so it drops digits past the 10th decimal place of the
fractional year.
An earlier draft of this module applied :func:`_round20` here by mistake;
the step-0 Fortran cross-check (tests/ref/smoke_ref.F90) caught the ~4e-11
divergence. It now uses :func:`_round10`. This is a transcription fix
(Python was not matching the Fortran), not one of the deferred numeric
defects.
"""

import math
from enum import IntEnum

from .constants import EPS, GRAV, RV, SDENSD, SDENSW, TREF
from .state import FasstState

__all__ = [
    "Phase",
    "dense",
    "spheats",
    "thconds",
    "met_date",
    "head",
    "soilhumid",
    "vap_press",
    "maxinfiltrate",
    "map_usda_soiltype_to_uscs",
]


class Phase(IntEnum):
    """Material phase selector -- Fortran `f1`, common to dense/spheats/thconds.

    SNOW is only a valid selector for :func:`dense` (Fortran source has no
    `case(4)` in `spheats` or `thconds`); passing it there falls through to
    the Fortran default-uninitialized behavior of 0.0, replicated below.
    """

    WATER_VAPOR = 0
    WATER = 1
    ICE = 2
    DRY_AIR = 3
    SNOW = 4


def _anint(x: float, p: int) -> float:
    """Fortran ``ANINT(x * 10**p) * 10**-p`` -- round half **away from zero**.

    NOT ``np.rint`` / ``round`` (both round half to *even*). At the 1e10 /
    1e20 scales used here the difference is reachable: e.g. a raw output of
    1.4589574548636986e-05 scaled by 1e20 is exactly ...698.5, which
    ``anint`` sends to ...699 and ``rint`` to ...698. The step-6 sp_humid
    reference cross-check (2026-09-03) exposed this as a systematic 1-ULP
    divergence on small outputs. NaN/inf pass through unchanged, as
    ``ANINT`` does. (Mirrors ``fasst/snow.py::_anint``.)
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def _round20(x: float) -> float:
    """Fortran idiom `anint(x*1d20)*1d-20` -- see module docstring."""
    return _anint(x, 20)


def _round10(x: float) -> float:
    """Fortran idiom `anint(x*1d10)*1d-10` -- used only by :func:`met_date`.

    For the values it sees (decimal years ~1900-2100) the `anint` step is
    active -- `year*1e10` keeps a fractional part -- so it genuinely rounds
    the fractional year to 10 decimal places, on top of the same
    double-rounding :func:`_round20` also has. See module docstring.
    """
    return _anint(x, 10)


# ******************************************************************************
def dense(temp1: float, wind1: float, phase: Phase) -> float:
    """Density (kg/m^3) as a function of temperature.

    Fortran: fasst_functions.F90 lines 3-136.
    Sources: air & water vapor -- http://users.wpi.edu/~icardi/PDF (100-1600K,
    standard pressure); water -- Hillel (1998) "Environmental Soil Physics",
    Noborio et al. (1996) Soil Sci Soc Am J 60: 1010-1021; ice --
    http://www.engineeringtoolbox.com.

    `temp1` is Kelvin; `wind1` is wind speed (m/s), used only for phase SNOW.
    """
    t1 = temp1 - TREF  # Celsius

    if phase == Phase.WATER_VAPOR:  # Kelvin (temp1)
        d = (
            4.192e-12 * temp1 * temp1 * temp1 * temp1
            - 1.25128e-8 * temp1 * temp1 * temp1
            + 1.45079e-5 * temp1 * temp1
            - 8.12253e-3 * temp1
            + 2.17634
        )
        d = min(1.14, max(6e-1, d))

    elif phase == Phase.WATER:  # Celsius (t1)
        d = (
            -3e-07 * t1 * t1 * t1 * t1
            + 7e-05 * t1 * t1 * t1
            - 9.92e-3 * t1 * t1
            + 8.666e-2 * t1
            + 999.81
        )
        d = min(1e3, max(962.0, d))

    elif phase == Phase.ICE:  # Celsius (t1)
        d = (
            -2e-10 * t1 * t1 * t1 * t1 * t1 * t1
            - 7e-8 * t1 * t1 * t1 * t1 * t1
            - 1e-5 * t1 * t1 * t1 * t1
            - 7e-4 * t1 * t1 * t1
            - 2.37e-2 * t1 * t1
            - 4.36e-1 * t1
            + 9.1612e2
        )
        d = min(9.257e2, max(9.162e2, d))

    elif phase == Phase.DRY_AIR:  # Kelvin (temp1)
        if temp1 > 0.0:
            d = 3.6077819e2 * temp1 ** (-1.00336)
        else:
            d = 1.0
        d = min(2.05, max(0.948, d))

    elif phase == Phase.SNOW:  # Kelvin (temp1); wet-snow-density fraction, then scaled
        t3 = max(0.0, min(wind1, 2.0))
        if 258.16 < temp1 < 278.15:
            t2 = 1.4 * (278.15 - temp1) ** (-1.15) + (8e-3 * t3**1.7 if t3 > 0.0 else 0.0)
            d = 1.0 - 9.51e-1 * math.exp(-t2) if t2 <= 5e1 else 1.0
        elif temp1 <= 258.16:
            t2 = 8e-3 * t3**1.7 if t3 > 0.0 else 0.0
            d = 1.0 - 9.04e-1 * math.exp(-t2) if t2 <= 5e1 else 1.0
        else:
            t2 = 1.4 + (8e-3 * t3**1.7 if t3 > 0.0 else 0.0)
            d = 1.0 - 9.51e-1 * math.exp(-t2) if t2 <= 5e1 else 1.0
        # Fortran carries two commented-out scale factors here (*1d3, *5.5d2)
        # alongside the active one below -- evidence the factor has changed
        # at least twice (doc point 6). SDENSW is the active one.
        d = max(SDENSD, min(SDENSW, d * SDENSW))

    else:
        # No Fortran `case default`: an out-of-range phase silently keeps
        # the function's 0d0 initializer. Replicated, not "fixed".
        d = 0.0

    return _round20(d)


# ******************************************************************************
def spheats(temp1: float, phase: Phase) -> float:
    """Specific heat (J/(kg*K)) as a function of temperature.

    Fortran: fasst_functions.F90 lines 193-246.
    Sources: air & water vapor -- http://users.wpi.edu/~ierardi/PDF; water --
    http://www.engineeringtoolbox.com; ice -- Jordan (1991) CRREL SR 91-16,
    p.17.

    `temp1` is Kelvin. `phase` SNOW is not handled (no Fortran `case(4)`)
    and falls through to 0.0, matching the Fortran source.
    """
    t1 = temp1 - TREF  # Celsius

    if phase == Phase.WATER_VAPOR:  # Kelvin (temp1)
        sh = (
            2.3888e-8 * temp1 * temp1 * temp1 * temp1
            - 6.5129e-5 * temp1 * temp1 * temp1
            + 6.6178e-2 * temp1 * temp1
            - 2.9086e1 * temp1
            + 6.6256e3
        )
        sh = min(3.26e3, max(2e3, sh))

    elif phase == Phase.WATER:  # Celsius (t1)
        sh = (
            -1e-9 * t1 * t1 * t1 * t1 * t1 * t1
            + 4e-7 * t1 * t1 * t1 * t1 * t1
            - 4e-5 * t1 * t1 * t1 * t1
            + 1.6e-3 * t1 * t1 * t1
            + 4.5e-3 * t1 * t1
            - 1.8731 * t1
            + 4210.0
        )
        sh = min(4219.0, max(4178.0, sh))

    elif phase == Phase.ICE:  # Kelvin (temp1) -- Jordan (1991) is linear in Kelvin;
        sh = -13.3 + 7.8 * temp1  # do not "fix" this to t1, see module docstring
        sh = min(2050.0, max(1389.0, sh))

    elif phase == Phase.DRY_AIR:  # Kelvin (temp1)
        sh = (
            1.9327e-10 * temp1 * temp1 * temp1 * temp1
            - 7.9999e-7 * temp1 * temp1 * temp1
            + 1.1407e-3 * temp1 * temp1
            - 4.489e-1 * temp1
            + 1.0575e3
        )
        sh = min(1.25e3, max(1e3, sh))

    else:
        sh = 0.0

    return _round20(sh)


# ******************************************************************************
def thconds(temp1: float, phase: Phase) -> float:
    """Thermal conductivity (W/(m*K)) as a function of temperature.

    Fortran: fasst_functions.F90 lines 249-300.
    Sources: air & water vapor -- http://users.wpi.edu/~icardi/PDF; water --
    Farouki (1981); ice -- http://www.engineeringtoolbox.com.

    `temp1` is Kelvin. `phase` SNOW is not handled (no Fortran `case(4)`)
    and falls through to 0.0, matching the Fortran source.

    Note (doc point 6): the water branch saturates to its floor (0.58)
    exactly at the freezing point (273.15 K) -- the regime that matters
    most for a freeze/thaw soil model. At T=273.15K the *raw* polynomial
    value (0.57037) is NOT what gets returned; the clamp is. Do not
    "simplify away" the clamp thinking it is inactive there.
    """
    t1 = temp1 - TREF  # Celsius

    if phase == Phase.WATER_VAPOR:  # Kelvin (temp1)
        tc = 8.3154e-5 * temp1 - 7.4556e-3
        tc = min(2.36e-2, max(6.95e-3, tc))

    elif phase == Phase.WATER:  # Kelvin (temp1) -- see module docstring
        tc = 1.8e-3 * temp1 + 0.0787
        tc = min(7.5e-1, max(5.8e-1, tc))

    elif phase == Phase.ICE:  # Celsius (t1)
        tc = 4e-7 * t1 * t1 * t1 + 1e-4 * t1 * t1 - 6.9e-3 * t1 + 2.2174
        tc = min(3.48, max(2.2174, tc))

    elif phase == Phase.DRY_AIR:  # Kelvin (temp1)
        tc = (
            1.5207e-11 * temp1 * temp1 * temp1
            - 4.8574e-8 * temp1 * temp1
            + 1.0184e-4 * temp1
            - 3.9333e-4
        )
        tc = min(1e-1, max(1.59e-2, tc))

    else:
        tc = 0.0

    return _round20(tc)


# ******************************************************************************
def met_date(year: float, doy: float, hr: float, minute: float) -> float:
    """Decimal calendar date.

    Fortran: fasst_functions.F90 lines 303-332.

    KNOWN DEFECTS (doc point 8) -- transcribed literally, NOT fixed:
      * The leap-year test (`year mod 4 == 0`) is incomplete: it ignores
        the century/400 exceptions, so 1900 and 2100 are (wrongly)
        classified as leap years.
      * The day-of-year divisor is off by one in both branches (367 for
        leap years which have 366 days; 366 for common years which have
        365) -- this looks deliberate (keeps the fraction strictly < 1)
        but means the decimal date is not linear in elapsed time, and
        dates in different years are not directly comparable by subtraction.
      * fasst_driver.F90 (line 1110) inverts this exact transform
        (multiplying by 367d0*tstps) to count missing met-file steps. It
        is internally consistent with this function's error. If `met_date`
        is ever corrected, the driver's inverse transform must be fixed in
        the same change, or gap detection will drift by ~0.3% (~24 hourly
        steps per year).

    Rounding: the Fortran source ends with `anint(met_date*1d10)*1d-10`
    (1e-10, not the usual 1e-20), so the fractional year is quantized to
    10 decimal places. Uses :func:`_round10`, not :func:`_round20` -- see
    module docstring.
    """
    mody = year - math.trunc(year * 2.5e-1) * 4e0

    if abs(mody) <= EPS:
        md = year + doy / 367e0 + hr / (24e0 * 367e0) + minute / (6e1 * 24e0 * 367e0)
    else:
        md = year + doy / 366e0 + hr / (24e0 * 366e0) + minute / (6e1 * 24e0 * 366e0)

    return _round10(md)


# ******************************************************************************
def head(state: FasstState, i: int, smt: float) -> float:
    """Pressure head (m) from soil water content, van Genuchten (1980).

    Fortran: fasst_functions.F90 lines 335-375.

    Not pure in the Fortran source -- reads `nsoilp(i,8..24)` and
    `pheadmin(i)` from global state; see module docstring.

    Returns `state.pheadmin[i]` (doc point 12: indistinguishable from a
    genuinely very dry soil) whenever the van Genuchten parameters
    (`nsoilp[i,10]`, `[i,11]`, `[i,12]`) are missing or zero.
    """
    nsoilp = state.nsoilp
    hd = 0.0

    if (nsoilp[i, 11] > EPS and nsoilp[i, 12] > EPS) and nsoilp[i, 10] > EPS:
        if nsoilp[i, 15] < smt < nsoilp[i, 24]:
            w1 = (smt - nsoilp[i, 8]) / (nsoilp[i, 9] - nsoilp[i, 8])  # unitless
            e1 = -1.0 / nsoilp[i, 12]
            e2 = 1.0 / nsoilp[i, 11]
            if w1 > EPS and abs(w1**e1 - 1.0) > EPS:
                hd = -(1.0 / nsoilp[i, 10]) * (w1**e1 - 1.0) ** e2 * 1e-2
        elif smt <= nsoilp[i, 15]:
            hd = state.pheadmin[i]
        elif smt >= nsoilp[i, 24]:
            hd = 0.0
    else:
        hd = state.pheadmin[i]

    return _round20(hd)  # m


# ******************************************************************************
def soilhumid(state: FasstState, i: int, ph: float, sms: float, st: float) -> float:
    """Soil relative humidity, Campbell (1985).

    Fortran: fasst_functions.F90 lines 378-422.

    Not pure in the Fortran source -- reads `ntype(i)`, `nsoilp(i,9)`,
    `ice(i)` from global state; see module docstring.

    Note: `ph <= 0` for unsaturated soil. For `ntype[i] >= 20` (concrete,
    asphalt, rock, snow, water, air) without ice, the Fortran source
    applies an unexplained factor of 5 (doc point 12) -- saturated at 1.0,
    so in practice any `sms > 0.2 * porosity` already reads as RH 1.0.
    """
    t1 = 0.0

    if state.ntype[i] < 20:
        if ph >= EPS:
            sh = 1.0
        else:
            if abs(st) > EPS:
                t1 = GRAV * ph / (RV * st)

            if abs(t1) > 7.0 or sms <= EPS:
                sh = sms / state.nsoilp[i, 9]
            else:
                sh = math.exp(t1)
    else:
        if state.ice[i] <= EPS:
            sh = 5.0 * sms / state.nsoilp[i, 9]  # unexplained factor of 5, see doc point 12
        else:
            sh = sms / state.nsoilp[i, 9]

    sh = min(1.0, max(0.0, sh))

    return _round20(sh)


# ******************************************************************************
def vap_press(state: FasstState, i: int, rh: float, ap: float) -> float:
    """Vapor pressure (Pa), Tetens' formula.

    Fortran: fasst_functions.F90 lines 425-462.

    Not pure in the Fortran source -- reads `stt(i)`, `ice(i)`, `nnodes`,
    `hm` from global state; see module docstring. Despite its declared
    arguments, this function is temperature-dependent through `state.stt[i]`
    -- invisible at the call site (doc "Consequência prática").

    The 610.78 constant duplicates `constants.VPSAT0` (same value, two
    independent declarations in the Fortran source -- doc point 9); kept
    literal here for transcription fidelity, not imported, so the
    duplication stays visible until the phase-2 Tetens unification
    (doc PENDENCIAS) resolves it. The same (a, b) coefficient pairs are
    also duplicated in sp_humid.F90.
    """
    if (state.ice[i] <= 0.0 and i <= state.nnodes) or (i > state.nnodes and state.hm <= EPS):
        a, b = 17.269, 35.86  # over water; no snow
    else:
        a, b = 21.8745, 7.66  # over ice/snow

    t1 = state.stt[i] - b  # K
    if abs(a * (state.stt[i] - TREF) / t1) > 5e1 or abs(t1) <= EPS:
        vp = ap * 1e2 * rh  # Pa
    else:
        vp = 610.78 * rh * math.exp(a * (state.stt[i] - TREF) / t1)  # Pa

    return _round20(vp)


# ******************************************************************************
def maxinfiltrate(state: FasstState, i: int, smi: float) -> float:
    """Maximum infiltration rate / sortivity (m^2/s).

    Fortran: fasst_functions.F90 lines 465-542.

    Not pure in the Fortran source -- reads `nsoilp(i,7..12)` from global
    state; see module docstring.

    One verified, deliberately UNCORRECTED defect (doc point 5):
      * Quadrature bias: the Riemann sum's first slice uses `dw = w(1) - 0`
        instead of `w(1) - w(0)`, over-weighting a spurious initial slice
        proportional to `w(0)`; and `do while(j < n)` stops at j=n-1,
        omitting the j=n endpoint (theta = porosity). Both are systematic
        biases, not noise -- they do not shrink with larger n. Replicating
        them (not fixing them) is the phase-1 choice, because this
        function sets the maximum infiltration rate in new_profile.F90 and
        a "fix" changes the surface water balance in a way that is hard to
        audit mid-transcription.

    NOT a defect (doc point 4 revisited, step-4 validation 2026-09-03): the
    Fortran `theta = smi + (nsoilp(i,9) - smi)*real(j)/real(n)` was thought
    to compute the index ratio in single precision. It does not. Operator
    precedence binds it as `((nsoilp(i,9)-smi) * real(j)) / real(n)`, so
    the multiply and divide run in float64; `real(j)` / `real(n)` only
    convert `j <= 149` and `n = 150` to float32, and integers below 2**24
    are exact in float32. So there is no precision loss. An earlier draft
    wrapped this in `np.float32(...)`, which *introduced* a ~2.7e-7
    relative error the Fortran never had -- removed below. `gen_functions_ref.F90`
    now matches to <= 2 ULP on maxinfiltrate.
    """
    nsoilp = state.nsoilp
    n = 150

    e1 = 1.0 / nsoilp[i, 12] if nsoilp[i, 12] > EPS else 0.0
    e2 = 1.0 / nsoilp[i, 11] if nsoilp[i, 11] > EPS else 0.0

    j = 1
    wold = 0.0
    s2 = 0.0
    while j < n:
        # left-associative float64, exact integer operands -- doc point 4 (see docstring)
        theta = smi + (nsoilp[i, 9] - smi) * j / n
        w = (theta - nsoilp[i, 8]) / (nsoilp[i, 9] - nsoilp[i, 8])
        w = min(max(EPS, w), 1.0)

        dw = w - wold
        if w > EPS and e1 > EPS:
            i1 = w ** (-e1 - 5e-1) if abs(-e1 - 5e-1) > EPS else 0.0

            if (nsoilp[i, 12] > EPS and abs(1.0 - w**e1) > EPS) and abs(
                1.0 - (1.0 - w**e1) ** nsoilp[i, 12]
            ) > EPS:
                i2 = (1.0 - (1.0 - w**e1) ** nsoilp[i, 12]) ** 2.0
            else:
                i2 = 1.0

            if abs(w ** (-e1) - 1.0) > EPS and abs(e2 - 1.0) > EPS:
                i3 = (w ** (-e1) - 1.0) ** (e2 - 1.0)
            else:
                i3 = 0.0
        else:
            i1 = 1.0
            i2 = 1.0
            i3 = 0.0

        integral = i1 * i2 * i3 * dw
        wold = w
        s2 = s2 + integral

        j = j + 1

    integral = s2
    if nsoilp[i, 10] * nsoilp[i, 11] * nsoilp[i, 12] > EPS:
        s2 = 2.0 * nsoilp[i, 7] * integral / (nsoilp[i, 10] * nsoilp[i, 11] * nsoilp[i, 12])  # cm^2/s

    return _round20(s2 * 1e-4)  # m^2/s


# ******************************************************************************
# LIS/STATSGO soil types:
#     1 = sand                2 = loamy sand          3 = sandy loam
#     4 = silt loam           5 = silt                6 = loam
#     7 = sandy clay loam     8 = silty clay loam     9 = clay loam
#    10 = sandy clay         11 = silty clay         12 = clay
#    13 = peat               14 = open water         15 = bedrock
#
# FASST/USCS soil types (SEDRIS EDCS_AC_SOIL_TYPES, stype(maxl)):
#   0 = unknown     5 = SW     9 = ML     12 = CH       15 = PT
#   1 = GW          6 = SP    10 = CL     13 = MH       16 = MC (SMSC) nonSEDRIS
#   2 = GP          7 = SM    11 = OL     14 = OH       17 = CM (CLML)
#   3 = GM          8 = SC                              18 = EVaporites (not used)
#  20 = COncrete   21 = ASphalt  Note: both of these are nonSEDRIS
#  25 = ROck       30 = SNow     Note: both of these are nonSEDRIS
#  26 = WAter                    Note: this is nonSEDRIS
#  27 = AIr                      Note: this is nonSEDRIS
#
# "Me" table -- active in the Fortran source (data statement, uncommented).
_USDA_TO_USCS = (6, 7, 7, 9, 9, 9, 8, 10, 10, 8, 12, 12, 15, 26, 25)

# "GSL" table -- present in the Fortran source only as a commented-out
# alternate `data` statement, differing in 7 of 15 entries. Which of the
# two is authoritative is an open, UNRESOLVED scientific question (doc
# point 10 / PENDENCIAS: "pergunte ao autor qual é a vigente") -- kept here
# only for traceability, never wired up.
# _USDA_TO_USCS_GSL = (6, 7, 8, 7, 9, 10, 8, 10, 10, 8, 10, 12, 14, 26, 3)


def map_usda_soiltype_to_uscs(lis: int) -> int:
    """Map a LIS/STATSGO USDA soil type index to a FASST/USCS soil type.

    Fortran: fasst_functions.F90 lines 545-586, `map_USDA_SoilType_to_USCS`
    (renamed to snake_case here, per this codebase's Python naming
    convention -- see fasst/state.py module docstring on case sensitivity).

    Deliberate behavior DIVERGENCE from the Fortran source (doc point 10):
    the Fortran `convert(lis)` lookup has no upper-bound check on `lis`; with
    `-fcheck=all` (active in the project makefile) an out-of-range `lis`
    aborts the program, but without it the read is out-of-bounds memory.
    Here it raises `ValueError` instead of silently indexing past the table
    or reading garbage -- an explicit, documented change in behavior, not a
    silent one.
    """
    if lis <= 0:
        return lis
    if lis > len(_USDA_TO_USCS):
        raise ValueError(f"USDA soil type {lis} out of range (expected 1..{len(_USDA_TO_USCS)})")
    return _USDA_TO_USCS[lis - 1]
