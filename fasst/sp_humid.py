"""Moist-air thermodynamics, transcribed from sp_humid.F90.

Fortran source: sp_humid.F90 (7 361 bytes, 189 lines, 1 subroutine). See
docs_transcricao/doc-transcricao-sp_humid.docx for the full
transcription rationale, numeric verification, and pending validation work
referenced throughout this module as "the transcription doc".

Given atmospheric pressure, air temperature, relative humidity, and a
pressure head, this computes vapor and saturation vapor pressure, vapor
and dry-air density, mixing ratio, wet-bulb temperature, and the
volumetric vapor content `thvc` used by new_profile.F90's vapor transport,
together with several temperature/head derivatives. This is the
thermodynamic closure of the surface boundary conditions; it is called
from module_canopy.F90, module_lowveg.F90, new_profile.F90, open_water.F90,
and surfenergy.F90.

Purity: unlike the four impure functions in fasst/functions.py, this
routine reads no `FasstState` fields -- its only external dependencies are
:func:`fasst.functions.dense` and seven module-global *constants*
(`eps, grav, Rv, Rd, Tref, vpsat0`, all already in fasst/constants.py).
:func:`sp_humid` below is therefore a pure function of its five arguments.

Returning 23 quantities, not 10 -- doc point 1
-------------------------------------------------
The Fortran subroutine declares 10 `intent(out)` arguments but computes and
discards 13 more physically-named quantities internally, 6 of which are
still routed through the closing `anint(x*1d20)*1d-20` block (doc point 9)
despite never being returned -- strong evidence they were `intent(out)` in
an earlier version whose signature was trimmed without cleaning the body.
Per the doc's recommendation (PENDENCIAS: "Recomendo o dataclass -- custo
nulo, ganho de diagnóstico"), :func:`sp_humid` returns a single
:class:`MoistAirState` with every one of those 23 fields, grouped and
commented below by which of the three tiers they fall into. This is a
value-preserving decision (nothing here changes the arithmetic), not a
speculative addition -- every discarded field already existed in the
Fortran local-variable list.

Known, UNCORRECTED defects (all doc-verified; do not silently fix)
----------------------------------------------------------------------
* doc point 2 -- when the "boiling point" branch is taken (temperature so
  high that `a*tc/t1` exceeds 50, or `t1 == 0`), `desdt`/`dedt`/`dedh` are
  never assigned and keep their zero initializers, which silently zeroes
  every downstream temperature derivative (`drvdt`, `dradt`, `dmrdt`,
  `dthvdt`). Since `dmrdt` feeds `surfenergy.F90`'s latent-heat
  linearization, an identically-zero Jacobian there would stall an
  implicit iteration rather than diverge -- a plausible-looking wrong
  answer with no warning. The trigger requires very high temperatures;
  whether this branch is ever hit in practice is unverified (doc
  PENDENCIAS: "Instrumentar: contador...").
* doc point 3 -- `t2` (`= ap - vpress`) is only assigned inside the
  `abs(temp1) > eps` branch; if `temp1` is (corrupted) zero, `t2` keeps its
  zero initializer, the `abs(t2) > eps` guard below fails, and `mixr`,
  `sphumid`, `dqdt`, `dmrdt`, `dmrdh` all silently stay zero.
* doc point 4 -- `rhoda` is clamped to `[0.95, 2.8]`, but `dradt`/`dradh`
  are the derivatives of the *unclamped* value. Wherever `rhoda` actually
  saturates, the returned derivative is inconsistent with the returned
  value -- and `dradt` feeds `dthvdt`, an output new_profile.F90 consumes.
* doc point 5 -- QUANTIFIED wet-bulb scale error. `delta`/`deltad` use
  `4.099e6`; the correct Clausius-Clapeyron-via-Tetens coefficient is
  `a*(Tref-b) = 17.269*237.29 = 4.098e3` -- a 1000x error. The psychrometric
  term `6.6*apres` (with `apres` in mbar) evaluates to roughly 100x the
  correct psychrometric constant `gamma ~= 6.65e-4*p` (p in Pa, ~66 Pa/K at
  1013 hPa). Additionally `dewpt = max(temp1, ...)` should physically be
  `min` (Td <= T always) -- verified numerically to make `dewpt` identical
  to `temp1` in every case tested, i.e. the entire dew-point calculation is
  inert. Net effect measured at 1013.25 hPa: the returned wet-bulb
  depression is 300x-500x too small (e.g. at 20C/RH50%, -0.010 K instead
  of the correct -5.511 K) -- `wetbulb` is, in practice, just `temp1`. This
  propagates into module_canopy.F90 (precipitation heat flux driving
  temperature difference, line 1237; `dense`/`spheats` evaluated at
  `wetbulba`, lines 987/992/1167). TRANSCRIBED EXACTLY, not corrected: doc
  point 5 explicitly recommends filing this with the measured numbers
  rather than fixing it during transcription, since a fix changes the
  canopy energy balance and requires field-data revalidation, not just a
  Fortran comparison.
* doc point 6 -- `thvc`'s normalizing density `t5` was changed from a
  commented-out `rhov + mixr*rhoda` (~1e-2 scale, giving `thvc ~ 0.5`) to
  the active `rhow + mixr*rhoda` where `rhow = dense(temp1, water)` is
  liquid-water density (~1e3 scale, giving `thvc ~ 1e-5`) -- five orders of
  magnitude apart, undocumented in the Fortran source. Evaluating liquid
  water density at the *air* temperature to normalize a vapor content is
  physically questionable. `thvc`/`dthvdt`/`dthvdh` feed new_profile.F90's
  vapor transport. Transcribed as the active version; doc PENDENCIAS flags
  confirming with the author whether this was intentional before trusting
  it.
* doc point 7 -- `0.622` (`= Rd/Rv`, commented as such) is a literal that
  differs from the imported `Rd/Rv` by a relative 5.5e-5 -- two sources of
  truth for the same ratio (appears 4x). The Tetens (a, b) pairs
  (17.269/35.86 water; 21.8745/7.66 ice) duplicate the identical pairs in
  fasst/functions.py's `vap_press`. `vpsat0` is correctly imported and used
  once (line 90) but `610.78` is hand-written again at line 132 -- an
  inconsistency within this same Fortran file. All kept literal; unifying
  is phase-2 work (doc PENDENCIAS).
* doc point 8 -- `dqdt` (specific humidity derivative) uses denominator
  `ap`, the usual `q ~= 0.622*e/p` approximation; the *exact* form for
  `q = r/(1+r)` would use a different denominator. `dqdt` is a discarded
  output today so this doesn't affect any live result, but
  module_lowveg.F90 line 136 receives this subroutine's `dmrdt` position
  into a caller-side variable literally named `dqdtf` -- a naming
  mismatch between what the caller calls it and what is actually passed,
  worth auditing (doc PENDENCIAS).
* doc point 10 -- `rh2` is a FRACTION (0-1) here, but `fasst_global`'s
  `dmet1(iw,5)` (relative humidity) is stored as a PERCENTAGE
  (fasst_main.F90 line 164 clamps it to `<= 100`). Every one of the 5
  call sites needs auditing individually to confirm which scale it passes
  -- a factor-of-100 error here would survive an entire transcription
  undetected (doc PENDENCIAS).

A Python-specific translation note (not a Fortran defect)
--------------------------------------------------------------
The Fortran guard `if(dabs(a*tc/t1) > 5d1 .or. dabs(t1) <= eps)` computes
`a*tc/t1` *before* checking whether `t1` is zero. Fortran does not raise on
this (IEEE 754 float division by zero yields +-Infinity, and
`dabs(Infinity) > 50` is true, routing to the same branch the `dabs(t1) <=
eps` clause would have chosen anyway) -- but Python's `/` operator raises
`ZeroDivisionError` on `x / 0.0`. `t1 == 0` (`temp1` exactly `35.86` K over
water or `7.66` K over ice) is a MANDATORY boundary case in the doc's
validation grid, so this is not a corner that can be left to crash. Below,
the two clauses are evaluated in the opposite order
(`abs(t1) <= EPS or abs(a * tc / t1) > 5e1`) so Python's left-to-right
short-circuiting skips the division exactly when Fortran's IEEE arithmetic
would have produced Infinity anyway -- same boolean result on every input,
no crash.

Four OTHER raw divisions are unguarded in the Fortran source itself
(`dRHdt`, `dRHdh` at the very top of the routine -- unconditional, before
even the `ic` phase branch; `dqdt`'s `0.622/ap`; and `wetbulb`'s final
`.../(deltap + 6.6*apres)`). Fortran silently turns their zero-denominator
case into +-Infinity or NaN and keeps going; `temp1 == 0` is in fact one of
the doc's MANDATORY boundary cases (doc point 3), so `dRHdt`/`dRHdh` must
not crash on it. `apres == 0` is outside the doc's validation grid, but
`dqdt` and `wetbulb`'s divisions are given the same treatment for
consistency -- there is no reason for one unguarded division to crash and
another not to. All four use the :func:`_fdiv` helper below instead of
`/`, which reproduces plain IEEE-754 float division (silent `inf`/`nan`,
no exception) rather than Python's `ZeroDivisionError`-raising `/`. Every
other division in this function sits inside its own `if abs(x) > eps`
branch already, exactly as in the Fortran source, and is left as plain `/`
since that guard already makes a zero denominator provably unreachable
there.
"""

import math
from dataclasses import dataclass

import numpy as np

from .constants import EPS, GRAV, RD, RV, TREF, VPSAT0
from .functions import Phase, dense

__all__ = ["MoistAirState", "sp_humid"]


def _round20(x: float) -> float:
    """Fortran idiom `anint(x*1d20)*1d-20`. Full analysis in
    fasst/functions.py's module docstring ("The `anint(x*1d20)*1d-20`
    idiom"); in brief:

      * never the identity -- the two multiplies are inexact, so the round
        trip costs ~1-2 ULP at any magnitude (`_round20(962.0)` ==
        961.9999999999999);
      * for `|x| >~ 4.5e-5` (`x*1e20 >= 2^52`) the `anint` is a no-op on
        the scaled value and only that double-rounding remains;
      * for `|x| <~ 4.5e-5` the `anint` actively rounds, quantizing the
        result with an absolute ~1e-20 step -- this is where sp_humid's
        small (~1e-5) `mixr` / `dmrdt` / `rhov` outputs sit.

    Duplicated rather than imported so each transcribed file stays a
    self-contained record of its own source file's idioms.

    `anint` rounds half AWAY FROM ZERO -- not `np.rint` (half to even). The
    step-6 sp_humid reference cross-check (2026-09-03) showed the two
    disagree by 1 ULP whenever `x*1e20` lands exactly on a `.5` tie, which
    is common for the small (~1e-5) `mixr` / `dmrdt` / `rhov` outputs.
    """
    if not math.isfinite(x):
        return x
    xs = x * 1e20
    return math.copysign(math.floor(abs(xs) + 0.5), xs) * 1e-20


def _fdiv(x: float, y: float) -> float:
    """Plain IEEE-754 float division: `x / 0.0` yields +-inf or nan instead
    of raising, matching what Fortran's arithmetic silently does at the
    three call sites below that have no `eps` guard in the source. See
    module docstring, "A Python-specific translation note".
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


@dataclass
class MoistAirState:
    """Return value of :func:`sp_humid` -- all 23 quantities the Fortran
    subroutine computes, not just its 10 `intent(out)` arguments. See
    module docstring, "Returning 23 quantities, not 10".
    """

    # -- primary outputs (Fortran intent(out), same order) ------------------
    mixr: float      # mixing ratio [kg/kg]
    dmrdt: float      # d(mixr)/dtemp1 [1/K]
    vpress: float      # vapor pressure [Pa]
    wetbulb: float      # wet-bulb temperature [K] -- see doc point 5, effectively == temp1
    vpsat: float      # saturation vapor pressure [Pa]
    rhov: float      # water vapor density [kg/m^3]
    rhoda: float      # dry air density [kg/m^3], clamped to [0.95, 2.8]
    thvc: float      # volumetric vapor content, normalized by liquid water density -- see doc point 6
    dthvdt: float      # d(thvc)/dtemp1
    dthvdh: float      # d(thvc)/dh

    # -- computed and discarded in the Fortran source (doc point 1), but
    # still routed through the closing anint(x*1d20)*1d-20 block (doc point 9,
    # lines 172-187) despite never being returned --------------------------
    sphumid: float      # specific humidity [kg/kg] -- NOT the same quantity as mixr
    dqdt: float      # d(sphumid)/dtemp1 [1/K] -- inexact denominator, see doc point 8
    drvdt: float      # d(rhov)/dtemp1 [kg/(m^3*K)]
    rhovs: float      # saturated water vapor density [kg/m^3]
    drvsdt: float      # d(rhovs)/dtemp1 [kg/(m^3*K)]
    dradt: float      # d(rhoda)/dtemp1 -- inconsistent with saturated rhoda, doc point 4

    # -- computed and discarded, never reach the anint block either --------
    desdt: float      # d(vpsat)/dtemp1 [Pa/K]
    dedt: float      # d(vpress)/dtemp1 [Pa/K]
    dedh: float      # d(vpress)/dh [Pa/m]
    drvdh: float      # d(rhov)/dh [kg/m^4]
    dradh: float      # d(rhoda)/dh [kg/m^4]
    dmrdh: float      # d(mixr)/dh [1/m] -- consumed internally to compute dthvdh
    dewpt: float      # dew point [K] -- inert, always == temp1, see doc point 5


def sp_humid(ic: int, apres: float, temp1: float, rh2: float, ph: float) -> MoistAirState:
    """Moist-air thermodynamic closure.

    Fortran: sp_humid.F90, subroutine `sp_humid`. Formulation for `vpress`
    comes from the WES program `teten.c`. No subroutines called.

    Args (doc point 10 -- units and the argument convention `ic`):
        ic: 0 = over water; anything else = over ice/snow. Callers pass
            an already-computed phase flag (`d0i`/`f1i` at the call sites).
        apres: air pressure (mbar). Converted internally to Pa; used raw
            (still in mbar) in the wet-bulb psychrometric term -- see doc
            point 5 -- an inconsistency preserved from the Fortran source.
        temp1: air temperature (K).
        rh2: relative humidity as a FRACTION in [0, 1], NOT a percentage --
            see doc point 10 and its warning about `dmet1(iw,5)` being
            stored as a percentage elsewhere in the model. Callers must be
            individually audited.
        ph: pressure head (m); negative in unsaturated soil.

    Returns:
        A :class:`MoistAirState` with all 23 computed quantities.
    """
    if ic == 0:  # over water
        a, b = 17.269, 35.86
    else:  # over ice/snow
        a, b = 21.8745, 7.66

    ap = apres * 1e2  # Pa
    t1 = temp1 - b  # K
    tc = temp1 - TREF  # C

    drhdt = _fdiv(-ph * GRAV, RV * temp1 * temp1)  # 1/K (dRH/dtemp1) -- unguarded in Fortran, see docstring
    drhdh = _fdiv(GRAV, RV * temp1)  # 1/m (dRH/dh) -- unguarded in Fortran, see docstring

    # -- vapor pressures and their derivatives -------------------------------
    desdt = 0.0
    dedt = 0.0
    dedh = 0.0
    # Clause order swapped vs. the Fortran source to avoid a Python
    # ZeroDivisionError at t1 == 0 -- see module docstring, "A
    # Python-specific translation note".
    if abs(t1) <= EPS or abs(a * tc / t1) > 5e1:
        vpress = ap  # boiling point
        vpsat = vpress / rh2 if rh2 > EPS else vpress
        # desdt, dedt, dedh stay 0.0 here -- DEFECT, preserved, see doc point 2
    else:
        vpsat = min(ap, VPSAT0 * math.exp(a * tc / t1))  # Pa (saturation vapor pressure)
        vpress = min(vpsat, rh2 * vpsat)  # Pa (vapor pressure)

        desdt = vpsat * a * (1.0 / t1 - tc / (t1 * t1))  # Pa/K (dvpsat/dtemp1)
        dedt = vpress * (drhdt + a * (1.0 / t1 - tc / (t1 * t1)))  # Pa/K (dvpress/dtemp1)
        dedh = drhdh * vpress  # Pa/m (dvpress/dh)

    # -- dry air and vapor densities and their derivatives -------------------
    t2 = 0.0
    rhovs = 0.0
    rhov = 0.0
    drvsdt = 0.0
    drvdt = 0.0
    drvdh = 0.0
    dradt = 0.0
    dradh = 0.0
    if abs(temp1) > EPS:
        c1 = 1.0 / (RV * temp1)  # kg/J
        c2 = 1.0 / (RD * temp1)  # kg/J
        t2 = ap - vpress  # Pa

        rhovs = c1 * vpsat  # kg/m^3 (saturated water vapor density)
        rhov = c1 * vpress  # kg/m^3 (water vapor density)

        drvsdt = c1 * (desdt - vpsat / temp1)  # kg/(m^3*K) (drhovs/dtemp1)
        drvdt = c1 * (dedt - vpress / temp1)  # kg/(m^3*K) (drhov/dtemp1)
        drvdh = c1 * dedh  # kg/m^4 (drhov/dh)

        rhoda = min(2.8, max(0.95, t2 * c2))  # kg/m^3 (dry air density)
        dradt = -c2 * (dedt + t2 / temp1)  # kg/(m^3*K) (drhoa/dtemp1)
        dradh = -c2 * dedh  # kg/m^4 (drhoa/dh)
    else:
        rhoda = 0.95
        # t2 stays 0.0 here -- DEFECT, preserved, see doc point 3

    # -- specific humidity and mixing ratio and their derivatives -----------
    # NOTE: 0.622 = Rd/Rv (literal, not the imported ratio) -- see doc point 7
    mixr = 0.0
    sphumid = 0.0
    dqdt = 0.0
    dmrdt = 0.0
    dmrdh = 0.0
    if abs(t2) > EPS:
        mixr = 0.622 * vpress / t2  # unitless [kg/kg] (mixing ratio)
        sphumid = mixr / (1.0 + mixr)  # unitless [kg/kg] (specific humidity)

        dqdt = _fdiv(0.622, ap) * dedt  # 1/K (dsphumid/dtemp1) -- unguarded in Fortran; see doc point 8
        dmrdt = dedt * 0.622 * (1.0 / t2 + vpress / (t2 * t2))  # 1/K (dmr/dtemp1)
        dmrdh = dedh * 0.622 * (1.0 / t2 + vpress / (t2 * t2))  # 1/m (dmr/dh)

    # -- dew point ------------------------------------------------------------
    if vpress <= 1e-15:
        vpress = 0.0  # prevent program crashes
    if abs(vpress) > EPS and rh2 > EPS:
        t3 = math.log(vpress / (rh2 * 610.78))  # unitless -- 610.78 duplicates VPSAT0, see doc point 7
    else:
        t3 = 0.0

    if abs(t3 - a) > EPS:
        # dmax1 here should physically be dmin1 (Td <= T always) -- inert
        # bug, preserved exactly. See doc point 5.
        dewpt = max(temp1, (t3 * b - TREF * a) / (t3 - a))  # K (dewpt temperature)
    else:
        dewpt = temp1

    # -- wet bulb temperature -- see doc point 5 for the verified scale error
    delta = 4.099e6 * vpsat / (t1 * t1) if abs(t1) > EPS else 0.0  # Pa
    t4 = dewpt - b  # K
    deltad = 4.099e6 * vpress / (t4 * t4) if abs(t4) > EPS else 0.0  # Pa
    deltap = (delta + deltad) * 5e-1  # Pa

    wetbulb = temp1 - _fdiv(vpsat - vpress, deltap + 6.6 * apres)  # K (wetbulb temperature) -- unguarded in Fortran, see docstring

    # -- dthetav/dT and dthetav/dh -- see doc point 6 on the t5 rescale ------
    rhow = dense(temp1, 0.0, Phase.WATER)  # density of LIQUID WATER (~1000 kg/m^3)
    t5 = rhow + mixr * rhoda
    thvc = 0.0
    dthvdt = 0.0
    dthvdh = 0.0
    if abs(t5) > EPS:
        thvc = mixr * rhoda / t5
        dthvdt = (dmrdt * rhoda + dradt * mixr - thvc * (dmrdt * rhoda + dradt * mixr)) / t5
        dthvdh = (dmrdh * rhoda + dradh * mixr - thvc * (dmrdh * rhoda + dradh * mixr)) / t5

    return MoistAirState(
        mixr=_round20(mixr),
        dmrdt=_round20(dmrdt),
        vpress=_round20(vpress),
        wetbulb=_round20(wetbulb),
        vpsat=_round20(vpsat),
        rhov=_round20(rhov),
        rhoda=_round20(rhoda),
        thvc=_round20(thvc),
        dthvdt=_round20(dthvdt),
        dthvdh=_round20(dthvdh),
        sphumid=_round20(sphumid),
        dqdt=_round20(dqdt),
        drvdt=_round20(drvdt),
        rhovs=_round20(rhovs),
        drvsdt=_round20(drvsdt),
        dradt=_round20(dradt),
        # not rounded in the Fortran source -- preserved as-is, see docstring
        desdt=desdt,
        dedt=dedt,
        dedh=dedh,
        drvdh=drvdh,
        dradh=dradh,
        dmrdh=dmrdh,
        dewpt=dewpt,
    )
