"""Soil/snow temperature profile (Crank-Nicholson), transcribed from
soil_tmp.F90.

Fortran source: soil_tmp.F90 (29 998 bytes, 993 lines, single subroutine).
Called once per inner-convergence iteration from ``new_profile.F90:1209``
(``Call soil_tmp(...)`` -- see :mod:`fasst.new_profile` docstring,
"Injected untranscribed dependencies", ``soil_tmp_fn``). This module
implements that contract directly: no subroutines are called from here
(header comment, verified -- there genuinely are none), and every function
it uses (``dense``, ``spheats``, ``head``, ``soilhumid``, ``vap_press``) is
already transcribed in :mod:`fasst.functions`. **This is therefore the
first file in the ``new_profile`` dependency chain whose own numeric
validation is not blocked by a missing dependency** -- see VALIDATION
below.

What :func:`soil_tmp` does
----------------------------
For the active node range ``1..ntemp`` (soil, snow, mixed and vegetation
nodes; air nodes, ``ntype == 27``, are explicitly skipped) it (1) builds a
Crank-Nicholson tridiagonal system for the temperature INCREMENT
``delstt(i)`` at each node, folding in conductive exchange with neighbors,
phase-change energy (freeze/thaw, evaporation/condensation), advective
heat from water/vapor flow, and the surface energy-balance linearization
(``isurfg``/``isurff`` and their ``T``-derivatives) at the top node(s);
(2) solves that system with the Thomas algorithm, in up to TWO independent
segments when an air node splits the active range (see point 1); (3) walks
back down from ``ntemp`` to ``1``, applying the increment with several
independent magnitude clamps (a per-iteration cap ``dtopi``, a
"continuity in time" cap ``maxdel``, a "continuity in space" cap against a
reference node); (4) at soil/water nodes (``ntype(i) <= 25``), partitions
the resulting temperature change into ice/water mass exchange (conserving
``bftm(i)``, the node's fixed exchangeable water+ice budget for the
current timestep) and updates the node's water-vapor content
(``rhov``/``rhoda``/``wvc``) via :func:`~fasst.functions.vap_press`.

Fortran line ranges::

    1-58       declarations
    59-186     zero-out; phase constants; `sn` (top node) by icase; `fs`,
               `tf2` (fusion vs. sublimation latent heat)
    188-576    MAIN LOOP over i=1..ntemp (skipping ntype==27): phase-change
               energy terms, tridiagonal coefficients A/B/C/D, surface
               linearization by icase, `rhs_errort`/`sert` tracking
    579-661    Thomas algorithm, up to 2 segments (is1..if1, is2..if2)
    663-991    walk-back loop i=ntemp..1: apply + clamp `delstt`,
               continuity checks, ice/water mass partition, water-vapor
               update

Uses the functions :func:`~fasst.functions.dense`,
:func:`~fasst.functions.spheats`, :func:`~fasst.functions.head`,
:func:`~fasst.functions.soilhumid`, :func:`~fasst.functions.vap_press`
(all already transcribed; header comment verified accurate here, unlike
new_profile's and surfenergy's).

Not pure -- reads AND writes ``state``
---------------------------------------
Reads: ``ntot, nnodes, icase, step, mflag, elev, nz, nsoilp, ntype,
node_type, mstflag, ntemp, stt, too, wvc, soil_moist, ice, phead, bftm,
hfol, hm, fv1, flowu, flowl, grthcond, grspheat, dsmdh, source, sink,
timstep, tmelt, istart, isurfoldf, isurfoldg, iw, meltfl, deltat_fasst,
refreeze`` plus the constants ``EPS, GRAV, LHFUS, LHSUB, RD, RV, TREF``.
**Writes**: ``state.stt, state.wvc, state.ice, state.soil_moist,
state.phead, state.bftm`` (all per-node, indices ``1..ntemp`` or
``1..nnodes`` depending on the block) and ``state.tmelt[state.iw]``. Also
mutates its own ``rhov``/``rhoda`` array PARAMETERS in place
(``intent(inout)`` in the Fortran) -- these are the caller's transient
per-timestep locals (:mod:`fasst.new_profile`'s own ``rhov``/``rhoda``),
not :class:`FasstState` fields; ``FasstState`` has no ``rhov``/``rhoda`` of
its own.

Point 1 -- up to two independent Thomas-algorithm segments around an
air-node gap
------------------------------------------------------------------------
``sn2`` (index of the FIRST ``ntype(i) == 27`` node found scanning
``i=1..ntemp``) and ``sn3`` (index of the LAST such node) partition the
active range into up to two solve segments (``is1..if1``, ``is2..if2``),
skipping the air block in between. The 5-way branch (lines 582-599)
assumes any air nodes present form a single contiguous block between
``sn2`` and ``sn3`` -- not verified against ``node_type`` assignment
elsewhere, just replicated. Inside the sweep, a node with ``ntype(j) ==
27`` gets ``delstt(j) = 0`` and neither ``gam(j)`` nor ``bet`` are updated
for it; since ``A(j) = B(j) = C(j) = 0`` for any such node (they are never
re-assigned away from the initial zero-out for ``ntype(i) == 27``, only
``D(i)`` is redundantly re-zeroed at line 572), the NEXT real node's
``gam = C(air)/bet = 0`` correctly decouples it from the segment on the
other side of the gap. Consistent, but easy to break by "simplifying" the
skip logic -- preserved literally, including the harmless duplicate
``is2 = 0; if2 = 0`` (already zeroed once, lines 69-70; re-zeroed again at
lines 580-581 right before the branch tree).

Point 2 -- ``y = 60.0`` unconditionally overwrites a 4-branch computation
------------------------------------------------------------------------
Lines 743-753: ``y`` is computed via a 4-way branch tree (elevation-based,
snow-depth-based, foliage-height-based, or ``1 - deep``), each a genuine
piece of physics -- and then the very next line unconditionally reassigns
``y = 6d1`` (60.0), discarding all four branches. The intended formula is
still visible, commented out, immediately after:
``!dmin1(dmax1(maxdel,4d1*y),6d1)``. ``deep`` (the quartic polynomial at
lines 733-738, unrelated coefficients 4.388/-11.498/10.947/-4.7845/0.9546,
unnamed) is NOT dead -- it is the live multiplier for ``maxdel`` at line
740 (``maxdel = dmax1(maxdel*deep, mxd1*deep)``), which IS used by the
"continuity in time" clamp at line 756. Only the 4-branch ``y`` computation
is dead; ``maxdel``/``deep`` survive. Replicated literally: the 4 branches
are still computed (so a future re-activation of the commented formula
needs no rewiring), and then thrown away by the same unconditional
override. ``y`` is REUSED for an unrelated purpose (ice/water mass
exchange, kg) later in the same node's processing (point 5) -- one Python
local, reused twice, matching the Fortran's single ``y``.

Point 3 -- water/air RH exclusion: real slip, narrow consequence
--------------------------------------------------------------------
Line 959: ``if((ntype(i) /= 26.or.ntype(i) /= 27).and.node_type(i) /=
'VG') then``. The ``.or.`` clause is a tautology (always true for any
integer), so the condition reduces to plain ``node_type(i) /= 'VG'`` --
the intended exclusion of water (26) and air (27) nodes from the
``soilhumid``-based RH path never actually applies. This is a genuine
slip (almost certainly meant to be ``.and.``, mirroring the
``node_type(i) /= 'WA'.and.node_type(i) /= 'AI'`` idiom used everywhere
else in the project), but its observable damage is nil here: this whole
block already sits inside ``if(i <= nnodes.and.ntype(i) <= 25)`` (line
832), which excludes ``ntype`` 26 and 27 before the tautology is ever
evaluated. So the water/air branch (``rh = dmet(5)*1d-2``, ``if(ntype(i)
== 26) rh = 1d0``) is reachable only through the OUTER guard already doing
the exclusion job -- meaning the inner ``ntype(i) == 26`` special case
inside that branch is the part that is actually dead code (a water node
can never reach it, since the outer guard already excludes ``ntype <=
25`` -- wait, 26 is NOT <= 25, so a water node never enters this block at
all). Replicated literally (the ``.or.``, not ``.and.``) -- do not "fix"
it to ``.and.`` believing it changes behavior; it does not, here.

Point 4 -- third hand-typed copy of the WES latent-heat-of-vaporization fit
--------------------------------------------------------------------------
Line 206: ``lhe = 2500775.6d0 - 2369.729d0*(stt(i) - Tref)``. Same fit
already flagged as duplicated in ``sflux.py`` (``_LHVAP_LITERAL``) and
``surfenergy.py`` (``_LHEVAP_INTERCEPT``/``_LHEVAP_SLOPE``) -- this is now
a THIRD independent hand-written copy in the project. Kept literal here
too (``_LHEVAP_INTERCEPT``/``_LHEVAP_SLOPE`` re-declared, not imported from
either sibling module -- see those modules' own docstrings on why each
copy stays self-contained).

Point 5 -- ice/water mass partition conserves ``bftm``, computed once per
timestep
------------------------------------------------------------------------
``bftm(i)`` (a :class:`FasstState` field) is the node's total exchangeable
water+ice mass, in kg -- but it is only SET on the very first
inner-iteration of a timestep (``if(ii == 1.and.iter == 0) bftm(i) = smm +
icm``, line 837); every subsequent call within the same timestep reuses
the value fixed on that first call as the conservation target for the
ice/water repartition. Getting the ``ii``/``iter`` guard wrong (e.g.
re-deriving ``bftm`` every call) would silently change the mass-balance
target mid-timestep.

Point 6 -- ``ftempt``, ``dfwi``, ``dfiw``: computed, never read (live code)
--------------------------------------------
Lines 807-811 compute ``ftempt`` (vegetation-surface vs. air temperature,
by ``node_type``); lines 828-829 compute ``dfwi = rhow/rhoi`` and ``dfiw =
rhoi/rhow`` (water/ice density ratios). None of the three is read anywhere
in the ACTIVE code afterward, nor is any an ``intent(out)`` argument --
``dfwi``/``dfiw``'s only readers are the commented-out moisture
-reconciliation block at lines 942-954 (point 7). Same class as
new_profile's ``full`` -- computed here for structural fidelity, with no
live downstream consumer.

Point 7 -- the ``mstflag``/measured-soil-temperature override, and its
commented-out moisture reconciliation
------------------------------------------------------------------------
When ``mstflag == 1`` and node ``nnodes`` has a valid (non-missing)
measured temperature in ``dmet(13)``, ``stt(nnodes)`` is forced to that
value and ``fflag = 1`` is set (lines 698-703); ``fflag`` then re-forces
the same override at TWO further points (784, 799) after the continuity
clamps would otherwise have adjusted ``stt(nnodes))``. The block that
would reconcile ``soil_moist``/``ice`` with this forced temperature
(lines 942-954, commented out) is inactive -- so forcing the measured
temperature can leave the node's ice/water partition inconsistent with
its (now overridden) temperature. Not reconciled here either, per the
"replicate silence" decision already made for :mod:`fasst.new_profile`
(same class of decision, doc point 3 there).

Point 8 -- two unrelated fields named ``hfol`` / ``hfol_tot``
------------------------------------------------------------------
Line 730: ``d1 = dmax1(1d0, hfol, hm)``. ``state.hfol`` here is NOT
``state.hfol_tot`` (the "total foliage height" field new_profile.F90's
``icase`` topology logic uses) -- they are two distinct
:class:`FasstState` fields. Easy to conflate when cross-referencing this
file against ``new_profile.py``; kept as ``state.hfol`` literally.

Point 9 -- ``dhdT``/``dhdTo`` alternate formula, commented out
--------------------------------------------------------------
Lines 243-246 carry a commented-out alternate closed-form fit for
``dhdT``/``dhdTo`` (``(7/71.89)*(-0.1425 - 9.52d-4*(T - Tref))*phead``);
the ACTIVE formula (lines 240, 242) is the logarithmic one,
``-(Rv/grav)*log(rh)``. Not transcribed (inactive), noted for anyone
comparing against an older version of the model.

Point 10 -- two declared-but-unused locals
-------------------------------------------
``f2`` and ``tmp`` (declared line 44-56, ``f2`` zeroed at line 129) are
never read anywhere else in the file -- verified with a whole-word search,
not just a substring scan (this project has a real ``f2`` in other files
that IS used; this one, in this file, is not). Neither has a Python
counterpart below.

Point 11 -- the same rounding idiom at three different precisions
-------------------------------------------------------------------
``anint(x*1d15)*1d-15`` (most quantities: the tridiagonal coefficients,
``delstt``, ``bftm``, the ice/water mass terms, ``wvc``), ``anint(x*1d10)
*1d-10`` (``rhs_errort``, ``stt``, ``ice``, ``soil_moist``), and
``anint(x*1d20)*1d-20`` (``lhs4``/``lhs5``, ``rhov``, ``rhoda``, the mixing
-ratio fraction). Preserved per-site, not unified -- same discipline as
``new_profile.py``'s doc point 5.

Local-array sizing deviation (same as ``new_profile.py``)
-------------------------------------------------------------
The Fortran declares ``A(ntot)``, ``delstt(ntot)``, etc. as Fortran 90
AUTOMATIC arrays (sized fresh at each call from the current value of the
module variable ``ntot`` -- legal because ``ntot`` is in scope via
``use``, unlike ``new_profile.F90``'s ``x(14,maxn)``, which used a
compile-time ``parameter`` bound that ``ntot`` could exceed). There is no
latent out-of-bounds risk to work around here. They are still sized
``MAXN + EXTRAN`` below rather than dynamically at ``state.ntot``, purely
for consistency with every other node-indexed local array in this
project's Python transcription (and because ``state.ntot <= MAXN +
EXTRAN`` always holds) -- a fixed upper bound, not a fidelity concession.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

from .constants import EPS, EXTRAN, GRAV, LHFUS, LHSUB, MAXN, RD, RV, TREF
from .functions import Phase, dense, head, soilhumid, spheats, vap_press
from .state import FasstState

__all__ = ["SoilTmpResult", "soil_tmp"]

_NODESZ = MAXN + EXTRAN

# Hand-written WES latent-heat-of-vaporization fit -- third independent
# copy in the project (see module docstring, point 4).
_LHEVAP_INTERCEPT = 2500775.6  # J/kg
_LHEVAP_SLOPE = 2369.729       # J/(kg*K)


def _reals(n: int) -> np.ndarray:
    return np.zeros(n + 1, dtype=np.float64)


def _anint(x: float, p: int) -> float:
    """Fortran ``ANINT(x * 10**p) * 10**-p`` -- round half away from zero."""
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    return math.copysign(math.floor(abs(xs) + 0.5), xs) * (10.0**-p)


def _round10(x: float) -> float:
    return _anint(x, 10)


def _round15(x: float) -> float:
    return _anint(x, 15)


def _round20(x: float) -> float:
    return _anint(x, 20)


def _aint(x: float) -> float:
    """Fortran ``AINT``/``DINT`` with no scaling -- truncate toward zero."""
    if not math.isfinite(x):
        return x
    return float(math.trunc(x))


class SoilTmpResult(NamedTuple):
    """The 4 ``intent(out)`` values of ``soil_tmp.F90``.

    ``iflag``: 0/1, "is there relevant ice at any active node after this
    update" (line 955) -- unrelated to ``new_profile``'s own reuse of a
    same-named local for a 1-7 convergence diagnostic (see
    :mod:`fasst.new_profile` docstring, point 4).
    ``rhs_errort``: largest-magnitude tridiagonal residual ``D(i)``, signed
    (``sert`` carries its value; ``rhs_errort`` its magnitude).
    ``errort``: largest ``abs(delstt(i))`` over soil nodes (``i <=
    nnodes``) -- the convergence driver for ``new_profile``'s ``while ...
    errort > allowerrort`` loop.

    A :class:`NamedTuple`, not a plain ``@dataclass``, so it unpacks
    positionally -- confirmed against the only real caller,
    ``fasst.new_profile.new_profile``'s ``soil_tmp_fn`` contract (see its
    docstring), which does
    ``iceflag, rhs_errort, errort, sert = soil_tmp_fn(...)``.
    """

    iflag: int
    rhs_errort: float
    errort: float
    sert: float


# ===========================================================================
def soil_tmp(
    state: FasstState,
    isurfg: float,
    isurff: float,
    disurfg: float,
    disurff: float,
    disurffg: float,
    disurfgf: float,
    sigflo: float,
    rhotot: float,
    airo: float,
    dmet: np.ndarray,
    grthcondo: np.ndarray,
    grspheato: np.ndarray,
    zt: np.ndarray,
    delzt: np.ndarray,
    told: np.ndarray,
    iceo: np.ndarray,
    wvco: np.ndarray,
    flowuo: np.ndarray,
    flowlo: np.ndarray,
    fv1o: np.ndarray,
    sinko: np.ndarray,
    sourceo: np.ndarray,
    ph_old: np.ndarray,
    sm_old: np.ndarray,
    thvc: np.ndarray,
    dthvdt: np.ndarray,
    rhov: np.ndarray,
    rhoda: np.ndarray,
    ii: int,
    iter_: int,
) -> SoilTmpResult:
    """Advance the soil/snow temperature profile by one Crank-Nicholson step.

    Fortran: soil_tmp.F90 lines 1-993, ``subroutine soil_tmp``. See module
    docstring for structure and the numbered fidelity notes.

    Args:
        state: model state; see module docstring, "Not pure".
        isurfg, isurff, disurfg, disurff, disurffg, disurfgf: surface
            energy-balance sum and its ``T``-derivatives, from the LAST
            :func:`~fasst.surfenergy.surfenergy` call (Fortran
            ``intent(in)``) -- see :mod:`fasst.new_profile`'s
            ``soil_tmp_fn`` contract for exactly which local feeds each.
        sigflo, rhotot, airo: low-veg fraction, blended air density, and
            previous-timestep air temperature (all read-only).
        dmet: 13-element (phantom slot 0) derived-meteorology array.
        grthcondo, grspheato, zt, delzt, told, iceo, wvco, flowuo, flowlo,
            fv1o, sinko, sourceo, ph_old, sm_old, thvc, dthvdt: previous
            -timestep / previous-iteration per-node arrays (read-only).
        rhov, rhoda: per-node water-vapor / dry-air density arrays,
            ``intent(inout)`` -- mutated in place for nodes with ``i <=
            nnodes and ntype[i] <= 25`` (point 3); left untouched
            elsewhere, same as the Fortran.
        ii, iter_: pass-through sub-step / inner-iteration indices, used
            only for the ``bftm`` first-call guard (point 5).

    Returns:
        A :class:`SoilTmpResult`. Also mutates ``state`` in place (see
        module docstring, "Not pure").
    """
    ntot = state.ntot
    ntemp = state.ntemp
    nnodes = state.nnodes
    icase = state.icase
    step = state.step
    iw = state.iw

    A = _reals(_NODESZ)
    B = _reals(_NODESZ)
    C = _reals(_NODESZ)
    D = _reals(_NODESZ)
    gam = _reals(_NODESZ)
    delstt = _reals(_NODESZ)
    tcheck = _reals(_NODESZ)
    icheck = _reals(_NODESZ)
    wcheck = _reals(_NODESZ)
    scheck = _reals(_NODESZ)
    pcheck = _reals(_NODESZ)
    lhs4 = _reals(_NODESZ)
    lhs5 = _reals(_NODESZ)

    d0i = Phase.WATER_VAPOR
    d1i = Phase.WATER
    d2i = Phase.ICE
    d3i = Phase.DRY_AIR
    d4i = Phase.SNOW

    sn2 = 0
    sn3 = 0
    fflag = 0
    iflag = 0
    cflag = 0

    sn = nnodes
    if icase == 1 or icase == 4:
        sn = ntemp
    elif icase == 3:
        sn = ntemp - 1

    fs = 1.0 / float(step)
    tf2 = LHFUS
    if state.meltfl == "s":
        tf2 = LHSUB

    # -- solve for parameters needed in the tridiagonal matrix ---------------
    rhs_errort = -abs(state.mflag)
    errort = -abs(state.mflag)
    # `sert` is intent(out) in the Fortran and is only ever assigned inside
    # `if abs(D[i]) > rhs_errort`, which is false for every i if every D(i)
    # comes out NaN (e.g. a degenerate all-zero-temperature startup state) --
    # IEEE754 comparisons against NaN are always false. The Fortran would
    # then return `sert` genuinely uninitialized; Python cannot replicate
    # that, so it is seeded to 0.0 here (same precedent as
    # fasst.surfenergy's `smr`/`taf`, fasst.new_profile's `hatm` etc.).
    sert = 0.0

    for i in range(1, ntemp + 1):
        if state.ntype[i] != 27:
            q1lhsa = 0.0
            q1lhsb = 0.0
            q1lhsc = 0.0
            q1rhs = 0.0
            q2lhsa = 0.0
            q2lhsb = 0.0
            q2lhsc = 0.0
            q2rhs = 0.0
            qb = 0.0
            lhs1 = 0.0
            lhs2 = 0.0
            lhs3 = 0.0

            lhe = _LHEVAP_INTERCEPT - _LHEVAP_SLOPE * (state.stt[i] - TREF)

            f1 = state.deltat_fasst / (2.0 * delzt[i])

            rhow = dense(state.stt[i], 0.0, d1i)
            rhowo = dense(told[i], 0.0, d1i)
            rhoi = dense(state.stt[i], 0.0, d2i)

            out = 5e-1 * (
                spheats(state.stt[i], d1i) * rhow * state.sink[i]
                + spheats(told[i], d1i) * rhowo * sinko[i]
            )
            if abs(out) < EPS:
                out = 0.0

            in_ = 5e-1 * (
                spheats(state.stt[i], d1i) * rhow * state.source[i]
                + spheats(told[i], d1i) * rhowo * sourceo[i]
            )
            if abs(in_) < EPS:
                in_ = 0.0

            rhsDi = 0.0
            rhsi = 0.0
            rhsDl = 0.0
            rhsw = 0.0
            rhsDv = 0.0
            rhsv = 0.0

            sph = 5e-1 * (state.grspheat[i] + grspheato[i])

            # -- phase-change energies -------------------------------------
            if i <= nnodes and state.ntype[i] != 26:
                rh = soilhumid(state, i, state.phead[i], state.soil_moist[i], state.stt[i])
                rho = soilhumid(state, i, ph_old[i], sm_old[i], told[i])
                pres = dmet[11] + 1e-2 * ((state.elev - state.nz[i]) + state.phead[i]) * rhow * GRAV
                po = dmet[11] + 1e-2 * ((state.elev - state.nz[i]) + ph_old[i]) * rhowo * GRAV
                dhdT = 0.0
                if abs(rh) > EPS:
                    dhdT = -(RV / GRAV) * math.log(rh)
                dhdTo = 0.0
                if abs(rho) > EPS:
                    dhdTo = -(RV / GRAV) * math.log(rho)

                # freezing releases heat -> warms soil/water around it;
                # ice - iceo > 0 -> warming; ice - iceo < 0 -> cooling
                if state.stt[i] > TREF:
                    didT = -(rhow * LHFUS / (rhoi * GRAV * TREF)) * state.dsmdh[i]
                else:
                    didT = -(rhow * LHFUS / (rhoi * GRAV * state.stt[i])) * state.dsmdh[i]

                lhs5[i] = didT

                rhsDi = 0.0
                rhsi = 0.0
                if state.ice[i] + iceo[i] > EPS:
                    rhsDi = 5e-1 * (
                        (LHFUS + spheats(state.stt[i], d2i) * abs(state.stt[i] - TREF)) * rhoi
                        + (LHFUS + spheats(told[i], d2i) * abs(told[i] - TREF)) * dense(told[i], 0.0, d2i)
                    )
                    rhsi = -rhsDi * (state.ice[i] - iceo[i])
                    lhs1 = lhs5[i]

                # water warming/cooling
                if state.soil_moist[i] + sm_old[i] > EPS:
                    t1 = 5e-1 * (
                        spheats(state.stt[i], d1i) * abs(state.stt[i] - TREF) * rhow
                        + spheats(told[i], d1i) * abs(told[i] - TREF) * rhowo
                    )
                    t1 = _round15(t1)

                    Wo = GRAV * (ph_old[i] - told[i] * dhdTo)
                    Wt = GRAV * (state.phead[i] - state.stt[i] * dhdT)
                    t2 = 5e-1 * (Wt * rhow + Wo * rhowo)
                    rhsDl = _round15(t1 + t2)
                    rhsw = rhsDl * (state.soil_moist[i] - sm_old[i])
                    lhs3 = dhdT * state.dsmdh[i]

                # evaporation/condensation energy
                rhsDv = 0.0
                rhsv = 0.0
                if state.wvc[i] + wvco[i] > EPS:
                    rhovo = 0.622 * vap_press(state, i, rho, po) / (RV * told[i])
                    rhovi = 0.622 * vap_press(state, i, rh, pres) / (RV * state.stt[i])
                    rhsDv = (
                        5e-1 * (lhe + spheats(state.stt[i], d0i) * abs(state.stt[i] - TREF)) * rhovi
                        + (lhe + spheats(told[i], d0i) * abs(told[i] - TREF)) * rhovo
                    )
                    rhsv = rhsDv * (state.wvc[i] - wvco[i])

                t1 = dthvdt[i] * (state.nsoilp[i, 2] - (state.soil_moist[i] + state.ice[i]))
                t2 = -thvc[i] * (dhdT + lhs1) * state.dsmdh[i]

                lhs4[i] = t1 + t2
                if state.wvc[i] + wvco[i] > EPS:
                    lhs2 = lhs4[i]

            elif i <= nnodes and state.node_type[i] == "WA":
                # open water
                sph = 5e-1 * (state.grspheat[i] + grspheato[i])

                t1 = 5e-1 * (
                    spheats(state.stt[i], d1i) * abs(state.stt[i] - TREF) * rhow
                    + spheats(told[i], d1i) * abs(told[i] - TREF) * rhowo
                )
                rhsDl = t1
                rhsw = rhsDl * (state.soil_moist[i] - sm_old[i])
                lhs3 = 0.0

                didT = -(rhow / rhoi) * LHFUS / (state.timstep * 3.6e3)
                lhs5[i] = didT
                if state.ice[i] > EPS or iceo[i] > EPS:
                    rhsDi = 5e-1 * (
                        (LHFUS + spheats(state.stt[i], d2i) * abs(state.stt[i] - TREF)) * rhoi
                        + (LHFUS + spheats(told[i], d2i) * abs(told[i] - TREF)) * dense(told[i], 0.0, d2i)
                    )
                    rhsi = rhsDi * (state.ice[i] - iceo[i])
                    lhs1 = lhs5[i]

            elif i > nnodes:
                # snow
                if state.node_type[i] == "HM" or state.node_type[i] == "MX":
                    hs = delzt[i]

                    rhsDi = LHFUS * rhotot + state.grspheat[i] * abs(state.stt[i] - TREF)
                    rhsi = rhsDi * (-state.refreeze / hs)
                    didT = tf2 * (-state.refreeze / hs) / (hs * GRAV * TREF)
                    didT = 0.0

                    rhsi = rhsi * fs
                    didT = didT * fs
                    lhs1 = didT

                    # snow evaporation/condensation energy
                    if state.wvc[i] > EPS or wvco[i] > EPS:
                        rhsDv = 0.0
                        rhsv = 0.0

                        rhsDv = lhe * rhotot + state.grspheat[i] * abs(state.stt[i] - TREF)
                        rhsv = rhsDv * (state.wvc[i] - wvco[i]) * fs

                        lhs2 = -(rhov[i] / rhow) * didT

            lhs1 = _round20(lhs1)
            lhs2 = _round20(lhs2)
            lhs3 = _round20(lhs3)
            lhs4[i] = _round20(lhs4[i])
            lhs5[i] = _round20(lhs5[i])

            # -- interior and bottom nodes ------------------------------------
            if i <= nnodes - 1:
                if i == 1:
                    kthl = 0.0
                    kthlo = 0.0
                    if state.node_type[1] != "SN":
                        qb = 2.0 * 7.5e-2
                else:
                    kthl = (state.grthcond[i - 1] * delzt[i - 1] + state.grthcond[i] * delzt[i]) / (
                        delzt[i - 1] + delzt[i]
                    )
                    kthlo = (grthcondo[i - 1] * delzt[i - 1] + grthcondo[i] * delzt[i]) / (
                        delzt[i - 1] + delzt[i]
                    )

                kthu = (state.grthcond[i + 1] * delzt[i + 1] + state.grthcond[i] * delzt[i]) / (
                    delzt[i + 1] + delzt[i]
                )
                kthuo = (grthcondo[i + 1] * delzt[i + 1] + grthcondo[i] * delzt[i]) / (
                    delzt[i + 1] + delzt[i]
                )

            else:
                # top soil node, veg & snow nodes if present
                kthl = (state.grthcond[i - 1] * delzt[i - 1] + state.grthcond[i] * delzt[i]) / (
                    delzt[i - 1] + delzt[i]
                )
                kthlo = (grthcondo[i - 1] * delzt[i - 1] + grthcondo[i] * delzt[i]) / (
                    delzt[i - 1] + delzt[i]
                )

                if i == ntemp:
                    kthu = 0.0
                    kthuo = 0.0
                else:
                    kthu = (state.grthcond[i + 1] * delzt[i + 1] + state.grthcond[i] * delzt[i]) / (
                        delzt[i + 1] + delzt[i]
                    )
                    kthuo = (grthcondo[i + 1] * delzt[i + 1] + grthcondo[i] * delzt[i]) / (
                        delzt[i + 1] + delzt[i]
                    )

                if icase == 0:
                    # no veg, no snow
                    q1lhsb = disurfg
                    q1rhs = state.isurfoldg + isurfg
                elif icase == 1 or icase == 4:
                    # snow/no veg, or veg buried by snow
                    if i == ntemp:
                        out = 0.0
                        in_ = 0.0
                        q1lhsb = disurfg
                        q1rhs = state.isurfoldg + isurfg
                    else:
                        q1lhsb = 0.0
                        q1rhs = 0.0
                elif icase == 2:
                    # vegetation, no snow
                    if i == ntemp:
                        out = 0.0
                        in_ = 0.0
                        q1lhsb = disurff
                        q2lhsb = -disurfgf
                        q1lhsa = disurffg
                        q2lhsa = -disurfg
                        q1rhs = state.isurfoldf + isurff
                        q2rhs = -(state.isurfoldg + isurfg)
                    else:
                        q1lhsb = disurfg
                        q1lhsc = disurfgf
                        q1rhs = state.isurfoldg + isurfg
                elif icase == 3:
                    # vegetation, snow
                    if i == ntemp:
                        out = 0.0
                        in_ = 0.0
                        q1lhsb = disurff
                        q2lhsb = -disurfgf
                        q1lhsa = disurffg
                        q2lhsa = -disurfg
                        q1rhs = state.isurfoldf + isurff
                        q2rhs = -(state.isurfoldg + isurfg)
                    elif ntemp - 1 >= i >= nnodes + 1:
                        out = 0.0
                        in_ = 0.0
                        if i == ntemp - 1:
                            q1lhsb = disurfg
                            q1lhsc = disurfgf
                            q1rhs = state.isurfoldg + isurfg

            if i == 1:
                dzl = 0.0
                dzu = 1.0 / (zt[i + 1] - zt[i])
                tl = state.stt[i]
                tlo = told[i]
                tu = 5e-1 * (state.stt[i] + state.stt[i + 1])
                tuo = 5e-1 * (told[i] + told[i + 1])
                rhsl = qb
                rhsu = kthu * (state.stt[i + 1] - state.stt[i]) * dzu + kthuo * (
                    told[i + 1] - told[i]
                ) * dzu
            elif 1 < i <= nnodes - 1:
                dzl = 1.0 / (zt[i] - zt[i - 1])
                dzu = 1.0 / (zt[i + 1] - zt[i])
                tl = 5e-1 * (state.stt[i] + state.stt[i - 1])
                tlo = 5e-1 * (told[i] + told[i - 1])
                tu = 5e-1 * (state.stt[i] + state.stt[i + 1])
                tuo = 5e-1 * (told[i] + told[i + 1])
                rhsl = kthl * (state.stt[i] - state.stt[i - 1]) * dzl + kthlo * (
                    told[i] - told[i - 1]
                ) * dzl
                rhsu = kthu * (state.stt[i + 1] - state.stt[i]) * dzu + kthuo * (
                    told[i + 1] - told[i]
                ) * dzu
            else:
                dzl = 1.0 / (zt[i] - zt[i - 1])
                dzu = 0.0
                tl = 5e-1 * (state.stt[i] + state.stt[i - 1])
                tlo = 5e-1 * (told[i] + told[i - 1])
                rhsl = kthl * (state.stt[i] - state.stt[i - 1]) * dzl + kthlo * (
                    told[i] - told[i - 1]
                ) * dzl + q2rhs
                if i == ntemp:
                    tu = 5e-1 * (state.stt[i] + dmet[4])
                    tuo = 5e-1 * (told[i] + airo)
                    rhsu = q1rhs
                else:
                    tu = 5e-1 * (state.stt[i] + state.stt[i + 1])
                    tuo = 5e-1 * (told[i] + told[i + 1])
                    dzu = 1.0 / (zt[i + 1] - zt[i])
                    rhsu = kthu * (state.stt[i + 1] - state.stt[i]) * dzu + kthuo * (
                        told[i + 1] - told[i]
                    ) * dzu + q1rhs

            # A(i) = coeff. for delstt(i-1); B(i) = coeff. for delstt(i);
            # C(i) = coeff. for delstt(i+1); D(i) = rhs. z positive upwards
            # from the bottom.
            if i == 1:
                B[i] = (
                    f1 * kthu / (zt[i + 1] - zt[i])
                    - f1 * (in_ - out)
                    + sph
                    + f1 * 5e-1 * (state.flowu[i] - state.flowl[i])
                    - rhsDi * lhs1
                    + rhsDv * lhs2
                    + rhsDl * lhs3
                )
                C[i] = -f1 * kthu / (zt[i + 1] - zt[i]) - f1 * 5e-1 * state.flowu[i]
            elif i == ntemp:
                A[i] = (
                    -f1 * kthl / (zt[i] - zt[i - 1])
                    - f1 * 5e-1 * state.flowl[i]
                    + f1 * (q2lhsa - q1lhsa)
                )
                B[i] = (
                    f1 * kthl / (zt[i] - zt[i - 1])
                    - f1 * (in_ - out)
                    + sph
                    + f1 * 5e-1 * (state.flowu[i] - state.flowl[i])
                    - rhsDi * lhs1
                    + rhsDv * lhs2
                    + rhsDl * lhs3
                    + f1 * (q1lhsb - q2lhsb)
                )
            else:
                A[i] = -f1 * kthl / (zt[i] - zt[i - 1]) - f1 * 5e-1 * state.flowl[i]
                B[i] = (
                    f1 * (kthl / (zt[i] - zt[i - 1]) + kthu / (zt[i + 1] - zt[i]))
                    + f1 * 5e-1 * (state.flowu[i] - state.flowl[i])
                    - f1 * (in_ - out)
                    - rhsDi * lhs1
                    + rhsDv * lhs2
                    + rhsDl * lhs3
                    - f1 * q1lhsb
                    + sph
                )
                C[i] = (
                    -f1 * kthu / (zt[i + 1] - zt[i])
                    - f1 * 5e-1 * state.flowu[i]
                    + f1 * (q1lhsc - q2lhsc)
                )

            D[i] = (
                -sph * (state.stt[i] - told[i])
                + rhsi
                - rhsv
                - rhsw
                + f1 * (rhsu - rhsl)
                - f1 * (state.fv1[i] + fv1o[i])
                + f1 * (in_ - out) * (state.stt[i] + told[i])
                - f1 * (state.flowu[i] * abs(tu - TREF) - state.flowl[i] * abs(tl - TREF))
                - f1 * (flowuo[i] * abs(tuo - TREF) - flowlo[i] * abs(tlo - TREF))
            )

            A[i] = _round15(A[i])
            B[i] = _round15(B[i])
            C[i] = _round15(C[i])
            D[i] = _round15(D[i])

            if abs(D[i]) > rhs_errort:
                rhs_errort = abs(D[i])
                sert = D[i]

        elif state.ntype[i] == 27:
            D[i] = 0.0
            if sn2 == 0:
                sn2 = i
            sn3 = i
    rhs_errort = _round10(rhs_errort)

    # -- solve the matrix eqn. for delta sst(i) (Thomas algorithm) ----------
    is2 = 0
    if2 = 0
    if sn2 == 0 and sn3 == 0:  # no air nodes
        is1 = 1
        if1 = ntemp
    elif sn2 == 1 and sn3 == ntemp:  # all air nodes
        is1 = 0
        if1 = 0
    elif sn2 > 1 and sn3 == ntemp:  # air on top
        is1 = 1
        if1 = sn2 - 1
    elif sn2 == 1 and sn3 < ntemp:  # air on bottom
        is1 = sn3 + 1
        if1 = ntemp
    elif sn2 > 1 and sn3 < ntemp:  # air layer (interior)
        is1 = 1
        if1 = sn2 - 1
        is2 = sn3 + 1
        if2 = ntemp
    else:
        is1 = 0
        if1 = 0

    if is1 > 0 and if1 > 0:
        bet = B[is1]
        if abs(bet) <= 1e-10:
            tsign = 1.0
            if bet < 0.0:
                tsign = -1.0
            delstt[is1] = tsign * abs(D[is1])
        else:
            delstt[is1] = D[is1] / bet

        for j in range(is1 + 1, if1 + 1):
            if state.ntype[j] != 27:
                gam[j] = C[j - 1] / bet
                bet = B[j] - A[j] * gam[j]
                if abs(bet) <= 1e-10 or abs(A[j]) <= EPS:
                    tsign = 1.0
                    if bet < 0.0:
                        tsign = -1.0
                    delstt[j] = tsign * abs(D[j])
                else:
                    delstt[j] = (D[j] - A[j] * delstt[j - 1]) / bet
            else:
                delstt[j] = 0.0

        for j in range(if1 - 1, is1 - 1, -1):
            delstt[j] = delstt[j] - gam[j + 1] * delstt[j + 1]

    if is2 > 0 and if2 > 0:
        bet = B[is2]
        if abs(bet) <= 1e-10:
            tsign = 1.0
            if bet < 0.0:
                tsign = -1.0
            delstt[is2] = tsign * abs(D[is2])
        else:
            delstt[is2] = D[is2] / bet

        for j in range(is2 + 1, if2 + 1):
            if state.ntype[j] != 27:
                gam[j] = C[j - 1] / bet
                bet = B[j] - A[j] * gam[j]
                if abs(bet) <= 1e-10 or abs(A[j]) <= EPS:
                    tsign = 1.0
                    if bet < 0.0:
                        tsign = -1.0
                    delstt[j] = tsign * abs(D[j])
                else:
                    delstt[j] = (D[j] - A[j] * delstt[j - 1]) / bet
            else:
                delstt[j] = 0.0

        for j in range(if2 - 1, is2 - 1, -1):
            delstt[j] = delstt[j] - gam[j + 1] * delstt[j + 1]

    # -- determine new temps, ice and vapor contents; adjust water ----------
    # contents accordingly
    errort = -99999.9
    dtop = 2.5
    if iter_ == 0:
        dtop = 1e-1 * dtop
    mxd1 = 4e1 * (state.nz[nnodes] - state.nz[1])

    for i in range(ntemp, 0, -1):
        delstt[i] = _round15(delstt[i])

        tcheck[i] = state.stt[i]
        wcheck[i] = state.wvc[i]
        icheck[i] = state.ice[i]
        scheck[i] = state.soil_moist[i]
        pcheck[i] = state.phead[i]

        if abs(delstt[i]) > errort and i <= nnodes:
            errort = abs(delstt[i])

        if iter_ > 0 and (abs(state.ice[i] - iceo[i]) > EPS and i <= nnodes):
            dtopi = 1e-1 * dtop
        else:
            dtopi = dtop

        if abs(delstt[i]) > dtopi:
            tsign = 0.0
            if delstt[i] != 0.0:
                tsign = abs(delstt[i]) / delstt[i]
            delstt[i] = tsign * dtopi
            delstt[i] = _round15(delstt[i])

        state.stt[i] = state.stt[i] + delstt[i]

        if i == nnodes and state.mstflag == 1:
            # use measured soil temp -- Fortran uses AINT (truncate), not
            # ANINT, here.
            if _aint(abs(dmet[13] - state.mflag) * 1e5) * 1e-5 > EPS:
                state.stt[nnodes] = dmet[13]
                fflag = 1

        if state.ntype[i] == 27:
            if i > nnodes or (state.ntype[nnodes] == 27 or state.ntype[1] == 27):
                state.stt[i] = dmet[4]
            else:
                if 1 < i < ntemp:
                    state.stt[i] = 5e-1 * (state.stt[sn3 + 1] + state.stt[sn2 - 1])
                else:
                    state.stt[i] = dmet[4]

        maxdel = 1e1 * state.timstep

        if iw != 1:
            if abs(dmet[4] - airo) > EPS:
                maxdel = min(mxd1, 1e1 * abs(dmet[4] - airo))

        d1 = 0.0
        if i < nnodes:
            d1 = min(1.0, state.elev - state.nz[i])
        elif i >= nnodes and ntemp != nnodes:
            d1 = max(1.0, state.hfol, state.hm)

        deep = 1.0
        if d1 > EPS:
            deep = 4.388 * (d1**4e0) - 11.498 * (d1**3e0) + 10.947 * d1 * d1 - 4.7845 * d1 + 9.546e-1
            deep = min(1.0, max(1e-1, deep))

        maxdel = max(maxdel * deep, mxd1 * deep)

        # `y` computed via 4 branches then unconditionally overwritten --
        # see module docstring, point 2. Computed literally, then discarded.
        if state.elev - state.nz[i] >= 1.0:
            y = 1.0 + min(1.0, 1e-2 * ((state.elev - state.nz[i]) - 1.0))
        elif i >= nnodes and state.hm >= 1.0:
            y = 1.0 + min(1.0, 1e-2 * (state.hm - 1.0))
        elif i >= nnodes and state.hfol >= 1.0:
            y = 1.0 + min(1.0, 1e-2 * (state.hfol - 1.0))
        else:
            y = 1.0 - deep

        y = 6e1

        # -- check for continuity in time ------------------------------------
        if abs(state.stt[i] - state.too[i]) > maxdel:
            if state.stt[i] > state.too[i]:
                state.stt[i] = state.too[i] + maxdel
            else:
                state.stt[i] = state.too[i] - maxdel

        # -- check for vertical continuity ------------------------------------
        if state.ntype[i] != 27:
            if i == ntemp or (i != ntemp and cflag == 0):
                cflag = 1
                if abs(dmet[4] - state.stt[i]) > 4e1:
                    if state.stt[i] > dmet[4]:
                        j = 1
                        while state.stt[i] > dmet[4] + 4e1:
                            state.stt[i] = tcheck[i] - j * maxdel
                            j = j + 1
                    else:
                        j = 1
                        while state.stt[i] < dmet[4] - 4e1:
                            state.stt[i] = tcheck[i] + j * maxdel
                            j = j + 1

                if i == nnodes and fflag == 1:
                    state.stt[i] = dmet[13]
                state.tmelt[iw] = state.stt[i]
            else:
                k = sn
                if sn2 > 0:
                    k = sn2 - 1

                if abs(state.stt[k] - state.stt[i]) > y and k > 1:
                    if state.stt[i] > state.stt[k]:
                        state.stt[i] = state.stt[k] + y
                    else:
                        state.stt[i] = state.stt[k] - y

                if i == nnodes and fflag == 1:
                    state.stt[i] = dmet[13]
                if icase == 3 and i == sn:
                    state.tmelt[iw] = state.stt[i]

        state.stt[i] = _round10(state.stt[i])

        if state.node_type[i] == "VG":
            ftempt = state.stt[i]
        else:
            ftempt = dmet[4]

        if state.stt[i] <= TREF:
            if told[i] <= TREF:
                rhow = dense(273.16, 0.0, d1i)
            else:
                rhow = dense(told[i], 0.0, d1i)
            rhoi = dense(state.stt[i], 0.0, d2i)
        else:
            rhow = dense(state.stt[i], 0.0, d1i)
            if told[i] <= TREF:
                rhoi = dense(told[i], 0.0, d2i)
            else:
                rhoi = dense(TREF, 0.0, d2i)
        dfwi = rhow / rhoi
        dfiw = rhoi / rhow

        # -- update ice content -----------------------------------------------
        if i <= nnodes and state.ntype[i] <= 25:
            smm = state.soil_moist[i] * rhow
            icm = state.ice[i] * rhoi
            stc = _round15(state.stt[i] - TREF)

            if ii == 1 and iter_ == 0:
                state.bftm[i] = smm + icm
            state.bftm[i] = _round15(state.bftm[i])

            if stc < EPS or state.ice[i] > 0.0:
                minice = 0.0
                maxice = max(0.0, min(icm + smm, state.bftm[i]))

                if abs(stc) > EPS:
                    y = abs(lhs5[i]) * max(abs(delstt[i]), abs(stc))
                else:
                    y = abs(lhs5[i] * delstt[i])

                if abs(lhs5[i]) <= EPS and state.ice[i] > EPS:
                    didT = -(rhow / rhoi) * LHFUS / (GRAV * TREF * abs(delzt[i]))
                    if abs(stc) > EPS:
                        y = abs(didT) * max(abs(delstt[i]), abs(stc))
                    else:
                        y = abs(didT * delstt[i])

                if state.node_type[i] == "WA":
                    if abs(stc) > EPS:
                        y = abs(lhs5[i]) * min(abs(delstt[i]), abs(stc))
                    else:
                        y = abs(lhs5[i] * delstt[i])

                if y < EPS:
                    y = 0.0
                y = _round15(y)

                if stc >= EPS:
                    # melting
                    icm = (icheck[i] - y) * rhoi
                    if icm < minice:
                        y = icheck[i] * rhoi - minice
                        icm = minice
                    smm = scheck[i] * rhow + y * rhoi
                elif stc < EPS and state.soil_moist[i] > EPS:
                    # freezing
                    icm = (icheck[i] + y) * rhoi
                    if icm > maxice:
                        y = icm - maxice
                        icm = maxice
                    smm = max(0.0, scheck[i] * rhow - y * rhoi)
                elif stc < EPS and state.node_type[i] == "WA":
                    icm = (icheck[i] + y) * rhoi
                    if icm > maxice:
                        y = icm - maxice
                        icm = maxice
                    smm = max(0.0, scheck[i] * rhow - y * rhoi)

                icm = _round15(icm)
                state.ice[i] = max(0.0, icm / rhoi)
                state.ice[i] = _round10(state.ice[i])
                smm = _round15(smm)
                state.soil_moist[i] = max(0.0, smm / rhow)
                state.soil_moist[i] = _round10(state.soil_moist[i])

                if (icheck[i] <= EPS and state.ice[i] > EPS) or (
                    icheck[i] > EPS and state.ice[i] <= EPS
                ):
                    if iw != state.istart:
                        state.stt[i] = TREF

                if abs((smm + icm) - state.bftm[i]) > 1e-5 and state.node_type[i] != "WA":
                    smm = state.bftm[i] - icm
                    state.soil_moist[i] = _round10(smm / rhow)

                    if state.ice[i] > EPS:
                        if state.soil_moist[i] <= 1e-10:
                            state.soil_moist[i] = 0.0
                    else:
                        if state.soil_moist[i] < state.nsoilp[i, 15]:
                            state.soil_moist[i] = state.nsoilp[i, 15]
                        elif state.soil_moist[i] > state.nsoilp[i, 24]:
                            state.soil_moist[i] = state.nsoilp[i, 24]
                    if state.node_type[i] != "WA" and state.node_type[i] != "AI":
                        state.phead[i] = head(state, i, state.soil_moist[i])

                elif state.ice[i] <= EPS:
                    if state.soil_moist[i] < state.nsoilp[i, 15]:
                        state.soil_moist[i] = state.nsoilp[i, 15]
                        if state.node_type[i] != "WA" and state.node_type[i] != "AI":
                            state.phead[i] = head(state, i, state.soil_moist[i])
                    elif state.soil_moist[i] > state.nsoilp[i, 24]:
                        state.soil_moist[i] = state.nsoilp[i, 24]
                        if state.node_type[i] != "WA" and state.node_type[i] != "AI":
                            state.phead[i] = head(state, i, state.soil_moist[i])
                    elif state.node_type[i] == "WA":
                        state.soil_moist[i] = 1.0

                if state.ice[i] > EPS:
                    iflag = 1

            # -- update water vapor content -----------------------------------
            if (state.ntype[i] != 26 or state.ntype[i] != 27) and state.node_type[i] != "VG":
                # `.or.` is a tautology -- see module docstring, point 3.
                # Reduces to `node_type[i] != 'VG'`. Preserved literally.
                rh = soilhumid(state, i, state.phead[i], state.soil_moist[i], state.stt[i])
                pres = dmet[11] + 1e-2 * ((state.elev - state.nz[i]) + state.phead[i]) * rhow * GRAV
            else:
                rh = dmet[5] * 1e-2
                if state.ntype[i] == 26:
                    rh = 1.0
                pres = dmet[11]

            vpress = vap_press(state, i, rh, pres)
            mixr = 0.0
            if pres * 1e2 - vpress > EPS:
                mixr = 0.622 * vpress / (pres * 1e2 - vpress)

            rhov[i] = 0.622 * vpress / (RV * state.stt[i])
            rhov[i] = _round20(rhov[i])
            rhoda[i] = max(0.95, min(2.8, (pres * 1e2 - vpress) / (RD * state.stt[i])))
            rhoda[i] = _round20(rhoda[i])

            t1 = min(1.0, max(mixr * rhoda[i] / (rhow + mixr * rhoda[i]), 0.0))
            t1 = _round20(t1)
            state.wvc[i] = t1 * (state.nsoilp[i, 2] - (state.soil_moist[i] + state.ice[i]))
            if state.ntype[i] == 27:
                state.wvc[i] = t1

            state.wvc[i] = max(
                0.0,
                min(state.wvc[i], state.nsoilp[i, 2] - (state.soil_moist[i] + state.ice[i])),
            )
            state.wvc[i] = _round15(state.wvc[i])

    return SoilTmpResult(iflag=iflag, rhs_errort=rhs_errort, errort=errort, sert=sert)
