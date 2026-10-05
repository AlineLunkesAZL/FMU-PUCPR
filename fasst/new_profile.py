"""Time integration of the soil/snow/vegetation profile, transcribed from
new_profile.F90.

Fortran source: new_profile.F90 (66 693 bytes, 1 860 lines) -- the largest
file in the project after module_canopy and module_snow, and the only one
consisting of a single subroutine (no attached helpers). See
docs_transcricao/doc-transcricao-new_profile.docx for the full transcription
rationale ("the transcription doc" below); this docstring restates its 15
numbered points plus the extra findings from the pre-transcription review.

What :func:`new_profile` does
------------------------------
Called once per meteorological timestep (fasst_main.F90:251, inside the main
loop over met records -- fan-in = 1). It (1) picks the internal sub-step
count ``step`` -- cut to 1/4 of the standard value when snow or type-2
precipitation is present (``lflag``); (2) sets the timestep's node topology
(``node_type``/``ntemp``) from ``icase in {0,1,2,3,4}`` -- bare soil, snow
over bare soil, vegetation without snow, vegetation with partial snow,
vegetation buried by snow; (3) solves the coupled heat/moisture balance in
an inner ``while`` loop (``soil_moisture`` + ``soil_tmp`` + a mid-loop
``surfenergy`` call, iterated to convergence or ``maxiter = 10``); (4)
repeats that loop ``step`` times, once per internal sub-step, in
``for ii in range(1, step + 1)``; (5) afterwards does a freeze/thaw check,
adjusts soil-node thickness for ice expansion/contraction, and computes the
output energy terms via :func:`~fasst.sflux.sflux` (once for bare
soil/snow, once more for foliage if present, once more for melt energy if
applicable).

Fortran line ranges (transcription-doc Table 1)::

    1-236      declarations, SAVE, zeroing
    238-284    internal (sub-)timestep selection
    286-432    first_time==0 initialization (infer_test branch restores a
               partial checkpoint)
    434-498    met interpolation for the timestep (delvar/delvar1/dmet)
    499-859    profile topology (node_type, ntemp, zt) per icase 0-4
    861-877    node thicknesses (delzt)
    879-999    initial surface flux (1st surfenergy) + max infiltration
    960-1042   surface water balance (hpond, qtop, qtopv, overland)
    1070-1614  MAIN LOOP over ii=1..step: sub_divide, albedo_emis, inner
               convergence loop (soil_moisture / soil_tmp / 2nd surfenergy
               / flow_param), sub-step closeout (3rd surfenergy at ii==step)
    1619-1705  final freezing check; delzs/nz adjustment for ice volume
    1707-1858  energy closure and outputs (sflux x1-3; melt; tot_moist;
               sumrunoff; infl_cum; commented debug diagnostics)

Calls: :func:`~fasst.sp_humid.sp_humid` (x5 call sites),
:func:`~fasst.surfenergy.surfenergy` (x3), :func:`~fasst.sflux.sflux` (x3),
plus six subroutines that are NOT YET transcribed to Python (see next
section). Uses the functions :func:`~fasst.functions.dense`,
:func:`~fasst.functions.soilhumid`, :func:`~fasst.functions.maxinfiltrate`,
:func:`~fasst.functions.head`, :func:`~fasst.functions.vap_press` (all
already transcribed; ``maxinfiltrate`` is called ONLY here in the whole
project).

Point 1 (header defect): the Fortran header comment lists ``smooth`` as a
called subroutine and omits ``sp_humid``. ``smooth`` is never called (its
one call site, line 1520, is commented out, and there is no ``smooth.F90``
among the 30 source files); ``sp_humid`` IS called 3 times (426, 607,
1267). Not transcribed / not a defect here -- the header comment is simply
not authoritative, same class as surfenergy point 13 and module_snow point
2.

Injected untranscribed dependencies (point 15 -- BLOCKS numeric validation)
----------------------------------------------------------------------------
Of the 9 subroutines ``new_profile`` calls, only ``surfenergy``, ``sflux``
and ``sp_humid`` have a Python module today; they are imported directly
below. The other six -- ``soil_tmp``, ``soil_moisture``, ``th_param``,
``flow_param``, ``sub_divide``, ``albedo_emis`` -- have no ``fasst/*.py``
counterpart yet. Unlike ``surfenergy``'s ``lowveg_met`` (reached only on
the vegetated path), all six sit on :func:`new_profile`'s UNCONDITIONAL
main path: this function cannot run to completion without a real
implementation injected for each. They are accepted as required
keyword-only callables so the wiring can be built and unit-tested with
stubs today, per the contract below (verified against each ``.F90``
header's actual ``intent`` declarations, not guessed):

* ``albedo_emis_fn(state, oldsd, sd, su) -> None`` -- Fortran
  ``albedo_emis(oldsd,sd,su)``, called at new_profile.F90:1075. ALL THREE
  arguments are ``intent(in)`` in the Fortran; there is no ``intent(out)``
  argument at all. Everything it computes (``albedo_fasst``, ``sgralbedo``,
  ``emis``, ``sgremis``, plus its own internal ``save:: albedoo``) is
  written straight into module-global state. Must mutate ``state`` in
  place; returns nothing.
* ``th_param_fn(state, n, pdens, pdensnew, rhov, rhoda) -> float`` --
  Fortran ``th_param(n,pdens,pdensnew,rhov,rhoda,rhotot)``, called at
  new_profile.F90:566, 1343, 1555 (with ``n = ntemp`` the first time, ``n =
  ntot`` afterwards). Only ``rhotot`` is ``intent(out)``; ``rhov``/``rhoda``
  are ``intent(in)`` (read-only, despite superficially looking like they
  might get updated). Writes ``state.grthcond``/``state.grspheat`` (indices
  ``1..n``) in place; returns ``rhotot``.
* ``flow_param_fn(state, isn, qtop, qtopv, simeltsm, zt, delz) ->
  tuple[float, float, float, float]`` -- Fortran
  ``flow_param(isn,qtop,qtopv,simeltsm,zt,delz,qbot,klhtop,kvhtop,kvttop)``,
  called at new_profile.F90:937 (``isn = nnodes``), 1041, 1240, 1562 (``isn
  = 1``). Returns ``(qbot, klhtop, kvhtop, kvttop)`` -- the four declared
  ``intent(out)`` scalars, in that order. Its header comment documents that
  it also computes ``dsmdh`` (a :class:`~fasst.state.FasstState` field) and
  reads/writes several more module globals (``khu``, ``khl``, ...) that
  never appear in its formal argument list -- treat it as free to mutate
  ``state`` beyond the four returned scalars.
* ``sub_divide_fn(state, ii, delvar, delvar1, dmet) -> None`` -- Fortran
  ``sub_divide(ii,delvar,delvar1,dmet)``, called once per sub-step at
  new_profile.F90:1073. ``dmet`` is the sole ``intent(out)`` -- a 13-element
  (phantom-indexed, slots 1..13) array that must be fully overwritten in
  place; nothing is returned.
* ``soil_moisture_fn(state, sm_old, wvco, iceo, sourceo, sinkro, vino,
  told, thvc, dthvdh, runoff) -> tuple[float, float, float, float]`` --
  Fortran ``soil_moisture(sm_old,wvco,iceo,sourceo,sinkro,vino,told,thvc,
  dthvdh,runoff,rhs_errorm,errorm,smerror,sert1)``, called at
  new_profile.F90:1196, guarded by ``errorm > allowerrorm``. Returns
  ``(rhs_errorm, errorm, smerror, sert1)`` in that order; ``runoff`` is
  ``intent(inout)`` and must be mutated in place. Per its own header
  comment ("Newton-Raphson technique"), it is expected to update
  ``state.soil_moist``, ``state.phead``, ``state.sink``, ``state.sinkr``,
  ``state.source``, ``state.ice``, ``state.wvc`` and ``state.dsmdh`` in
  place (all module globals it imports but that are not in its formal
  argument list). ``errorm`` re-enters the ``while errorm > allowerrorm``
  loop condition -- this is the convergence driver for the moisture half of
  the coupled solve.
* ``soil_tmp_fn(state, isurfg, isurff, disurfg, disurff, disurffg,
  disurfgf, sigflo, rhotot, airo, dmet, grthcondo, grspheato, zt, delzt,
  told, iceo, wvco, flowuo, flowlo, fv1o, sinko, sourceo, ph_old, sm_old,
  thvc, dthvdt, rhov, rhoda, ii, iter_) -> tuple[int, float, float,
  float]`` -- Fortran ``soil_tmp(iflag,isurfg,isurff,disurfg,disurff,
  disurffg,disurfgf,sigflo,rhotot,airo,dmet,grthcondo,grspheato,zt,delzt,
  told,iceo,wvco,flowuo,flowlo,fv1o,sinko,sourceo,ph_old,sm_old,thvc,
  dthvdt,rhov,rhoda,rhs_errort,errort,sert,ii,iter)``, called at
  new_profile.F90:1209. Returns ``(iflag, rhs_errort, errort, sert)`` --
  the caller's local ``iceflag`` receives the ``iflag`` output positionally
  (a genuine, if uninformative, ``intent(out)``: it is overwritten on every
  call and never read before the NEXT call inside the ii/iter loop, only
  read again after the whole ``for ii`` loop finishes, where it has already
  been reset -- see point 4). ``rhov``/``rhoda`` are ``intent(inout)`` and
  must be mutated in place (unlike ``th_param_fn``'s, which are read-only).
  Per its header comment ("Crank-Nicholson"), it updates ``state.stt``,
  ``state.wvc``, ``state.ice``, ``state.phead``, ``state.soil_moist``,
  ``state.source``, ``state.sink`` in place. ``errort`` re-enters the
  ``while ... errort > allowerrort`` loop condition -- the temperature half
  of the coupled solve.

``lowveg_met`` (see :mod:`fasst.surfenergy`) is a SEVENTH injected
dependency, reachable transitively through the three :func:`surfenergy`
calls whenever ``icase in (2, 3)``; forwarded here as a passthrough
keyword-only argument.

Not pure -- reads AND writes ``state``
---------------------------------------
Reads essentially the whole profile-related subset of :class:`FasstState`
(``iw, ntot, nnodes, ntemp, node_type, maxlines, met, stt, phead, hm, wvc,
delzs, dmet1, albedo_fasst, deltat_fasst, deltati, elev, emis, ftemp,
hfol_tot, hpond, icase, icaseo, iend, istart, infer_test, isurfoldg, mflag,
newsd, pi(via constants), Rd/Rv(constants), refreeze, refreezei, rough,
sgralbedo, sigfl, sloper, step, stepi, stll, stls, storll, storls, timstep,
toptemp, veg_flagl, veg_flagh, flowu, flowl, vin, isurfoldf, fv1, nsoilp,
ice, soil_moist, sink, source, sinkr, too, woo, ioo, smoo, phoo, nz, sdens,
grthcond, grspheat, ntype, tot_moist, lhes, bftm, delzsi, nzi, lheatf,
sdown, sup, tmelt, evap_heat, ponding, irdown, irup, pheat1, sheat, lheat,
cheat, cheat1``) plus the constants ``EPS, GRAV, HT_MIN, HT_MINM, KVEG, PI,
RD, RV, SPFLAG, TREF``. **Writes** essentially the same set back (every
per-node array above except the read-only inputs ``nsoilp``/``bftm``) plus
``state.icase, state.icaseo, state.ntemp, state.ftemp, state.toptemp,
state.hpond, state.sgralbedo, state.albedo_fasst, state.emis,
state.deltat_fasst, state.step, state.delzs, state.nz``. A caller reading
only the returned :class:`NewProfileResult` misses nearly everything --
this function is the textbook case for point 2 below.

Point 2 -- ~21 SAVE variables: the routine itself is a stateful object
------------------------------------------------------------------------
The ``save::`` block (Fortran lines 45-46) lists ``iwm, iwt, iwmr, iwtr,
iwsm, ntempo, infl_cum, maxrt, maxt, maxrm, maxm, tcum, klhtop, s2i,
thetai, klhtopi, airo, x1, x, totice, stcalc`` -- all of them persist from
one call of ``new_profile`` (one met timestep) to the next. This is the
central point of the transcription: these become fields of
:class:`NewProfileState`, instantiated once per model run/point and passed
in on every call (the ``save`` parameter below), exactly the same pattern
:class:`FasstState` uses for the module-level ``save::`` blocks, but scoped
to a single subroutine instead of a whole module. A naive
``def new_profile(...)`` with ordinary Python locals would lose all memory
between timesteps. ``save.stcalc`` in particular is the converged
temperature from the previous timestep, read at line 580 to seed ``stt``
for snow/mixed nodes before anything else is computed this timestep.

Of these 21, nine (``maxrt, maxt, maxrm, maxm, iwm, iwt, iwmr, iwtr, iwsm``)
are initialized at ``first_time == 0`` and then NEVER reassigned anywhere
in the active (uncommented) code -- the only code that would update them
(new_profile.F90:1460-1479) is entirely commented out. They are kept as
:class:`NewProfileState` fields for structural fidelity to the ``save::``
list, but no line of :func:`new_profile` below writes them after
initialization; do not expect them to change.

Point 3 -- convergence: 7 exit codes computed, none reported
----------------------------------------------------------------
The inner ``while errorm > allowerrorm or errort > allowerrort`` loop
(Fortran line 1091) can exit via 7 distinct paths (``iflag = 1..7``, lines
1350-1404), each a different stagnation/oscillation heuristic -- including
``iflag == 6`` (line 1371), which UNDOES the current iteration and restores
state from the ``x(...)`` snapshot when the sign of total ice changes
between ``save.totice`` and its previous value while the error is
diverging. After the loop, ``iflag`` is never read for diagnostics. The
ENTIRE non-convergence logging block (Fortran lines 1481-1513, which would
``write(10, ...)`` "No convergence in new_profile; day..., hour...,
minute...", set ``error_code = 1`` and ``code5 = 1``) is commented out in
the source. So today: if the profile fails to converge within
``maxiter = 10`` iterations, the Fortran silently accepts the last computed
state and moves on -- no exception, no warning, no counter.

Per the user's explicit instruction for this transcription pass ("mantenha
a fidelidade do código original, deixando futuras melhorias para outra
etapa"), this is decided explicitly here, not by omission: **the silence
is replicated.** No ``logging.warning``, no counter, no exception is added
for non-convergence at ``iter >= maxiter``. Reintroducing that log (or
adding a Python-native diagnostic) is left as an open, explicitly-deferred
decision for a later phase (transcription-doc PENDENCIAS).

Point 4 -- ``iflag``: two unrelated meanings, one name
------------------------------------------------------------
Directly tied to point 3: ``iflag`` is reused for the final freezing check
(Fortran line 1620 onward -- there it is just 0/1, "is there relevant ice
at any node?"), with no relation to the 1-7 convergence diagnostic values
the same variable carried a few lines earlier (and, inside the loop, is a
DIFFERENT variable from ``soil_tmp``'s returned ``iflag``/``iceflag`` --
see the injected-dependency section above). Transcribed as two independent
Python locals (``iflag`` reused only within the convergence loop scope
proper, per its original Fortran role, and a second use for the freezing
check) to avoid conflating them, even though the source itself conflates
the name. Same class as surfenergy's ``beta``, doc point 12 there.

Point 5 -- inconsistent rounding of the same physical quantity
----------------------------------------------------------------
``ice(i) = (refreeze + refreezei)*f1`` is computed identically in three
``icase`` branches (1, 3, 4) but rounded at different precisions:
``anint(ice(i)*1d20)*1d-20`` on the ``icase == 1`` branch (line 689) versus
``anint(ice(i)*1d10)*1d-10`` on ``icase == 3`` (line 781) and ``icase ==
4`` (line 850). Replicated literally, per branch -- NOT unified.

Point 6 -- ``tot_moist``: inverted rounding factors structurally zero it
---------------------------------------------------------------------------
Line 1834: ``tot_moist(iw) = anint(summoist*1d-10)*1d10``. Every other
occurrence of this idiom in the file (and the project) is
``anint(x*1dN)*1d-N`` -- multiply by a LARGE factor before rounding. Here
the exponents are swapped: ``summoist`` (order 0.01-1 m) times ``1d-10`` is
~1e-12; ``anint`` of that is 0; ``0 * 1d10 == 0``. So ``tot_moist[iw]`` is
structurally zero. Project-wide search: ``tot_moist`` is read back only at
new_profile.F90:1065 (``sumo = tot_moist(iw-1)``), and ``sumo`` feeds only
the fully-commented ``write(88,...)``/``write(99,...)`` debug mass-balance
block at the end of the file. No live output is affected today.
Replicated literally in phase 1 (including the resulting zero) -- not
"corrected". Line 1839's ``extra = anint(extra*1d-10)*1d10`` uses the exact
same inverted idiom and is replicated the same way, for the same reason
(also feeds only the same dead debug block).

Point 7 -- ``full``: counted, never read
---------------------------------------------
Line 1184, inside ``do i=1,ntot``: ``if(nsoilp(i,24)-soil_moist(i) < eps)
full = full + 1``. Its only would-be consumer (Fortran lines 1191-1207,
``if(full < nnodes.or. ...)``) is entirely commented out; the live branch
(line 1195, ``if(errorm > allowerrorm)``) never inspects it. Dead
computation, same class as module_snow's ``overc`` and sflux's
``dchw``/``dchw1``. ``full`` is computed here (as a plain Python local, for
literal fidelity) and never used afterward, matching the source.

Point 8 -- pond-albedo block duplicated, unnamed 1.33 refraction index
--------------------------------------------------------------------------
Fortran lines 885-898 and 1022-1036 are the identical 14-line block
verbatim: when ``hpond > eps``, compute the Fresnel albedo of a water
surface from the solar zenith angle and water's refractive index, written
as the bare literal ``1.33d0`` (both sites), unnamed. Per the
transcription doc (fase futura), this is left duplicated here -- both call
sites reproduced separately with the bare ``1.33`` -- rather than factored
into a shared helper; that refactor is explicitly out of scope for this
pass.

Point 9 -- ``sp_humid``'s ``rh2`` at the 3 call sites here is already a
fraction, not a percentage (contribution to the open ``sp_humid.py`` audit
of "5 callers x rh2 scale"): all three sites (426, 607, 1267) pass
``soilhumid(...)`` (already 0-1) for soil nodes, ``1.0`` for ``'WA'`` nodes,
or ``1d-2*dmet(5)`` for ``'AI'`` nodes -- the percent-to-fraction
conversion is explicit in that last case. These three sites can be marked
audited-and-correct for scale.

Point 10 -- ``d0i``/``d1i``/``d2i``/``d3i``: bare-phase selector -> ``Phase``
--------------------------------------------------------------------------------
Zeroed in the initial zero-out block (lines 91-94); ``d1i=1, d2i=2, d3i=3``
are assigned once (lines 234-236) and never change again; ``d0i`` is never
reassigned at all, staying ``0`` for the whole call, used correctly as
"water vapor" in the ``sp_humid`` calls. Same ``Phase`` convention already
catalogued in :mod:`fasst.functions` (``Phase.WATER_VAPOR = 0``,
``Phase.WATER = 1``, ``Phase.ICE = 2``, ``Phase.DRY_AIR = 3``) -- used here
as :class:`~fasst.functions.Phase` members, not bare ints, for both their
role as a ``phase``/``ic`` argument (:func:`dense`, :func:`sp_humid`) AND
their reuse as the ``vflag`` argument of :func:`~fasst.sflux.sflux`
(``d1i`` = foliage call, ``d0i`` = ground/melt calls) -- the Fortran
literally reuses the same two integer locals for both roles, so
``Phase.WATER``/``Phase.WATER_VAPOR`` (``IntEnum``, so plain-``int``
compatible) are passed as ``vflag`` too, preserving that reuse rather than
introducing a second pair of "flag" locals the source does not have.

Point 11 -- 3rd ``surfenergy`` call writes persistent ``lhes(iw)``; the
first two discard into scratch ``lh1``
--------------------------------------------------------------------------
The 3 ``surfenergy`` calls (lines 902, 1233, 1543) are positionally
near-identical, but the last argument differs: the first two pass the
local ``lh1`` (never read afterward); the third (line 1547, only when
``ii == step`` -- the sub-step closeout) passes ``lhes(iw)``, the
allocatable array indexed by timestep (fasst_global.F90:183). Only the
LAST call's ``lh`` output is retained for use outside this routine.
Preserved here by discarding ``result.lh`` on the first two calls and
writing ``state.lhes[state.iw] = result.lh`` only on the third.

Point 12 -- commented I/O and diagnostics: debug mass-balance vestiges
--------------------------------------------------------------------------
``write(*,*)`` (lines 1745-1746, energy-closure ``fcheck``/
``sheatf-lheatf(iw)``/etc.), ``write(88,...)`` and ``write(99,101/102,...)``
(lines 1841-1858, water balance, hardcoded units 88/99, named formats
101/102) are all commented out, while the calculations feeding them
(``summoist``, ``sumin``, ``extra``) remain active. Same class as
surfenergy's hardcoded unit 10 and this file's own point-3 convergence
log. None of these are transcribed as output statements; the feeding
calculations are computed (points 6/12) but never printed/logged.

Point 13 -- ``maxinfiltrate``: sole caller in the project
---------------------------------------------------------------
Line 918: ``s2i = maxinfiltrate(nnodes,thetai)``. ``fasst/functions.py``
already documents (its own "only new_profile" note) that this function has
a structural quadrature bias and an incidental single-precision-looking
(but actually not lossy) term -- neither was "fixed" there, and neither is
touched here. Both affect ``in_rate``/``infl_max`` here, i.e. the maximum
soil water infiltration rate this routine computes. If ``maxinfiltrate``
is ever corrected, the water balance computed by ``new_profile`` changes
with it.

Point 14 -- ``x``/``x1``: SAVE, but writes always precede any same-call
read -- unverified inertness hypothesis, NOT assumed
--------------------------------------------------------------------------
``x``, ``x1`` (see point 2) are only READ when ``ii > 1`` (line 1101) and
are unconditionally rewritten at the end of every inner ``while`` iteration
(line 1166 onward), before any subsequent read. This suggests persisting
``x``/``x1`` across CALLS of ``new_profile`` (i.e. across timesteps, via
SAVE) has no observable effect: the first sub-step (``ii == 1``) of a
timestep never reads ``x``, and every later sub-step will already have
been written within the SAME call. If that reasoning is right, ``x``/``x1``
could be ordinary Python locals without loss of fidelity -- but unlike
point 2's ``stcalc`` (real cross-timestep memory, already traced), this
has NOT been tested. **Not assumed here**: ``x``/``x1`` remain
:class:`NewProfileState` fields (true SAVE semantics, matching the Fortran
declaration) pending a dedicated test forcing ``step == 1`` across
successive calls with and without the persistence, per the transcription
doc's VALIDACAO section.

Point 5 (sizing) / extra finding -- ``x``/``stcalc`` array bound
--------------------------------------------------------------------
Fortran declares ``x(14,maxn)`` and ``stcalc(maxn)`` (``maxn = 100``), but
both are written in ``do i=1,ntot`` loops (lines 1168-1185, 1373-1388) and
``stcalc`` in ``do i=1,ntemp`` (1622). ``ntot`` can exceed ``maxn`` at
runtime (``ntot = nnodes + extran``, and ``nnodes`` can itself reach
``maxn``): the Fortran would then silently write past the declared bound.
Sized here as ``MAXN + EXTRAN`` (matching how ``fasst.state`` sizes
``stt``/``ice``/``phead``) instead of reproducing a latent Fortran
out-of-bounds write as a Python ``IndexError``. Every other ``ntot``-sized
local array below (``zt, delzt, told, wvco, iceo, ph_old, sm_old, ...``)
uses the same ``MAXN + EXTRAN`` sizing for the same reason -- the Fortran
declares them ``(ntot)``, which is not a compile-time constant there
either. Row 14 of ``x``/``y`` (used for ``flowl``) is never explicitly
zeroed by the Fortran's own init loops (which only run ``do i=1,13``); this
makes no observable difference here since the phantom-slot array
constructor below already zero-fills everything.

Extra finding -- ``delvar`` re-rounded twice, ``delvar1`` never on its own
------------------------------------------------------------------------------
Not one of the transcription doc's 15 points. Inside ``do ll=1,12`` (lines
463-478, filling ``delvar1``), the closing line is
``delvar(ll) = anint(delvar(ll)*1d20)*1d-20`` -- it re-rounds ``delvar``
(already rounded once in the PRIOR loop, lines 447-461), not ``delvar1``
(the array this loop actually fills). ``delvar1`` itself is never
independently quantized. Reads exactly like a copy-paste slip in the
Fortran. Replicated literally below -- not "corrected" to round
``delvar1`` instead.

Locals read before any assignment reaches them on some path
-----------------------------------------------------------------
Fortran locals are not implicitly zeroed; a handful here are read on a
path that never assigns them first in that call:

* ``hatm`` is not in the SAVE list and the zero-out block does not cover
  it; it is read at line 940 (``qv_max`` on the ``hm <= eps`` branch)
  while only ever assigned inside the ``i > ntemp`` (air-layer) branch of
  the big node loop (line 631). If the profile has no air layers above
  ``ntemp`` that assignment never runs before line 940 is reached.
* ``qtop``, ``qbot``, ``sigflo`` are set at lines 407-409 only on the
  ``infer_test == 0`` path of the ``first_time == 0`` branch;
  ``infer_test == 1`` skips straight past them to ``flow_param`` at line
  937, which reads ``qtop``.
* the initial-loop reads of ``dmet(11)`` for ``'WA'``/``'AI'`` nodes
  (lines 420, 423) happen BEFORE ``dmet`` is first populated at line
  482 onward -- ``dmet`` is a plain (non-SAVE) local, so this is reading
  whatever the compiler happened to leave on the stack.

Fortran's undefined-memory behavior in these cases cannot be replicated in
Python (there is no leftover stack garbage to read). Following the
precedent set in :mod:`fasst.surfenergy` (``smr``/``taf``, module
docstring / doc point there), every such local is initialized to ``0.0``
at the top of this function, with the same one-line rationale repeated at
each site: Fortran would silently use whatever was already there
(typically ``0d0`` in practice, but not guaranteed).

A known, NOT modeled, cross-procedure hazard: ``pdens`` aliasing
-----------------------------------------------------------------------
``new_profile``'s own ``pdens`` argument is ``intent(in)``, but
``surfenergy.F90``'s ``pdens`` dummy is ``intent(inout)`` (see
:class:`~fasst.surfenergy.SurfaceEnergy`'s ``pdens`` field docstring: it is
unconditionally zeroed inside ``surfenergy`` and only written afterward).
Because Fortran subroutine calls without an explicit interface pass
arguments by reference regardless of declared intent, it is possible the
Fortran runtime silently mutates ``new_profile``'s own ``pdens`` through
this aliasing on every ``surfenergy`` call. Python has no equivalent
by-reference float aliasing, so this is NOT modeled: ``pdens`` is treated
here as genuinely read-only, matching its DOCUMENTED ``intent(in)``
contract rather than its possible runtime behavior. Flagged, not resolved
-- outside this pass's scope.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Tuple

import numpy as np

from .constants import (
    EPS,
    EXTRAN,
    GRAV,
    HT_MIN,
    HT_MINM,
    KVEG,
    MAXCOL,
    MAXN,
    PI,
    RD,
    RV,
    SPFLAG,
    TREF,
)
from .functions import Phase, dense, head, maxinfiltrate, soilhumid, vap_press
from .indices import DerivedMet as D
from .indices import MetCol
from .sflux import sflux
from .sp_humid import sp_humid
from .state import FasstState
from .surfenergy import LowVegMetResult, surfenergy

__all__ = ["NewProfileState", "NewProfileResult", "new_profile"]

_MAXITER = 10  # Fortran integer(ip),parameter:: maxiter = 10
_ALLOWERRORT = 2e-2  # K, acceptable error for temperature iterations
_ALLOWERRORM = 2e-3  # m, acceptable error for pressure-head iterations


# --- phantom-slot array builders (index 0 unused; see fasst.state) --------
# Node-indexed locals are sized MAXN + EXTRAN, not the Fortran ``ntot`` --
# see module docstring, "x/stcalc array bound".
_NODESZ = MAXN + EXTRAN


def _reals(n: int) -> np.ndarray:
    return np.zeros(n + 1, dtype=np.float64)


def _reals2(n: int, m: int) -> np.ndarray:
    return np.zeros((n + 1, m + 1), dtype=np.float64)


# --- Fortran rounding idioms (mirrors fasst/functions.py::_anint) ---------
def _anint(x: float, p: int) -> float:
    """Fortran ``ANINT(x * 10**p) * 10**-p`` -- round half away from zero."""
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def _round20(x: float) -> float:
    return _anint(x, 20)


def _round18(x: float) -> float:
    return _anint(x, 18)


def _round15(x: float) -> float:
    return _anint(x, 15)


def _round10(x: float) -> float:
    return _anint(x, 10)


def _round10_inv(x: float) -> float:
    """Fortran idiom ``anint(x*1d-10)*1d10`` -- see module docstring, point 6."""
    return _anint(x, -10)


def _trunc10(x: float) -> float:
    """Fortran ``DINT(x*1d10)*1d-10`` -- truncate toward zero at 1e-10."""
    if not math.isfinite(x):
        return x
    return math.trunc(x * 1e10) * 1e-10


def _aint(x: float) -> float:
    """Fortran ``AINT``/``DINT`` with no scaling -- truncate toward zero."""
    if not math.isfinite(x):
        return x
    return float(math.trunc(x))


@dataclass
class NewProfileState:
    """The ~21-field ``save::`` block of ``new_profile.F90`` (lines 40-46),
    promoted to instance state -- see module docstring, point 2.

    Instantiate one ``NewProfileState`` per run/grid-point (alongside its
    :class:`~fasst.state.FasstState`) and pass it into every call of
    :func:`new_profile` for that point; it is mutated in place.
    """

    iwm: int = 0
    iwt: int = 0
    iwmr: int = 0
    iwtr: int = 0
    iwsm: int = 0
    ntempo: int = 0
    infl_cum: float = 0.0
    tcum: float = 0.0
    klhtop: float = 0.0
    s2i: float = 0.0
    thetai: float = 0.0
    klhtopi: float = 0.0
    airo: float = 0.0
    totice: float = 0.0
    maxrt: float = 0.0
    maxt: float = 0.0
    maxrm: float = 0.0
    maxm: float = 0.0
    x1: np.ndarray = field(default_factory=lambda: _reals(4))
    x: np.ndarray = field(default_factory=lambda: _reals2(14, _NODESZ))
    stcalc: np.ndarray = field(default_factory=lambda: _reals(_NODESZ))


@dataclass
class NewProfileResult:
    """The 5 ``intent(inout)`` scalar arguments of ``new_profile.F90``.

    ``pdens``, ``oldsd``, ``simelt`` are ``intent(in)`` in the Fortran and
    are not returned. See module docstring for why ``code5`` never actually
    changes (its only writer is commented out) and why ``phie``/``fsup``/
    ``firup`` are conditionally-updated pass-throughs, not pure outputs.
    """

    code5: int
    first_time: int
    phie: float
    fsup: float
    firup: float


# Type aliases for the six untranscribed dependencies -- see module
# docstring, "Injected untranscribed dependencies", for the exact contract
# (argument order, what is returned vs. mutated in place) of each.
AlbedoEmisFn = Callable[[FasstState, float, float, float], None]
ThParamFn = Callable[[FasstState, int, float, float, np.ndarray, np.ndarray], float]
FlowParamFn = Callable[
    [FasstState, int, float, float, float, np.ndarray, np.ndarray],
    Tuple[float, float, float, float],
]
SubDivideFn = Callable[[FasstState, int, np.ndarray, np.ndarray, np.ndarray], None]
SoilMoistureFn = Callable[..., Tuple[float, float, float, float]]
SoilTmpFn = Callable[..., Tuple[int, float, float, float]]


# ===========================================================================
def new_profile(
    state: FasstState,
    save: NewProfileState,
    code5: int,
    first_time: int,
    pdens: float,
    oldsd: float,
    simelt: float,
    phie: float,
    fsup: float,
    firup: float,
    *,
    albedo_emis_fn: AlbedoEmisFn,
    th_param_fn: ThParamFn,
    flow_param_fn: FlowParamFn,
    sub_divide_fn: SubDivideFn,
    soil_moisture_fn: SoilMoistureFn,
    soil_tmp_fn: SoilTmpFn,
    lowveg_met: Optional[Callable[..., LowVegMetResult]] = None,
) -> NewProfileResult:
    """Advance the soil/snow/vegetation profile by one meteorological
    timestep.

    Fortran: new_profile.F90 lines 1-1860, ``subroutine new_profile``. See
    module docstring for structure, the six injected dependencies'
    contracts, and the numbered fidelity notes (points 1-14 plus the extra
    findings).

    Args:
        state: model state; see module docstring, "Not pure", for the
            (very long) list of fields read and written.
        save: this routine's own persistent (Fortran ``save::``) state --
            see :class:`NewProfileState` and module docstring point 2.
        code5, first_time, pdens, oldsd, simelt, phie, fsup, firup: the
            Fortran subroutine's own arguments, in declaration order.
            ``pdens``, ``oldsd``, ``simelt`` are read-only (``intent(in)``).
            ``fsup``/``firup`` are accepted for interface parity with the
            Fortran ``intent(inout)`` declaration (same precedent as
            :func:`~fasst.surfenergy.surfenergy`'s ``pdens`` argument) but
            their incoming values are unconditionally discarded -- the
            Fortran zeroes them (lines 137, 140) before anything reads them.
            All five of ``code5, first_time, phie, fsup, firup`` are
            threaded through :class:`NewProfileResult`.
        albedo_emis_fn, th_param_fn, flow_param_fn, sub_divide_fn,
            soil_moisture_fn, soil_tmp_fn: required injected
            implementations of the six not-yet-transcribed dependencies --
            see module docstring for each one's exact contract. All six sit
            on the unconditional main path; this function cannot complete
            without real implementations.
        lowveg_met: forwarded to every :func:`~fasst.surfenergy.surfenergy`
            call; required only if ``sigfl > eps`` and ``icase != 4`` is
            ever reached (see :mod:`fasst.surfenergy`).

    Returns:
        A :class:`NewProfileResult`. Also mutates ``state`` and ``save`` in
        place (see module docstring).
    """
    iw = state.iw
    nnodes = state.nnodes

    # -- zero-out block (Fortran lines 87-236) ------------------------------
    ii = 0
    iter_ = 0
    sn = nnodes
    d0i = Phase.WATER_VAPOR
    d1i = Phase.WATER_VAPOR  # reassigned to Phase.WATER just below (d1i=1)
    d2i = Phase.WATER_VAPOR
    d3i = Phase.WATER_VAPOR
    iflag = 0
    iceflag = 0
    full = 0
    lflag = 0
    isn = 0
    errorm = 0.0
    rhs_errorm = 0.0
    smerror = 0.0
    errort = 0.0
    rhs_errort = 0.0
    evaprate = 0.0
    infl_max = 0.0
    a = 0.0
    b = 0.0
    t1 = 0.0
    vp = 0.0
    stempt = 0.0
    precip = 0.0
    taf = 0.0
    simeltsm = 0.0
    sheatf = 0.0
    sheatg = 0.0
    lhtf = 0.0
    lheatg = 0.0
    pheatf = 0.0
    pheatg = 0.0
    rhsurf = 0.0
    cc1 = 0.0
    d1 = 0.0
    rpp = 0.0
    dqdtf = 0.0
    dqdtg = 0.0
    sfac = 0.0
    sumrunoff = 0.0
    mixrf = 0.0
    mixrg = 0.0
    mixra = 0.0
    mixraf = 0.0
    kave = 0.0
    temp2 = 0.0
    sthick = 0.0
    fsdown = 0.0
    fsup = 0.0
    firdown = 0.0
    fvir = 0.0
    firup = 0.0
    fpheat = 0.0
    pdensnew = 0.0
    excess = 0.0
    isf = 0.0
    isg = 0.0
    fcheat = 0.0
    fcheat1 = 0.0
    cheati = 0.0
    cheatv = 0.0
    cheatw = 0.0
    fheat = 0.0
    rhotot = 0.0
    iru = 0.0
    pht1 = 0.0
    sh = 0.0
    lh = 0.0
    ch = 0.0
    ch1 = 0.0
    eh = 0.0
    mixrgs = 0.0
    fcheck = 0.0
    ftempo = 0.0
    dedt = 0.0
    t2 = 0.0
    swd = 0.0
    swu = 0.0
    ird = 0.0
    dnet1 = 0.0
    dnet2 = 0.0
    disf = 0.0
    disg = 0.0
    disfg = 0.0
    disgf = 0.0
    in_rate = 0.0
    plimit = 0.0
    qtopi = 0.0
    hpondi = 0.0
    ph = 0.0
    sert = 0.0
    qtopv = 0.0
    sumthick = 0.0
    state.melt[iw] = 0.0
    albedoi = 0.0
    emisi = 0.0
    Z = 0.0
    r = 0.0
    sert1 = 0.0
    mat = 0.0
    mat1 = 0.0
    mat2 = 0.0
    mat3 = 0.0
    mat4 = 0.0
    toticeo = 0.0
    rstep = 0.0
    qv_max = 0.0
    D = 0.0
    kvhtop = 0.0
    kvttop = 0.0
    tsign = 0.0
    lh1 = 0.0
    shs = 0.0
    lhs = 0.0
    mixas = 0.0
    mixgs = 0.0

    # locals read on a path that never assigns them first in this call --
    # Fortran leaves these as undefined stack memory; 0.0 is the honest
    # Python substitute (see module docstring, "Locals read before...").
    hatm = 0.0
    qtop = 0.0
    qbot = 0.0
    sigflo = 0.0
    pres = 0.0
    rh = 0.0

    t0 = _reals(4)
    delvar = _reals(MAXCOL)
    delvar1 = _reals(12)
    dmet = _reals(13)
    sinko = _reals(_NODESZ)
    zt = _reals(_NODESZ)
    delzt = _reals(_NODESZ)
    iceo = _reals(_NODESZ)
    wvco = _reals(_NODESZ)
    told = _reals(_NODESZ)
    fv1o = _reals(_NODESZ)
    ph_old = _reals(_NODESZ)
    sourceo = _reals(_NODESZ)
    vino = _reals(_NODESZ)
    flowuo = _reals(_NODESZ)
    flowlo = _reals(_NODESZ)
    rhov = _reals(_NODESZ)
    rhoda = _reals(_NODESZ)
    runoff = _reals(_NODESZ)
    sm_old = _reals(_NODESZ)
    grspheato = _reals(_NODESZ)
    grthcondo = _reals(_NODESZ)
    sinkro = _reals(_NODESZ)
    thvc = _reals(_NODESZ)
    dthvdt = _reals(_NODESZ)
    dthvdh = _reals(_NODESZ)
    temp = _reals(_NODESZ)
    y1 = _reals(4)
    y = _reals2(14, _NODESZ)

    ntot = state.ntot

    for i in range(1, ntot + 1):
        delzt[i] = 0.0
        zt[i] = 0.0
        state.flowu[i] = 0.0
        state.flowl[i] = 0.0
        state.fv1[i] = 0.0
        told[i] = state.stt[i]
        ph_old[i] = state.phead[i]
        vino[i] = state.vin[i]
        rhov[i] = 0.0
        rhoda[i] = 0.0
        runoff[i] = 0.0

    d1i = Phase.WATER
    d2i = Phase.ICE
    d3i = Phase.DRY_AIR

    # -- determine the new (internal) timestep (Fortran lines 238-284) -----
    state.deltat_fasst = state.deltati
    state.step = state.stepi

    lflag = 0
    if iw >= 2:
        if _aint(state.met[iw - 1, MetCol.PT]) == 2 or _aint(state.met[iw, MetCol.PT]) == 2:
            lflag = 1
    else:
        if _aint(state.met[iw, MetCol.PT]) == 2:
            lflag = 1
    if state.hm > EPS:
        lflag = 2

    if lflag == 1:
        state.deltat_fasst = state.deltati * 2.5e-1
        rstep = state.timstep * 3.6e3 / state.deltat_fasst
        state.step = int(_aint(rstep))
    elif lflag == 2:
        state.deltat_fasst = state.deltati * 2.5e-1
        rstep = state.timstep * 3.6e3 / state.deltat_fasst
        state.step = int(_aint(rstep))

    f1 = 1.0 / float(state.step)
    f1 = _round20(f1)
    f2 = 1.0 / state.deltat_fasst
    f2 = _round20(f2)

    # -- maximum fluid infiltration ------------------------------------------
    albedoi = state.albedo_fasst
    emisi = state.emis

    if abs(state.nsoilp[nnodes, 7] - SPFLAG) > EPS:
        infl_max = 1e-2 * state.nsoilp[nnodes, 7]
    else:
        infl_max = 0.0
    overland = 0.0
    plimit = 1e-1 * state.rough * max(0.0, min(1.0, math.cos(state.sloper)))

    # -- inflow due to snow and ice melting ----------------------------------
    simeltsm = simelt * f1
    simeltsm = _anint(simeltsm, 15)

    # -- initialize variables if first time through (lines 286-432) --------
    if first_time == 0:
        state.ntemp = 0
        save.ntempo = 0
        state.ntemp = nnodes
        save.ntempo = state.ntemp
        state.icase = 0
        state.icaseo = 0
        state.hpond = 0.0
        save.infl_cum = 0.0
        save.tcum = 0.0
        save.klhtop = 0.0
        save.klhtop = state.nsoilp[nnodes, 7] * 1e-2
        save.s2i = 0.0
        save.thetai = 0.0
        save.klhtopi = 0.0
        save.airo = 0.0
        save.totice = 0.0
        phie = 6.5e-1
        state.ftemp = state.dmet1[iw, 4]

        save.maxrt = -abs(state.mflag)
        save.maxt = -abs(state.mflag)
        save.maxrm = -abs(state.mflag)
        save.maxm = -abs(state.mflag)
        maxsm = -abs(state.mflag)
        save.iwm = 0
        save.iwt = 0
        save.iwmr = 0
        save.iwtr = 0
        save.iwsm = 0

        if state.infer_test == 1:
            for i in range(1, ntot + 1):
                told[i] = 0.0
                wvco[i] = 0.0
                iceo[i] = 0.0
                ph_old[i] = 0.0
                sm_old[i] = 0.0
                sourceo[i] = 0.0
                sinko[i] = 0.0
                sinkro[i] = 0.0
                flowuo[i] = 0.0
                flowlo[i] = 0.0
                fv1o[i] = 0.0
                vino[i] = 0.0
                save.stcalc[i] = 0.0

                told[i] = state.stt[i]
                wvco[i] = state.wvc[i]
                iceo[i] = state.ice[i]
                ph_old[i] = state.phead[i]
                sm_old[i] = state.soil_moist[i]
                sourceo[i] = state.source[i]
                sinko[i] = state.sink[i]
                sinkro[i] = state.sinkr[i]
                flowuo[i] = state.flowu[i]
                flowlo[i] = state.flowl[i]
                fv1o[i] = state.fv1[i]
                vino[i] = state.vin[i]
                save.stcalc[i] = state.stt[i]

            if iw >= 2:
                save.airo = state.met[iw - 1, MetCol.TMP] + TREF
            else:
                save.airo = state.met[iw, MetCol.TMP] + TREF
        else:
            for i in range(1, ntot + 1):
                wvco[i] = 0.0
                iceo[i] = 0.0
                ph_old[i] = 0.0
                sourceo[i] = 0.0
                sinko[i] = 0.0
                sinkro[i] = 0.0
                flowuo[i] = 0.0
                flowlo[i] = 0.0
                fv1o[i] = 0.0
                vino[i] = 0.0

                state.too[i] = 0.0
                state.woo[i] = 0.0
                state.ioo[i] = 0.0
                state.smoo[i] = 0.0
                state.phoo[i] = 0.0

                if iw == state.istart:
                    wvco[i] = state.wvc[i]
                    iceo[i] = state.ice[i]
                    ph_old[i] = state.phead[i]
                    sourceo[i] = state.source[i]
                    sinko[i] = state.sink[i]
                    sinkro[i] = state.sinkr[i]
                    flowuo[i] = state.flowu[i]
                    flowlo[i] = state.flowl[i]
                    fv1o[i] = state.fv1[i]
                    vino[i] = state.vin[i]
                    state.too[i] = state.stt[i]
                    state.smoo[i] = state.soil_moist[i]
                    state.woo[i] = state.wvc[i]
                    state.ioo[i] = state.ice[i]
                    state.phoo[i] = state.phead[i]
                    save.stcalc[i] = 0.0
                    save.stcalc[i] = state.stt[i]

            if iw == state.istart:
                for i in range(1, 5):
                    save.x1[i] = 0.0
                for i in range(1, 14):
                    for ii in range(1, MAXN + 1):
                        save.x[i, ii] = 0.0

            qtop = 0.0
            qbot = 0.0
            sigflo = 0.0
        pdensnew = pdens

        for i in range(1, nnodes + 1):
            pres = (
                state.dmet1[state.istart, 11]
                + 1e-2
                * (state.elev - state.nz[i] + abs(state.phead[i]))
                * dense(state.stt[i], 0.0, d1i)
                * GRAV
            )
            if state.node_type[i] != "WA" and state.node_type[i] != "AI":
                rh = soilhumid(state, i, state.phead[i], state.soil_moist[i], state.stt[i])
            elif state.node_type[i] == "WA":
                rh = 1.0
                pres = dmet[11]
            elif state.node_type[i] == "AI":
                rh = 1e-2 * dmet[5]
                pres = dmet[11]

            m = sp_humid(d0i, pres, state.stt[i], rh, state.phead[i])
            t0[1], t0[2], t0[3], t0[4] = m.mixr, m.dmrdt, m.vpress, m.wetbulb
            rhov[i], rhoda[i] = m.rhov, m.rhoda
            t1 = m.vpsat
            thvc[i], dthvdt[i], dthvdh[i] = m.thvc, m.dthvdt, m.dthvdh
    else:
        state.sdens[iw] = state.sdens[iw - 1]

    # -- initialize "old" variables (lines 434-445) -------------------------
    sigflo = state.sigfl
    ftempo = state.ftemp

    icode = 0
    for i in range(1, ntot + 1):
        flowuo[i] = state.flowu[i]
        flowlo[i] = state.flowl[i]
        fv1o[i] = state.fv1[i]
        if state.ice[i] > 0.0:
            icode = 1

    # -- met interpolation variables (lines 447-478) -------------------------
    for ll in range(1, MAXCOL + 1):
        delvar[ll] = 0.0
        if iw >= 2:
            if (
                _aint(abs(state.met[iw, ll] - state.mflag) * 1e5) * 1e-5 > EPS
                and _aint(abs(state.met[iw - 1, ll] - state.mflag) * 1e5) * 1e-5 > EPS
            ):
                delvar[ll] = (state.met[iw, ll] - state.met[iw - 1, ll]) * f1
                if ll == 10 or ll == 12:
                    delvar[ll] = max(0.0, state.met[iw, ll] / (state.timstep * 3.6e3))
            else:
                delvar[ll] = state.mflag
        delvar[ll] = _round20(delvar[ll])

    for ll in range(1, 13):
        delvar1[ll] = 0.0
        dmet[ll] = 0.0
        if iw >= 2:
            if (
                _aint(abs(state.dmet1[iw, ll] - state.mflag) * 1e5) * 1e-5 > EPS
                and _aint(abs(state.dmet1[iw - 1, ll] - state.mflag) * 1e5) * 1e-5 > EPS
            ):
                delvar1[ll] = (state.dmet1[iw, ll] - state.dmet1[iw - 1, ll]) * f1
                if ll == 6 or ll == 7:
                    delvar1[ll] = max(0.0, state.dmet1[iw, ll] / (state.timstep * 3.6e3))
            else:
                delvar1[ll] = state.mflag
        # NOT a typo in this transcription -- the Fortran re-rounds `delvar`,
        # not `delvar1`, here. See module docstring, "delvar re-rounded twice".
        delvar[ll] = _round20(delvar[ll])

    isn = iw - 1
    if iw == 1:
        isn = iw
    dmet[1] = state.dmet1[isn, 1]
    dmet[2] = state.dmet1[isn, 2]
    dmet[3] = state.dmet1[isn, 3]
    dmet[4] = state.dmet1[isn, 4]
    dmet[5] = state.dmet1[isn, 5]
    dmet[6] = delvar1[6] * math.cos(state.sloper)
    dmet[6] = _anint(dmet[6], 15)
    if state.met[iw, MetCol.PT] == 2:
        runoff[:] = delvar1[6] * math.sin(state.sloper)
    dmet[7] = delvar1[7] * math.cos(state.sloper)
    dmet[7] = _anint(dmet[7], 15)
    dmet[8] = state.dmet1[isn, 8]
    dmet[9] = state.dmet1[isn, 9]
    dmet[10] = state.dmet1[isn, 10]
    dmet[11] = state.dmet1[isn, 11]
    dmet[12] = state.dmet1[isn, 12]
    dmet[13] = state.dmet1[isn, 13]

    # -- adjust node count for snow / low veg (lines 499-564) ---------------
    state.ntemp = nnodes
    state.icaseo = state.icase
    state.icase = 0
    if state.veg_flagl == 0 or state.sigfl <= EPS:
        if state.hm > EPS:
            state.icase = 1
            if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                state.ntemp = nnodes + 2
            else:
                state.ntemp = nnodes + 1
            for i in range(nnodes + 1, state.ntemp + 1):
                state.node_type[i] = "HM"
    else:
        if state.hm <= EPS:
            state.icase = 2
            state.ntemp = nnodes + 1
            state.node_type[state.ntemp] = "VG"
        else:
            if state.hfol_tot - state.hm <= EPS:
                state.icase = 4
                if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                    if abs(state.hfol_tot - state.hm) <= EPS or abs(state.hfol_tot - HT_MIN) <= EPS:
                        state.ntemp = nnodes + 2
                        if abs(state.hfol_tot - state.hm) <= EPS:
                            state.node_type[state.ntemp] = "MX"
                        else:
                            state.node_type[state.ntemp] = "HM"
                    else:
                        state.ntemp = nnodes + 3
                        state.node_type[state.ntemp] = "HM"
                        if state.hfol_tot > HT_MIN:
                            state.node_type[nnodes + 2] = "MX"
                        else:
                            state.node_type[nnodes + 2] = "HM"
                    state.node_type[nnodes + 1] = "MX"
                elif state.hm <= HT_MIN:
                    if abs(state.hfol_tot - state.hm) <= EPS:
                        state.ntemp = nnodes + 1
                    else:
                        state.ntemp = nnodes + 2
                    state.node_type[nnodes + 1] = "MX"
                    if state.ntemp == nnodes + 2:
                        state.node_type[nnodes + 2] = "MX"
            elif state.hfol_tot - state.hm > EPS:
                state.icase = 3
                if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                    state.ntemp = nnodes + 3
                else:
                    state.ntemp = nnodes + 2
                state.node_type[state.ntemp] = "VG"
                state.node_type[nnodes + 1] = "MX"
                if state.ntemp == nnodes + 3:
                    state.node_type[nnodes + 2] = "MX"

    rhotot = th_param_fn(state, state.ntemp, pdens, pdensnew, rhov, rhoda)

    state.ftemp = dmet[4]
    sn = nnodes

    ntemp = state.ntemp
    icase = state.icase
    icaseo = state.icaseo

    for i in range(1, ntot + 1):
        told[i] = state.too[i]
        iceo[i] = state.ioo[i]
        wvco[i] = state.woo[i]
        ph_old[i] = state.phoo[i]
        grthcondo[i] = state.grthcond[i]
        grspheato[i] = state.grspheat[i]

        if state.node_type[i] in ("HM", "SN", "MX"):
            state.stt[i] = save.stcalc[i]

        if abs(told[i]) <= EPS:
            told[i] = 1.0

        if i <= nnodes:
            sm_old[i] = state.smoo[i]

            zt[i] = state.nz[i]
            sourceo[i] = state.source[i]
            sinko[i] = state.sink[i]
            sinkro[i] = state.sinkr[i]

            if abs(state.stt[i]) <= EPS:
                state.stt[i] = 1.0
            pres = dmet[11] - 1e-2 * state.phead[i] * dense(state.stt[i], 0.0, d1i) * GRAV
            if state.node_type[i] != "WA" and state.node_type[i] != "AI":
                rh = soilhumid(state, i, state.phead[i], state.soil_moist[i], state.stt[i])
            elif state.node_type[i] == "WA":
                rh = 1.0
                pres = dmet[11]
            elif state.node_type[i] == "AI":
                rh = 1e-2 * dmet[5]
                pres = dmet[11]

            m = sp_humid(d0i, pres, state.stt[i], rh, state.phead[i])
            t0[1], t0[2], t0[3], t0[4] = m.mixr, m.dmrdt, m.vpress, m.wetbulb
            rhov[i], rhoda[i] = m.rhov, m.rhoda
            t1 = m.vpsat
            thvc[i], dthvdt[i], dthvdh[i] = m.thvc, m.dthvdt, m.dthvdh
            if state.ntype[i] == 27:
                state.wvc[i] = rhov[i] / rhoda[i]

        elif i > ntemp:
            if ntemp < save.ntempo:
                # Fortran: told(i) = too(i+1); ph_old(i) = phoo(i+1) -- reads
                # the STATE arrays at i+1, not the local told/ph_old. Guarded
                # against i+1 running past the array (see module docstring,
                # "x/stcalc array bound" -- same class of latent OOB).
                told[i] = state.too[i + 1] if i + 1 < state.too.shape[0] else 0.0
                iceo[i] = 0.0
                wvco[i] = 0.0
                ph_old[i] = state.phoo[i + 1] if i + 1 < state.phoo.shape[0] else 0.0
            pres = dmet[11]
            sourceo[i] = 0.0
            sinko[i] = 0.0
            sinkro[i] = 0.0

            rh = 1e-2 * dmet[5]
            vp = vap_press(state, i, rh, pres)

            state.phead[i] = -vp / (RV * state.stt[i] * GRAV) + qtopv * state.deltat_fasst
            if state.hpond > EPS:
                state.phead[i] = 0.0
            if state.phead[i] > 0.0:
                state.phead[i] = 0.0
            state.phead[i] = _round20(state.phead[i])
            if i == ntot:
                hatm = state.phead[i]

            state.ice[i] = 0.0
            state.wvc[i] = 0.0
            state.stt[i] = dmet[4]
            state.node_type[i] = "AI"
            zt[i] = state.elev + max(state.hm, state.hfol_tot) + 10.0 * float(i - nnodes)

        elif nnodes < i <= ntemp:
            if icase == 1:
                sn = ntemp
                if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                    if i == nnodes + 1:
                        zt[i] = state.elev + (state.hm - HT_MIN)
                    else:
                        zt[i] = state.elev + state.hm
                else:
                    zt[i] = state.elev + state.hm
                state.nz[i] = zt[i]

                if icaseo == 0 or ntemp != save.ntempo:
                    told[i] = min(dmet[4], TREF)
                    state.too[i] = min(dmet[4], TREF)
                    state.stt[i] = min(dmet[4], TREF)
                    y[1, i] = state.stt[i]
                    save.x[1, i] = state.stt[i]
                    state.toptemp = state.stt[i]
                    iceo[i] = 0.0
                    wvco[i] = phie
                    ph_old[i] = 0.0
                elif icaseo == 1 and ntemp > save.ntempo and i == ntemp:
                    told[i] = told[i - 1]
                    state.too[i] = state.too[i - 1]
                    state.stt[i] = state.stt[i - 1]
                    y[1, i] = state.stt[i]
                    save.x[1, i] = state.stt[i]
                    state.toptemp = state.stt[i]
                    iceo[i] = iceo[i - 1]
                    wvco[i] = wvco[i - 1]
                    ph_old[i] = 0.0
                elif icaseo == 1 and ntemp < save.ntempo and i == ntemp:
                    state.toptemp = state.stt[i]

                rh = 1.0
                vp = vap_press(state, i, rh, pres)
                rhoda[i] = (pres - vp) / (RD * state.stt[i])
                rhov[i] = vp / (RV * state.stt[i])

                state.wvc[i] = phie
                state.ice[i] = (state.refreeze + state.refreezei) * f1
                state.ice[i] = _round20(state.ice[i])

                state.phead[i] = 0.0

            elif icase == 2:
                zt[i] = state.elev + state.hfol_tot
                state.nz[i] = zt[i]

                state.wvc[i] = 0.0
                wvco[i] = 0.0
                state.ice[i] = 0.0
                iceo[i] = 0.0
                if icaseo == 4 or icaseo == 3 or (first_time == 0 and state.infer_test == 0):
                    state.stt[i] = ftempo
                    told[i] = ftempo
                    state.too[i] = ftempo
                    # Fortran: ph_old(i) = phoo(i+1) -- state array at i+1.
                    ph_old[i] = state.phoo[i + 1] if i + 1 < state.phoo.shape[0] else 0.0

                rh = 1e-2 * dmet[5]
                vp = vap_press(state, i, rh, pres)

                state.phead[i] = -vp / (RV * state.stt[i] * GRAV) + qtopv * state.deltat_fasst
                if state.hpond > EPS:
                    state.phead[i] = 0.0
                if state.phead[i] > 0.0:
                    state.phead[i] = 0.0
                state.phead[i] = _round20(state.phead[i])

            elif icase == 3:
                if i == ntemp:
                    zt[i] = state.elev + state.hfol_tot
                    state.nz[i] = zt[i]

                    if icaseo == 2:
                        told[i] = ftempo
                        state.too[i] = ftempo
                        state.stt[i] = ftempo
                    state.wvc[i] = 0.0
                    wvco[i] = 0.0
                    state.ice[i] = 0.0
                    iceo[i] = 0.0

                    rh = 1e-2 * dmet[5]
                    vp = vap_press(state, i, rh, pres)

                    state.phead[i] = -vp / (RV * state.stt[i] * GRAV) + qtopv * state.deltat_fasst
                    if state.hpond > EPS:
                        state.phead[i] = 0.0
                    if state.phead[i] > 0.0:
                        state.phead[i] = 0.0
                    state.phead[i] = _round20(state.phead[i])
                else:
                    sn = i
                    if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                        if i == nnodes + 1:
                            zt[i] = state.elev + (state.hm - HT_MIN)
                        else:
                            zt[i] = state.elev + state.hm
                    else:
                        zt[i] = state.elev + state.hm
                    state.nz[i] = zt[i]
                    state.node_type[i] = "MX"

                    if icaseo == 2:
                        told[i] = min(dmet[4], TREF)
                        state.too[i] = min(dmet[4], TREF)
                        state.stt[i] = min(dmet[4], TREF)
                        wvco[i] = phie
                        state.toptemp = state.stt[i]
                        iceo[i] = 0.0
                    elif icaseo == 3 and ntemp > save.ntempo and i == ntemp - 1:
                        told[i] = state.toptemp
                        state.too[i] = state.toptemp
                        state.stt[i] = state.toptemp
                        wvco[i] = wvco[i - 1]
                        iceo[i] = iceo[i - 1]

                    rh = 1.0
                    vp = vap_press(state, i, rh, pres)
                    rhoda[i] = (pres - vp) / (RD * state.stt[i])
                    rhov[i] = vp / (RV * state.stt[i])

                    state.wvc[i] = phie
                    state.ice[i] = (state.refreeze + state.refreezei) * f1
                    state.ice[i] = _round10(state.ice[i])

                    state.phead[i] = 0.0

            elif icase == 4:
                sn = ntemp
                if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                    if ntemp == nnodes + 3:
                        if i == nnodes + 1:
                            zt[i] = state.elev + min(state.hfol_tot, HT_MIN)
                            state.ftemp = state.stt[i]
                        elif i == nnodes + 2:
                            zt[i] = state.elev + (state.hm - max(state.hfol_tot, HT_MIN))
                        elif i == ntemp:
                            zt[i] = state.elev + state.hm
                    elif ntemp == nnodes + 2:
                        if i == nnodes + 1:
                            zt[i] = state.elev + state.hfol_tot
                        else:
                            zt[i] = state.elev + state.hm
                else:
                    if i == nnodes + 1:
                        zt[i] = state.elev + state.hfol_tot
                        state.ftemp = state.stt[i]
                    else:
                        zt[i] = state.elev + state.hm
                state.nz[i] = zt[i]

                if icaseo == 2:
                    told[i] = min(dmet[4], TREF)
                    state.too[i] = min(dmet[4], TREF)
                    state.stt[i] = min(dmet[4], TREF)
                    state.toptemp = state.stt[i]
                    iceo[i] = 0.0
                    wvco[i] = phie
                elif icaseo == 3 and i == ntemp:
                    state.stt[i] = state.stt[i - 1]
                    told[i] = told[i - 1]
                    state.too[i] = state.too[i - 1]
                    state.toptemp = state.stt[i]
                    iceo[i] = iceo[i - 1]
                    wvco[i] = wvco[i - 1]
                elif icaseo == 4 and ntemp > save.ntempo and i == ntemp:
                    state.stt[i] = state.stt[i - 1]
                    told[i] = told[i - 1]
                    state.too[i] = state.too[i - 1]
                    state.toptemp = state.stt[i]
                    iceo[i] = iceo[i - 1]
                    wvco[i] = wvco[i - 1]

                rh = 1.0
                vp = vap_press(state, i, rh, pres)
                rhoda[i] = (pres - vp) / (RD * state.stt[i])
                rhov[i] = vp / (RV * state.stt[i])

                state.wvc[i] = phie
                state.ice[i] = (state.refreeze + state.refreezei) * f1
                state.ice[i] = _round10(state.ice[i])

                state.phead[i] = 0.0

        state.source[i] = 0.0
        state.sink[i] = 0.0
        state.sinkr[i] = 0.0

    # -- node thicknesses (lines 861-877) ------------------------------------
    for i in range(1, ntot + 1):
        zt[i] = _round20(zt[i])
        if i > ntemp:
            delzt[i] = 0.0
        elif nnodes <= i <= ntemp:
            if i == ntemp:
                delzt[i] = zt[i] - zt[i - 1]
            else:
                delzt[i] = 5e-1 * (zt[i + 1] - zt[i - 1])
            delzt[i] = _round18(delzt[i])
        else:
            delzt[i] = _round18(state.delzs[i])

    # -- initial surface energy flux (lines 879-905) -------------------------
    stempt = state.stt[sn]
    kave = 5e-1 * (state.grthcond[sn] + state.grthcond[sn - 1])
    if state.hm > EPS:
        kave = state.grthcond[sn]
    sthick = _trunc10(zt[sn] - zt[sn - 1])

    if state.hpond > EPS:
        if state.met[iw, MetCol.ZEN] >= 9e1:
            state.sgralbedo = 0.0
        else:
            Z = state.met[iw, MetCol.ZEN] * PI / 180.0
            r = math.asin(math.sin(Z) / 1.33)
            state.albedo_fasst = 5e-1 * (
                math.sin(Z - r) ** 2 / math.sin(Z + r) ** 2
                + math.tan(Z - r) ** 2 / math.tan(Z + r) ** 2
            )
        state.emis = state.nsoilp[nnodes, 4]
    else:
        state.albedo_fasst = albedoi
        state.emis = emisi

    ii = 0
    iter_ = 0
    fe = surfenergy(state, ii, iter_, sn, kave, sthick, pdens, dmet, lowveg_met=lowveg_met)
    evaprate = fe.evapcm
    rhsurf = fe.rhsurf
    pheatf, pheatg = fe.pheatf, fe.pheatg
    mixrf, mixraf, mixrg, mixra = fe.mixrf, fe.mixraf, fe.mixrgrs, fe.mixra
    lhtf, lheatg = fe.lhtf, fe.lheatg
    sheatf, sheatg = fe.sheatf, fe.sheatg
    dqdtf, dqdtg = fe.dqdtf, fe.dqdtg
    sfac, d1, rpp, cc1 = fe.sfac, fe.d1, fe.rpp, fe.cc1
    isg, isf = fe.isurfg, fe.isurff
    disg, disf, disfg, disgf = fe.disurfg, fe.disurff, fe.disurffg, fe.disurfgf
    lh1 = fe.lh  # discarded -- see module docstring, point 11

    rain = dmet[6]
    snow = dmet[7]

    if iw == state.istart and state.infer_test == 0:
        state.isurfoldg = isg
        state.isurfoldf = isf

    # -- maximum infiltration rate (lines 915-933) ---------------------------
    if save.tcum <= EPS or save.infl_cum <= EPS:
        save.thetai = state.soil_moist[nnodes]
        save.s2i = maxinfiltrate(state, nnodes, save.thetai)
        save.klhtopi = save.klhtop
        in_rate = infl_max
    elif save.tcum > EPS and save.infl_cum > EPS and save.s2i > EPS:
        t1 = 2.0 * (save.infl_cum - save.klhtopi * save.tcum) * infl_max / save.s2i
        if abs(t1) < 5e1:
            t1 = math.exp(2.0 * (save.infl_cum - save.klhtopi * save.tcum) * infl_max / save.s2i)
            if math.cos(state.sloper) - 1.0 / t1 >= EPS:
                in_rate = infl_max * min(1.0, math.cos(state.sloper))
            else:
                in_rate = infl_max * min(1.0, math.cos(state.sloper) - 1.0 / t1 + 1.0)
            in_rate = max(0.0, in_rate)

    # -- maximum allowed evaporation rate (lines 935-958) --------------------
    isn = nnodes
    qbot, save.klhtop, kvhtop, kvttop = flow_param_fn(
        state, isn, qtop, qtopv, simeltsm, zt, delzt
    )
    if state.hm <= EPS:
        qv_max = kvhtop * max(
            0.0,
            (hatm - state.phead[nnodes] + 1.0) + kvttop * (dmet[4] - state.stt[nnodes]),
        )
    else:
        if phie > EPS and state.stt[sn] > EPS:
            D = max(0.0, (phie ** (5.0 / 3.0)) * (2.12e-5 * (state.stt[sn] / TREF) ** 2.0))
        kvhtop = max(
            0.0,
            D
            * (vap_press(state, sn, 1.0, dmet[11]) / (RV * state.stt[sn]))
            * GRAV
            / (dense(state.stt[sn], 0.0, d1i) * RV * state.stt[sn]),
        )
        qv_max = kvhtop * max(0.0, (hatm - state.phead[sn] + 1.0))

    tsign = 0.0
    if evaprate != 0.0:
        tsign = abs(evaprate) / evaprate
    elif qv_max != 0.0:
        tsign = abs(qv_max) / qv_max
    qtopv = 0.0
    qtopv = -tsign * min(abs(evaprate), abs(qv_max))
    if abs(evaprate) <= EPS:
        qtopv = qv_max

    # -- surface water balance (lines 960-999) --------------------------------
    extra = 0.0
    hpondi = state.hpond
    if state.hm <= EPS and state.node_type[nnodes] != "WA":
        precip = rain
        if _aint(state.met[iw, MetCol.PT]) != 2:
            precip = 0.0

        if precip > EPS and qtopv > EPS:
            qtopv = 0.0

        if state.hpond > EPS:
            state.hpond = hpondi + (qtopv + precip) * state.deltat_fasst + simeltsm
            extra = qtopv
            qtopv = 0.0
            if state.hpond < EPS:
                extra = 0.0
                qtopv = state.hpond * f2
                state.hpond = 0.0

            qtopi = state.hpond * f2
            if state.hpond < EPS:
                state.hpond = 0.0
        else:
            qtopi = precip + state.hpond * f2 + simeltsm * f2
    else:
        qtopi = state.hpond * f2 + simeltsm * f2
        if state.hm > EPS:
            extra = qtopv
            qtopv = 0.0
    qtopi = _round20(qtopi)
    qtopv = _round20(qtopv)

    qtop = qtopi

    if state.node_type[nnodes] != "WA":
        if state.hpond > EPS or state.nsoilp[nnodes, 24] - (
            state.soil_moist[nnodes] + state.ice[nnodes]
        ) < EPS:
            qtop = max(
                0.0,
                (state.nsoilp[nnodes, 24] - (state.soil_moist[nnodes] + state.ice[nnodes]))
                * state.delzs[nnodes]
                * f2,
            )
            state.hpond = max(0.0, state.hpond - qtop * state.deltat_fasst)
            if state.hpond > plimit:
                overland = overland + (state.hpond - plimit)
                state.hpond = plimit
    elif state.node_type[nnodes] == "WA":
        qtop = 0.0
        state.hpond = 0.0
        overland = 0.0
    state.hpond = _round20(state.hpond)

    if state.hpond > EPS:
        if state.met[iw, MetCol.ZEN] >= 9e1:
            state.sgralbedo = 0.0
        else:
            Z = state.met[iw, MetCol.ZEN] * PI / 180.0
            r = math.asin(math.sin(Z) / 1.33)
            state.albedo_fasst = 5e-1 * (
                math.sin(Z - r) ** 2 / math.sin(Z + r) ** 2
                + math.tan(Z - r) ** 2 / math.tan(Z + r) ** 2
            )
        state.albedo_fasst = _round20(state.albedo_fasst)
        state.emis = state.nsoilp[nnodes, 4]
    else:
        state.albedo_fasst = albedoi
        state.emis = emisi

    hatm = hatm + qtopv * state.deltat_fasst

    isn = 1
    qbot, save.klhtop, kvhtop, kvttop = flow_param_fn(
        state, isn, qtop, qtopv, simeltsm, zt, delzt
    )

    if first_time == 0:
        first_time = 1

    if iw == 1:
        sumo = 0.0
        for i in range(1, ntot + 1):
            temp[i] = 0.0
            if i <= nnodes:
                rhow = dense(state.stt[i], 0.0, d1i)
                if state.ice[i] > EPS:
                    rhoi = dense(state.stt[i], 0.0, d2i)
                    temp[i] = state.soil_moist[i] + state.ice[i] * rhoi / rhow
                else:
                    temp[i] = state.soil_moist[i]
                if state.wvc[i] > EPS:
                    temp[i] = temp[i] + state.wvc[i] * rhov[i] / rhow

                sumo = sumo + temp[i] * state.delzs[i]
    elif iw >= 2:
        sumo = state.tot_moist[iw - 1]

    # *************************************************************************
    # MAIN CALCULATION LOOP (lines 1070-1614)
    for ii in range(1, state.step + 1):
        sub_divide_fn(state, ii, delvar, delvar1, dmet)

        albedo_emis_fn(state, oldsd, dmet[1], dmet[8])

        iter_ = 0
        iflag = 0
        errorm = 1.1e0 * _ALLOWERRORM
        errort = 1.1e0 * _ALLOWERRORT
        sert = rhs_errort
        sert1 = rhs_errorm
        # ---------------------------------------------------------------
        # INNER (Newton) ITERATION LOOP
        while errorm > _ALLOWERRORM or errort > _ALLOWERRORT:
            lflag = 0
            if iter_ > 0:
                mat1 = mat
                mat4 = mat3
                mat3 = mat2
            mat = sert
            mat2 = sert1

            if ii > 1:
                if iter_ == 0:
                    for i in range(1, ntot + 1):
                        if (y[3, i] > EPS and save.x[3, i] < EPS) or (
                            y[3, i] < EPS and save.x[3, i] > EPS
                        ):
                            told[i] = y[1, i]
                            wvco[i] = y[2, i]
                            iceo[i] = y[3, i]
                            sm_old[i] = y[4, i]
                            ph_old[i] = y[5, i]
                        else:
                            told[i] = save.x[1, i]
                            wvco[i] = save.x[2, i]
                            iceo[i] = save.x[3, i]
                            sm_old[i] = save.x[4, i]
                            ph_old[i] = save.x[5, i]

                        sinko[i] = save.x[6, i]
                        sinkro[i] = save.x[7, i]
                        sourceo[i] = save.x[8, i]

                        if i <= nnodes:
                            flowuo[i] = y[9, i]
                            flowlo[i] = y[14, i]
                            fv1o[i] = y[10, i]
                            vino[i] = y[11, i]

                        grthcondo[i] = y[12, i]
                        grspheato[i] = y[13, i]
                    state.isurfoldg = save.x1[1]
                    state.isurfoldf = save.x1[2]
                    sigflo = save.x1[3]
                    mixrgo = save.x1[4]
                else:
                    for i in range(1, ntot + 1):
                        told[i] = y[1, i]
                        wvco[i] = y[2, i]
                        iceo[i] = y[3, i]
                        sm_old[i] = y[4, i]
                        ph_old[i] = y[5, i]

                        sinko[i] = y[6, i]
                        sinkro[i] = y[7, i]
                        sourceo[i] = y[8, i]

                        if i <= nnodes:
                            flowuo[i] = y[9, i]
                            flowlo[i] = y[14, i]
                            fv1o[i] = y[10, i]
                            vino[i] = y[11, i]

                        grthcondo[i] = y[12, i]
                        grspheato[i] = y[13, i]
                    state.isurfoldg = y1[1]
                    state.isurfoldf = y1[2]
                    sigflo = y1[3]
                    mixrgo = y1[4]

            full = 0
            for i in range(1, ntot + 1):
                save.x[1, i] = state.stt[i]
                save.x[2, i] = state.wvc[i]
                save.x[3, i] = state.ice[i]
                save.x[4, i] = state.soil_moist[i]
                save.x[5, i] = state.phead[i]
                save.x[6, i] = state.sink[i]
                save.x[7, i] = state.sinkr[i]
                save.x[8, i] = state.source[i]
                if i <= nnodes:
                    save.x[9, i] = state.flowu[i]
                    save.x[14, i] = state.flowl[i]
                    save.x[10, i] = state.fv1[i]
                    save.x[11, i] = state.vin[i]
                save.x[12, i] = state.grthcond[i]
                save.x[13, i] = state.grspheat[i]
                if state.nsoilp[i, 24] - state.soil_moist[i] < EPS:
                    full = full + 1
            save.x1[1] = isg
            save.x1[2] = isf
            save.x1[3] = state.sigfl
            save.x1[4] = mixrg

            if errorm > _ALLOWERRORM:
                rhs_errorm, errorm, smerror, sert1 = soil_moisture_fn(
                    state, sm_old, wvco, iceo, sourceo, sinkro, vino, told,
                    thvc, dthvdh, runoff,
                )
            else:
                lflag = 2

            # `iceflag` is soil_tmp's `iflag`, which is intent(out) ONLY in
            # the Fortran (no input value is needed) -- see module
            # docstring's soil_tmp_fn contract.
            iceflag, rhs_errort, errort, sert = soil_tmp_fn(
                state, isg, isf, disg, disf, disfg, disgf, sigflo,
                rhotot, save.airo, dmet, grthcondo, grspheato, zt, delzt,
                told, iceo, wvco, flowuo, flowlo, fv1o, sinko, sourceo,
                ph_old, sm_old, thvc, dthvdt, rhov, rhoda, ii, iter_,
            )

            if iw == state.istart and ii == 1 and iter_ == 0:
                for i in range(1, nnodes + 1):
                    if state.ice[i] > EPS and iceo[i] <= EPS:
                        iceo[i] = state.ice[i]
                        state.ioo[i] = state.ice[i]
                        sm_old[i] = state.soil_moist[i]
                        state.smoo[i] = state.soil_moist[i]
                        wvco[i] = state.wvc[i]
                        state.woo[i] = state.wvc[i]
                        ph_old[i] = state.phead[i]
                        state.phoo[i] = state.phead[i]

            # -- surface energy flux, mid-iteration --------------------------
            stempt = state.toptemp

            fe = surfenergy(state, ii, iter_, sn, kave, sthick, pdens, dmet, lowveg_met=lowveg_met)
            evaprate = fe.evapcm
            rhsurf = fe.rhsurf
            pheatf, pheatg = fe.pheatf, fe.pheatg
            mixrf, mixraf, mixrg, mixra = fe.mixrf, fe.mixraf, fe.mixrgrs, fe.mixra
            lhtf, lheatg = fe.lhtf, fe.lheatg
            sheatf, sheatg = fe.sheatf, fe.sheatg
            dqdtf, dqdtg = fe.dqdtf, fe.dqdtg
            sfac, d1, rpp, cc1 = fe.sfac, fe.d1, fe.rpp, fe.cc1
            isg, isf = fe.isurfg, fe.isurff
            disg, disf, disfg, disgf = fe.disurfg, fe.disurff, fe.disurffg, fe.disurfgf
            lh1 = fe.lh  # discarded -- see module docstring, point 11

            # -- hydraulic parameters -----------------------------------------
            isn = 1
            qbot, save.klhtop, kvhtop, kvttop = flow_param_fn(
                state, isn, qtop, qtopv, simeltsm, zt, delzt
            )

            state.toptemp = state.stt[nnodes]
            state.ftemp = dmet[4]

            toticeo = save.totice
            save.totice = 0.0
            for i in range(1, ntot + 1):
                save.totice = save.totice + state.ice[i]
                if abs(state.stt[i]) <= EPS:
                    state.stt[i] = 1.0

                if i <= nnodes:
                    pres = dmet[11] - 1e-2 * state.phead[i] * dense(state.stt[i], 0.0, d1i) * GRAV
                    if state.node_type[i] != "WA" and state.node_type[i] != "AI":
                        rh = soilhumid(state, i, state.phead[i], state.soil_moist[i], state.stt[i])
                    elif state.node_type[i] == "WA":
                        rh = 1.0
                        pres = dmet[11]
                    elif state.node_type[i] == "AI":
                        rh = 1e-2 * dmet[5]
                        pres = dmet[11]

                    m = sp_humid(d0i, pres, state.stt[i], rh, state.phead[i])
                    t0[1], t0[2], t0[3], t0[4] = m.mixr, m.dmrdt, m.vpress, m.wetbulb
                    rhov[i], rhoda[i] = m.rhov, m.rhoda
                    t1 = m.vpsat
                    thvc[i], dthvdt[i], dthvdh[i] = m.thvc, m.dthvdt, m.dthvdh

                    if state.ntype[i] == 27:
                        state.wvc[i] = rhov[i] / rhoda[i]
                elif i > ntemp:
                    state.stt[i] = dmet[4]
                    state.ice[i] = 0.0
                    state.wvc[i] = 0.0
                    pres = dmet[11]
                    rh = 1e-2 * dmet[5]
                    vp = vap_press(state, i, rh, pres)

                    state.phead[i] = -vp / (RV * state.stt[i] * GRAV) + qtopv * state.deltat_fasst
                    if state.hpond > EPS:
                        state.phead[i] = 0.0
                    if state.phead[i] > 0.0:
                        state.phead[i] = 0.0
                    state.phead[i] = _round20(state.phead[i])
                    if i == ntot:
                        hatm = state.phead[i]

                elif nnodes < i <= ntemp:
                    if icase == 1:
                        state.toptemp = state.stt[ntemp]
                        rh = 1.0
                        vp = vap_press(state, i, rh, pres)
                        rhoda[i] = (pres - vp) / (RD * state.stt[i])
                        rhov[i] = vp / (RV * state.stt[i])
                        state.phead[i] = 0.0
                    elif icase == 2:
                        state.ftemp = state.stt[ntemp]

                        pres = dmet[11]
                        rh = 1e-2 * dmet[5]
                        vp = vap_press(state, i, rh, pres)

                        state.phead[i] = -vp / (RV * state.stt[i] * GRAV) + qtopv * state.deltat_fasst
                        if state.hpond > EPS:
                            state.phead[i] = 0.0
                        if state.phead[i] > 0.0:
                            state.phead[i] = 0.0
                        state.phead[i] = _round20(state.phead[i])
                    elif icase == 3:
                        if state.node_type[i] == "VG":
                            state.ftemp = state.stt[i]

                            pres = dmet[11]
                            rh = 1e-2 * dmet[5]
                            vp = vap_press(state, i, rh, pres)

                            state.phead[i] = -vp / (RV * state.stt[i] * GRAV) + qtopv * state.deltat_fasst
                            if state.hpond > EPS:
                                state.phead[i] = 0.0
                            if state.phead[i] > 0.0:
                                state.phead[i] = 0.0
                            state.phead[i] = _round20(state.phead[i])
                        elif state.node_type[i] == "MX":
                            state.toptemp = state.stt[i]

                            rh = 1.0
                            vp = vap_press(state, i, rh, pres)
                            rhoda[i] = (pres - vp) / (RD * state.stt[i])
                            rhov[i] = vp / (RV * state.stt[i])
                            state.phead[i] = 0.0
                    elif icase == 4:
                        state.toptemp = state.stt[ntemp]
                        if state.node_type[i] == "MX":
                            state.ftemp = state.stt[i]

                        rh = 1.0
                        vp = vap_press(state, i, rh, pres)
                        rhoda[i] = (pres - vp) / (RD * state.stt[i])
                        rhov[i] = vp / (RV * state.stt[i])
                        state.phead[i] = 0.0

            # -- thermal parameters --------------------------------------------
            rhotot = th_param_fn(state, ntot, pdens, pdensnew, rhov, rhoda)

            kave = 5e-1 * (state.grthcond[sn] + state.grthcond[sn - 1])
            sthick = _trunc10(zt[sn] - zt[sn - 1])

            iter_ = iter_ + 1
            if iter_ >= _MAXITER:
                iflag = 1
                break
            elif iter_ >= 3:
                if abs(mat + sert) < 1.0:
                    iflag = 2
                    break
                elif abs(mat1 - sert) < 5.0 and (mat1 * mat < EPS and mat * sert < EPS):
                    iflag = 3
                    break
                elif abs(mat1 + sert) < 5.0 and (mat1 * mat < EPS and mat * sert < EPS):
                    iflag = 4
                    break
                elif abs(mat - sert) < 1.0:
                    iflag = 5
                    break
                elif abs(sert) >= abs(mat) and (
                    (save.totice <= 0.0 and toticeo > EPS) or (save.totice > 0.0 and toticeo <= EPS)
                ):
                    iflag = 6
                    for i in range(1, ntot + 1):
                        state.stt[i] = save.x[1, i]
                        state.wvc[i] = save.x[2, i]
                        state.ice[i] = save.x[3, i]
                        state.soil_moist[i] = save.x[4, i]
                        state.phead[i] = save.x[5, i]
                        state.sink[i] = save.x[6, i]
                        state.sinkr[i] = save.x[7, i]
                        state.source[i] = save.x[8, i]
                        state.flowu[i] = save.x[9, i]
                        state.flowl[i] = save.x[14, i]
                        state.fv1[i] = save.x[10, i]
                        state.vin[i] = save.x[11, i]
                        state.grthcond[i] = save.x[12, i]
                        state.grspheat[i] = save.x[13, i]
                    isg = save.x1[1]
                    isf = save.x1[2]
                    state.sigfl = save.x1[3]
                    mixrg = save.x1[4]
                    break
                elif (abs(mat3 - sert1) < 3.0 and abs(mat4 - mat2) < 3.0) and mat4 * mat3 * mat2 * sert1 > EPS:
                    if lflag == 0:
                        iflag = 7
                        break
        # END OF INNER ITERATION LOOP -----------------------------------------

        # -- soil profile final check for this sub-step (lines 1518-1564) -----
        if ii == state.step:
            for i in range(nnodes + 1, ntemp + 1):
                if (state.node_type[i] == "HM" or state.node_type[i] == "MX") and state.hm > EPS:
                    if state.stt[i] > TREF:
                        state.stt[i] = TREF
                    if told[i] > TREF:
                        told[i] = TREF
                    shs = sheatg
                    lhs = lheatg
                    mixas = mixraf
                    mixgs = mixrg

            for i in range(1, nnodes + 1):
                if state.node_type[i] == "SN" and state.stt[i] > TREF:
                    state.stt[i] = TREF
                if state.node_type[i] == "SN" and told[i] > TREF:
                    told[i] = TREF

            state.toptemp = state.stt[sn]
            temp2 = state.stt[sn - 1]
            stempt = state.toptemp

            fe = surfenergy(state, ii, iter_, sn, kave, sthick, pdens, dmet, lowveg_met=lowveg_met)
            evaprate = fe.evapcm
            rhsurf = fe.rhsurf
            pheatf, pheatg = fe.pheatf, fe.pheatg
            mixrf, mixraf, mixrg, mixra = fe.mixrf, fe.mixraf, fe.mixrgrs, fe.mixra
            lhtf, lheatg = fe.lhtf, fe.lheatg
            sheatf, sheatg = fe.sheatf, fe.sheatg
            dqdtf, dqdtg = fe.dqdtf, fe.dqdtg
            sfac, d1, rpp, cc1 = fe.sfac, fe.d1, fe.rpp, fe.cc1
            isg, isf = fe.isurfg, fe.isurff
            disg, disf, disfg, disgf = fe.disurfg, fe.disurff, fe.disurffg, fe.disurfgf
            state.lhes[iw] = fe.lh  # ONLY this call retains lh -- point 11

            if state.node_type[nnodes] == "WA":
                state.hpond = 0.0
                overland = 0.0
            state.hpond = _round20(state.hpond)

            rhotot = th_param_fn(state, ntot, pdens, pdensnew, rhov, rhoda)

            kave = 5e-1 * (state.grthcond[sn] + state.grthcond[sn - 1])
            if state.hm > EPS:
                kave = state.grthcond[sn]
            sthick = _trunc10(zt[sn] - zt[sn - 1])

            isn = 1
            qbot, save.klhtop, kvhtop, kvttop = flow_param_fn(
                state, isn, qtop, qtopv, simeltsm, zt, delzt
            )

        for i in range(1, ntot + 1):
            y[1, i] = state.stt[i]
            y[2, i] = state.wvc[i]
            y[3, i] = state.ice[i]
            y[4, i] = state.soil_moist[i]
            y[5, i] = state.phead[i]
            y[6, i] = state.sink[i]
            y[7, i] = state.sinkr[i]
            y[8, i] = state.source[i]
            y[9, i] = state.flowu[i]
            y[14, i] = state.flowl[i]
            y[10, i] = state.fv1[i]
            y[11, i] = state.vin[i]
            y[12, i] = state.grthcond[i]
            y[13, i] = state.grspheat[i]

        y1[1] = isg
        y1[2] = isf
        y1[3] = state.sigfl
        y1[4] = mixrg
    # END OF MAIN LOOP -------------------------------------------------------
    # *************************************************************************

    # -- final freezing check (lines 1619-1644) -------------------------------
    iflag = 0
    for i in range(1, ntemp + 1):
        save.stcalc[i] = state.stt[i]

        if state.ntype[i] < 20:
            if state.soil_moist[i] > 1e-7 and state.stt[i] - TREF < 0.0:
                state.ice[i] = state.bftm[i] / dense(state.stt[i], 0.0, d2i)
                state.soil_moist[i] = 0.0
                state.phead[i] = head(state, i, state.soil_moist[i])
            elif state.ice[i] > 1e-7 and state.stt[i] - TREF >= 0.0:
                state.soil_moist[i] = state.bftm[i] / dense(state.stt[i], 0.0, d1i)
                state.phead[i] = head(state, i, state.soil_moist[i])
                state.ice[i] = 0.0

        if (i <= nnodes and (state.ice[i] > EPS or state.ioo[i] > EPS)) and (
            state.node_type[i] != "WA" and state.node_type[i] != "AI"
        ):
            iflag = 1

    state.toptemp = state.stt[sn]
    temp2 = state.stt[sn - 1]
    stempt = state.toptemp

    # -- adjust node thickness for freezing/thawing (lines 1650-1705) --------
    if iflag == 1:
        iceflag = 0
        for i in range(nnodes, 0, -1):
            if (state.ice[i] > EPS or state.ioo[i] > EPS) and (
                state.node_type[i] != "WA" and state.node_type[i] != "AI"
            ):
                # Fortran REUSES the single subroutine-scope `d1` local here
                # (line 1656) -- the SAME variable surfenergy's latent-term
                # denominator was stored into earlier. This deliberately
                # overwrites it with an unrelated ice-expansion ratio; the
                # leftover value (whichever assignment ran last in this
                # `for i` loop, which counts down to i=1) is what flows into
                # `d1` at the final sflux calls below (lines 1739/1749/1769),
                # NOT surfenergy's value, whenever iflag == 1 here. Faithful
                # to the source, not "fixed" -- do not give this its own
                # variable.
                d1 = (dense(TREF, 0.0, d1i) - dense(state.stt[i], 0.0, d2i)) / dense(
                    state.stt[i], 0.0, d2i
                )
                if i == nnodes:
                    state.delzs[i] = state.delzs[i] + d1 * (state.ice[i] - state.ioo[i]) * state.delzs[
                        i
                    ] / state.nsoilp[i, 2]
                    if state.delzs[i] < state.delzsi[i]:
                        state.delzs[i] = state.delzsi[i]
                    if state.ice[i] <= EPS and state.delzs[i] > state.delzsi[i]:
                        state.delzs[i] = state.delzsi[i]
                    state.delzs[i] = _round20(state.delzs[i])

                    sumthick = state.delzs[i]
                else:
                    state.delzs[i] = state.delzs[i] + d1 * (state.ice[i] - state.ioo[i]) * state.delzs[
                        i
                    ] / state.nsoilp[i, 2]
                    if state.delzs[i] < state.delzsi[i]:
                        state.delzs[i] = state.delzsi[i]
                    if state.ice[i] <= EPS and state.delzs[i] > state.delzsi[i]:
                        state.delzs[i] = state.delzsi[i]
                    state.delzs[i] = _round20(state.delzs[i])

                    state.nz[i] = sumthick + 5e-1 * (state.nz[i + 1] - state.nz[i])
                    state.nz[i] = _round20(state.nz[i])
                    state.nz[i] = state.elev - state.nz[i]
                    state.nz[i] = _round20(state.nz[i])

                    sumthick = sumthick + state.delzs[i]
                    sumthick = _round20(sumthick)
            else:
                if i == nnodes:
                    sumthick = state.delzs[i]
                else:
                    state.nz[i] = sumthick + 5e-1 * (state.nz[i + 1] - state.nz[i])
                    state.nz[i] = _round20(state.nz[i])
                    state.nz[i] = state.elev - state.nz[i]
                    state.nz[i] = _round20(state.nz[i])

                    sumthick = sumthick + state.delzs[i]
                    sumthick = _round20(sumthick)
            if abs(state.delzs[i] - state.delzsi[i]) > EPS:
                iceflag = 1

        if iceflag == 0:
            for i in range(1, nnodes + 1):
                state.delzs[i] = state.delzsi[i]
                state.nz[i] = state.nzi[i]

    for i in range(1, ntot + 1):
        state.too[i] = told[i]
        state.woo[i] = wvco[i]
        state.ioo[i] = iceo[i]
        state.phoo[i] = ph_old[i]
        state.smoo[i] = sm_old[i]

    state.storll = state.stll
    state.storls = state.stls
    save.airo = dmet[4]

    if state.newsd > EPS:
        pdensnew = pdens
    else:
        pdensnew = 0.0

    # -- energy closure and outputs (lines 1725-1858) -------------------------
    taf = (
        (1.0 - 7e-1 * state.sigfl) * state.dmet1[iw, 4]
        + 6e-1 * state.sigfl * state.ftemp
        + 1e-1 * state.sigfl * state.toptemp
    )

    if state.veg_flagl != 0 and (state.sigfl > EPS and state.hfol_tot > 0.0):
        t1 = 0.0
        t2 = 1.0
        t0[1] = 0.0
        t0[2] = 0.0
        lh = lhtf
        sh = sheatf
        ph = pheatf

        ff = sflux(
            state, int(d1i), state.step, sn,
            state.dmet1[iw, 8], state.dmet1[iw, 1], state.met[iw, MetCol.IRUP],
            state.dmet1[iw, 2], t2, state.ftemp, cc1, taf, ph, sh, lh, mixraf,
            mixrf, KVEG, state.toptemp, state.hfol_tot, rpp, rhsurf, t1,
            t0[1], t0[2], dqdtf, dqdtg, d1,
        )
        fsdown, fsup = ff.radswd, ff.radswu
        firdown, firup = ff.radlwd, ff.radlwu
        pheatf, sheatf = ff.phw, ff.shw
        state.lheatf[iw] = ff.lhw
        fcheat, fcheat1 = ff.chw, ff.chw1
        fcheck = ff.net
        dnet1, dnet2 = ff.dnet1, ff.dnet2

    # -- check for closure; evap_heat should be 0.0 unless snow is present ---
    ff = sflux(
        state, int(d0i), state.step, sn,
        state.dmet1[iw, 8], state.dmet1[iw, 1], state.met[iw, MetCol.IRUP],
        state.dmet1[iw, 2], sfac, stempt, cc1, taf, pheatg, sheatg, lheatg,
        mixraf, mixrg, kave, temp2, sthick, rhsurf, rpp, state.ioo[sn],
        state.woo[sn], state.smoo[sn], dqdtg, dqdtf, d1,
    )
    state.sdown[iw], state.sup[iw] = ff.radswd, ff.radswu
    state.irdown[iw], state.irup[iw] = ff.radlwd, ff.radlwu
    state.pheat1[iw], state.sheat[iw], state.lheat[iw] = ff.phw, ff.shw, ff.lhw
    state.cheat[iw], state.cheat1[iw] = ff.chw, ff.chw1
    state.evap_heat[iw] = ff.net
    dnet1, dnet2 = ff.dnet1, ff.dnet2

    if iw + 1 < state.iend:
        # Fortran uses ANINT (round) here, unlike the AINT (truncate) used by
        # the delvar/delvar1 tolerance checks above -- `_anint(x, 5)` already
        # implements `anint(x*1d5)*1d-5`, no separate *1e5/*1e-5 needed.
        if (
            _anint(abs(state.dmet1[iw, 8] - state.mflag), 5) <= EPS
            and _anint(abs(state.dmet1[iw + 1, 8] - state.mflag), 5) > EPS
        ):
            state.dmet1[iw, 8] = state.albedo_fasst * state.dmet1[iw, 1]

    # -- melt energy -----------------------------------------------------------
    if state.hm > EPS or state.node_type[nnodes] == "SN":
        if state.tmelt[iw] > TREF:
            taf = (
                (1.0 - 7e-1 * state.sigfl) * state.dmet1[iw, 4]
                + 6e-1 * state.sigfl * state.ftemp
                + 1e-1 * state.sigfl * state.tmelt[iw]
            )

            ff = sflux(
                state, int(d0i), state.step, sn,
                state.dmet1[iw, 8], state.dmet1[iw, 1], state.met[iw, MetCol.IRUP],
                state.dmet1[iw, 2], sfac, state.tmelt[iw], cc1, taf, pheatg,
                shs, lhs, mixas, mixgs, kave, temp2, sthick, rhsurf, rpp,
                state.ioo[ntemp], state.woo[ntemp], state.smoo[ntemp], dqdtg,
                dqdtf, d1,
            )
            swd, swu = ff.radswd, ff.radswu
            ird, iru = ff.radlwd, ff.radlwu
            pht1, sh, lh, ch, ch1 = ff.phw, ff.shw, ff.lhw, ff.chw, ff.chw1
            eh = ff.net
            dnet1, dnet2 = ff.dnet1, ff.dnet2
            state.melt[iw] = eh
        else:
            state.melt[iw] = state.evap_heat[iw]
    else:
        state.melt[iw] = 0.0
    state.melt[iw] = state.melt[iw] * float(state.stepi)
    state.melt[iw] = _round10(state.melt[iw])

    overland = _round10(overland)
    state.hpond = _round10(state.hpond)
    state.ponding[iw] = state.hpond
    save.ntempo = state.ntemp

    sumrunoff = 0.0
    for i in range(ntemp, 0, -1):
        sumrunoff = sumrunoff + runoff[i]
    sumrunoff = _round20(sumrunoff)

    if precip > EPS and state.hm <= EPS:
        save.infl_cum = save.infl_cum + qtop * state.deltat_fasst
        save.tcum = save.tcum + state.timstep * 3.6e3
    elif qtop > EPS and state.hm > EPS:
        save.infl_cum = save.infl_cum + qtop * state.deltat_fasst
        save.tcum = save.tcum + state.timstep * 3.6e3
    else:
        save.infl_cum = 0.0
        save.tcum = 0.0
    if save.infl_cum < EPS:
        save.infl_cum = 0.0

    # -- moisture self-check (debug information; write() calls all commented
    # out in the Fortran source -- see module docstring, point 12) ----------
    summoist = 0.0
    sumin = 0.0
    for i in range(1, ntot + 1):
        temp[i] = 0.0
        if i <= nnodes:
            rhow = dense(state.stt[i], 0.0, d1i)
            if state.ice[i] > EPS:
                rhoi = dense(state.stt[i], 0.0, d2i)
                temp[i] = state.soil_moist[i] + state.ice[i] * rhoi / rhow
            else:
                temp[i] = state.soil_moist[i]

            if state.wvc[i] > EPS:
                temp[i] = temp[i] + state.wvc[i] * rhov[i] / rhow

            summoist = summoist + temp[i] * state.delzs[i]
            sumin = sumin + state.vin[i]
    # inverted-exponent idiom, structurally zero -- see module docstring,
    # point 6. Replicated literally, not "corrected".
    state.tot_moist[iw] = _round10_inv(summoist)

    if (state.hpond > EPS or hpondi > EPS) or overland > EPS:
        extra = extra + (max(0.0, state.hpond - hpondi) + overland) * f2
    extra = _round10_inv(extra)

    return NewProfileResult(code5=code5, first_time=first_time, phie=phie, fsup=fsup, firup=firup)
