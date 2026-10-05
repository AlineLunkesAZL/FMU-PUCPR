"""Mutable model state, transcribed from fasst_global.F90.

Fortran source: fasst_global.F90, lines 50-191 (everything that is not a
``parameter``). The ``save::`` statements scattered across the Fortran module
(8 of them, covering practically every non-allocatable scalar and array) are
redundant in Fortran 90/2003 -- module variables already have implicit SAVE
-- but they are informative: they signal that the original author knew that
100% of this state survives across calls.

Architectural note (docs_transcricao/2026-08-24-doc-transcricao-fasst_global.docx,
"Papel arquitetural" / point 1 / PENDENCIAS): FASST does not pass state
through arguments -- it passes state through module-shared memory. This is
the direct modernization of a Fortran 77 COMMON block, and reproducing it
literally as Python module-level globals would carry the same bug class
into Python: state leaking between successive runs in multi-point mode
(which is exactly why module_zerovars.F90's zero_parameters /
zero_mxl_params exist upstream). Instead, every symbol below is a field of
:class:`FasstState`, instantiated once per run. Passing a `FasstState`
around replaces the old COMMON-by-module-inclusion; call sites en masse
still need their ~30 signatures rewritten to accept it, which is out of
scope for this file (fasst_global.F90 itself has zero executable lines).

Indexing: 1-based Fortran arrays are represented with one unused "phantom"
slot at index 0 (see :mod:`fasst.indices` docstring and doc point 2b) so
that ``arr(j)`` in Fortran reads as ``arr[j]`` in Python without an offset.
The sole exception is ``avect``, which is declared ``avect(0:4, nclayers)``
in Fortran -- i.e. it is *already* 0-based in the original source, so it is
given no phantom offset here (see doc point 8).

Fixed-width character arrays (``stype``, ``node_type``, ``sclass``,
``nclass``, and the scalar ``meltfl``) are transcribed as NumPy fixed-width
unicode arrays / plain str. Fortran string comparison pads with trailing
spaces before comparing (``'SN' == 'SN  '`` is true); Python does not. Any
code ported from a Fortran callsite that compares these fields MUST use
``.strip()`` (or an Enum) on both sides -- never compare raw literals (see
doc point 10).

Case sensitivity: Fortran is case-insensitive, Python is not. The Fortran
source itself is inconsistent (e.g. this module declares ``Tref`` while
fasst_main.F90 line 455 writes ``tref``, and this module declares ``Zd``
next to ``zh``). All fields here use canonical lowercase names; audit each
call site being transcribed for spelling variants before wiring it up (see
doc point 9).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .constants import EXTRAN, MAXL, MAXN, MAXP, MT, NCLAYERS, SPFLAG


# --- phantom-slot array builders (index 0 unused, see module docstring) ---
def _ints(n: int) -> np.ndarray:
    return np.zeros(n + 1, dtype=np.int64)


def _reals(n: int) -> np.ndarray:
    return np.zeros(n + 1, dtype=np.float64)


def _reals2(n: int, m: int) -> np.ndarray:
    return np.zeros((n + 1, m + 1), dtype=np.float64)


def _chars(n: int, width: int) -> np.ndarray:
    return np.full(n + 1, "", dtype=f"<U{width}")


@dataclass
class FasstState:
    """One FASST run's worth of module-global state (fasst_global.F90).

    Instantiate one `FasstState` per run/grid-point instead of relying on
    module-level globals, and call :meth:`reset` between successive points
    in multi-point mode. See module docstring for the rationale.
    """

    # -- met file pointer information (fasst_global.F90 lines 50-59) -------
    maxlines: int = 0
    ncols: int = 0
    istart: int = 0
    iend: int = 0
    met_count: int = 0
    single_multi_flag: int = 0
    water_flag: int = 0
    moverlapr: int = 0
    mflag: float = 0.0          # missing-data sentinel read from the met file header (typically -999)
    timeoffset: float = 0.0
    elev: float = 0.0
    lat: float = 0.0
    mlong: float = 0.0
    timstep: float = 0.0
    slope_fasst: float = 0.0
    albedo_fasst: float = 0.0
    emis: float = 0.0
    iheight: float = 0.0
    aspect: float = 0.0
    sloper: float = 0.0

    # -- initial conditions, soil structure (lines 62-82) -------------------
    nlayers: int = 0
    nnodes: int = 0
    refn: int = 0
    soiltype: np.ndarray = field(default_factory=lambda: _ints(MAXL))
    ntype: np.ndarray = field(default_factory=lambda: _ints(MAXN))
    icourse: np.ndarray = field(default_factory=lambda: _ints(MAXN))
    tm: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    zti: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    zm: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    sm: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    soil_moist: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    nz: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    delzs: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    grthcond: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    grspheat: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    stt: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    ice: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    phead: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    wvc: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    sgralbedo: float = 0.0
    sgremis: float = 0.0
    gwl: float = 0.0
    lthick: np.ndarray = field(default_factory=lambda: _reals(MAXL))
    rho_fac: np.ndarray = field(default_factory=lambda: _reals(MAXL))
    rough: float = 0.0
    isoilp: np.ndarray = field(default_factory=lambda: _reals2(MAXL, MAXP))
    soilp: np.ndarray = field(default_factory=lambda: _reals2(MAXL, MAXP))
    nsoilp: np.ndarray = field(default_factory=lambda: _reals2(MAXN, MAXP))
    nzi: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    delzsi: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    # fixed-width strings -- see module docstring on Fortran padding semantics
    stype: np.ndarray = field(default_factory=lambda: _chars(MAXL, 2))              # 'SN','WA','AI','SM'
    node_type: np.ndarray = field(default_factory=lambda: _chars(MAXN + EXTRAN, 2))  # 'SN','WA','AI','SM'
    sclass: np.ndarray = field(default_factory=lambda: _chars(MAXL, 4))             # 'USCS','USDA'
    nclass: np.ndarray = field(default_factory=lambda: _chars(MAXN + EXTRAN, 4))    # 'USCS','USDA'

    # -- snow and ice parameters (lines 85-88) -------------------------------
    sn_istat: np.ndarray = field(default_factory=lambda: _ints(2))
    hsaccum: float = 0.0
    hi: float = 0.0
    dsnow: float = 0.0
    refreeze: float = 0.0
    newsd: float = 0.0
    atopf: float = 0.0
    km: np.ndarray = field(default_factory=lambda: _reals(EXTRAN))
    sphm: np.ndarray = field(default_factory=lambda: _reals(EXTRAN))
    refreezei: float = 0.0
    sn_stat: np.ndarray = field(default_factory=lambda: _reals(15))  # positions undocumented, see doc point 8
    vsmelt: float = 0.0
    vimelt: float = 0.0
    iswe: float = 0.0
    hm: float = 0.0

    # -- foliage variables (lines 113-125) -----------------------------------
    vegl_type: int = 0
    vegh_type: int = 0
    veg_flagl: int = 0
    veg_flagh: int = 0
    iseason: int = 0
    biome_source: int = 0
    new_vtl: int = 0
    new_vth: int = 0
    isigfl: float = 0.0
    iepf: float = 0.0
    ifola: float = 0.0
    ihfol: float = 0.0
    isigfh: float = 0.0
    izh: float = 0.0
    iheightn: float = 0.0
    zh: float = 0.0
    zd: float = 0.0             # Fortran: Zd (mixed case -- see module docstring, doc point 9)
    albf: float = 0.0
    ftemp: float = 0.0
    z0l: float = 0.0
    chnf: float = 0.0
    trmlm: float = 0.0
    trmhm: float = 0.0
    rsl: float = 0.0
    storll: float = 0.0
    storls: float = 0.0
    sqrt_chnf: float = 0.0
    sigfh: float = 0.0
    uaf: float = 0.0
    hfol_tot: float = 0.0
    lail: float = 0.0
    sigfl: float = 0.0
    state: float = 0.0          # Fortran field literally named `state` -- unrelated to the FasstState class itself
    fola: float = 0.0
    ilail: float = 0.0
    epf: float = 0.0
    hfol: float = 0.0
    stll: float = 0.0
    stls: float = 0.0
    veg_prp: np.ndarray = field(default_factory=lambda: _reals2(18, 17))  # 18 = number of FASST vegetation types (unnamed magic number in Fortran)
    rk: np.ndarray = field(default_factory=lambda: _reals2(18, MAXN))
    frl: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    frh: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    sinkr: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    ifoliage_type: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    ilai: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    iclump: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    irho: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    laif: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    itau: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    ialp: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    ieps: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    dzveg: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    storcl: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    storcs: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    idzveg: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS))
    stor: np.ndarray = field(default_factory=lambda: _reals(NCLAYERS + NCLAYERS))
    # avect(0:4, nclayers): already 0-based in Fortran -- NO phantom offset (see module docstring)
    avect: np.ndarray = field(default_factory=lambda: np.zeros((5, NCLAYERS), dtype=np.float64))

    # -- calculation variables (lines 153-158) -------------------------------
    iw: int = 0
    stepi: int = 0
    error_code: int = 0
    error_type: int = 0
    infer_test: int = 0
    oldpos: int = 0
    freq_id: int = 0
    vitd_index: int = 0
    ecount: int = 0
    mstflag: int = 0
    step: int = 0
    icase: int = 0
    icaseo: int = 0
    dstart: int = 0             # partial-maintenance vestige: also appears commented as parameter=14 (doc point 11)
    mpos: int = 0
    deltat_fasst: float = 0.0
    deltati: float = 0.0
    pheadmin: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    hpond: float = 0.0
    toptemp: float = 0.0
    ptemp: float = 0.0
    mgap: float = 0.0
    isurfoldg: float = 0.0
    isurfoldf: float = 0.0
    meltfl: str = ""             # character(len=1) -- see module docstring on padded comparisons

    # -- profile parameters (lines 166-173) ----------------------------------
    ntot: int = 0
    ntemp: int = 0
    khu: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    khl: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    vin: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    flowu: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    flowl: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    fv1: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    source: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    sink: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    too: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    smoo: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    woo: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    ioo: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    phoo: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))
    bftm: np.ndarray = field(default_factory=lambda: _reals(MAXN))
    dsmdh: np.ndarray = field(default_factory=lambda: _reals(MAXN + EXTRAN))

    # -- allocatable arrays (lines 179-191) ----------------------------------
    # Unlike every block above, these carry no `save::` in the Fortran source
    # (doc point 11) and are sized at runtime from the met file line count
    # (met_count), not from a fixed `parameter`. They are intentionally left
    # unallocated (None) here -- allocation is the responsibility of whatever
    # transcribes the driver code that currently does `allocate(...)`
    # (fasst_driver.F90), not of this module.
    slushy: Optional[np.ndarray] = None
    dmet1: Optional[np.ndarray] = None          # dmet1(:, DerivedMet.*), 13 cols -- see fasst.indices
    sdens: Optional[np.ndarray] = None
    canopy_temp: Optional[np.ndarray] = None
    airt: Optional[np.ndarray] = None
    ft: Optional[np.ndarray] = None
    tt: Optional[np.ndarray] = None
    surfice: Optional[np.ndarray] = None
    surfmoist: Optional[np.ndarray] = None
    surficep: Optional[np.ndarray] = None
    surfmoistp: Optional[np.ndarray] = None
    surfci: Optional[np.ndarray] = None
    surfrci: Optional[np.ndarray] = None
    surfcbr: Optional[np.ndarray] = None
    lhes: Optional[np.ndarray] = None
    frthick: Optional[np.ndarray] = None
    twthick: Optional[np.ndarray] = None
    cheat: Optional[np.ndarray] = None
    met: Optional[np.ndarray] = None            # met(:, MetCol.*) -- see fasst.indices
    cheat1: Optional[np.ndarray] = None
    sdown: Optional[np.ndarray] = None
    sup: Optional[np.ndarray] = None
    irdown: Optional[np.ndarray] = None
    irup: Optional[np.ndarray] = None
    pheat1: Optional[np.ndarray] = None
    lheat: Optional[np.ndarray] = None
    evap_heat: Optional[np.ndarray] = None
    sheat: Optional[np.ndarray] = None
    melt: Optional[np.ndarray] = None
    tmelt: Optional[np.ndarray] = None
    surfd: Optional[np.ndarray] = None
    surfemis: Optional[np.ndarray] = None
    surfemisf: Optional[np.ndarray] = None
    surfemisc: Optional[np.ndarray] = None
    lheatf: Optional[np.ndarray] = None
    cevap: Optional[np.ndarray] = None
    ponding: Optional[np.ndarray] = None
    tot_moist: Optional[np.ndarray] = None
    tot_thick: Optional[np.ndarray] = None
    sstate: Optional[np.ndarray] = None          # character(len=1), allocatable

    def reset(self) -> None:
        """Reinitialize every fixed-size field to the state that
        ``zero_parameters`` (module_zerovars.F90) leaves them in, in place.

        Cross-check performed 2026-09-03 (Passo 3,
        docs_transcricao/2026-09-03-testes-sem-fortran.md): most of
        ``zero_parameters`` is ``x = 0d0`` / ``x = 0`` / ``x = ' '``, which
        equals the dataclass defaults. The exceptions -- the values it sets
        to a *non-zero* sentinel -- are listed in
        :data:`_ZERO_PARAMETERS_SENTINELS` and applied here after the
        defaults are restored:

          * ``freq_id = 1``
          * ``iswe = ftemp = isigfl = iepf = ifola = ihfol = SPFLAG`` (999.0)
          * ``albf = 1.0``            (``1d0 - fola`` with ``fola`` just zeroed)
          * ``rho_fac[1..MAXL] = 1.0``
          * ``isoilp[1..MAXL, 1..MAXP] = soilp[...] = SPFLAG``
          * ``ifoliage_type / ilai / iclump / irho / itau / ialp / ieps /
            idzveg [1..NCLAYERS] = SPFLAG``

        Two Fortran zeroing GAPS are deliberately NOT reproduced here (they
        are reported to the author instead, step-10 report; a genuinely
        clean Python reset is the safer default):

          * ``zero_parameters`` clears ``sn_stat(1..13)`` but the array is
            ``sn_stat(15)`` -- elements 14-15 keep stale values in Fortran.
            Here the whole ``sn_stat`` is zeroed.
          * ``zero_mxl_params`` clears ``dmet1(:,1..10)`` but ``dmet1`` has
            13 columns -- 11-13 keep stale values in Fortran. ``dmet1`` is
            allocatable here (left to the driver), so this method does not
            touch it at all.

        The scalar flags ``ttest`` / ``mtest`` and the water params
        ``wtype`` / ``wdepth`` / ``wvel`` that ``zero_parameters`` also
        writes are subroutine ``intent(out)`` arguments in the Fortran, not
        module state, so they have no ``FasstState`` field and are out of
        scope here.

        Allocatable arrays (sized by ``met_count``) are left untouched:
        their lifecycle belongs to the driver code (``zero_mxl_params`` runs
        after ``allocate`` in the Fortran), not to this method.
        """
        fresh = FasstState()
        for name in self.__dataclass_fields__:
            if name not in _ALLOCATABLE_FIELDS:
                setattr(self, name, getattr(fresh, name))

        # non-zero sentinels from zero_parameters (module_zerovars.F90)
        self.freq_id = 1
        for name in _ZERO_PARAMETERS_SPFLAG_SCALARS:
            setattr(self, name, SPFLAG)
        self.albf = 1.0
        self.rho_fac[1:MAXL + 1] = 1.0
        self.isoilp[1:MAXL + 1, 1:MAXP + 1] = SPFLAG
        self.soilp[1:MAXL + 1, 1:MAXP + 1] = SPFLAG
        for name in _ZERO_PARAMETERS_SPFLAG_NCLAYERS:
            getattr(self, name)[1:NCLAYERS + 1] = SPFLAG


# Scalars that zero_parameters (module_zerovars.F90) sets to `spflag` (999.0)
# rather than to 0. `iswe` is written twice there (`0d0` then `spflag`).
_ZERO_PARAMETERS_SPFLAG_SCALARS = (
    "iswe", "isigfl", "iepf", "ifola", "ihfol", "ftemp",
)

# nclayers-sized arrays that zero_parameters fills with `spflag`.
_ZERO_PARAMETERS_SPFLAG_NCLAYERS = (
    "ifoliage_type", "ilai", "iclump", "irho", "itau", "ialp", "ieps", "idzveg",
)


_ALLOCATABLE_FIELDS = frozenset(
    {
        "slushy", "dmet1", "sdens", "canopy_temp", "airt", "ft", "tt",
        "surfice", "surfmoist", "surficep", "surfmoistp", "surfci",
        "surfrci", "surfcbr", "lhes", "frthick", "twthick", "cheat", "met",
        "cheat1", "sdown", "sup", "irdown", "irup", "pheat1", "lheat",
        "evap_heat", "sheat", "melt", "tmelt", "surfd", "surfemis",
        "surfemisf", "surfemisc", "lheatf", "cevap", "ponding",
        "tot_moist", "tot_thick", "sstate",
    }
)
