"""Column-index enums for the FASST data-interface arrays.

Three column layouts, all transcribed straight from the Fortran source --
this module is pure cataloguing of the data interface, it contains no
executable model logic:

* :class:`MetCol`    -- ``met(timestep, :)``, the raw met file
                        (fasst_global.F90 lines 9-45, ``ip_*`` parameters).
* :class:`DerivedMet`-- ``dmet1(timestep, :)``, 13 columns of *derived*
                        meteorology. The Fortran has no ``ip_*``-style
                        constant for these; the mapping lives only in
                        fasst_main.F90 lines 158-191. Catalogued here per
                        the fasst_global transcription doc
                        (PENDENCIAS / "catalogar dmet1"), Passo 2 of
                        docs_transcricao/2026-09-03-*.md.
* :class:`SoilProp`  -- column ``j`` of ``isoilp(layer, j)`` /
                        ``nsoilp(node, j)`` / ``soilp(layer, j)``, the 25
                        soil properties. In the Fortran this is only a
                        comment block (fasst_main.F90 lines 68-94);
                        promoted to an enum here so ``fasst/functions.py``'s
                        bare ``nsoilp[i, 8]`` / ``[i, 24]`` indices can be
                        read against a name.

Kept 1-based deliberately ("phantom slot" strategy -- see
docs_transcricao/2026-08-24-doc-transcricao-fasst_global.docx, point 2b,
option chosen for phase 1). ``met``, ``dmet1`` and ``*soilp`` arrays in
:mod:`fasst.state` are allocated with an unused row/column 0, so a Fortran
index such as ``met(iw, ip_tsol)`` maps directly onto
``met[iw, MetCol.TSOL]`` and ``nsoilp(i, 12)`` onto
``nsoilp[i, SoilProp.VG_M]`` without an off-by-one shift. This preserves
line-by-line correspondence with the Fortran source for the first
transcription pass, at the cost of one wasted slot per array dimension.
A later pass may switch to 0-based, shifted indices (option 2a) once the
whole codebase has been transcribed and cross-checked.
"""

from enum import IntEnum

from .constants import MAXCOL, MAXP


class MetCol(IntEnum):
    """Column layout of ``met(timestep, MetCol.*)``."""

    YEAR = 1     # year
    DOY = 2      # day of year
    HR = 3       # hour (local time)
    MIN = 4      # minute (local time)
    AP = 5       # air pressure (mbar)
    TMP = 6      # air temperature (C)
    RH = 7       # relative humidity (%)
    WS = 8       # wind speed (m/s)
    WDIR = 9     # wind direction (+ clockwise from N)
    PREC = 10    # precipitation rate (mm/hr)
    PT = 11      # precipitation type
    PREC2 = 12   # precipitation snow rate (mm/hr)
    PT2 = 13     # precipitation snow type
    LCD = 14     # low cloud amount (fractional amount)
    LHGT = 15    # low cloud hgt (km)
    LCT = 16     # low cloud type
    MCD = 17     # middle cloud amount (fractional amount)
    MHGT = 18    # middle cloud hgt (km)
    MCT = 19     # middle cloud type
    HCD = 20     # high cloud amount (fractional amount)
    HHGT = 21    # high cloud hgt (km)
    HCT = 22     # high cloud type
    TSOL = 23    # total solar flux (W/m^2)
    DIR = 24     # direct solar flux (W/m^2)
    DIF = 25     # diffuse solar flux (W/m^2)
    UPSOL = 26   # reflected solar flux (W/m^2)
    IR = 27      # downwelling IR flux (W/m^2)
    IRUP = 28    # upwelling IR flux (W/m^2)
    ZEN = 29     # solar zenith angle (degrees)
    AZ = 30      # solar azimuth angle from North clockwise
    SD = 31      # snow depth (m)
    TSOIL = 32   # soil temperature (C)
    HI = 33      # surface ice thickness (m). This is for roads, etc.
    VIS = 34     # visibility (km). This is a pass-through param.
    AER = 35     # aerosol type. This is a pass-through param.


assert len(MetCol) == MAXCOL, "MetCol must cover YEAR..AER with no gaps or duplicates"


class DerivedMet(IntEnum):
    """Column layout of ``dmet1(timestep, DerivedMet.*)`` -- derived meteorology.

    Fortran source: fasst_main.F90 lines 158-191. There is NO ``ip_*``
    parameter for these columns in fasst_global.F90; ``dmet1`` is filled
    element by element at the top of the per-timestep loop in
    ``fasst_main``. Each column is assembled from one raw ``met`` column,
    with the small fix-ups noted below, then all 13 are quantised with
    ``anint(x*1d15)*1d-15`` (fasst_main.F90:189-191).

    Column -> source (fasst_main.F90):
      1 TSOL_NET   met(:,ip_tsol); forced to 0 if < 1e-2 W/m^2   (:158-159)
      2 IR_DOWN    met(:,ip_ir)                                    (:160)
      3 WIND       met(:,ip_ws)                                    (:161)
      4 TAIR_K     met(:,ip_tmp) + Tref  -- air temperature in K   (:162)
      5 RH_PCT     min(met(:,ip_rh), 100)  -- still a PERCENT here (:164)
      6 PRECIP_M   met(:,ip_prec)*1e-3  -- m per timestep (rain and/or snow),
                   with met(:,ip_prec)/ip_pt zeroed when <= eps    (:166-172)
      7 SNOW_M     met(:,ip_prec2)*1e-3 -- m per timestep (snow only),
                   with met(:,ip_prec2)/ip_pt2 zeroed when <= eps  (:174-179)
      8 UPSOL      met(:,ip_upsol); forced to 0 when |TSOL_NET| <= eps (:181-182)
      9 SOLAR_DIR  met(:,ip_dir)                                   (:183)
     10 SOLAR_DIF  met(:,ip_dif)                                   (:184)
     11 AIR_PRESS  met(:,ip_ap)  -- mbar                           (:185)
     12 IR_UP      met(:,ip_irup)                                  (:186)
     13 TSOIL      met(:,ip_tsoil)                                 (:187)

    KNOWN DEFECT (see Passo 3 / step-10 report): ``zero_mxl_params``
    (module_zerovars.F90:282-284) zeroes only columns 1..10 of ``dmet1``;
    columns 11..13 (AIR_PRESS, IR_UP, TSOIL) are never cleared between
    multi-point runs. Preserve literally when transcribing; do not "fix".

    UNIT TRAP: column 5 (RH_PCT) is a percentage (0..100). ``sp_humid``
    treats its ``rh2`` argument as a fraction (0..1). See Passo 6's
    complementary audit -- a factor-100 error can hide across a whole
    transcription here.
    """

    TSOL_NET = 1
    IR_DOWN = 2
    WIND = 3
    TAIR_K = 4
    RH_PCT = 5
    PRECIP_M = 6
    SNOW_M = 7
    UPSOL = 8
    SOLAR_DIR = 9
    SOLAR_DIF = 10
    AIR_PRESS = 11
    IR_UP = 12
    TSOIL = 13


MAXDMET = 13  # dmet1's active column count (fasst_main.F90 `do i=1,13`)
assert len(DerivedMet) == MAXDMET, "DerivedMet must cover columns 1..13 with no gaps"


class SoilProp(IntEnum):
    """Column ``j`` of ``isoilp(layer,j)`` / ``nsoilp(node,j)`` / ``soilp(layer,j)``.

    Verbatim transcription of the comment block in fasst_main.F90 lines
    68-94 ("SOIL PROPERTIES ... where i=layers, k=nodes, j=soil property"):

       1  bulk dry density (g/cm^3)
       2  porosity
       3  albedo_fasst
       4  emissivity
       5  quartz content
       6  dry thermal conductivity (W/m*K)
       7  saturated hydraulic conductivity (cm/s)
       8  minimum volumetric water content        [sic "minumum" in source]
       9  maximum volumetric water content
      10  van Genuchten's alpha (1/cm)
      11  van Genuchten's n
      12  van Genuchten's m  (m = 1 - 1/n)
      13  specific heat of dry soil (J/kg*K)
      14  organic fraction (vol/vol)
      15  minimum allowed volumetric water content
      16  wilting point volumetric water content
      17  field capacity volumetric water content
      18  % sand
      19  % silt
      20  % clay
      21  % carbon
      22  plastic limit                           [sic "plactic" in source]
      23  percent fines passing #200 sieve
      24  0.999d0*soilp(i,9)  => max water content, capped to avoid
          numerical instability
      25  % gravel

    Cross-check with ``fasst/functions.py`` (bare integer indices today):
    ``head`` reads 8, 9, 10, 11, 12, 15, 24; ``soilhumid`` reads 9;
    ``maxinfiltrate`` reads 7, 8, 9, 10, 11, 12.
    """

    BULK_DRY_DENSITY = 1
    POROSITY = 2
    ALBEDO = 3
    EMISSIVITY = 4
    QUARTZ = 5
    DRY_THCOND = 6
    KSAT = 7
    THETA_MIN = 8          # theta_r in the van Genuchten sense (functions.head)
    THETA_MAX = 9          # theta_s / porosity used by head, soilhumid, maxinfiltrate
    VG_ALPHA = 10
    VG_N = 11
    VG_M = 12
    DRY_SPECIFIC_HEAT = 13
    ORGANIC_FRACTION = 14
    THETA_MIN_ALLOWED = 15
    WILTING_POINT = 16
    FIELD_CAPACITY = 17
    PCT_SAND = 18
    PCT_SILT = 19
    PCT_CLAY = 20
    PCT_CARBON = 21
    PLASTIC_LIMIT = 22
    PCT_FINES_200 = 23
    THETA_MAX_CAPPED = 24  # 0.999 * soilp(i,9); functions.head upper bound
    PCT_GRAVEL = 25


assert len(SoilProp) == MAXP, "SoilProp must cover the 25 columns of maxp"
