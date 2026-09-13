"""Physical and numerical constants, transcribed from fasst_global.F90.

Fortran source: fasst_global.F90, lines 3-6 (kinds), 48-49 and 62-65, 109 and
136 (dimensioning parameters), 90-103 (snow/ice), 109-111 (foliage), 136-151
(calculation constants).

The Fortran kinds `ip`/`dp` (selected_int_kind(8), selected_real_kind(15,307))
are discarded: Python's native `int` is arbitrary precision and NumPy's
float64 already matches `dp`. The commented-out `ip1`/`dp1` (single-precision
kinds) are likewise dropped -- they are vestiges of an abandoned mixed
precision build and were never active.

Every literal below is copied verbatim from the Fortran source, including the
ones that disagree with current CODATA/reference values (sigma, pi, vK). Do
NOT "fix" them to math.pi / scipy.constants: a bit-exact comparison against
the Fortran output is the only reliable way to catch a transcription error,
and it breaks at the first changed digit. Divergences are tracked as separate
technical debt (see docs_transcricao/2026-08-24-doc-transcricao-fasst_global.docx,
point 6) and must be resolved -- if ever -- as an isolated, validated change.
"""

import sys

# --- Dimensioning / sizing parameters --------------------------------------
# fasst_global.F90 lines 48-49, 62-65, 109, 136
MAXCOL = 35      # maximum number of columns of met data
MOVERLAP = 15    # maximum number of lines of metfile overlap
MAXL = 10        # maximum number of soil layers
MAXN = 100       # maximum number of nodes
MAXP = 25        # maximum number of soil parameters
EXTRAN = 3       # extra nodes for snow, vegetation, air
NCLAYERS = 3     # number of canopy layers
MT = 30          # maximum number of soil types

# --- Calculation constants (fasst_global.F90 lines 138-151) ---------------
PI = 3.141592654               # NOT math.pi -- rel. err vs. CODATA ~4.1e-10
SIGMA = 5.669e-8                # Stefan-Boltzmann const. (W/m^2*K^4); CODATA 5.670374419e-8, rel. err ~2.4e-4
SPFLAG = 999.0                  # missing soil, veg parameter flag
ILIM = 0.998                    # closeness of ice to maximum water content
WLIM = 0.998                    # closeness of watervapor to porosity
EPS = sys.float_info.epsilon    # tolerance limit for equality (== epsilon(1d0) in Fortran)
# NOTE: EPS is used throughout the model as "is this exactly zero, up to the
# last bit" (if(dabs(x) <= eps)), NOT as "is this physically negligible" --
# preserve this usage literally when transcribing callers; do not replace
# with a physically-motivated tolerance (see doc point 4).
VK = 4e-1                       # von Karman's constant (unitless); commented alt. value in Fortran: 0.35
TREF = 273.15                   # reference temperature (K)
GRAV = 9.81                     # gravitational acceleration (m/s^2)
RV = 8.3143 / 1.8015e-2         # specific gas constant for water vapor (J/kg*K)
RD = 8.3143 / 2.8964e-2         # specific gas constant for dry air (J/kg*K)
VPSAT0 = 610.78                 # saturation vapor pressure at 273.15 K (Pa)
HT_MIN = 2.5e-1                 # minimum top thickness of snow (m) [SNTHERM = 1.67e-2] -- see doc point 7
HT_MINM = 1e-3                  # absolute minimum snow depth (m) [SNTHERM = 2e-3]

# --- Snow and ice parameters (fasst_global.F90 lines 90-103) --------------
SNALBEDO = 0.8       # new snow albedo (unitless); commented alt. value in Fortran: 0.9
SOALBEDO = 0.5       # old snow albedo (unitless)
SEMIS = 0.97         # snow emissivity (unitless); commented alt. value in Fortran: 0.92
STHDIFF = 2e-07      # snow thermal diffusivity (m^2/s)
STHCOND = 0.3492     # snow thermal conductivity (W/m*K)
SDENSW = 550.0       # wet snow density (kg/m^3)
SDENSD = 50.0        # dry snow density (kg/m^3)
IALBEDO = 0.7        # ice albedo (unitless)
IEMIS = 0.9          # ice emissivity (unitless)
ITHDIFF = 1.167e-06  # ice thermal diffusivity (m^2/s)
ITHCOND = 2.1648     # ice thermal conductivity (W/m*K)
IDENS = 916.5        # ice density (kg/m^3)
LHFUS = 3.335e05     # latent heat of fusion (J/kg)
LHSUB = 2.838e06     # latent heat of sublimation (J/kg = (m/s)^2)

# --- Foliage constants (fasst_global.F90 lines 110-111) --------------------
KVEG = 0.38      # average vegetation thermal conductivity (W/m*K)
SPHVEG = 3.5e3   # average vegetation specific heat (J/kg*K)
