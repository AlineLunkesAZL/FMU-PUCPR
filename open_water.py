"""
open_water.py -- traducao de open_water.F90

Fonte Fortran: open_water.F90 (362 linhas). Duas sub-rotinas:
`open_water` (publica) e `find_index` (auxiliar, "attached to this
subroutine" no cabecalho do original).

`open_water` calcula, para o passo `state.iw`, as condicoes de uma
superficie de agua aberta (lago, rio, canal): temperatura da agua a partir
apenas da latitude e do mes (tabela `scorr`, rotina recebida de Guy Seelye,
AER -- o proprio autor marca como "stop-gap"), variaveis meteorologicas
derivadas (`dmet1`), temperatura/umidade do dossel e da vegetacao baixa
(quando existem), neve nova, crescimento de gelo, modelo de neve,
albedo/emissividade da superficie, e a gravacao do arquivo de reinicio
(`read_old_data`, modo 2) nos ultimos `moverlap` passos.

Assinatura original:
    subroutine open_water(wtype,sprint,wvel,wdepth,sdensi,first_loop,
                          code2,code4,oldsd)
    integer(ip),intent(in):: wtype,sprint
    real(dp),intent(in):: wvel,wdepth,sdensi
    integer(ip),intent(inout):: first_loop,code2,code4
    real(dp),intent(inout):: oldsd

Arquitetura
-----------
`open_water(state, ow, wtype, sprint, wvel, wdepth, sdensi, first_loop,
code2, code4, oldsd)` recebe `state: FasstState` primeiro e um
`OpenWaterState` (`ow`) em seguida. Os argumentos `intent(inout)` voltam
numa tupla: `first_loop, code2, code4, oldsd = open_water(...)`.

`OpenWaterState` reune as variaveis `save` que persistem entre chamadas:
  * `first_snow` (unica `save::` explicita de `open_water`);
  * os objetos de estado salvo dos modulos chamados, ja definidos no
    projeto: `CanopySavedState` (module_canopy), `AlbedoEmisState`
    (albedo_emis, usado por canopy_met), `LowVegState` (module_lowveg) e
    `SnowState` (module_snow).
Instancie um `OpenWaterState` por execucao e reutilize-o em todas as
chamadas -- o driver cria UM antes do STEP 4, nunca um por ponto, porque
um `save` do Fortran dura o programa inteiro.

No Fortran os `save` de `canopy_met`, `albedo_emis`, `lowveg_met` e
`module_snow` pertencem a essas rotinas, e sao os MESMOS para
`open_water` e para `fasst_main`: o driver cria o `FasstMainState` com
`canopy`/`alb`/`lowveg`/`snow` apontando para os deste objeto.

Constantes (`EPS`, `PI`, `SPFLAG`, `TREF`, `SIGMA`, `IALBEDO`, `IEMIS`,
`IDENS`, `ITHCOND`, `STHCOND`, `LHFUS`, `MOVERLAP`, `NCLAYERS`) vem de
`fasst.constants`; colunas de `met` sao `MetCol.*`. As colunas de `dmet1`
(13) nao tem enum e sao usadas por numero, como nas demais traducoes.
Indexacao 1-based identica ao Fortran (posicao 0 sem uso), inclusive
`scorr(1..12)`.

Dependencias externas (todas fora deste arquivo)
-------------------------------------------------
  canopy_met        fasst.module_canopy   (real)
  low_veg_prop      fasst.module_lowveg   (real)
  lowveg_met        fasst.module_lowveg   (real)
  sp_humid, dense   fasst.sp_humid / fasst.functions (reais; so a
                    interface usada em module_canopy/module_lowveg)
  read_old_data     fasst.read_old_data   (real)
  snow_model        fasst.snow            (real)

`snow_model(state, snow, sprint, water_flag, first_snow, sdensi, sdf)
-> (first_snow, sdf, phie1)` -- o `call snow_model(sprint,water_flag,
first_snow,sdensi,hsaccum,phie)` do Fortran passa `hsaccum` em `sdf` e
recebe `phie` em `phie1`. Isolado no adaptador `_snow_model`.

Pontos de atencao (comportamento do original preservado)
---------------------------------------------------------
1. `lowveg_met(...,dmet1,...)`: o original passa o array INTEIRO
   `dmet1(maxlines,13)` para o dummy `real(dp),intent(inout):: dmet(13)`
   (module_lowveg.F90:26, tamanho explicito). Isso e associacao de
   sequencia: `dmet(k)` e o k-esimo elemento de `dmet1` em ordem de
   COLUNA -- com `maxlines >= 13`, `dmet1(k,1)`, a radiacao solar dos
   passos 1..13, e nao a linha do passo `iw`. As escritas de
   `lowveg_met` em `dmet(6)`/`dmet(7)` caem nas mesmas posicoes. Fisica
   sem sentido, mas e o que o Fortran faz: conferido com gfortran
   (TDM-GCC), que compila sem aviso e associa exatamente assim. A
   traducao anterior passava `dmet1(iw,1:13)`; foi trocada pela forma
   literal, regra do repo (`indices.py:112`: "Preserve literally when
   transcribing; do not 'fix'"). Isolado em `_dmet_seq` /
   `_dmet_seq_grava`. Conferido contra o FASST compilado: com a forma
   literal o `_ro_soil.dat` do caso `multi` de `cmp_fortran_agua.py` sai
   identico; com a linha `iw` (mutante `DMET_LINHA=1`) diverge.
2. `watervel` e `waterdep` sao calculados (linhas 97-117) e nunca usados.
3. `phie` nao e inicializada na lista de zeragem: so recebe valor quando
   ha neve nova. Sem neve nova ela segue indefinida para `snow_model` e
   `read_old_data` (no Fortran, lixo de pilha). Aqui comeca em 0.0.
4. O `else` de `if(iw >= 2)` (linhas 237-239) e inalcancavel: iw == 1 ja
   foi tratado pelo ramo anterior (`iw == 1 .or. |hi| <= eps`).
5. `find_index` deixa `index1 = 0` se `tdate` nao cai em nenhum mes (doy
   366 em ano nao bissexto, doy <= 0 ou > 366). No Fortran `scorr(0)` le
   memoria fora do vetor; aqui usa a posicao 0 fantasma (0.0).
6. `dmet1(iw,4) = tmp + 273.15d0` usa literal, nao `Tref` (mesmo valor).
7. O cabecalho diz "calls read_old_data, missing_met, find_index,
   snow_model" e "uses the function met_date": `missing_met` e `met_date`
   NAO sao usadas no codigo; ja `canopy_met`, `low_veg_prop`,
   `lowveg_met`, `sp_humid` e `dense` sao chamadas e nao constam.
8. `read_old_data` (modo 2) usa `wstart = max(iend-(moverlap-1), istart)`
   com a constante `moverlap`, como initwater (initsurface usa
   `moverlapr`). O bloco anterior, comentado no original para os dados
   "imets" (linhas 301-310), nao foi traduzido.
9. O sdensi devolvido por `read_old_data` e descartado: `sdensi` e
   `intent(in)` em `open_water`.
10. `d1i = 4` em `dense(...)` (linha 178) corresponde a `Phase.SNOW`, como
    em module_canopy.py.
11. O ramo `else` da escrita no arquivo de neve usa `format 91`
    `(f10.2,8(f15.8),a5,f15.8)`; `meltfl` (len=1) sai alinhado a direita
    em `a5`.
"""

import math
from dataclasses import dataclass, field

from .albedo_emis import AlbedoEmisState
from .constants import (
    EPS, IALBEDO, IDENS, IEMIS, ITHCOND, LHFUS, MOVERLAP, NCLAYERS, PI,
    SIGMA, SPFLAG, STHCOND, TREF,
)
from .functions import Phase, dense
from .indices import MetCol
from .module_canopy import CanopySavedState, canopy_met
from .module_lowveg import LowVegState, low_veg_prop, lowveg_met
from .read_old_data import read_old_data
from .snow import SnowState, snow_model
from .sp_humid import sp_humid
from .state import FasstState

__all__ = ["OpenWaterState", "open_water", "find_index"]


@dataclass
class OpenWaterState:
    """Variaveis `save` de `open_water` e dos modulos que ela chama."""

    first_snow: int = 0
    canopy: CanopySavedState = field(default_factory=CanopySavedState)
    alb: AlbedoEmisState = field(default_factory=AlbedoEmisState)
    lowveg: LowVegState = field(default_factory=LowVegState)
    snow: SnowState = field(default_factory=SnowState)


# data scorr / data walbedo, wemis (open_water.F90, linhas 49-53); posicao 0 fantasma
_SCORR = [0.0, -12.5, -11.475, -10.45, -6.5, -2.55, -0.05, 2.45, 1.225, 0.0, -4.225, -8.45, -10.475]
_WALBEDO = 0.9
_WEMIS = 0.98

# data de find_index (linhas 336-340); posicao 0 fantasma
_DAYIMNT = [0, 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
_DAYEMNT = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
_INDEXM = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
_INDEXP = [0, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6]


def _aint(x, p=0):
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direcao a zero."""
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _anint(x, p=0):
    """Fortran ``ANINT(x*10**p)*10**-p`` -- arredonda para longe de zero.

    Diferente das copias de outros arquivos, usa `floor` e a parte
    fracionaria (exata) em vez de `floor(|x|+0.5)`, que erra para
    0.49999999999999994.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    a = abs(xs)
    r = math.floor(a)
    if a - r >= 0.5:
        r += 1.0
    return math.copysign(r, xs) * (10.0**-p)


def _fdiv(x, y):
    """Divisao IEEE-754: x/0 -> +-inf ou nan, sem levantar excecao."""
    x = float(x)
    y = float(y)
    if y == 0.0:
        if x == 0.0 or math.isnan(x):
            return math.nan
        neg = (math.copysign(1.0, x) < 0) != (math.copysign(1.0, y) < 0)
        return -math.inf if neg else math.inf
    return x / y


def _dsqrt(x):
    """DSQRT: raiz de negativo vira NaN (em vez de ValueError)."""
    return math.sqrt(x) if x >= 0.0 else math.nan


def _pow4(x):
    """``x**4d0`` com estouro -> +inf, como em IEEE."""
    try:
        return math.pow(x, 4.0)
    except OverflowError:
        return math.inf


def _dmet_pos(state, k):
    """Posicao (linha, coluna) do k-esimo elemento de `dmet1(maxlines,13)`
    em ordem de coluna -- a associacao de sequencia do ponto de atencao 1."""
    n = state.dmet1.shape[0] - 1                 # maxlines (posicao 0 fantasma)
    return (k - 1) % n + 1, (k - 1) // n + 1


def _dmet_seq(state):
    """`dmet(1:13)` como o dummy de `lowveg_met` o enxerga (ponto 1)."""
    return [float(state.dmet1[_dmet_pos(state, k)]) for k in range(1, 14)]


def _dmet_seq_grava(state, dmet):
    for k in range(1, 14):
        state.dmet1[_dmet_pos(state, k)] = dmet[k - 1]


def _snow_model(state, ow, sprint, sdensi):
    """Adaptador de `call snow_model(sprint,water_flag,first_snow,sdensi,hsaccum,phie)`.

    `phie` e `intent(out)` la (`phie1`): o valor de entrada nao e lido.
    """
    first_snow, hsaccum, phie = snow_model(
        state, ow.snow, sprint, state.water_flag, ow.first_snow, sdensi,
        state.hsaccum)
    ow.first_snow = first_snow
    state.hsaccum = hsaccum
    return phie


def open_water(state: FasstState, ow: OpenWaterState, wtype, sprint, wvel, wdepth,
               sdensi, first_loop, code2, code4, oldsd):
    iw = state.iw

    # ---------------- zero-out variables (open_water.F90, linhas 55-80) ----------------
    index1 = 0
    tdate = 0
    d1i = 0
    f1i = 0
    c1 = 0.0
    c2 = 0.0
    t1 = 0.0
    watervel = 0.0
    waterdep = 0.0
    dh = 0.0
    mixrgr = 0.0
    rhsurf = 0.0
    rhoaf = 0.0
    cf = 0.0
    pheatf = 0.0
    rpp = 0.0
    mixra = 0.0
    mixrf = 0.0
    wetbulba = 0.0
    pdens = 0.0
    smm = 0.0
    dates = 0.0
    xtemp = 0.0
    taf = 0.0
    stempt = 0.0
    # nao zeradas no original (ver ponto de atencao 3)
    phie = 0.0
    wstart = 0
    timeo = [[0.0] * 5 for _ in range(MOVERLAP + 1)]  # timeo(moverlap,4)

    if first_loop == 0:
        ow.first_snow = 0                       # first time through snow routines
        state.vsmelt = 0.0                      # snow melt (m)
        state.vimelt = 0.0                      # ice melt (m)
        first_loop = 1

    lflag = 1                                   # northern hemisphere
    if state.lat < 0.0:
        lflag = -1                              # southern hemisphere

    # ---------------- Main Calculation Loop (linhas 92-) ----------------
    tdate = int(state.met[iw, MetCol.DOY])

    index1 = find_index(tdate, int(state.met[iw, MetCol.YEAR]), lflag, index1)

    # water velocity (linhas 97-117) -- watervel/waterdep nao sao usados adiante
    if wtype == 1 or wtype == _aint(SPFLAG):    # open water
        watervel = 0.0                          # m/s
        waterdep = 10.0                         # m
    else:                                       # rivers, canals...
        if abs(wvel) <= EPS or abs(wvel - SPFLAG) <= EPS:     # unknown water velocity
            watervel = 1.0
        elif abs(wvel - 1) <= EPS:                            # >0 - <=1.5 m/s
            watervel = 0.75
        elif abs(wvel - 2) <= EPS:                            # >1.5 m/s
            watervel = 3.0

        if abs(wdepth) <= EPS or abs(wdepth - SPFLAG) <= EPS:  # unknown water depth
            waterdep = 0.5
        elif abs(wdepth - 1) <= EPS:                          # >0 - <=0.8 m
            waterdep = 0.4
        elif abs(wdepth - 2) <= EPS:                          # >0.8 m
            waterdep = 1.0

    # water surface temperature (linhas 119-123)
    lat = state.lat
    state.met[iw, MetCol.TSOIL] = (TREF + 56.182 - 0.95454 * abs(lat)
                                    + 0.0021307 * abs(lat) * abs(lat)
                                    + _SCORR[index1])                  # K
    state.tt[iw] = state.met[iw, MetCol.TSOIL]

    # initialize the met calc variable dmet1(iw,) (linhas 125-148)
    state.dmet1[iw, 1] = state.met[iw, MetCol.TSOL]
    if state.dmet1[iw, 1] < 1e-2:
        state.dmet1[iw, 1] = 0.0
    state.dmet1[iw, 2] = state.met[iw, MetCol.IR]
    state.dmet1[iw, 3] = state.met[iw, MetCol.WS]
    state.dmet1[iw, 4] = state.met[iw, MetCol.TMP] + 273.15                # K
    state.dmet1[iw, 5] = max(1e-1, state.met[iw, MetCol.RH])               # %
    if state.met[iw, MetCol.PREC] <= EPS:
        state.met[iw, MetCol.PREC] = 0.0
        state.met[iw, MetCol.PT] = 1
    state.dmet1[iw, 6] = state.met[iw, MetCol.PREC] * 1e-3                 # m/timstep
    if state.met[iw, MetCol.PREC2] <= EPS:
        state.met[iw, MetCol.PREC2] = 0.0
        state.met[iw, MetCol.PT2] = 1
    state.dmet1[iw, 7] = state.met[iw, MetCol.PREC2] * 1e-3                # m/timstep
    state.dmet1[iw, 8] = state.met[iw, MetCol.UPSOL]
    if abs(state.dmet1[iw, 1]) <= EPS:
        state.dmet1[iw, 8] = 0.0
    state.dmet1[iw, 9] = state.met[iw, MetCol.DIR]
    state.dmet1[iw, 10] = state.met[iw, MetCol.DIF]
    state.dmet1[iw, 11] = state.met[iw, MetCol.AP]
    state.dmet1[iw, 12] = state.met[iw, MetCol.IRUP]
    state.dmet1[iw, 13] = state.met[iw, MetCol.TSOIL]

    # calculate the canopy temperature and moisture (linhas 150-160)
    if state.veg_flagh == 1:
        canopy_met(state, ow.canopy, ow.alb, oldsd)               # correct met for canopy
        if state.dmet1[iw, 7] > 0.0:
            state.met[iw, MetCol.PT2] = 3.0
    else:
        state.surfemisc[iw] = state.mflag
        for i in range(1, NCLAYERS + 1):
            state.canopy_temp[i, iw] = state.mflag
    if state.error_type == 2:
        code2 = 1

    state.airt[1, iw] = state.met[iw, MetCol.TMP] + TREF
    state.airt[2, iw] = state.dmet1[iw, 4]

    # get the low vegetation properties (linhas 165-170)
    if state.veg_flagl == 1:
        low_veg_prop(state, state.dmet1[iw, 4], state.dzveg[NCLAYERS])
    else:
        state.sigfl = 0.0

    # determine density and amount of new snow if no existing snow on ground (linhas 172-196)
    state.newsd = 0.0
    pdens = 0.0
    smm = state.dmet1[iw, 6] + state.dmet1[iw, 7]
    if ((_aint(state.met[iw, MetCol.PT]) == 3 or _aint(state.met[iw, MetCol.PT2]) == 3)
            and smm > 0.0):
        d1i = 4
        pdens = dense(state.dmet1[iw, 4], state.dmet1[iw, 3], Phase.SNOW)      # kg/m^3 (d1i=4)
        if _aint(state.met[iw, MetCol.PT]) == 3:
            state.newsd = state.dmet1[iw, 6] * (pdens * 1e-3)                  # m
        else:
            state.newsd = state.dmet1[iw, 7] * (pdens * 1e-3)
        state.newsd = _anint(state.newsd, 10)
        state.dsnow = 2e-2 * ((3.0 / (4.0 * PI)) * 5e-5) ** (1.0 / 3.0)
        state.dsnow = _anint(state.dsnow, 10)
        phie = (1.0 - pdens * 1e-3) * 0.95
        phie = _anint(phie, 10)
        if abs(state.hsaccum) <= EPS and state.newsd < 2e-3:
            state.dsnow = 0.0
            phie = 0.0
            state.newsd = 0.0
            pdens = 0.0

    state.hm = state.hsaccum + state.newsd + state.hi                          # m

    # calculate the low vegetation temperature (linhas 200-224)
    if state.veg_flagl == 1:
        if state.hm > EPS:
            f1i = 1
        c2 = 1.0
        moist = sp_humid(f1i, state.dmet1[iw, 11], state.met[iw, MetCol.TSOIL], c2, c1)
        mixrgr = moist.mixr             # as outras 9 saidas (t1) sao descartadas

        d1i = 0
        rhsurf = 1.0
        if state.hm <= EPS:
            stempt = TREF
        else:
            stempt = state.dmet1[iw, 13]
        taf = ((1.0 - 7e-1 * state.sigfl) * state.dmet1[iw, 4]
               + state.sigfl * (6e-1 * state.ftemp + 1e-1 * stempt))           # foliage temp at the ground surface (K)
        dmet = _dmet_seq(state)                   # sic: dmet1 inteiro (ponto 1)
        rhoaf, cf, pheatf, rpp, mixra, mixrf, dqdtf = lowveg_met(
            state, ow.lowveg, d1i, d1i, int(state.met[iw, MetCol.PT]),
            int(state.met[iw, MetCol.PT2]), mixrgr, rhsurf, taf, dmet, wetbulba)
        _dmet_seq_grava(state, dmet)
        state.ft[iw] = state.ftemp
        state.surfemisf[iw] = state.epf
    else:
        state.ft[iw] = state.mflag
        state.surfemisf[iw] = state.mflag

    # compute any ice growth (linhas 226-244)
    if state.met[iw, MetCol.TSOIL] <= TREF:
        if iw == 1 or abs(state.met[iw, MetCol.HI]) <= EPS:
            state.met[iw, MetCol.HI] = (_dsqrt(2.0 * ITHCOND / (IDENS * LHFUS))
                                        * _dsqrt(abs(state.met[iw, MetCol.TMP]) * state.timstep * 3600.0))
            if state.met[iw, MetCol.TMP] > 0.0:
                state.met[iw, MetCol.HI] = 0.0
        else:
            dh = (1.0 / (IDENS * LHFUS)) * _fdiv(
                -(state.met[iw, MetCol.TMP]),
                _fdiv(state.met[iw, MetCol.HI], ITHCOND) + _fdiv(state.hsaccum, STHCOND))
            if iw >= 2:     # sempre verdadeiro aqui (ponto de atencao 4)
                state.met[iw, MetCol.HI] = state.met[iw - 1, MetCol.HI] + dh * state.timstep * 3.6e3
            else:
                state.met[iw, MetCol.HI] = dh * state.timstep * 3.6e3
            state.met[iw, MetCol.HI] = max(0.0, state.met[iw, MetCol.HI])
    if _aint(abs(state.met[iw, MetCol.HI] - state.mflag) * 1e5) * 1e-5 <= EPS:
        state.met[iw, MetCol.HI] = 0.0

    # compute snow depth (linhas 246-268)
    if (state.hsaccum > EPS
            or (_aint(state.met[iw, MetCol.PT]) == 3 or _aint(state.met[iw, MetCol.PT2]) == 3)):
        phie = _snow_model(state, ow, sprint, sdensi)

        state.vsmelt = state.vsmelt * 1e-2                                      # m
        if state.error_type == 4:
            code4 = 1
    else:
        state.vsmelt = 0.0

        if sprint == 1:
            dates = (state.met[iw, MetCol.DOY] + state.met[iw, MetCol.HR] / 24.0
                     + state.met[iw, MetCol.MIN] / (24.0 * 60.0))
            xtemp = 0.0
            state.meltfl = "0"
            # format 91: (f10.2,8(f15.8),a5,f15.8)
            state.unit55.write(
                "{:10.2f}".format(dates)
                + "".join("{:15.8f}".format(v) for v in
                          (xtemp, xtemp, xtemp, xtemp, xtemp, xtemp, state.hsaccum, xtemp))
                + "{:>5s}".format(state.meltfl)
                + "{:15.8f}".format(xtemp) + "\n")
    state.met[iw, MetCol.SD] = state.hsaccum
    oldsd = state.hsaccum

    # compute the met parameters (linhas 270-281)
    if abs(state.met[iw, MetCol.HI]) <= EPS:
        state.met[iw, MetCol.UPSOL] = state.met[iw, MetCol.TSOL] * _WALBEDO
        if state.met[iw, MetCol.TSOIL] > EPS:
            state.met[iw, MetCol.IRUP] = _WEMIS * SIGMA * _pow4(state.met[iw, MetCol.TSOIL])
        state.surfemis[iw] = _WEMIS
    else:
        state.met[iw, MetCol.UPSOL] = state.met[iw, MetCol.TSOL] * IALBEDO
        if state.met[iw, MetCol.TSOIL] > EPS:
            state.met[iw, MetCol.IRUP] = IEMIS * SIGMA * _pow4(state.met[iw, MetCol.TSOIL])
        state.surfemis[iw] = IEMIS

    state.surfci[iw] = state.mflag
    state.surfrci[iw] = state.mflag
    state.surfcbr[iw] = state.mflag
    state.surfd[iw] = state.mflag
    if abs(state.met[iw, MetCol.HI]) <= EPS:
        state.surfmoist[iw] = 1.0
        state.surfice[iw] = 0.0
        state.surfd[iw] = 1e3
    else:
        state.surfmoist[iw] = 0.0
        state.surfice[iw] = 1.0
        state.surfd[iw] = IDENS
    state.surfmoistp[iw] = state.surfmoist[iw]
    state.surficep[iw] = state.surfice[iw]

    if state.met[iw, MetCol.HI] > EPS:
        state.tt[iw] = TREF

    # write the new input file for the next model run (linhas 301-316)
    # (bloco "imets" comentado no original nao traduzido -- ponto de atencao 8)
    wstart = max(state.iend - (MOVERLAP - 1), state.istart)
    if iw >= wstart:
        d1i = 2
        _, phie = read_old_data(state, d1i, iw, wstart, sdensi, phie, timeo)

    return first_loop, code2, code4, oldsd


def find_index(tdate, year, lflag, index1=0):
    """Indice do mes (1-12) em `scorr` a partir do dia do ano (open_water.F90, linhas 321-361).

    `index1` (intent(out) no original) so e atribuido se `tdate` cai em
    algum mes; caso contrario mantem o valor recebido (por referencia, no
    Fortran) -- por isso o parametro com o valor atual e devolvido.
    """
    iday = 0
    eday = 0
    tflag = 0
    i = 1
    while i <= 12 and tflag == 0:
        mody = year - _aint(year * 2.5e-1) * 4.0
        if abs(mody) <= EPS and i > 2:
            iday = 1
        if abs(mody) <= EPS and i > 1:
            eday = 1

        if tdate >= _DAYIMNT[i] + iday and tdate <= _DAYEMNT[i] + eday:
            if lflag < 0:
                index1 = _INDEXP[i]
            if lflag > 0:
                index1 = _INDEXM[i]
            tflag = 1

        i += 1
    return index1
