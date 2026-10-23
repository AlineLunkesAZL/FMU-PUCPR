"""
fasst_main.py -- traducao de fasst_main.F90

Fonte Fortran: fasst_main.F90 (586 linhas). Uma unica sub-rotina,
`fasst_main`, chamada pelo driver uma vez por passo de tempo `iw` no
caminho de TERRA (`water_flag <= 2`). E a casca que integra, em ordem, os
modelos de fluxo solar/IR, temperatura e umidade do solo, resistencia do
solo, neve e gelo, organizada pelos marcadores STEP do proprio original:

  STEP 1   (148-198)  met do passo (`dmet1`), propriedades da vegetacao baixa
  STEP 2   (200-216)  temperatura/umidade do dossel (`canopy_met`)
  STEP 3   (218-271)  neve nova, `new_profile` (perfis de temperatura e
                      umidade do solo), emissividades
  STEP 3a  (274-281)  fluxos solar e IR emitidos pela superficie
  STEP 4   (284-287)  resistencia do solo (`soil_strength`)
  STEP 5   (290-348)  modelo de neve e gelo de superficie
  STEP 6   (351-581)  umidade/gelo medios de superficie, saida nodal
                      (unidade 3), profundidade de congelamento/degelo,
                      arquivo de reinicio (`read_old_data`, modo 2)

Assinatura original:
    subroutine fasst_main(nprint,sprint,frozen,nstr_flag,sdensi,surfdepth,
                          code2,code3,code4,code5,first_loop,phie,pdens,
                          oldsd,oldhi)
    integer(ip),intent(in):: nprint,sprint,frozen
    integer(ip),intent(in):: nstr_flag(maxn)
    real(dp),intent(in):: sdensi,surfdepth
    integer(ip),intent(inout):: code2,code3,code4,code5,first_loop
    real(dp),intent(inout):: phie,pdens,oldsd,oldhi

Estado desta transcricao
--------------------------
COMPLETA, e ligada no driver: `fasst_driver.main` chama `fasst_main` a
cada `iw` do caminho de terra. Os STEPS 1, 2 e 3 (linhas 148-271) sao a
traducao do Allan (`fasst_main_allan.py`, fundida aqui); os STEPS 3a, 4, 5
e 6 e o corpo que encadeia os sete, do Lucas. `new_profile` e a real,
ligada pelo adaptador `_new_profile`.

Arquitetura
-------------
Uma funcao por STEP (`step1`, `step2`, `step3`, `step3a`, `step4`,
`step5`, `step6`); o corpo de `fasst_main` so as chama em sequencia. Isto
se afasta da regra do driver ("um `program`, UM escopo", fasst_driver.py,
Arquitetura), e pode, porque as fronteiras foram conferidas: `d1i`/`d2i`
sao reatribuidas na 353 antes de qualquer leitura; `simelt` so vive no
STEP 3; `num2`, `dwi`, `sfdp`, `nsstate`, `j`, `k`, `wstart` e `timeo` so
no STEP 6, e todas comecam do mesmo valor que a zeragem do topo (115-124)
lhes da. As unicas locais que atravessam uma fronteira sao `fsup`/`firup`
(STEP 3 -> 3a), que viram retorno de `step3` e parametro de `step3a`.

Os `intent(inout)` que cada STEP toca voltam numa tupla, como no resto do
repo: `phie, oldsd, oldhi, code2, code4 = step5(...)`. Campos de
`fasst_global` (hsaccum, hi, vsmelt, met, surfci, ...) sao mutados em
`state`.

`FasstMainState` reune os `save` de `fasst_main` (`first_snow`,
`first_time`, linha 118) e os das rotinas que `fasst_main` e
`open_water` chamam (`snow`, `canopy`, `alb`, `lowveg`).
Esses quatro TEM de ser os mesmos objetos de `OpenWaterState`: no Fortran
os `save` pertencem as rotinas chamadas, nao a quem chama. O driver cria
um `FasstMainState` so, junto do `OpenWaterState`, com os quatro
compartilhados. Ja `first_snow` NAO e compartilhado: `open_water` tem o
seu proprio `save`. `profile` (os `save` de `new_profile`) e so daqui:
`open_water` nao chama `new_profile`.

Constantes (`EPS`, `PI`, `TREF`, `NCLAYERS`, `MOVERLAP`) vem de
`fasst.constants`; colunas de `met` sao `MetCol.*`. As 13 colunas de
`dmet1` nao tem enum e sao usadas por numero. Indexacao 1-based identica
ao Fortran (posicao 0 sem uso).

Dependencias externas (todas reais)
-------------------------------------
  low_veg_prop(state, fatemp, c3_hgt)  fasst.module_lowveg (STEP 1).
  canopy_met(state, canopy, alb, oldsd) fasst.module_canopy (STEP 2).
  new_profile(state, save, code5, first_time, pdens, oldsd, simelt, phie,
              fsup, firup, *, <7 dependencias injetadas>) -> NewProfileResult
                      fasst.new_profile (STEP 3), pelo adaptador
                      `_new_profile`. `th_param`, `flow_param` e `soil_tmp`
                      entram diretas; as outras quatro tem assinatura
                      diferente do contrato e ganham um invólucro:
                      `albedo_emis` recebe `fm.alb`; `sub_divide` devolve
                      `dmet` e o invólucro o copia no array de
                      `new_profile`; `soil_moisture` devolve `runoff` e o
                      invólucro o copia no lugar; `lowveg_met` recebe
                      `fm.lowveg`, le `dmet` 0-based (o de `surfenergy`
                      e 1-based: o invólucro converte e copia de volta,
                      e `intent(inout)`) e devolve 7-tupla, que vira
                      `LowVegMetResult` com o `wetbulba` de entrada
                      (module_lowveg nao o altera).
  soil_strength(state, nstr_flag) -> (ss_strengtho, ll_strengtho, scbro)
                      fasst.soil_strength; a ordem confere com o dummy de
                      soil_strength.F90 (`nstr_flag,ss_strengtho,
                      ll_strengtho,scbro`), que a linha 287 liga a
                      `surfci(iw),surfrci(iw),surfcbr(iw)`.
  snow_model(state, snow, sprint, water_flag, first_snow, sdensi, sdf)
      -> (first_snow, sdf, phie1)      fasst.snow. O `call snow_model(...,
                      hsaccum,phie)` passa `hsaccum` em `sdf` e recebe
                      `phie` em `phie1`. Isolado no adaptador `_snow_model`,
                      no molde de `open_water._snow_model`.
  icethick(state, hinit) -> hfinal     fasst.icethick; `call icethick(
                      oldhi,hi)` vira `state.hi = icethick(state, oldhi)`.
  read_old_data (modo 2), dense.

Pontos de atencao (comportamento do original preservado)
---------------------------------------------------------
1. `vsmelt = vsmelt*1d-2` (linha 298) converte cm -> m, mas
   `module_snow.F90:302` (e `snow.py:477`) JA devolve `vsmelt` em metros.
   No ramo agua-com-gelo o valor sai 100x menor; no ramo de terra (316) nao
   ha conversao nenhuma. Defeito do original, preservado. `open_water.F90`
   faz a mesma conversao dupla.
2. STEP 5 testa o tipo de precipitacao com `aint` (293, 342); o STEP 3 mistura
   `aint` (223) e `anint` (228) sobre o mesmo campo. Para `pt` inteiro e
   igual; para `pt = 3.7` o STEP 3 entra no bloco de neve (aint = 3) e
   escolhe `dmet1(iw,7)` (anint = 4). Transcrito literalmente.
3. `if(error_type == 2) code2 = 1` depois de `icethick` (344) nao pode ser
   disparado por ela: `icethick` tem fan-out zero e nunca escreve
   `error_type`. So repete o teste da linha 213 sobre o `error_type` que
   `canopy_met` (STEP 2) deixou -- e ainda roda quando `icethick` nao foi
   chamada (`hsaccum > eps`).
4. O comentario da linha 277 diz "- to reflect away from surface", mas o
   codigo SOMA `fsup`: `upsol = sup(iw) + fsup`.
5. `icethick` so roda com `|hsaccum| <= eps` (343), mas `oldhi = hi` (345)
   roda sempre que o bloco externo e atingido. Com neve, `hi` nao muda e
   `oldhi` recebe o `hi` antigo.
6. Os tres sub-ramos do STEP 5 parecem cobrir todo `node_type(nnodes)`, mas
   nao cobrem `ice(nnodes)` NaN num no `'WA'`: `NaN > eps` e `NaN <= eps`
   sao ambos falsos e o terceiro exige no /= 'WA'. Nenhum ramo dispara;
   `vsmelt` e `hsaccum` ficam como estavam. Mesmo comportamento aqui.
7. STEP 4 so roda com `water_flag == 0`. Com `water_flag` 1 ou 2,
   `surfci/surfrci/surfcbr(iw)` ficam com o que `zero_mxl_params` deixou.
8. O `format 91` daqui (335) NAO e o de `open_water.F90`: e
   `(i10,2x,f10.3,8(f15.8),a5,3(f15.8))`, com o ANO na frente e 14 valores.
   `meltfl` (len=1) sai alinhado a direita em `a5`. A linha dos ramos COM
   neve e a do `format 91` de `module_snow`, escrita por `snow.py`.
9. O bloco de gelo de superficie (342) testa so `ip_pt` (2 = chuva, 4 =
   chuva congelante), nunca `ip_pt2`.
10. STEP 6, media de superficie (366-370): o terceiro ramo do `do while`
   testa `aint((elev-nzi(j))*1d5)*1d-5` -- TRUNCA, onde os outros dois e a
   guarda do laco arredondam com `anint`. E quando nenhum dos tres ramos
   pega, `num2` fica com o valor da volta anterior. Os dois preservados.
11. STEP 6, congelamento/degelo: o `do while` interno de cada ramo compara
   `ice(j)` com `pfrozen` SEM o fator `dwi` (gelo -> agua) que o teste
   de fora usa. (No ramo nao congelado, os dois testes leem
   `sstate(iw-1)`, nao o recem-escrito: f -> t -> u, um passo por vez --
   e o pretendido, nao defeito.)
12. `timeo` e `timeo(2,4)` aqui (122) e `timeo(moverlap,4)` em
   read_old_data. So o modo 3 escreve nele; no modo 2, o unico que
   `fasst_main` usa (580), a declaracao menor e inofensiva.
13. A mensagem `' Start of run'` (139/143) sai pelo list-directed do
   gfortran: branco de inicio de registro, `I12` e um branco antes de
   texto que segue numero -- `"  Start of run         100 of ..."`,
   conferido no TDM-GCC.
14. `dmet1(iw,5) = dmin1(met(iw,ip_rh),1d2)` so limita a umidade a 100%. A
   linha comentada logo acima (165, `dmax1(1d-1, rh)`) e a usada em
   open_water; aqui nao ha piso, entao umidade 0 ou negativa passa.
15. `met` e alterado no proprio lugar: precipitacao <= eps vira 0 com tipo
   1 (166-179); depois de canopy_met, `met(iw,ip_pt2) = 3` se houver neve
   (204); e `met(iw,ip_tsoil)` e SOBRESCRITO com a temperatura do no
   superior (256). `dmet1(iw,13)` guarda a temperatura original.
16. O bloco 239-244 (zerar `dsnow`, `phie`, `newsd` e `pdens` sem neve
   acumulada e `newsd < 2d-3`) esta COMENTADO aqui; em open_water esta
   ativo. Nao transcrito.
17. `dsnow` nao depende de entrada nenhuma: vale sempre
   2d-2*((3/(4*pi))*5d-5)**(1/3) arredondado a 10 casas (~4,571d-4).
18. `canopy_met` recebe `oldsd`, a profundidade de neve do passo ANTERIOR
   (`oldsd = hsaccum` so no STEP 5).
19. O arredondamento a 15 casas de `dmet1` (189-191) acontece antes de
   `canopy_met`, que altera `dmet1` sem arredondar de novo.
20. `pdens` e zerado a cada chamada; so e diferente de zero com neve nova.
   Divisao por zero em `newsd` (pdens = 0) segue IEEE (inf ou NaN), via
   `_fdiv`.
21. O cabecalho (97-105) cita low_veg_prop, canopy_met, new_profile,
   soil_strength, snow_model, icethick, dense e met_date. `met_date` e
   declarada (122) e nunca usada; `read_old_data` (580) e chamada e nao
   consta.
"""

import math
from dataclasses import dataclass, field

from .albedo_emis import AlbedoEmisState, albedo_emis
from .constants import EPS, MOVERLAP, NCLAYERS, PI, TREF
from .flow_param import flow_param
from .functions import Phase, dense
from .icethick import icethick
from .indices import MetCol
from .module_canopy import CanopySavedState, canopy_met
from .module_lowveg import LowVegState, low_veg_prop, lowveg_met
from .new_profile import NewProfileState, new_profile
from .read_old_data import read_old_data
from .snow import SnowState, snow_model
from .soil_moisture import soil_moisture
from .soil_strength import soil_strength
from .soil_tmp import soil_tmp
from .state import FasstState
from .sub_divide import sub_divide
from .surfenergy import LowVegMetResult
from .th_param import th_param
from .write_outputs import _fmt_e, _fmt_f, _fmt_i

__all__ = ["FasstMainState", "fasst_main", "step1", "step2", "step3", "step3a",
           "step4", "step5", "step6"]


@dataclass
class FasstMainState:
    """Variaveis `save` de `fasst_main` (linha 118) e os `save` das rotinas
    que ela chama.

    `canopy`, `alb`, `lowveg` e `snow` tem de ser os MESMOS objetos de
    `OpenWaterState` (ver docstring): no Fortran sao `save` das rotinas
    chamadas, um so por programa.
    """

    first_snow: int = 0
    first_time: int = 0
    snow: SnowState = field(default_factory=SnowState)
    canopy: CanopySavedState = field(default_factory=CanopySavedState)
    alb: AlbedoEmisState = field(default_factory=AlbedoEmisState)
    lowveg: LowVegState = field(default_factory=LowVegState)
    profile: NewProfileState = field(default_factory=NewProfileState)   # so fasst_main chama new_profile


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


def _anint(x, p=0):
    """Fortran ``ANINT(x*10**p)*10**-p`` -- arredonda para longe de zero."""
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    return math.copysign(math.floor(abs(xs) + 0.5), xs) * (10.0**-p)


def _aint(x, p=0):
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direcao a zero."""
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _snow_model(state, fm, sprint, sdensi):
    """Adaptador de `call snow_model(sprint,water_flag,first_snow,sdensi,hsaccum,phie)`.

    `phie` e `intent(out)` la (`phie1`): o valor de entrada nao e lido.
    """
    first_snow, hsaccum, phie = snow_model(
        state, fm.snow, sprint, state.water_flag, fm.first_snow, sdensi,
        state.hsaccum)
    fm.first_snow = first_snow
    state.hsaccum = hsaccum
    return phie


def _write_snw_zero(state):
    """Linha de neve "vazia" na unidade 55 (linhas 306-312 / 326-332)."""
    iw = state.iw
    dates = (state.met[iw, MetCol.DOY] + state.met[iw, MetCol.HR] / 24.0
             + state.met[iw, MetCol.MIN] / (24.0 * 60.0))
    xtemp = 0.0
    state.meltfl = "0"
    # format 91: (i10,2x,f10.3,8(f15.8),a5,3(f15.8)) -- ponto de atencao 8
    state.unit55.write(
        "{:10d}".format(int(state.met[iw, MetCol.YEAR])) + "  "
        + "{:10.3f}".format(dates)
        + "".join("{:15.8f}".format(v) for v in
                  (xtemp, xtemp, xtemp, xtemp, xtemp, xtemp, state.hsaccum, xtemp))
        + "{:>5s}".format(state.meltfl)
        + "".join("{:15.8f}".format(v) for v in (xtemp, xtemp, xtemp)) + "\n")


def step3a(state: FasstState, fsup, firup):
    """STEP 3a (linhas 274-281): fluxos solar e IR emitidos pela superficie.

    So preenche `upsol`/`irup` quando faltam no met (= mflag). `fsup`/`firup`
    sao as contribuicoes da vegetacao baixa calculadas por `new_profile`.
    """
    iw = state.iw
    met = state.met
    if _aint(abs(met[iw, MetCol.UPSOL] - state.mflag) * 1e5) * 1e-5 <= EPS:
        met[iw, MetCol.UPSOL] = state.sup[iw] + fsup         # sic: soma (ponto de atencao 4)
    if abs(met[iw, MetCol.TSOL]) <= EPS:
        met[iw, MetCol.UPSOL] = 0.0

    if _aint(abs(met[iw, MetCol.IRUP] - state.mflag) * 1e5) * 1e-5 <= EPS:
        met[iw, MetCol.IRUP] = state.irup[iw] + firup


def step4(state: FasstState, nstr_flag):
    """STEP 4 (linhas 284-287): resistencia do solo, so em terra (ponto 7)."""
    iw = state.iw
    if state.water_flag == 0:
        state.surfci[iw], state.surfrci[iw], state.surfcbr[iw] = soil_strength(state, nstr_flag)


def step5(state: FasstState, fm: FasstMainState, sprint, sdensi, phie, oldsd,
          oldhi, code2, code4):
    """STEP 5 (linhas 290-348): modelo de neve e gelo de superficie.

    Devolve `(phie, oldsd, oldhi, code2, code4)`.
    """
    iw = state.iw
    nnodes = state.nnodes
    met = state.met
    top = state.node_type[nnodes]

    # snow_model calculates the snow depth (linhas 292-334)
    if ((state.hsaccum > EPS or top == "SN")
            or (_aint(met[iw, MetCol.PT]) == 3 or _aint(met[iw, MetCol.PT2]) == 3)):
        if top == "WA" and state.ice[nnodes] > EPS:                 # water with ice
            phie = _snow_model(state, fm, sprint, sdensi)

            state.vsmelt = state.vsmelt * 1e-2                      # m (sic: ja esta em m, ponto 1)
            if state.error_type == 4:
                code4 = 1
        elif (top == "WA" and state.ice[nnodes] <= EPS) or top == "AI":   # water, no ice
            state.vsmelt = 0.0
            state.hsaccum = 0.0

            if sprint == 1:
                _write_snw_zero(state)
        elif top != "WA" and top != "AI":                           # ground
            phie = _snow_model(state, fm, sprint, sdensi)

            if state.error_type == 4:
                code4 = 1
        # sem else: 'WA' com ice NaN nao cai em nenhum ramo (ponto 6)
    else:
        state.vsmelt = 0.0
        state.hsaccum = 0.0

        if sprint == 1:
            _write_snw_zero(state)
    met[iw, MetCol.SD] = state.hsaccum
    oldsd = state.hsaccum

    # surface ice thickness model (linhas 339-348)
    if top != "WA":
        if (state.hi > EPS
                or (_aint(met[iw, MetCol.PT]) == 2 or _aint(met[iw, MetCol.PT]) == 4)):   # so ip_pt (ponto 9)
            if abs(state.hsaccum) <= EPS:
                state.hi = icethick(state, oldhi)
            if state.error_type == 2:                               # nunca vem de icethick (ponto 3)
                code2 = 1
            oldhi = state.hi                                        # roda mesmo sem icethick (ponto 5)
    met[iw, MetCol.HI] = state.hi

    return phie, oldsd, oldhi, code2, code4


def _dwi(state, k):
    """`dense(stt(k),0d0,d2i)/dense(stt(k),0d0,d1i)` -- gelo / agua."""
    return dense(state.stt[k], 0.0, Phase.ICE) / dense(state.stt[k], 0.0, Phase.WATER)


def _q5(x):
    return _anint(x, 5)


def step1(state: FasstState):
    """STEP 1 (linhas 148-198): met do passo em `dmet1(iw,·)` e propriedades
    da vegetacao baixa."""
    iw = state.iw
    met = state.met
    dmet1 = state.dmet1
    # initialize the met calc variable dmet1(iw,)
    dmet1[iw, 1] = met[iw, MetCol.TSOL]
    if dmet1[iw, 1] < 1e-2:
        dmet1[iw, 1] = 0.0
    dmet1[iw, 2] = met[iw, MetCol.IR]
    dmet1[iw, 3] = met[iw, MetCol.WS]
    dmet1[iw, 4] = met[iw, MetCol.TMP] + TREF                  # K

    # dmin1 -> min: com NaN o resultado do Fortran depende do compilador
    # (gfortran 13: -O0 propaga o NaN, como o min do Python; -O2 devolve 100)
    dmet1[iw, 5] = min(met[iw, MetCol.RH], 1e2)                # % (ponto 14)
    if met[iw, MetCol.PREC] <= EPS:
        met[iw, MetCol.PREC] = 0.0
        met[iw, MetCol.PT] = 1

    dmet1[iw, 6] = met[iw, MetCol.PREC] * 1e-3                 # m/timstep (rain and/or snow)
    if dmet1[iw, 6] <= EPS:
        met[iw, MetCol.PT] = 1.0

    if met[iw, MetCol.PREC2] <= EPS:
        met[iw, MetCol.PREC2] = 0.0
        met[iw, MetCol.PT2] = 1
    dmet1[iw, 7] = met[iw, MetCol.PREC2] * 1e-3                # m/timstep (snow only)
    if dmet1[iw, 7] <= EPS:
        met[iw, MetCol.PT2] = 1.0

    dmet1[iw, 8] = met[iw, MetCol.UPSOL]
    if abs(dmet1[iw, 1]) <= EPS:
        dmet1[iw, 8] = 0.0
    dmet1[iw, 9] = met[iw, MetCol.DIR]
    dmet1[iw, 10] = met[iw, MetCol.DIF]
    dmet1[iw, 11] = met[iw, MetCol.AP]
    dmet1[iw, 12] = met[iw, MetCol.IRUP]
    dmet1[iw, 13] = met[iw, MetCol.TSOIL]                      # original (ponto 15)

    for i in range(1, 14):
        dmet1[iw, i] = _anint(dmet1[iw, i] * 1e15) * 1e-15    # antes de canopy_met (ponto 19)

    # get the low vegetation properties
    if state.veg_flagl == 1:
        low_veg_prop(state, dmet1[iw, 4], state.dzveg[NCLAYERS])
    else:
        state.sigfl = 0.0


def step2(state: FasstState, fm: FasstMainState, oldsd, code2):
    """STEP 2 (linhas 200-216): temperatura/umidade do dossel. Devolve `code2`."""
    iw = state.iw
    met = state.met
    if state.veg_flagh == 1:
        canopy_met(state, fm.canopy, fm.alb, oldsd)            # oldsd do passo ANTERIOR (ponto 18)
        if state.dmet1[iw, 7] > 0.0:
            met[iw, MetCol.PT2] = 3.0
        if state.dmet1[iw, 6] <= EPS:
            met[iw, MetCol.PT] = 1.0
        if state.dmet1[iw, 7] <= EPS:
            met[iw, MetCol.PT2] = 1.0
    else:
        state.surfemisc[iw] = state.mflag
        for i in range(1, NCLAYERS + 1):
            state.canopy_temp[i, iw] = state.mflag
    if state.error_type == 2:
        code2 = 1

    state.airt[1, iw] = met[iw, MetCol.TMP] + TREF
    state.airt[2, iw] = state.dmet1[iw, 4]
    return code2


def _new_profile(state, fm, code5, pdens, oldsd, simelt, phie):
    """Adaptador de `call new_profile(code5,first_time,pdens,oldsd,simelt,phie,fsup,firup)`.

    Liga as dependencias injetadas de `new_profile` as funcoes reais,
    ajustando as assinaturas que diferem do contrato (ver docstring).
    `fsup`/`firup` entram 0: `new_profile` os zera antes de ler.
    Devolve `(code5, phie, fsup, firup)`.
    """
    def _albedo_emis(st, o, sd, su):
        albedo_emis(st, fm.alb, o, sd, su)

    def _sub_divide(st, ii, delvar, delvar1, dmet):
        dmet[1:14] = sub_divide(st, ii, delvar, delvar1)[1:14]     # intent(out) no lugar

    def _soil_moisture(st, *args):
        runoff = args[-1]
        r, rhs_errorm, errorm, smerror, sert1 = soil_moisture(st, *args)
        runoff[:] = r                                               # intent(inout) no lugar
        return rhs_errorm, errorm, smerror, sert1

    def _lowveg_met(st, ii, iter_, pt1, pt2, mixrgr, wetness, taf, dmet, wetbulba):
        # `dmet` aqui e 1-based (posicao 0 fantasma); module_lowveg.lowveg_met
        # o le 0-based (dmet[0] == dmet(1)). intent(inout): volta copiado.
        seq = [float(v) for v in dmet[1:14]]
        rhoaf, cf, pheatf, rpp, mixra, mixrf, dqdtf = lowveg_met(
            st, fm.lowveg, ii, iter_, pt1, pt2, mixrgr, wetness, taf, seq, wetbulba)
        dmet[1:14] = seq
        # module_lowveg nao altera `wetbulba`: volta o valor de entrada
        return LowVegMetResult(rhoaf=rhoaf, cf=cf, pheatf=pheatf, rpp=rpp, mixra=mixra,
                               mixrf=mixrf, dqdtf=dqdtf, wetbulba=wetbulba)

    res = new_profile(
        state, fm.profile, code5, fm.first_time, pdens, oldsd, simelt, phie, 0.0, 0.0,
        albedo_emis_fn=_albedo_emis, th_param_fn=th_param, flow_param_fn=flow_param,
        sub_divide_fn=_sub_divide, soil_moisture_fn=_soil_moisture,
        soil_tmp_fn=soil_tmp, lowveg_met=_lowveg_met)
    fm.first_time = res.first_time
    return res.code5, res.phie, res.fsup, res.firup


def step3(state: FasstState, fm: FasstMainState, code3, code5, pdens, oldsd, phie):
    """STEP 3 (linhas 218-271): neve nova e `new_profile`.

    Devolve `(code3, code5, pdens, phie, fsup, firup)`.
    """
    iw = state.iw
    met = state.met
    dmet1 = state.dmet1
    # determine density and amount of new snow if no existing snow on ground
    state.newsd = 0.0
    pdens = 0.0                                                # zerado a cada passo (ponto 20)
    smm = dmet1[iw, 6] + dmet1[iw, 7]
    if (_aint(met[iw, MetCol.PT]) == 3 or _aint(met[iw, MetCol.PT2]) == 3) and smm > 0.0:
        pdens = dense(dmet1[iw, 4], dmet1[iw, 3], Phase.SNOW)  # kg/m^3 (d1i = 4)

        if _anint(met[iw, MetCol.PT]) == 3:                    # sic: anint (ponto 2)
            state.newsd = _fdiv(dmet1[iw, 6] * 1e3, pdens)     # m snow
        else:
            state.newsd = _fdiv(dmet1[iw, 7] * 1e3, pdens)
        state.newsd = _anint(state.newsd * 1e10) * 1e-10
        state.dsnow = 2e-2 * ((3.0 / (4.0 * PI)) * 5e-5) ** (1.0 / 3.0)   # constante (ponto 17)
        state.dsnow = _anint(state.dsnow * 1e10) * 1e-10
        phie = (1.0 - pdens * 1e-3) * 0.95
        phie = _anint(phie * 1e10) * 1e-10
        # linhas 239-244 comentadas no original (ponto 16)

    # calculate the new moisture and temperature profiles
    simelt = state.vsmelt + state.vimelt
    state.hm = state.hsaccum + state.newsd + state.hi          # m

    code5, phie, fsup, firup = _new_profile(state, fm, code5, pdens, oldsd, simelt, phie)

    if state.error_type == 3:
        code3 = 1

    met[iw, MetCol.TSOIL] = state.stt[state.nnodes]            # sobrescreve o met (ponto 15)
    state.surfemis[iw] = state.emis

    if state.veg_flagl == 1:
        state.ft[iw] = state.ftemp
        state.surfemisf[iw] = state.epf
    else:
        state.ft[iw] = state.mflag
        state.surfemisf[iw] = state.mflag
    state.tt[iw] = state.toptemp

    state.vsmelt = 0.0
    state.vimelt = 0.0
    state.refreeze = 0.0
    state.refreezei = 0.0
    return code3, code5, pdens, phie, fsup, firup


def _surface_averages(state, surfdepth):
    """STEP 6 (linhas 351-440): umidade/gelo medios da camada de superficie."""
    iw = state.iw
    nn = state.nnodes
    elev, nz, nzi, ns = state.elev, state.nz, state.nzi, state.nsoilp
    sm, ice = state.soil_moist, state.ice
    surfm = surfi = surfmp = surfip = surfdd = denom = 0.0
    num2 = 0.0
    sfdp = surfdepth + 1e-4

    if state.water_flag == 0:
        j = nn
        while j >= 2 and _q5(elev - nzi[j]) <= sfdp:
            if j == nn:
                num2 = 5e-1 * (elev - nzi[j - 1])
            elif _q5(elev - nzi[j]) <= sfdp and _q5(elev - nzi[j - 1]) <= sfdp:
                num2 = 5e-1 * (nz[j + 1] - nzi[j - 1])
            elif _aint(elev - nzi[j], 5) <= sfdp and _q5(elev - nzi[j - 1]) > sfdp:   # sic: aint
                if _q5(elev - nzi[j]) - sfdp <= EPS:
                    num2 = (elev - nzi[j]) - denom
                else:
                    num2 = 5e-1 * ((elev - nzi[j]) - surfdepth)
            # sem else: `num2` fica com o valor da volta anterior
            denom = denom + num2
            t1 = ns[j, 1]
            surfdd = surfdd + ns[j, 1] * num2

            surfm = surfm + sm[j] * num2
            if state.ntype[j] != 27 and abs(t1) > EPS:
                surfmp = surfmp + (sm[j] / t1) * num2
            else:
                surfmp = surfm

            if ice[j] > EPS:
                dwi = _dwi(state, j)
                surfi = surfi + ice[j] * num2 * dwi
                if state.ntype[j] != 27 and abs(t1) > EPS:
                    surfip = surfip + (ice[j] / t1) * num2 * dwi
                else:
                    surfip = surfi

            j = j - 1
    else:
        denom = elev - nz[nn - 1]
        surfdd = ns[nn, 1] * denom

        surfm = 5e-1 * (sm[nn] + sm[nn - 1]) * denom
        if state.ntype[nn] != 27:
            surfmp = 5e-1 * (sm[nn] + sm[nn - 1]) * denom / ns[nn, 1]
        else:
            surfmp = surfm

        if ice[nn] > EPS or ice[nn - 1] > EPS:
            dwi = _dwi(state, nn)
            dwi1 = _dwi(state, nn - 1)
            surfi = 5e-1 * (ice[nn] * dwi + ice[nn - 1] * dwi1) * denom
            if state.ntype[nn] != 27:
                surfip = 5e-1 * (ice[nn] * dwi + ice[nn - 1] * dwi1) * denom / ns[nn, 1]
            else:
                surfip = surfi

    if denom > EPS:
        state.surfmoist[iw] = surfm / denom
        state.surfice[iw] = surfi / denom
        state.surfmoistp[iw] = surfmp / denom
        state.surficep[iw] = surfip / denom
        state.surfd[iw] = (surfdd / denom) * 1e3                        # kg/m^3
    else:
        state.surfmoist[iw] = sm[nn]
        state.surfice[iw] = ice[nn]
        state.surfmoistp[iw] = sm[nn] / ns[nn, 2]
        state.surficep[iw] = ice[nn] / ns[nn, 2]
        state.surfd[iw] = ns[nn, 1] * 1e3                               # kg/m^3


def _pfrozen(state, frozen, k, dwi):
    return float(frozen) * 1e-2 * (state.ice[k] * dwi + state.soil_moist[k])


def _write_nodes(state, frozen):
    """STEP 6 (linhas 443-465): saida nodal na unidade 3, `format 1001`
    `(5(i5),3x,a2,1x,2(f9.3),4(e12.4),3x,a1)`."""
    iw = state.iw
    nn = state.nnodes
    met = state.met
    for j in range(1, nn + 1):
        k = nn - j + 1
        dwi = _dwi(state, k)
        pfrozen = _pfrozen(state, frozen, k, dwi)
        if state.ice[k] * dwi >= pfrozen and state.ice[k] > EPS:
            nsstate = "f"
        else:
            nsstate = "u"
            if state.ntype[k] > 19 and state.stt[k] <= TREF:
                nsstate = "f"
        state.unit3.write(
            _fmt_i(int(met[iw, MetCol.YEAR]), 5) + _fmt_i(int(met[iw, MetCol.DOY]), 5)
            + _fmt_i(int(met[iw, MetCol.HR]), 5) + _fmt_i(int(met[iw, MetCol.MIN]), 5)
            + _fmt_i(j, 5) + "   " + str(state.node_type[k]).ljust(2)[:2] + " "
            + _fmt_f(abs(state.elev - state.nz[k]), 9, 3) + _fmt_f(state.stt[k], 9, 3)
            + _fmt_e(state.soil_moist[k], 12, 4) + _fmt_e(state.ice[k], 12, 4)
            + _fmt_e(state.soil_moist[k] + state.ice[k], 12, 4) + _fmt_e(state.wvc[k], 12, 4)
            + "   " + nsstate + "\n")


def _freeze_thaw(state, frozen):
    """STEP 6 (linhas 467-562): profundidade de congelamento/degelo e estado."""
    iw = state.iw
    nn = state.nnodes
    elev, nz, ice = state.elev, state.nz, state.ice
    state.sstate[iw] = "u"
    state.frthick[iw] = elev - nz[nn]
    state.twthick[iw] = elev - nz[nn]
    dwi = _dwi(state, nn)
    pfrozen = _pfrozen(state, frozen, nn, dwi)

    fts = 0
    if ice[nn] * dwi >= pfrozen and ice[nn] > EPS:
        state.sstate[iw] = "f"
        for k in range(nn - 1, 0, -1):
            dwi = _dwi(state, k)
            pfrozen = _pfrozen(state, frozen, k, dwi)

            if (ice[k] * dwi >= pfrozen and ice[k] > EPS) and fts == 0:
                if k >= 2:
                    state.frthick[iw] = (nz[k] + nz[k - 1]) * 5e-1
                    state.frthick[iw] = elev - state.frthick[iw]
                elif k == 1:
                    state.twthick[iw] = 0.0
                    state.frthick[iw] = elev - nz[k]
            else:
                if k == nn - 1:
                    state.frthick[iw] = elev - (nz[k] + nz[k + 1]) * 5e-1
                fts = 1
                j = k
                while j >= 2 and (ice[j] < pfrozen or abs(ice[j]) <= EPS):   # sic: sem dwi
                    state.twthick[iw] = (nz[j] + nz[j - 1]) * 5e-1
                    state.twthick[iw] = elev - state.twthick[iw]
                    j = j - 1
                    dwi = _dwi(state, j)
                    pfrozen = _pfrozen(state, frozen, j, dwi)

                if j == 1 and (ice[j] * dwi < pfrozen or ice[j] <= EPS):
                    state.twthick[iw] = elev - nz[1]
                break
    else:
        if iw >= 2:
            if state.sstate[iw - 1] == "f":
                state.sstate[iw] = "t"
            if state.sstate[iw - 1] == "t":           # o anterior, nao o recem-escrito: f -> t -> u
                state.sstate[iw] = "u"

        state.twthick[iw] = (nz[nn] + nz[nn - 1]) * 5e-1
        state.twthick[iw] = elev - state.twthick[iw]

        for k in range(nn - 1, 0, -1):
            dwi = _dwi(state, k)
            pfrozen = _pfrozen(state, frozen, k, dwi)

            if (ice[k] * dwi < pfrozen or abs(ice[k]) <= EPS) and fts == 0:
                if k >= 2:
                    state.twthick[iw] = (nz[k] + nz[k - 1]) * 5e-1
                    state.twthick[iw] = elev - state.twthick[iw]
                elif k <= 1:
                    state.frthick[iw] = 0.0
                    state.twthick[iw] = elev - nz[k]
            else:
                if k == nn - 1:
                    state.twthick[iw] = elev - (nz[k] + nz[k + 1]) * 5e-1
                fts = 1
                j = k
                while j >= 2 and (ice[j] >= pfrozen and ice[j] > EPS):   # sic: sem dwi
                    state.frthick[iw] = (nz[j] + nz[j - 1]) * 5e-1
                    state.frthick[iw] = elev - state.frthick[iw]
                    j = j - 1
                    dwi = _dwi(state, j)
                    pfrozen = _pfrozen(state, frozen, j, dwi)
                if j == 1 and ice[j] * dwi >= pfrozen:
                    state.frthick[iw] = elev - nz[1]
                break
    if str(state.node_type[nn]).strip() == "WA":
        state.met[iw, MetCol.HI] = state.frthick[iw]


def step6(state: FasstState, nprint, frozen, sdensi, surfdepth, phie):
    """STEP 6 (linhas 351-581): medias de superficie, saida nodal,
    congelamento/degelo e arquivo de reinicio. Devolve `phie`."""
    iw = state.iw
    _surface_averages(state, surfdepth)
    if nprint == 1:
        _write_nodes(state, frozen)
    _freeze_thaw(state, frozen)

    state.tot_thick[iw] = _anint(state.elev - state.nz[1], 10)

    state.error_type = 0

    wstart = max(state.iend - (MOVERLAP - 1), state.istart)
    if iw >= wstart:
        # `timeo(2,4)` aqui contra `timeo(moverlap,4)` em read_old_data: so o
        # modo 3 escreve nele, entao a declaracao menor e inofensiva no modo 2.
        timeo = [[0.0] * 5 for _ in range(MOVERLAP + 1)]
        _, phie = read_old_data(state, 2, iw, wstart, sdensi, phie, timeo)
    return phie


def _print_start(iw, iend):
    # write(*,*)' Start of run',iw,'of',iend,'runs' -- list-directed do
    # gfortran: branco no inicio do registro, I12, e um branco antes de
    # texto que segue um numero (conferido no TDM-GCC).
    print("  Start of run%12d of%12d runs" % (iw, iend))


def fasst_main(state: FasstState, fm: FasstMainState, nprint, sprint, frozen,
               nstr_flag, sdensi, surfdepth, code2, code3, code4, code5,
               first_loop, phie, pdens, oldsd, oldhi):
    """`subroutine fasst_main` (fasst_main.F90, 586 linhas), um passo `iw`.

    `fm` leva os `save` (ver `FasstMainState`). Devolve os `intent(inout)`:
    `(code2, code3, code4, code5, first_loop, phie, pdens, oldsd, oldhi)`.
    """
    iw = state.iw
    if first_loop == 0:
        fm.first_time = 0                                   # first time through new profile routines
        fm.first_snow = 0                                   # first time through snow routines
        state.vsmelt = 0.0                                  # snow melt (m)
        state.vimelt = 0.0                                  # ice melt (m)
        first_loop = 1

    if iw - _anint(iw * 1e-2) * 1e2 == 0 and state.single_multi_flag == 0:
        _print_start(iw, state.iend)
    elif state.iend > 50 and (iw - _anint(iw * 2e-2) * 5e1 == 0
                              and state.single_multi_flag == 1):
        _print_start(iw, state.iend)

    step1(state)
    code2 = step2(state, fm, oldsd, code2)
    code3, code5, pdens, phie, fsup, firup = step3(state, fm, code3, code5, pdens, oldsd, phie)
    step3a(state, fsup, firup)
    step4(state, nstr_flag)
    phie, oldsd, oldhi, code2, code4 = step5(state, fm, sprint, sdensi, phie, oldsd,
                                             oldhi, code2, code4)
    phie = step6(state, nprint, frozen, sdensi, surfdepth, phie)
    return code2, code3, code4, code5, first_loop, phie, pdens, oldsd, oldhi


def demo():
    import io
    import sys

    from .fasst_driver import allocate_run_arrays

    mod = sys.modules[__name__]

    def novo(top="SO", ice=0.0):
        st = FasstState()
        st.maxlines = 3
        allocate_run_arrays(st)
        st.iw = 2
        st.nnodes = 4
        st.node_type[4] = top
        st.ice[4] = ice
        st.mflag = -999.0
        st.met[2, MetCol.PT] = 1.0
        st.met[2, MetCol.PT2] = 1.0
        st.unit55 = io.StringIO()
        return st

    # ---- STEP 1 ----
    st = novo()
    st.met[2, MetCol.TSOL] = 5e-3                        # < 1d-2 vira 0, e upsol junto
    st.met[2, MetCol.UPSOL] = 40.0
    st.met[2, MetCol.TMP] = 1.0
    st.met[2, MetCol.RH] = 130.0
    st.met[2, MetCol.PREC], st.met[2, MetCol.PT] = -0.0, 3.0
    st.met[2, MetCol.TSOIL] = 7.0
    step1(st)
    assert st.dmet1[2, 1] == 0.0 and st.dmet1[2, 8] == 0.0
    def q15(x):                                          # quantizado a 1d-15: 100 -> 100.00000000000001
        return _anint(x * 1e15) * 1e-15
    assert st.dmet1[2, 4] == q15(1.0 + TREF) and st.dmet1[2, 5] == q15(100.0)
    assert st.met[2, MetCol.PREC] == 0.0 and st.met[2, MetCol.PT] == 1
    assert st.dmet1[2, 13] == q15(7.0) and st.sigfl == 0.0

    # ---- STEP 3 (new_profile postico) ----
    real_np = mod._new_profile
    mod._new_profile = lambda s, f, c5, pd, od, sm, ph: (c5, ph, 1.0, 2.0)
    try:
        st = novo()
        st.nnodes = 1
        st.met[2, MetCol.PT], st.met[2, MetCol.PT2] = 3.7, 1.0   # aint 3 entra, anint 4 -> dmet1(7)
        st.dmet1[2, 3], st.dmet1[2, 4] = 2.0, 268.0
        st.dmet1[2, 6], st.dmet1[2, 7] = 0.0, 1e-2
        st.stt[1] = 270.0
        r = step3(st, FasstMainState(), 0, 0, 0.0, 0.0, 0.0)
        pd = dense(268.0, 2.0, Phase.SNOW)
        assert r[2] == pd and r[4:] == (1.0, 2.0)
        assert st.newsd == _anint(1e-2 * 1e3 / pd * 1e10) * 1e-10
        assert st.dsnow == _anint(2e-2 * ((3.0 / (4.0 * math.pi)) * 5e-5) ** (1.0 / 3.0) * 1e10) * 1e-10
        assert r[3] == _anint((1.0 - pd * 1e-3) * 0.95 * 1e10) * 1e-10
        assert st.met[2, MetCol.TSOIL] == 270.0 and st.hm == st.newsd
    finally:
        mod._new_profile = real_np

    # ---- STEP 3a ----
    st = novo()
    st.met[2, MetCol.UPSOL] = -999.0
    st.met[2, MetCol.IRUP] = -999.0000000000001          # sentinela quantizado ainda conta como faltante
    st.met[2, MetCol.TSOL] = 500.0
    st.sup[2], st.irup[2] = 100.0, 300.0
    step3a(st, 5.0, 7.0)
    assert st.met[2, MetCol.UPSOL] == 105.0              # SOMA fsup (ponto 4)
    assert st.met[2, MetCol.IRUP] == 307.0
    step3a(st, 5.0, 7.0)                                 # valor presente: nao mexe
    assert st.met[2, MetCol.UPSOL] == 105.0 and st.met[2, MetCol.IRUP] == 307.0
    st.met[2, MetCol.TSOL] = 0.0                         # noite: upsol zera mesmo presente
    step3a(st, 5.0, 7.0)
    assert st.met[2, MetCol.UPSOL] == 0.0 and st.met[2, MetCol.IRUP] == 307.0

    # ---- postiços para STEP 4/5 (conferem roteamento, nao a fisica) ----
    real = (mod.soil_strength, mod.snow_model)
    chamadas = []

    def _ss(state, nstr_flag):
        chamadas.append("ss")
        return 1.0, 2.0, 3.0

    def _sm(state, snow, sprint, water_flag, first_snow, sdensi, sdf):
        chamadas.append(("sm", first_snow, sdf))
        state.vsmelt = 0.5                               # m, como snow.py:477
        state.error_type = 4
        return 1, 0.25, 0.4

    mod.soil_strength, mod.snow_model = _ss, _sm
    try:
        # ---- STEP 4 ----
        st = novo()
        step4(st, [0] * 5)
        assert (st.surfci[2], st.surfrci[2], st.surfcbr[2]) == (1.0, 2.0, 3.0)
        st = novo()
        st.water_flag = 1
        chamadas.clear()
        step4(st, [0] * 5)
        assert chamadas == [] and st.surfci[2] == 0.0    # ponto 7

        fm = FasstMainState()

        # ---- STEP 5: terra com neve -> snow_model, vsmelt NAO convertido ----
        st = novo()
        st.hsaccum = 0.1
        r = step5(st, fm, 1, 200.0, 0.0, 9.0, 9.0, 0, 0)
        assert r == (0.4, 0.25, 9.0, 0, 1)               # hi=0, pt=1: bloco de gelo nao roda
        assert st.vsmelt == 0.5 and fm.first_snow == 1
        assert st.met[2, MetCol.SD] == 0.25 and st.met[2, MetCol.HI] == 0.0
        assert st.unit55.getvalue() == ""

        # ---- agua com gelo -> snow_model, vsmelt*1d-2 (ponto 1) ----
        st = novo("WA", 0.5)
        st.met[2, MetCol.PT] = 3.0
        r = step5(st, fm, 1, 200.0, 0.0, 9.0, 9.0, 0, 0)
        assert st.vsmelt == 0.005 and r[0] == 0.4 and r[4] == 1

        # ---- agua sem gelo / ar -> zera e escreve o format 91 local ----
        for top, ice in (("WA", 0.0), ("AI", 0.3)):
            st = novo(top, ice)
            st.hsaccum = 0.1
            st.vsmelt = 0.7
            st.met[2, MetCol.YEAR] = 2024.0
            st.met[2, MetCol.DOY] = 10.0
            st.met[2, MetCol.HR] = 12.0
            st.met[2, MetCol.MIN] = 0.0
            r = step5(st, fm, 1, 200.0, 0.123, 9.0, 9.0, 0, 0)
            assert r == (0.123, 0.0, 9.0, 0, 0)
            assert st.vsmelt == 0.0 and st.hsaccum == 0.0 and st.meltfl == "0"
            assert st.unit55.getvalue() == (
                "      2024" + "  " + "    10.500" + "     0.00000000" * 8 + "    0"
                + "     0.00000000" * 3 + "\n")

        # ---- sem neve: else externo; sprint=0 nao escreve ----
        st = novo()
        st.vsmelt = 0.7
        step5(st, fm, 0, 200.0, 0.0, 9.0, 9.0, 0, 0)
        assert st.vsmelt == 0.0 and st.unit55.getvalue() == ""

        # ---- 'WA' com ice NaN: nenhum ramo (ponto 6) ----
        st = novo("WA", math.nan)
        st.hsaccum = 0.1
        st.vsmelt = 0.7
        chamadas.clear()
        r = step5(st, fm, 1, 200.0, 0.0, 9.0, 9.0, 0, 0)
        assert chamadas == [] and st.vsmelt == 0.7 and r[1] == 0.1
    finally:
        mod.soil_strength, mod.snow_model = real

    # ---- gelo de superficie (icethick REAL) ----
    fm = FasstMainState()
    st = novo()
    st.met[2, MetCol.PT] = 2.0                           # chuva
    st.met[2, MetCol.PREC] = 1.0
    st.stt[4] = 260.0
    st.toptemp = 260.0
    st.melt[2] = -50.0
    st.lheat[2] = -50.0
    st.timstep = 1.0
    st.error_type = 2                                    # deixado pelo STEP 2 (ponto 3)
    r = step5(st, fm, 0, 200.0, 0.0, 0.0, 0.0, 0, 0)
    assert st.hi > 0.0 and r[2] == st.hi and r[3] == 1
    assert st.met[2, MetCol.HI] == st.hi

    # com neve: icethick nao roda, mas oldhi = hi roda (ponto 5)
    st = novo()
    st.hi = 0.02
    st.met[2, MetCol.PT2] = 3.0                          # neve so em pt2 -> ramo de neve
    mod.snow_model = lambda s, sn, sp, wf, fs, sd, sdf: (1, 0.3, 0.0)
    try:
        r = step5(st, fm, 0, 200.0, 0.0, 0.0, 0.0, 0, 0)
    finally:
        mod.snow_model = real[1]
    assert st.hsaccum == 0.3 and st.hi == 0.02 and r[2] == 0.02

    # pt2 = 2 sozinho nao liga o bloco de gelo (ponto 9); 'WA' tambem nao
    for top, pt, pt2 in (("SO", 1.0, 2.0), ("WA", 2.0, 1.0)):
        st = novo(top)
        st.met[2, MetCol.PT] = pt
        st.met[2, MetCol.PT2] = pt2
        st.met[2, MetCol.PREC] = 1.0
        st.timstep = 1.0
        r = step5(st, fm, 0, 200.0, 0.0, 0.0, 0.07, 0, 0)
        assert st.hi == 0.0 and r[2] == 0.07

    # ---- STEP 6: format 1001 da unidade 3 ----
    # gabarito: linha que o FASST compilado escreveu no caso `grama` de
    # cmp_fortran_agua.py (saida.grd), refeita aqui com os mesmos valores
    st = novo()
    st.nnodes = 2
    st.met[2, MetCol.YEAR], st.met[2, MetCol.DOY], st.met[2, MetCol.HR] = 2010.0, 100.0, 11.0
    st.elev = 1.0
    st.nz[1], st.nz[2] = 0.985, 1.0
    st.node_type[1], st.node_type[2] = "SM", "SM"
    st.stt[1], st.stt[2] = 38.696, 38.0
    st.soil_moist[1], st.ice[1], st.wvc[1] = 0.0, 0.1337, 0.0
    st.ntype[1] = st.ntype[2] = 7
    st.unit3 = io.StringIO()
    _write_nodes(st, 100)
    assert st.unit3.getvalue().splitlines()[1] == (
        " 2010  100   11    0    2   SM     0.015   38.696  0.0000E+00  0.1337E+00"
        "  0.1337E+00  0.0000E+00   f")

    print("demo OK")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
