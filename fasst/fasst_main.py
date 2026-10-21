"""
fasst_main.py -- traducao de fasst_main.F90 (PARCIAL: STEPS 3a, 4 e 5)

Fonte Fortran: fasst_main.F90 (586 linhas). Uma unica sub-rotina,
`fasst_main`, chamada pelo driver uma vez por passo de tempo `iw` no
caminho de TERRA (`water_flag <= 2`). E a casca que integra, em ordem, os
modelos de fluxo solar/IR, temperatura e umidade do solo, resistencia do
solo, neve e gelo, organizada pelos marcadores STEP do proprio original:

  STEP 1   (148-198)  met do passo (`dmet1`), propriedades da vegetacao baixa
  STEP 2   (200-216)  temperatura/umidade do dossel (`canopy_met`)
  STEP 3   (218-271)  neve nova, `new_profile` (perfis de temperatura e
                      umidade do solo), emissividades
  STEP 3a  (274-281)  fluxos solar e IR emitidos pela superficie  <-- aqui
  STEP 4   (284-287)  resistencia do solo (`soil_strength`)       <-- aqui
  STEP 5   (290-348)  modelo de neve e gelo de superficie         <-- aqui
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
PARCIAL. Transcritos os STEPS 3a, 4 e 5 (linhas 274-348); faltam os
STEPS 1, 2, 3 e 6 e o corpo de `fasst_main` que os encadeia. Ate la o
driver continua chamando o stub `fasst_driver.fasst_main` (que levanta
`NotImplementedError`), e o caminho de terra ainda para nele.

Arquitetura
-------------
Uma funcao por STEP: `step3a`, `step4`, `step5`. Quando os outros STEPS
chegarem, o corpo de `fasst_main` so as chama em sequencia. Isto se
afasta da regra do driver ("um `program`, UM escopo", fasst_driver.py,
Arquitetura), e pode, porque a fronteira foi conferida: nenhuma local
atribuida em 274-348 (`dates`, `xtemp`) e lida em 351-586, e `d1i`/`d2i`
sao reatribuidas na 353 antes de qualquer leitura. As unicas locais que
ENTRAM vindas de fora sao `fsup`/`firup`, saidas de `new_profile` no
STEP 3 (`NewProfileResult.fsup/.firup`), e viram parametros de `step3a`.

Os `intent(inout)` que cada STEP toca voltam numa tupla, como no resto do
repo: `phie, oldsd, oldhi, code2, code4 = step5(...)`. Campos de
`fasst_global` (hsaccum, hi, vsmelt, met, surfci, ...) sao mutados em
`state`.

`FasstMainState` reune os `save` de `fasst_main` (`first_snow`,
`first_time`, linha 118) e o `SnowState` de module_snow. Esse `snow`
TEM de ser o mesmo objeto de `OpenWaterState.snow`: no Fortran os `save`
de `module_snow` pertencem ao modulo, nao a quem chama, e sao os mesmos
para `open_water` e `fasst_main` (pendencia registrada em open_water.py).
Ja `first_snow` NAO e compartilhado: `open_water` tem o seu proprio
`save`. Os objetos de canopy/lowveg/albedo entram aqui quando os STEPS 1
e 2 forem transcritos.

Constantes (`EPS`) vem de `fasst.constants`; colunas de `met` sao
`MetCol.*`. Indexacao 1-based identica ao Fortran (posicao 0 sem uso).

Dependencias externas (todas reais)
-------------------------------------
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

Pontos de atencao (comportamento do original preservado)
---------------------------------------------------------
1. `vsmelt = vsmelt*1d-2` (linha 298) converte cm -> m, mas
   `module_snow.F90:302` (e `snow.py:477`) JA devolve `vsmelt` em metros.
   No ramo agua-com-gelo o valor sai 100x menor; no ramo de terra (316) nao
   ha conversao nenhuma. Defeito do original, preservado. `open_water.F90`
   faz a mesma conversao dupla.
2. STEP 5 testa o tipo de precipitacao com `aint` (293, 342); o STEP 3 mistura
   `aint` (223) e `anint` (228) sobre o mesmo campo. Para `pt` inteiro e
   igual; transcrito literalmente com `_aint`.
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
   `meltfl` (len=1) sai alinhado a direita em `a5`. Como `snow.py` omite a
   escrita da unidade 55 (snow.py:302), o `.snw` so recebe as linhas dos
   ramos sem neve -- pendencia de snow.py.
9. O bloco de gelo de superficie (342) testa so `ip_pt` (2 = chuva, 4 =
   chuva congelante), nunca `ip_pt2`.
"""

import math
from dataclasses import dataclass, field

from .constants import EPS
from .icethick import icethick
from .indices import MetCol
from .snow import SnowState, snow_model
from .soil_strength import soil_strength
from .state import FasstState

__all__ = ["FasstMainState", "step3a", "step4", "step5"]


@dataclass
class FasstMainState:
    """Variaveis `save` de `fasst_main` (linha 118) e o `save` de module_snow.

    `snow` deve ser o MESMO objeto de `OpenWaterState.snow` (ver docstring).
    """

    first_snow: int = 0
    first_time: int = 0
    snow: SnowState = field(default_factory=SnowState)


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

    print("demo OK")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
