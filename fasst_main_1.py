"""
fasst_main.py -- traducao PARCIAL de fasst_main.F90 (ate o STEP 3)

Fonte Fortran: fasst_main.F90 (586 linhas), sub-rotina `fasst_main`.

Escopo
------
A pedido, so as linhas 1-271 foram traduzidas: primeira chamada
(`first_loop`), STEP 1 (meteorologia derivada `dmet1` e propriedades da
vegetacao baixa), STEP 2 (temperatura e umidade do dossel) e STEP 3 (neve
nova, `new_profile` e variaveis de superficie). NAO foram traduzidos o
STEP 3a (fluxos solar e infravermelho ascendentes, linhas 274-281), o
STEP 4 (resistencia do solo), o STEP 5 (neve e gelo), o STEP 6 (umidade e
gelo medios da superficie, saida dos nos, espessuras de congelamento e
degelo, reinicio) e o `end subroutine`. A funcao termina logo depois do
STEP 3 e devolve os argumentos `intent(inout)` no estado em que estao
nesse ponto.

`fasst_main` e o "casco" de um passo de tempo: integra os modelos de fluxo
solar e infravermelho, derretimento de neve, temperatura do solo, umidade e
resistencia do solo e congelamento/degelo (cabecalho, linhas 22-43).

Assinatura original:
    subroutine fasst_main(nprint,sprint,frozen,nstr_flag,sdensi,surfdepth,
                          code2,code3,code4,code5,first_loop,phie,pdens,
                          oldsd,oldhi)
    integer(ip),intent(in):: nprint,sprint,frozen
    integer(ip),intent(in):: nstr_flag(maxn)
    real(dp),intent(in):: sdensi,surfdepth
    integer(ip),intent(inout):: code2,code3,code4,code5,first_loop
    real(dp),intent(inout):: phie,pdens,oldsd,oldhi

Arquitetura
-----------
`fasst_main(state, fm, nprint, sprint, frozen, nstr_flag, sdensi,
surfdepth, code2, code3, code4, code5, first_loop, phie, pdens, oldsd,
oldhi)` recebe `state: FasstState` primeiro e um `FasstMainState` (`fm`) em
seguida. Os argumentos `intent(inout)` voltam numa tupla, na ordem da lista
de argumentos do Fortran:
    code2, code3, code4, code5, first_loop, phie, pdens, oldsd, oldhi = \\
        fasst_main(...)
(`code4` e `oldhi` so mudam no STEP 5, ainda nao traduzido; seguem na
tupla para a assinatura nao mudar quando o resto for traduzido.)

`FasstMainState` reune as variaveis `save` que persistem entre chamadas:
  * `first_snow` e `first_time` (as duas `save::` explicitas de
    `fasst_main`, linha 118); `first_snow` so e zerada aqui e vai a
    `snow_model` no STEP 5; `first_time` vai a `new_profile`;
  * os objetos de estado salvo do dossel chamado no STEP 2:
    `CanopySavedState` (module_canopy) e `AlbedoEmisState` (albedo_emis).
Instancie um `FasstMainState` por execucao e reutilize-o em todas as
chamadas. No Fortran as variaveis `save` de module_canopy existem em uma
unica copia, compartilhada com `open_water`; aqui cada rotina-raiz tem a
sua. Para reproduzir o compartilhamento, crie os objetos uma vez no driver
e injete os mesmos `canopy`/`alb` em `OpenWaterState` e `FasstMainState`.

Constantes (`EPS`, `PI`, `TREF`, `NCLAYERS`) vem de `fasst.constants`;
colunas de `met` sao `MetCol.*`. As 13 colunas de `dmet1` nao tem enum e
sao usadas por numero, como nas demais traducoes (significado de cada uma
na secao "dmet1" do documento de transcricao). Indexacao 1-based identica
ao Fortran (posicao 0 sem uso).

Dependencias externas
---------------------
  canopy_met        fasst.module_canopy   (real)
  low_veg_prop      fasst.module_lowveg   (real)
  dense, Phase      fasst.functions       (reais; d1i = 4 -> Phase.SNOW)
  new_profile       STUB neste arquivo    (INTERFACE ASSUMIDA -- abaixo)

`new_profile` (new_profile.F90) nao foi fornecido. Como em `sort2`
(read_complex.py), ha um stub que levanta `NotImplementedError`, e o
adaptador `_new_profile` (unico ponto de contato) assume
    new_profile(state, code5, first_time, pdens, oldsd, simelt, phie)
        -> (code5, first_time, pdens, oldsd, phie, fsup, firup)
isto e, os seis primeiros argumentos do `call new_profile(...)` do Fortran
(com `state` na frente) e as duas saidas `fsup`/`firup`. Como o Fortran
passa tudo por referencia e nao se sabe qual `intent` cada um tem, todos os
argumentos que o original poderia alterar voltam na tupla (`simelt` nao:
nao e lido depois da chamada). Quando new_profile for traduzida, apague o
stub, importe a funcao real e ajuste so o adaptador.

Pontos de atencao (comportamento do original preservado)
---------------------------------------------------------
1. `dmet1(iw,5) = dmin1(met(iw,ip_rh),1d2)` so limita a umidade a 100%. A
   linha comentada logo acima (linha 165, `dmax1(1d-1, rh)`) e a usada em
   open_water; aqui nao ha piso, entao umidade 0 ou negativa passa.
2. `met` e alterado no proprio lugar: precipitacao <= eps vira 0 com tipo
   1 (linhas 166-179); depois de canopy_met, `met(iw,ip_pt2) = 3` se
   houver neve (linha 204); e `met(iw,ip_tsoil)` e SOBRESCRITO com a
   temperatura do no superior (linha 256). `dmet1(iw,13)` guarda a
   temperatura original.
3. A condicao de neve nova usa `aint(met(iw,ip_pt)) == 3` (linha 223), mas
   a escolha entre chuva+neve (`dmet1(iw,6)`) e so neve (`dmet1(iw,7)`)
   usa `anint` (linha 228). Para tipo 3,7, por exemplo, entra no bloco
   (aint = 3) e escolhe `dmet1(iw,7)` (anint = 4). Preservado.
4. O bloco das linhas 239-244 (zerar `dsnow`, `phie`, `newsd` e `pdens`
   quando nao ha neve acumulada e `newsd < 2d-3`) esta COMENTADO em
   fasst_main; em open_water esta ativo. Nao traduzido.
5. `dsnow` nao depende de nenhuma entrada: vale sempre
   2d-2*((3/(4*pi))*5d-5)**(1/3) arredondado a 10 casas (~4,571d-4).
6. `canopy_met` recebe `oldsd`, a profundidade de neve do passo ANTERIOR
   (`oldsd = hsaccum` so na linha 337, STEP 5).
7. O arredondamento a 15 casas de `dmet1` (linhas 189-191) acontece antes
   de `canopy_met`, que altera `dmet1` sem arredondar de novo.
8. `pdens` e zerado a cada chamada; so e diferente de zero com neve nova.
9. O cabecalho (linhas 97-105) cita low_veg_prop, canopy_met, new_profile,
   soil_strength, snow_model e icethick e as funcoes dense e met_date.
   `met_date` e declarada (linha 122) e nunca usada; `read_old_data`
   (linha 580) e chamada e nao consta.
10. Argumentos nao usados na parte traduzida: nprint, sprint, frozen,
   nstr_flag, sdensi, surfdepth, code4, oldhi (entram a partir do STEP 3a).
   Dos locais declarados nas linhas 121-126 so i, d1i, simelt, smm, fsup e
   firup sao usados aqui; os demais (j, k, fts, d2i, wstart, surfm, surfi,
   surfmp, surfip, denom, surfdd, num2, dates, xtemp, pfrozen, timeo, t1,
   sfdp, dwi, dwi1, nsstate e met_date) pertencem aos STEPs 3a-6 e nao foram
   declarados nem inicializados (linhas 130-138).
11. `write(*,*)' Start of run',iw,'of',iend,'runs'` virou `print`; o
   espacamento de campos da saida dirigida por lista do gfortran (inteiro em
   12 colunas) nao e reproduzido.
12. Divisao por zero em `newsd` (pdens = 0) segue IEEE (inf ou NaN), como no
   Fortran, via `_fdiv`.
"""

import math
from dataclasses import dataclass, field

from .albedo_emis import AlbedoEmisState
from .constants import EPS, NCLAYERS, PI, TREF
from .functions import Phase, dense
from .indices import MetCol
from .module_canopy import CanopySavedState, canopy_met
from .module_lowveg import low_veg_prop
from .state import FasstState

__all__ = ["FasstMainState", "fasst_main"]


@dataclass
class FasstMainState:
    """Variaveis `save` de `fasst_main` e dos modulos que ela chama."""

    first_snow: int = 0
    first_time: int = 0
    canopy: CanopySavedState = field(default_factory=CanopySavedState)
    alb: AlbedoEmisState = field(default_factory=AlbedoEmisState)


def _aint(x):
    """Fortran ``AINT(x)`` -- trunca em direcao a zero (NaN/inf passam)."""
    if not math.isfinite(x):
        return x
    return math.copysign(float(math.trunc(x)), x)


def _anint(x):
    """Fortran ``ANINT(x)`` -- arredonda para longe de zero.

    Usa `floor` e a parte fracionaria (exata) em vez de `floor(|x|+0.5)`,
    que erra para 0.49999999999999994.
    """
    if not math.isfinite(x):
        return x
    a = abs(x)
    r = math.floor(a)
    if a - r >= 0.5:
        r += 1.0
    return math.copysign(r, x)


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


# ---------------------------------------------------------------------------
# Dependencia externa nao fornecida
# ---------------------------------------------------------------------------
def new_profile(state, code5, first_time, pdens, oldsd, simelt, phie):
    """TODO: `new_profile` nao esta em nenhum arquivo recebido (ver docstring)."""
    raise NotImplementedError(
        "new_profile nao foi fornecida (fasst_main.F90, linha 251)")


def _new_profile(state, fm, code5, pdens, oldsd, simelt, phie):
    """Adaptador de ``call new_profile(code5,first_time,pdens,oldsd,simelt,phie,fsup,firup)``.

    Interface assumida -- ver docstring do modulo.
    """
    code5, first_time, pdens, oldsd, phie, fsup, firup = new_profile(
        state, code5, fm.first_time, pdens, oldsd, simelt, phie)
    fm.first_time = first_time
    return code5, pdens, oldsd, phie, fsup, firup


def fasst_main(state: FasstState, fm: FasstMainState, nprint, sprint, frozen, nstr_flag,
               sdensi, surfdepth, code2, code3, code4, code5, first_loop, phie, pdens,
               oldsd, oldhi):
    iw = state.iw

    # initialize certain check/test variables (linhas 129-138): as de uso
    # exclusivo do STEP 3a em diante (nao traduzido) nao foram declaradas;
    # d1i e simelt recebem valor antes de serem lidas (linhas 225 e 248).

    if first_loop == 0:
        fm.first_time = 0                       # first time through new profile routines
        fm.first_snow = 0                       # first time through snow routines
        state.vsmelt = 0.0                      # snow melt (m)
        state.vimelt = 0.0                      # ice melt (m)
        first_loop = 1

    # ---------------- Main Calculation Loop (linhas 148-155) ----------------
    if iw - _anint(iw * 1e-2) * 1e2 == 0 and state.single_multi_flag == 0:
        print(f" Start of run {iw} of {state.iend} runs")
    elif state.iend > 50 and (iw - _anint(iw * 2e-2) * 5e1 == 0
                              and state.single_multi_flag == 1):
        print(f" Start of run {iw} of {state.iend} runs")

    # ---------------- STEP 1 (linhas 157-198; sem marcador no codigo) ----------------
    # initialize the met calc variable dmet1(iw,)
    state.dmet1[iw, 1] = state.met[iw, MetCol.TSOL]
    if state.dmet1[iw, 1] < 1e-2:
        state.dmet1[iw, 1] = 0.0
    state.dmet1[iw, 2] = state.met[iw, MetCol.IR]
    state.dmet1[iw, 3] = state.met[iw, MetCol.WS]
    state.dmet1[iw, 4] = state.met[iw, MetCol.TMP] + TREF                      # K

    # dmin1 -> min: com NaN o resultado do Fortran depende do compilador (gfortran 13:
    # -O0 propaga o NaN, como o min do Python; -O2 devolve 100)
    state.dmet1[iw, 5] = min(state.met[iw, MetCol.RH], 1e2)                    # % (ponto de atencao 1)
    if state.met[iw, MetCol.PREC] <= EPS:
        state.met[iw, MetCol.PREC] = 0.0
        state.met[iw, MetCol.PT] = 1

    state.dmet1[iw, 6] = state.met[iw, MetCol.PREC] * 1e-3                     # m/timstep (rain and/or snow)
    if state.dmet1[iw, 6] <= EPS:
        state.met[iw, MetCol.PT] = 1.0

    if state.met[iw, MetCol.PREC2] <= EPS:
        state.met[iw, MetCol.PREC2] = 0.0
        state.met[iw, MetCol.PT2] = 1
    state.dmet1[iw, 7] = state.met[iw, MetCol.PREC2] * 1e-3                    # m/timstep (snow only)
    if state.dmet1[iw, 7] <= EPS:
        state.met[iw, MetCol.PT2] = 1.0

    state.dmet1[iw, 8] = state.met[iw, MetCol.UPSOL]
    if abs(state.dmet1[iw, 1]) <= EPS:
        state.dmet1[iw, 8] = 0.0
    state.dmet1[iw, 9] = state.met[iw, MetCol.DIR]
    state.dmet1[iw, 10] = state.met[iw, MetCol.DIF]
    state.dmet1[iw, 11] = state.met[iw, MetCol.AP]
    state.dmet1[iw, 12] = state.met[iw, MetCol.IRUP]
    state.dmet1[iw, 13] = state.met[iw, MetCol.TSOIL]

    for i in range(1, 14):
        state.dmet1[iw, i] = _anint(state.dmet1[iw, i] * 1e15) * 1e-15

    # get the low vegetation properties
    if state.veg_flagl == 1:
        low_veg_prop(state, state.dmet1[iw, 4], state.dzveg[NCLAYERS])
    else:
        state.sigfl = 0.0

    # ---------------- STEP 2 (linhas 200-216) ----------------
    # calculate the canopy temperature and moisture
    if state.veg_flagh == 1:
        canopy_met(state, fm.canopy, fm.alb, oldsd)                 # correct met for canopy
        if state.dmet1[iw, 7] > 0.0:
            state.met[iw, MetCol.PT2] = 3.0
        if state.dmet1[iw, 6] <= EPS:
            state.met[iw, MetCol.PT] = 1.0
        if state.dmet1[iw, 7] <= EPS:
            state.met[iw, MetCol.PT2] = 1.0
    else:
        state.surfemisc[iw] = state.mflag
        for i in range(1, NCLAYERS + 1):
            state.canopy_temp[i, iw] = state.mflag
    if state.error_type == 2:
        code2 = 1

    state.airt[1, iw] = state.met[iw, MetCol.TMP] + TREF
    state.airt[2, iw] = state.dmet1[iw, 4]

    # ---------------- STEP 3 (linhas 218-271) ----------------
    # determine density and amount of new snow if no existing snow on ground
    state.newsd = 0.0
    pdens = 0.0
    smm = state.dmet1[iw, 6] + state.dmet1[iw, 7]
    if ((_aint(state.met[iw, MetCol.PT]) == 3 or _aint(state.met[iw, MetCol.PT2]) == 3)
            and smm > 0.0):
        d1i = 4
        pdens = dense(state.dmet1[iw, 4], state.dmet1[iw, 3], Phase.SNOW)      # kg/m^3 (d1i = 4)

        if _anint(state.met[iw, MetCol.PT]) == 3:                  # ponto de atencao 3
            state.newsd = _fdiv(state.dmet1[iw, 6] * 1e3, pdens)   # m snow
        else:
            state.newsd = _fdiv(state.dmet1[iw, 7] * 1e3, pdens)
        state.newsd = _anint(state.newsd * 1e10) * 1e-10
        state.dsnow = 2e-2 * ((3.0 / (4.0 * PI)) * 5e-5) ** (1.0 / 3.0)
        state.dsnow = _anint(state.dsnow * 1e10) * 1e-10
        phie = (1.0 - pdens * 1e-3) * 0.95
        phie = _anint(phie * 1e10) * 1e-10
        # linhas 239-244 comentadas no original (ponto de atencao 4)

    # calculate the new moisture and temperature profiles
    simelt = state.vsmelt + state.vimelt
    state.hm = state.hsaccum + state.newsd + state.hi                          # m

    code5, pdens, oldsd, phie, fsup, firup = _new_profile(
        state, fm, code5, pdens, oldsd, simelt, phie)

    if state.error_type == 3:
        code3 = 1

    state.met[iw, MetCol.TSOIL] = state.stt[state.nnodes]
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

    # ---------------- fim da parte traduzida (fasst_main.F90, linha 271) ----------------
    # STEP 3a em diante (linhas 274-584) NAO traduzido, a pedido. Quando for
    # traduzido, continua aqui: `fsup` e `firup` (saidas de new_profile) sao
    # usados no STEP 3a; os argumentos nprint, sprint, frozen, nstr_flag,
    # sdensi, surfdepth, code4 e oldhi, a partir dos STEPs 4-6.
    return code2, code3, code4, code5, first_loop, phie, pdens, oldsd, oldhi
