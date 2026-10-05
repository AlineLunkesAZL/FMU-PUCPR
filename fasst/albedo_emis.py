"""
albedo_emis.py -- traducao de albedo_emis.F90

Fonte Fortran: albedo_emis.F90 (113 linhas). Uma unica sub-rotina
publica (`albedo_emis`), sem nenhuma sub-rotina chamada. Calcula o
albedo (albedo_fasst) e a emissividade (emis) da superficie no nó mais
superficial do perfil (nnodes) para o passo de tempo atual, a partir do
tipo de solo/água, da umidade, do ângulo zênite solar (para água), e do
estado de neve/gelo acumulado (sem neve/gelo, gelo sem neve, neve nova,
ou neve envelhecida -- com dois modelos empíricos de decaimento de
albedo combinados: Douville et al. 1995 e Roesch 2000). Quando há
radiação solar incidente/refletida medida (sd/su), o albedo medido tem
prioridade sobre todo o resto.

Assinatura original:
    subroutine albedo_emis(oldsd,sd,su)
    real(dp),intent(in):: oldsd,sd,su

Arquitetura
-------------
`albedo_emis` recebe `state: FasstState` como argumento explicito.
`EPS` e `TREF` vem de `fasst.constants`; `mflag` e `sigfl` sao campos
reais de `FasstState`. `ip_zen` (indice de coluna do meteorologico, no
Fortran) e `MetCol.ZEN`. Nao ha dependencias externas.

As duas saidas reais desta sub-rotina (`albedo_fasst`, `emis`) sao
campos de `FasstState`, atualizadas como efeito colateral -- nao ha
retorno em tupla.

`albedoo` e uma variavel local com o atributo `save` no Fortran (seu
valor persiste entre chamadas sucessivas da sub-rotina, diferente de
uma variavel local comum). Aqui e o campo `albedoo` do dataclass
:class:`AlbedoEmisState`, passado explicitamente e mutado in-place --
mesmo tratamento usado para as variaveis `save` de outros arquivos
deste projeto (`LowVegState`, `SnowState`, `CanopySavedState`). O valor
inicial de `albedoo` nao importa: pela logica do proprio Fortran,
`albedoo` so e lido no ramo de "neve envelhecida", e antes de qualquer
chamada chegar la pela primeira vez, ou o ramo de "neve nova" (que
sempre define `albedoo`) ou o proprio ramo de "neve envelhecida" com
`iw==1` (que tambem redefine `albedoo` antes de usa-lo) ja terao
rodado.

Pontos de atencao
-------------------
* `sgremis1` e inicializada com `nsoilp(nnodes,4)` logo no topo da
  sub-rotina (linha 30 do .F90), mas essa atribuicao e sempre
  sobrescrita mais adiante -- os tres ramos do bloco `if(ntype(nnodes)
  <= 19)/else if(...26)/else` cobrem 100% dos casos e cada um atribui
  `sgremis1` de novo, incondicionalmente. Preservado por fidelidade,
  mesmo sendo codigo morto (nunca influencia o resultado).
* `albedo1` (usada apenas quando `aflag==1`) so e de fato atribuida
  quando `aflag` tambem e atribuido para 1, no mesmo bloco -- ou seja,
  nunca e lida sem ter sido escrita antes nessa mesma chamada.
* PROVAVEL BUG DO FORTRAN, preservado sem correcao (fiel ao original):
  no ramo de agua (`ntype(nnodes)==26`), o albedo de Fresnel usa o
  fator `5d1` (50.0) -- `sgralbedo = 5d1*(sin²(Z-r)/sin²(Z+r) +
  tan²(Z-r)/tan²(Z+r))`. A formula classica de Fresnel (reflectancia
  nao polarizada, media das duas polarizacoes) usa o fator `5d-1`
  (0.5), nao `5d1`. Com o fator do Fortran, `sgralbedo` fica ~100x
  maior que o fisicamente esperado (Z=30°: 2.11 em vez de 0.0211) --
  ou seja, quase sempre bem maior que 1.0. Como `albedo_fasst` e
  inicializada em 1.0 e depois recebe `dmin1(albedo_fasst,sgralbedo)`,
  na pratica esse ramo nunca reduz o albedo da agua abaixo de 1.0 (o
  valor calculado e descartado pelo min). Mantido tal como no Fortran.
"""

import math
from dataclasses import dataclass

from .constants import EPS, IALBEDO, IEMIS, SEMIS, SNALBEDO, SOALBEDO, TREF
from .indices import MetCol
from .state import FasstState


@dataclass
class AlbedoEmisState:
    """Variavel `save` de `albedo_emis` -- ver docstring do modulo."""

    albedoo: float = 0.0


def _aint(x, p=0):
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direção a zero."""
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _anint(x, p=0):
    """Fortran ``ANINT(x*10**p)*10**-p`` -- arredonda para longe de zero."""
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def albedo_emis(state: FasstState, alb: AlbedoEmisState, oldsd, sd, su):
    # ---------------- inicializacao de variaveis locais (linhas 22-30) ----------------
    aflag = 0
    albedo1 = 0.0  # so lida quando aflag==1 (ver docstring)
    t1 = 0.0
    Z = 0.0
    r = 0.0
    f1 = 0.0
    salbedoo = 0.0
    salbedon = 0.0
    sgremis1 = 0.0
    sgremis1 = state.nsoilp[state.nnodes][4]  # sobrescrita nos 3 ramos abaixo -- ver docstring

    albedoo = alb.albedoo

    # ---------------- fator de correcao a partir de radiacao solar medida (linhas 32-42) ----------------
    state.albedo_fasst = 1.0
    if sd > EPS and _aint(abs(sd - state.mflag) * 1e5) * 1e-5 > EPS:
        if su > EPS and _aint(abs(su - state.mflag) * 1e5) * 1e-5 > EPS:
            albedo1 = (1.0 - state.sigfl) * su / sd
            if albedo1 > 1.0:
                albedo1 = 0.99
            aflag = 1
    elif sd < EPS and _aint(abs(sd - state.mflag) * 1e5) * 1e-5 > EPS:
        albedo1 = 0.0
        aflag = 0

    # ---------------- albedo/emissividade da superficie, sem neve/gelo (linhas 44-71) ----------------
    if state.ntype[state.nnodes] <= 19:
        # solo -- albedo funcao da umidade (decaimento exponencial)
        a = 0.0
        a = math.exp(-1.0) + 0.23865 * state.soil_moist[state.nnodes] / state.nsoilp[state.nnodes][9]
        state.sgralbedo = -state.nsoilp[state.nnodes][3] * math.log(a)
        state.sgralbedo = max(
            0.5 * state.nsoilp[state.nnodes][3],
            min(state.nsoilp[state.nnodes][3], state.sgralbedo),
        )

        # emissividade -- interpolacao linear com a umidade
        a = 0.0
        a = 0.99 - state.nsoilp[state.nnodes][4]
        sgremis1 = state.nsoilp[state.nnodes][4] + a * state.soil_moist[state.nnodes] / state.nsoilp[state.nnodes][9]
        sgremis1 = min(0.99, max(state.nsoilp[state.nnodes][4], sgremis1))
    elif state.ntype[state.nnodes] == 26:  # water
        # reflectancia de Fresnel (media das polarizacoes) a partir do
        # angulo zenital solar e do angulo de refracao (Snell, n=1.33)
        Z = state.met[state.iw][MetCol.ZEN] * math.pi / 180.0
        r = math.asin(math.sin(Z) / 1.33)
        if Z + r > 0.0:
            state.sgralbedo = 50.0 * (
                math.sin(Z - r) * math.sin(Z - r) / (math.sin(Z + r) * math.sin(Z + r))
                + math.tan(Z - r) * math.tan(Z - r) / (math.tan(Z + r) * math.tan(Z + r))
            )
        else:
            state.sgralbedo = state.nsoilp[state.nnodes][3]
        sgremis1 = state.nsoilp[state.nnodes][4]
    else:
        state.sgralbedo = state.nsoilp[state.nnodes][3]
        sgremis1 = state.nsoilp[state.nnodes][4]

    # ---------------- ajuste por neve/gelo acumulados (linhas 73-106) ----------------
    if abs(state.hsaccum + state.hi + state.newsd) <= EPS:
        # sem neve nem gelo
        state.emis = max(sgremis1, state.sgremis)
        state.albedo_fasst = min(state.albedo_fasst, state.sgralbedo)
    elif state.hi > 0.0 and abs(state.hsaccum + state.newsd) <= EPS:
        # gelo, sem neve
        state.emis = IEMIS
        state.albedo_fasst = IALBEDO
    elif state.hsaccum > oldsd or state.newsd > EPS:
        # neve nova
        state.emis = SEMIS
        state.albedo_fasst = SNALBEDO
        albedoo = state.albedo_fasst
    else:
        # neve envelhecida
        state.emis = SEMIS
        if state.iw == 1:
            albedoo = SNALBEDO

        # Douville et al. (1995), Climate Dynamics 12(1) p.21-35
        f1 = 24.0
        if abs(state.atopf) <= EPS:
            # sem derretimento -- decaimento linear
            salbedoo = albedoo - 8e-3 * state.timstep / f1
        else:
            # com derretimento -- decaimento exponencial em direcao a snalbedo
            salbedoo = SNALBEDO + (albedoo - SNALBEDO) * math.exp(-0.24 * state.timstep / f1)
        salbedoo = max(SOALBEDO, min(SNALBEDO, salbedoo))

        # Roesch (2000) -- decaimento polinomial com |temperatura - Tref|
        t1 = abs(state.toptemp - TREF)
        salbedon = (
            SNALBEDO
            - 0.07582627 * t1
            - 5.5360168e-3 * t1 * t1
            - 5.2966269e-5 * t1 * t1 * t1
            + 4.2372742e-6 * t1 * t1 * t1 * t1
        )
        salbedon = max(SOALBEDO, min(SNALBEDO, salbedon))

        state.albedo_fasst = min(state.albedo_fasst, max(salbedon, salbedoo))
        albedoo = state.albedo_fasst

    # ---------------- override por radiacao medida + arredondamento (linhas 108-111) ----------------
    if aflag == 1:
        state.albedo_fasst = albedo1

    state.albedo_fasst = _anint(state.albedo_fasst, 20)
    state.emis = _anint(state.emis, 20)

    # persiste a variavel `save` do Fortran para a proxima chamada
    alb.albedoo = albedoo
