"""Balanço de energia e água de um dossel vegetal de n camadas (FASST).

Fonte Fortran: module_canopy.F90 (106 693 bytes, 3 012 linhas, 22
sub-rotinas: 2 pontos de entrada públicos — canopy_met e veg_proph — + 20
auxiliares dentro do módulo). Ver doc-transcricao-module_canopy.docx para o
inventário das 22 sub-rotinas e as pendências.

O que faz
---------
A cada passo de tempo, :func:`canopy_met` (1) recalcula as propriedades
óticas/térmicas de cada camada — refletância, transmitância, absorção,
emissividade, fator de agrupamento — a partir do tipo de vegetação e da
estação do ano (:func:`canopy_prop`); (2) monta a matriz de fatores de
forma/visão Sij entre céu, camadas foliares e solo, por traçado de raios
sobre uma distribuição estatística de ângulos foliares (:func:`scalc_2a` e
as seis sub-rotinas que orquestra: kernal, leaf_slope, mean_proj, gap_prob,
fmatrix, smatrix — Smith & Goltz, 1994); (3) calcula os perfis de vento,
umidade relativa e temperatura do ar dentro do dossel (:func:`roughness` —
Raupach, 1994; :func:`windprofile` — Yamazaki et al., 1992); (4) resolve a
interceptação de chuva/neve na folhagem, o derretimento de neve
interceptada e a evapotranspiração, limitada pela taxa de transpiração das
raízes (:func:`canopy_moist`); (5) resolve a temperatura de equilíbrio de
cada camada por um método de Newton generalizado sobre o balanço de
energia não linear (:func:`feval` + :func:`solve`, com inversão de matriz
via LU em :func:`matinv`/:func:`ludcmp`/:func:`lubksb` — Numerical
Recipes); e (6) atualiza as variáveis meteorológicas (solar, IR, vento,
temperatura do ar, UR) que chegam ao solo já atenuadas pelo dossel.
:func:`veg_proph` é uma segunda porta de entrada, independente de
canopy_met: inicializa veg_prp (tabela de propriedades por tipo de
vegetação), a fração de raiz por nó de solo (rk) e, para vegetação alta, a
geometria de camadas do dossel (zh, dzveg, laif).

Arquitetura
-------------
Cada uma das 22 funções recebe `state: FasstState` como argumento
explícito, para todo campo que é estado (a enorme maioria dos ~80 nomes
de fasst_global usados no arquivo). Constantes físicas fixas (EPS, TREF,
PI, SIGMA, SPFLAG, VK, RV, RD, LHFUS, SDENSW) vêm de `fasst.constants`;
índices de coluna do meteorológico (doy etc.), de `fasst.indices`
(MetCol). `NCLAYERS` (=3) também é constante de `fasst.constants`, não
campo de estado.

`_anint`/`_aint` (arredondamento e truncamento no estilo Fortran) e
`_fdiv` (divisão que devolve +-inf/nan em vez de levantar
ZeroDivisionError) são definidas neste arquivo — mesma convenção de
`fasst/functions.py`/`sp_humid.py`/`snow.py`/`module_lowveg.py`, que não
compartilham esses helpers entre arquivos.

`albedo_emis`, `dense`, `spheats`, `vap_press` e `sp_humid` são chamadas
reais, importadas de `fasst.albedo_emis`/`fasst.functions`/
`fasst.sp_humid`. A chamada a `albedo_emis` dentro de `canopy_met` usa os
três argumentos (oldsd, sd, su) do comentário Fortran original, sem
reatribuição — a função atualiza `state.albedo_fasst`/`state.emis` como
efeito colateral, lidos normalmente mais adiante (por exemplo, em
`canopy_prop`). `sp_humid` devolve um `MoistAirState` (23 campos
nomeados); os call sites em `canopy_moist` leem 4 desses campos por nome
(mixr, dqdt, vpress, wetbulb) em cada chamada. A leitura das tabelas de
vegetação (units 31/32/33) é um stub (`_read_veg_table_row`): arquivo de
abertura ainda não identificado/traduzido.

As 7 variáveis `save::` de `canopy_met` (x, epsc, rho, clump, tau, alp,
psi), que precisam persistir entre chamadas sucessivas, formam o
dataclass :class:`CanopySavedState`, passado explicitamente e mutado
in-place. `NC`, `NC2` e `NPSI` (tamanhos de alocação dos arrays de camada
e da matriz de fator de forma) são derivados do maior índice de fato
usado em cada família de array, com `NCLAYERS=3` fixo.

Indexação 1-based com folga no índice 0: todo array cujo índice Fortran é
usado literalmente no código (x(i), veg_prp(vegh_type,6), stt(refn)...) é
alocado aqui com uma posição a mais em cada dimensão 1-based — o índice 0
fica sempre sem uso, e os índices 1..N correspondem exatamente aos do
Fortran. Exceção: `avect`, cuja primeira dimensão já é 0-based no próprio
Fortran (`avect(0:4, nclayers)`) e em `fasst/state.py`
(`avect: ... np.zeros((5, NCLAYERS))`, sem a folga que as outras arrays
1-based recebem).

Todas as 22 funções mantêm o nome Fortran original e correspondência 1:1
com a sub-rotina fonte.

Três bugs do Fortran original, preservados de propósito
----------------------------------------------------------
Nenhum é corrigido silenciosamente, porque corrigir mudaria o resultado
numérico do modelo:

* :func:`ludcmp` — a checagem do último pivô da diagonal
  (``a(num_layers,num_layers) <= TINY1``) não usa ``dabs()``, ao contrário
  de todas as outras checagens de pivô na mesma sub-rotina. Qualquer pivô
  final negativo é substituído por TINY1 (positivo), corrompendo a
  decomposição LU nesses casos.
* :func:`gap_prob` — quando ``j=1`` e o primeiro ``t1`` estoura 50, o
  Fortran recalcula usando ``cosr(j-1) = cosr(0)``, fora dos limites
  declarados do array no original.
* :func:`veg_proph` — com fonte de bioma alternativa
  (``biome_source > 0``), a coluna 1 de ``veg_prp`` é lida da tabela mas
  descartada e zerada logo em seguida, por uma troca aparente de nome de
  variável (``veg_prp`` em vez de ``newveg_prp``) no Fortran de 1995.

Outros pontos de atenção (não são bugs): variáveis locais de
canopy_moist (rpf, rptr, ff) que não são reiniciadas a cada iteração do
laço interno; o caso de biome_source não reconhecido em veg_proph; o
idioma de arredondamento com expoente variável (k=20/15/10/5 conforme o
ponto do código); o dimensionamento dinâmico de Sij/Wir/Fijr em scalc_2a
a partir do parâmetro num_layers recebido (não de uma constante fixa); e
a coincidência de nome entre o parâmetro ``sigfh`` de qfunc e o estado
global ``state.sigfh`` (são variáveis diferentes).
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import EPS, LHFUS, NCLAYERS, PI, RD, RV, SDENSW, SIGMA, SPFLAG, TREF, VK
from .functions import Phase, dense, spheats, vap_press
from .indices import MetCol
from .sp_humid import sp_humid
from .state import FasstState
from .albedo_emis import albedo_emis, AlbedoEmisState

__all__ = [
    "CanopySavedState", "canopy_met", "canopy_prop", "canopy_moist",
    "roughness", "windprofile", "scalc_2a", "kernal", "leaf_slope",
    "mean_proj", "gap_prob", "fmatrix", "smatrix", "radiosity", "feval",
    "bfunc", "qfunc", "output", "solve", "matinv", "ludcmp", "lubksb",
    "veg_proph",
]

# Tamanhos de alocação (ver docstring do módulo, "NC, NC2, NPSI") -- maior
# índice de fato usado em cada família de array, com folga de índice 0.
NC = NCLAYERS + 3     # epsc, alp, rho, tau, clump, x, sabs, a, lai, ful, ... (até nclayers+1/+2)
NC2 = NCLAYERS + 3    # Sij, Wir (pré-alocação; scalc_2a devolve o tamanho exato)
NPSI = NCLAYERS + 3   # psi (usa até nclayers+1)

MAX_THETAR = 30   # parameter local de canopy_met
N_THETAR = 9
N_THETAK = 9
MAX_ITER = 10
TOL = 1e-2


@dataclass
class CanopySavedState:
    """Variáveis `save::` de `canopy_met` (module_canopy.F90).

    x, epsc, rho, clump, tau, alp, psi persistem entre chamadas sucessivas
    de `canopy_met` -- substituem a cláusula `save` do Fortran original.
    """

    x: np.ndarray = field(default_factory=lambda: np.zeros(NC))
    epsc: np.ndarray = field(default_factory=lambda: np.zeros(NC))
    rho: np.ndarray = field(default_factory=lambda: np.zeros(NC))
    clump: np.ndarray = field(default_factory=lambda: np.zeros(NC))
    tau: np.ndarray = field(default_factory=lambda: np.zeros(NC))
    alp: np.ndarray = field(default_factory=lambda: np.zeros(NC))
    psi: np.ndarray = field(default_factory=lambda: np.zeros(NPSI))


def _anint(x: float, p: int = 0) -> float:
    """Fortran ``ANINT(x*10**p)*10**-p`` -- round half *away* from zero.

    Duplicada por arquivo -- este projeto não tem um `fortran_compat`
    compartilhado (ver `fasst/functions.py`, `sp_humid.py`, `snow.py`,
    `module_lowveg.py`). NaN/inf passam inalterados, como o `ANINT` do
    Fortran.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def _aint(x: float, p: int = 0) -> float:
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direção a zero."""
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _fdiv(x: float, y: float) -> float:
    """Divisão IEEE-754 pura: ``x / 0.0`` -> +-inf / nan, sem levantar.

    Mesma função de `sp_humid.py`/`snow.py`/`module_lowveg.py`.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


# ---------------------------------------------------------------------------
# canopy_prop  (module_canopy.F90, linhas 561-860)
# ---------------------------------------------------------------------------
# Constantes de dados (Fortran DATA) -----------------------------------------
# ATENÇÃO: Fortran preenche arrays multidimensionais em ordem column-major.
# O `data df_alp /.../ ` do original preenche COLUNA por coluna, então cada
# linha de 3 valores do literal Fortran (com o comentário !needle / !broad /
# !mixed) corresponde a uma COLUNA de df_alp, não a uma linha.
NC1 = 3
NC1P = NC1 + 1  # folga no índice 0

_FOLIAGE_TYPE_CONST = 6  # `data foliagetype /6/` no Fortran (tipo de ângulo foliar, Dorman & Sellers)

_MCLUMP = np.zeros(NC1P)
_MCLUMP[1], _MCLUMP[2], _MCLUMP[3] = 0.5, 0.5, 1.0                    # topo, meio, fundo

_DF_ALP = np.zeros((NC1P, NC1P))       # df_alp[camada, vind]  (vind: 1=agulha,2=folha larga,3=mista)
_DF_ALP[1, 1], _DF_ALP[2, 1], _DF_ALP[3, 1] = 0.229, 0.214, 0.079     # agulha
_DF_ALP[1, 2], _DF_ALP[2, 2], _DF_ALP[3, 2] = 0.255, 0.046, 0.038     # folha larga
_DF_ALP[1, 3], _DF_ALP[2, 3], _DF_ALP[3, 3] = 0.242, 0.130, 0.058     # mista

_RHO_MIN = np.zeros(NC1P)
_RHO_MIN[1], _RHO_MIN[2], _RHO_MIN[3] = 0.21, 0.26, 0.24              # agulha, folha larga, mista
_RHO_MAX = 0.28

_TAU_MIN = 0.001
_TAU_MAX = np.zeros(NC1P)
_TAU_MAX[1], _TAU_MAX[2], _TAU_MAX[3] = 0.150, 0.150, 0.100           # agulha, folha larga, mista

_DZVEG1O = np.zeros(6)  # índices 1..5, folga em 0
_DZVEG1O[1], _DZVEG1O[2], _DZVEG1O[3], _DZVEG1O[4], _DZVEG1O[5] = 11.2, 9.0, 13.7, 16.3, 13.5
# vegh_type:                                                            3,    4,    5,    6,   18


# ---------------------------------------------------------------------------
# 2/22  canopy_prop -- propriedades oticas/termicas por camada (rho, tau, absorcao, emissividade, clumping)
# ---------------------------------------------------------------------------
def canopy_prop(state: FasstState, epsc, alp, rho, tau, clump, psi):
    """
    Define propriedades óticas/térmicas das camadas do dossel (refletância,
    transmitância, absorção de onda curta, emissividade de onda longa e
    fator de agrupamento) por tipo de vegetação e estação do ano.

    Tradução literal de canopy_prop (module_canopy.F90, linhas 561-860).
    `epsc, alp, rho, tau, clump, psi` entram como os valores SAVE da chamada
    anterior (intent(inout) no Fortran) e saem atualizados.
    """
    eps_old = np.zeros(NC)
    alp_old = np.zeros(NC)
    rho_old = np.zeros(NC)
    tau_old = np.zeros(NC)
    clump_old = np.zeros(NC)
    psi_old = np.zeros(NPSI)

    foliage_type = np.zeros(NC, dtype=int)
    season = 0
    vind = 0
    vind1 = 0
    ftg = 0.0
    m = 0.0

    epsc = epsc.copy()
    alp = alp.copy()
    rho = rho.copy()
    tau = tau.copy()
    clump = clump.copy()
    psi = psi.copy()

    for i in range(1, NCLAYERS + 1):
        eps_old[i] = epsc[i]
        alp_old[i] = alp[i]
        rho_old[i] = rho[i]
        tau_old[i] = tau[i]
        clump_old[i] = clump[i]

        foliage_type[i] = 0
        epsc[i] = 0.0
        alp[i] = 0.0
        rho[i] = 0.0
        tau[i] = 0.0
        clump[i] = 0.0

    for i in range(1, NCLAYERS + 2):
        psi_old[i] = psi[i]
        psi[i] = 0.0

    if state.albedo_fasst > EPS:
        psi[NCLAYERS + 1] = state.albedo_fasst
    else:
        if psi_old[NCLAYERS + 1] > EPS:
            psi[NCLAYERS + 1] = psi_old[NCLAYERS + 1]
        else:
            psi[NCLAYERS + 1] = state.nsoilp[state.nnodes, 3]

    # ---- determina a estação do ano ---------------------------------------
    doy = state.met[state.iw, MetCol.DOY]
    if state.lat >= 0.0:                                                        # hemisfério norte
        if doy >= 335 or doy <= 80:
            season = 1                                                       # inverno
        if 81 <= doy <= 151:
            season = 2                                                       # primavera
        if 152 <= doy <= 243:
            season = 3                                                       # verão
        if 244 <= doy <= 334:
            season = 4                                                       # outono
    else:                                                                    # hemisfério sul
        if doy >= 335 or doy <= 80:
            season = 3                                                       # verão
        if 81 <= doy <= 151:
            season = 4                                                       # outono
        if 152 <= doy <= 243:
            season = 1                                                       # inverno
        if 244 <= doy <= 334:
            season = 2                                                       # primavera

    if state.iw == state.istart and state.infer_test == 0:
        state.iseason = season

    ftg = 1.0 - 1.6e-3 * (298.0 - state.stt[state.refn]) * (298.0 - state.stt[state.refn])

    if state.vegh_type == 3 or state.vegh_type == 6:                               # perene, agulha e folha larga
        vind, vind1 = 1, 1
        if state.vegh_type == 6:
            vind, vind1 = 2, 4

        for j in range(1, NCLAYERS + 1):
            foliage_type[j] = _FOLIAGE_TYPE_CONST
            clump[j] = _MCLUMP[j]
            if j == NCLAYERS:
                rho[j] = _RHO_MAX
                tau[j] = _TAU_MIN
                epsc[j] = state.veg_prp[state.vegh_type, 12]
                alp[j] = state.veg_prp[state.vegh_type, 14]
            else:
                rho[j] = _RHO_MIN[vind]
                tau[j] = _TAU_MAX[vind]
                epsc[j] = state.veg_prp[state.vegh_type, 13]
                if j == 1:
                    alp[j] = state.veg_prp[state.vegh_type, 15]
                else:
                    m = (_DF_ALP[1, vind] - _DF_ALP[2, vind]) * state.dzveg[1] / _DZVEG1O[vind1]
                    alp[j] = _anint((state.veg_prp[state.vegh_type, 15] - m), 20)

    elif state.vegh_type == 4 or state.vegh_type == 5:                             # decíduo, agulha e folha larga
        vind, vind1 = 1, 2
        if state.vegh_type == 5:
            vind, vind1 = 2, 3

        for j in range(1, NCLAYERS + 1):
            foliage_type[j] = _FOLIAGE_TYPE_CONST
            clump[j] = _MCLUMP[j]
            if j == NCLAYERS:
                rho[j] = _RHO_MAX
                tau[j] = _TAU_MIN
                epsc[j] = state.veg_prp[state.vegh_type, 12]
                alp[j] = state.veg_prp[state.vegh_type, 14]
            else:
                rho[j] = _RHO_MIN[vind] + (1.0 - ftg) * (_RHO_MAX - _RHO_MIN[vind])
                if rho[j] < _RHO_MIN[1]:
                    rho[j] = _RHO_MIN[1]
                if rho[j] > _RHO_MAX:
                    rho[j] = _RHO_MAX

                tau[j] = _TAU_MAX[vind] - (1.0 - ftg) * (_TAU_MAX[vind] - _TAU_MIN)
                if tau[j] < _TAU_MIN:
                    tau[j] = _TAU_MIN
                if tau[j] > _TAU_MAX[vind]:
                    tau[j] = _TAU_MAX[vind]

                epsc[j] = state.veg_prp[state.vegh_type, 13] - (1.0 - ftg) * (
                    state.veg_prp[state.vegh_type, 13] - state.veg_prp[state.vegh_type, 12])
                if epsc[j] < state.veg_prp[state.vegh_type, 12]:
                    epsc[j] = state.veg_prp[state.vegh_type, 12]
                if epsc[j] > state.veg_prp[state.vegh_type, 13]:
                    epsc[j] = state.veg_prp[state.vegh_type, 13]

                if j == 1:
                    alp[j] = state.veg_prp[state.vegh_type, 15]
                else:
                    m = (_DF_ALP[1, vind] - _DF_ALP[2, vind]) * state.dzveg[1] / _DZVEG1O[vind1]
                    alp[j] = _anint((state.veg_prp[state.vegh_type, 15] - m), 20)

    elif state.vegh_type == 18:                                                 # mista
        vind, vind1 = 3, 5

        for j in range(1, NCLAYERS + 1):
            foliage_type[j] = _FOLIAGE_TYPE_CONST
            clump[j] = _MCLUMP[j]
            if j == NCLAYERS:
                rho[j] = _RHO_MAX
                tau[j] = _TAU_MIN
                epsc[j] = state.veg_prp[state.vegh_type, 12]
                alp[j] = state.veg_prp[state.vegh_type, 14]
            else:
                rho[j] = _RHO_MIN[3] + (1.0 - ftg) * (_RHO_MAX - _RHO_MIN[3])
                if rho[j] < _RHO_MIN[3]:
                    rho[j] = _RHO_MIN[3]
                if rho[j] > _RHO_MAX:
                    rho[j] = _RHO_MAX

                tau[j] = _TAU_MAX[3] - (1.0 - ftg) * (_TAU_MAX[3] - _TAU_MIN)
                if tau[j] < _TAU_MIN:
                    tau[j] = _TAU_MIN
                if tau[j] > _TAU_MAX[3]:
                    tau[j] = _TAU_MAX[3]

                epsc[j] = state.veg_prp[state.vegh_type, 13] - (1.0 - ftg) * (
                    state.veg_prp[state.vegh_type, 13] - state.veg_prp[state.vegh_type, 12])
                if epsc[j] < state.veg_prp[state.vegh_type, 12]:
                    epsc[j] = state.veg_prp[state.vegh_type, 12]
                if epsc[j] > state.veg_prp[state.vegh_type, 13]:
                    epsc[j] = state.veg_prp[state.vegh_type, 13]

                if j == 1:
                    alp[j] = state.veg_prp[state.vegh_type, 15]
                else:
                    m = (_DF_ALP[1, vind] - _DF_ALP[2, vind]) * state.dzveg[1] / _DZVEG1O[vind1]
                    alp[j] = _anint((state.veg_prp[state.vegh_type, 15] - m), 20)

    if season == state.iseason:
        for i in range(1, NCLAYERS + 1):
            if state.ifoliage_type[i] != SPFLAG:
                foliage_type[i] = int(state.ifoliage_type[i])
            if state.iclump[i] != SPFLAG:
                clump[i] = state.iclump[i]
            if state.irho[i] != SPFLAG:
                rho[i] = state.irho[i]
            if state.itau[i] != SPFLAG:
                tau[i] = state.itau[i]
            if state.ialp[i] != SPFLAG:
                alp[i] = state.ialp[i]
            if state.ieps[i] != SPFLAG:
                epsc[i] = state.ieps[i]

    # o espalhamento foliar é aproximado como 1/2 refletância + 1/2 transmitância
    for i in range(1, NCLAYERS + 1):
        if rho[i] + tau[i] > 1.0:
            if state.single_multi_flag == 0:
                print(f" rho + tau > 1 in layer {i}")
            tau[i] = 0.99 - rho[i]
        psi[i] = rho[i] + tau[i]

    if state.iw != state.istart:
        for i in range(1, NCLAYERS + 1):
            rho[i] = 0.5 * (rho[i] + rho_old[i])
            tau[i] = 0.5 * (tau[i] + tau_old[i])
            psi[i] = 0.5 * (psi[i] + psi_old[i])
            alp[i] = 0.5 * (alp[i] + alp_old[i])
            epsc[i] = 0.5 * (epsc[i] + eps_old[i])
            clump[i] = 0.5 * (clump[i] + clump_old[i])
        psi[NCLAYERS + 1] = 0.5 * (psi[NCLAYERS + 1] + psi_old[NCLAYERS + 1])

    for i in range(1, NCLAYERS + 1):
        if state.dzveg[i] > EPS:
            rho[i] = _anint(rho[i], 20)
            tau[i] = _anint(tau[i], 20)
            psi[i] = _anint(psi[i], 20)
            alp[i] = _anint(alp[i], 20)
            epsc[i] = _anint(epsc[i], 20)
            clump[i] = _anint(clump[i], 20)
        else:
            rho[i] = 0.0
            tau[i] = 0.0
            if state.albedo_fasst > EPS:
                psi[i] = state.albedo_fasst
            else:
                if psi_old[i] > EPS:
                    psi[i] = psi_old[i]
                else:
                    psi[i] = state.nsoilp[state.nnodes, 3]
            psi[i] = _anint(psi[i], 20)
            alp[i] = 1.0 - psi[i]
            alp[i] = _anint(alp[i], 20)
            epsc[i] = _anint(state.emis, 20)
            clump[i] = 0.0
    psi[NCLAYERS + 1] = _anint(psi[NCLAYERS + 1], 20)

    return season, foliage_type, epsc, alp, rho, tau, clump, psi


# ---------------------------------------------------------------------------
# 4/22  roughness -- altura de deslocamento zero-plano e comprimento de rugosidade do dossel
# ---------------------------------------------------------------------------
def roughness(lamda, beta, h, htemp):
    """
    Altura de deslocamento zero-plano e comprimento de rugosidade do dossel.
    Tradução literal de roughness (module_canopy.F90, linhas 1257-1326).
    Referência: Raupach (1994), Boundary Layer Meteorology, 71, p.211-216.
    """
    USUH_MAX = 0.3     # u*/Uh máximo (adimensional)
    CR = 0.3           # coef. de arrasto do elemento de rugosidade isolado
    CS = 0.003         # coef. de arrasto do substrato na altura h
    CW = 2.0           # constante, profundidade da subcamada de rugosidade
    CD1 = 7.5          # constante (adimensional)

    psih = math.log(CW) - 1.0 + 1.0 / CW                       # função de influência da subcamada de rugosidade
    temp1 = math.sqrt(CS + CR * (lamda + beta) * 0.5)          # adimensional
    usuh = min(temp1, USUH_MAX)                          # u*/Uh (adimensional)
    temp2 = math.sqrt(CD1 * (lamda + beta))
    h2 = h - htemp                                         # altura do dossel - (baixa vegetação e/ou neve+gelo) (m)

    if temp2 > 50.0:
        d = h2
    else:
        d = h2 * (1.0 - (1.0 - math.exp(-temp2)) / temp2)
    d = _anint(d, 20)

    if abs(-VK / usuh - psih) > 50.0:
        z0 = h2
    else:
        z0 = h2 * (1.0 - d / h2) * math.exp(-VK / usuh - psih)
    z0 = _anint(z0, 20)

    return d, z0


# ---------------------------------------------------------------------------
# 5/22  windprofile -- perfis normalizados de vento, UR e temperatura do ar dentro do dossel
# ---------------------------------------------------------------------------
def windprofile(num_layers, zh, glb, sigfh, dz):
    """
    Perfis normalizados de vento, umidade relativa e temperatura do ar num
    dossel de n camadas.
    Tradução literal de windprofile (module_canopy.F90, linhas 1329-1445).
    Referência: Yamazaki et al. (1992), J. Applied Meteorology, 31, p.86-103.
    """
    wind_prof = np.zeros(num_layers + 1)
    rh_prof = np.zeros(num_layers + 1)
    tai = np.zeros(num_layers + 1)

    cstar = 1.0                                            # densidade do dossel (adimensional)
    xc = math.log10(cstar)

    # eq. (32) da referência (versão exponencial, como no código da NASA)
    t1 = xc + 0.26
    cstar_prime = 0.0
    if abs((t1 + math.sqrt(t1 * t1 + 0.16)) * 0.5 - 0.3) < 50.0:
        cstar_prime = math.exp((t1 + math.sqrt(t1 * t1 + 0.16)) * 0.5 - 0.3)

    # função de ponderação para cstar, eq. (31) da referência
    f = 0.37 + (0.494 * (xc + 0.8)) / math.sqrt((xc + 0.8) * (xc - 0.5) + 1.1)

    zdcist = 0.0
    zb = 1.0 / zh
    for i in range(1, num_layers + 1):
        # ponto médio da camada, em % da altura total do dossel
        zhalf_layer = dz[i] * 0.5
        if i == 1:
            zdcist = zdcist
        else:
            zdcist = zdcist + dz[i - 1]                    # distância até o topo
        zpercent = (zh - (zhalf_layer + zdcist)) * zb

        if abs(1.1 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1) < 50.0:
            wind_prof[i] = math.exp(1.1 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1)
        if wind_prof[i] <= EPS:
            wind_prof[i] = 0.0

        if abs(5.5 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1) < 50.0:
            rh_prof[i] = math.exp(5.5 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1)
        if rh_prof[i] <= EPS:
            rh_prof[i] = 0.0

        if glb <= EPS:
            if abs(-0.05 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1) < 50.0:
                tai[i] = math.exp(-0.05 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1)
        else:
            if abs(0.05 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1) < 50.0:
                tai[i] = math.exp(0.05 * (zpercent - 1.0)) * (1.0 - sigfh * 0.1)
        if tai[i] <= EPS:
            tai[i] = 0.0

        wind_prof[i] = _anint(wind_prof[i], 20)
        rh_prof[i] = _anint(rh_prof[i], 20)
        tai[i] = _anint(tai[i], 20)

    return wind_prof, rh_prof, tai


# ---------------------------------------------------------------------------
# 7/22  kernal -- nucleo para projecao media da distribuicao de angulos foliares
# ---------------------------------------------------------------------------
def kernal(max_thetar, max_thetak, nthetar, nthetak):
    """
    Núcleo (kernel) usado na projeção média, por camada, da distribuição de
    ângulos foliares. Tradução literal de kernal (module_canopy.F90,
    linhas 1523-1610).
    """
    thetar = np.zeros(max_thetar + 1)
    sinr = np.zeros(max_thetar + 1)
    cosr = np.zeros(max_thetar + 1)
    thetak = np.zeros(max_thetak + 1)
    fk = np.zeros((max_thetar + 1, max_thetak + 1))

    f2 = PI / 180.0

    # seno e cosseno dos ângulos de integração sobre a elevação, distribuídos uniformemente
    f1 = 1.0 / float(nthetar)
    for i in range(2, nthetar + 1):
        fi = i - 1
        thetar[i] = (fi * f1) * PI * 0.5                # thetar, não 90.0
        sinr[i] = math.sin(thetar[i])
        cosr[i] = math.cos(thetar[i])

    # trata os pontos extremos (0 e 90 graus)
    thetar[1] = 0.1 * f2
    sinr[1] = math.sin(thetar[1])
    cosr[1] = math.cos(thetar[1])
    thetar[nthetar + 1] = 89.9 * f2                        # último thetar fixado em 89.9
    sinr[nthetar + 1] = math.sin(thetar[nthetar + 1])
    cosr[nthetar + 1] = math.cos(thetar[nthetar + 1])

    # mesmo procedimento para os ângulos de distribuição foliar
    f1 = 1.0 / float(nthetak)
    for j in range(2, nthetak + 1):
        fj = j - 1
        thetak[j] = (fj * f1) * PI * 0.5                # thetak, não 90.0
    thetak[1] = 0.1 * f2
    thetak[nthetak + 1] = 89.9 * f2

    # calcula o núcleo usado na projeção média (por camada) de thetak sobre thetar
    f1 = 2.0 / PI
    for j in range(1, nthetak + 2):
        for i in range(1, nthetar + 2):
            if thetak[j] <= (PI * 0.5 - thetar[i]):
                fk[i, j] = cosr[i] * math.cos(thetak[j])
            else:
                phi_prime = math.acos(-(1.0 / math.tan(thetak[j])) * (1.0 / math.tan(thetar[i])))
                fk[i, j] = f1 * cosr[i] * math.cos(thetak[j]) * (phi_prime - PI * 0.5 - math.tan(phi_prime))
            fk[i, j] = fk[i, j] * f1                        # ajuste p/ bater com o relatório SWOE
            fk[i, j] = _anint(fk[i, j], 20)

    return thetar, sinr, cosr, thetak, fk


# ---------------------------------------------------------------------------
# 8/22  leaf_slope -- PDF da distribuicao de inclinacao foliar (6 tipos)
# ---------------------------------------------------------------------------
def leaf_slope(num_layers, max_thetak, nthetak, foliage_type, thetak):
    """
    Função densidade de probabilidade (PDF) para as diferentes distribuições
    de inclinação foliar (planófila, erectófila, plagiófila, extremófila,
    uniforme, esférica).
    Tradução literal de leaf_slope (module_canopy.F90, linhas 1613-1669).
    """
    f = np.zeros((num_layers + 1, max_thetak + 1))
    f1 = 2.0 / PI

    for i in range(1, num_layers + 1):
        ft = foliage_type[i]
        if ft == 1 or ft == 2:
            for j in range(1, nthetak + 2):
                f[i, j] = f1 * (1.0 + math.cos(2.0 * thetak[j]))            # planófila, erectófila
        elif ft == 3 or ft == 4:
            for j in range(1, nthetak + 2):
                f[i, j] = f1 * (1.0 - math.cos(4.0 * thetak[j]))            # plagiófila, extremófila
        elif ft == 5:
            for j in range(1, nthetak + 2):
                f[i, j] = f1                                            # uniforme
        else:
            for j in range(1, nthetak + 2):
                f[i, j] = math.sin(thetak[j])                                # esférica

        for j in range(1, nthetak + 2):
            f[i, j] = _anint(f[i, j], 20)

    return f


# ---------------------------------------------------------------------------
# 9/22  mean_proj -- projecao media das folhas de cada camada nas direcoes thetar
# ---------------------------------------------------------------------------
def mean_proj(num_layers, max_thetak, max_thetar, nthetar, nthetak, fk, f):
    """
    Projeção média de todas as folhas de uma camada nas direções thetar,
    integrada pela regra do trapézio estendida (Numerical Recipes, fórmula
    4.1.11).
    Tradução literal de mean_proj (module_canopy.F90, linhas 1672-1755).
    """
    MAXTK = 30

    g_bar = np.zeros((num_layers + 1, max_thetar + 1))
    integrand = np.zeros((num_layers + 1, max_thetar + 1, MAXTK + 1))

    h = 1.0 / nthetak

    for i in range(1, num_layers + 1):
        for j in range(1, nthetar + 2):
            for k in range(1, nthetak + 2):
                integrand[i, j, k] = f[i, k] * fk[j, k]

    for i in range(1, num_layers + 1):
        for j in range(1, nthetar + 2):
            sumf = 0.0
            for k in range(2, nthetak + 1):
                sumf = sumf + integrand[i, j, k]
            g_bar[i, j] = h * (0.5 * (integrand[i, j, 1] + integrand[i, j, nthetak + 1]) + sumf)
            g_bar[i, j] = _anint(g_bar[i, j], 20)

    return g_bar


# ---------------------------------------------------------------------------
# 10/22  gap_prob -- probabilidade de intervalo (gap probability) por camada e direcao
# ---------------------------------------------------------------------------
def gap_prob(num_layers, max_thetar, nthetar, lai, sai, clump, g_bar, cosr):
    """
    Probabilidade de intervalo (gap probability) ao longo da direção do raio.
    `clump` é o fator de agrupamento de Markov entre 0 e 1 (1 = sem
    agrupamento; 0 = agrupamento total, ou seja po = 1).
    Tradução literal de gap_prob (module_canopy.F90, linhas 1758-1811).

    ATENÇÃO (doc ponto 2 — herdado do Fortran original, não é erro de tradução): quando
    j=1 e o primeiro cálculo de t1 estoura o limite de 50, o Fortran refaz
    o cálculo usando cosr(j-1) = cosr(0) — um índice fora dos limites
    declarados do array no original. Mantive o mesmo acesso aqui
    (cosr[0], que é sempre 0.0 nesta tradução por causa da folga do
    índice 0), reproduzindo o comportamento do Fortran nesse caso extremo.
    """
    po = np.zeros((num_layers + 1, max_thetar + 1))
    clmp = np.zeros(num_layers + 1)

    for i in range(1, num_layers + 1):
        clmp[i] = clump[i]
        if clump[i] < EPS:
            clmp[i] = 1e-2
        for j in range(1, nthetar + 2):
            t1 = _fdiv((lai[i] + sai[i]) * clmp[i] * g_bar[i, j], cosr[j])
            if abs(t1) > 50.0:
                t1 = _fdiv(1.1 * (lai[i] + sai[i]) * clmp[i] * g_bar[i, j], cosr[j - 1])
            t1 = min(1.0, t1)
            po[i, j] = math.exp(-t1)
            po[i, j] = _anint(po[i, j], 20)

    return po


# ---------------------------------------------------------------------------
# 11/22  fmatrix -- coeficientes de fator de forma Fijr entre sumidouro/fonte
# ---------------------------------------------------------------------------
def fmatrix(num_layers, max_thetar, nthetar, po):
    """
    Coeficientes de fator de forma/visão Fijr entre sumidouro/fonte, por
    camada, para a direção thetar.

    Índices de sumidouro/fonte:
        1       - Céu
        2       - Camada foliar 1
        3       - Camada foliar 2
        ...
        n       - Camada foliar n
        n+1     - Solo
    (os índices de LAI começam em 1, não em 2 — comentário do Fortran original)

    Tradução literal de fmatrix (module_canopy.F90, linhas 1814-1944).
    """
    Fijr = np.zeros((num_layers + 3, num_layers + 3, max_thetar + 1))

    for r in range(1, nthetar + 2):
        # ---- coeficientes de contribuição para o CÉU ------------------------
        Fijr[1, 1, r] = 0.0                                              # do céu para o céu

        for k in range(2, num_layers + 2):                               # das camadas foliares para o céu
            prod = 1.0
            if k > 2:
                for m in range(2, k):
                    prod = prod * po[m - 1, r]
            Fijr[1, k, r] = (1.0 - po[k - 1, r]) * prod

        prod = 1.0                                                       # do solo para o céu
        for k in range(2, num_layers + 2):
            prod = prod * po[k - 1, r]
        Fijr[1, num_layers + 2, r] = prod

        # ---- coeficientes de contribuição para as CAMADAS FOLIARES ----------
        for i in range(1, num_layers + 3):
            for j in range(2, num_layers + 2):
                if j < i:
                    prod = 1.0
                    if j + 1 < i:
                        for k in range(j + 1, i):
                            prod = prod * po[k - 1, r]
                    if i != num_layers + 2:
                        f = math.sqrt(po[j - 1, r]) * (1.0 - po[i - 1, r]) * prod
                    else:
                        f = math.sqrt(po[j - 1, r]) * prod
                elif i == j:
                    f = 2.0 * (1.0 - math.sqrt(po[j - 1, r]))
                else:                                                     # i_sumidouro > i_fonte
                    prod = 1.0
                    if j > i + 1:
                        for k in range(i + 1, j):
                            prod = prod * po[k - 1, r]
                    if i != 1:                                            # ou seja, se não for o céu
                        f = math.sqrt(po[j - 1, r]) * (1.0 - po[i - 1, r]) * prod
                    else:
                        f = math.sqrt(po[j - 1, r]) * prod                    # CÉU
                Fijr[j, i, r] = f

        # ---- coeficientes de contribuição para o SOLO -----------------------
        prod = 1.0                                                       # do céu para o solo
        for k in range(2, num_layers + 2):
            prod = prod * po[k - 1, r]
        Fijr[num_layers + 2, 1, r] = prod

        for k in range(2, num_layers + 2):                                # da camada foliar para o solo
            prod = 1.0
            if k + 1 <= num_layers + 1:
                for m in range(k + 1, num_layers + 2):
                    prod = prod * po[m - 1, r]
            Fijr[num_layers + 2, k, r] = (1.0 - po[k - 1, r]) * prod

        Fijr[num_layers + 2, num_layers + 2, r] = 0.0                     # do solo para o solo

        for i in range(1, num_layers + 3):
            for j in range(1, num_layers + 3):
                Fijr[i, j, r] = _anint(Fijr[i, j, r], 20)

    return Fijr


# ---------------------------------------------------------------------------
# 12/22  smatrix -- integra Fijr sobre os angulos de visao -> Sij
# ---------------------------------------------------------------------------
def smatrix(num_layers, max_thetar, nthetar, sinr, cosr, Fijr):
    """
    Integra os fatores de visão sobre todos os ângulos de visão para obter
    os coeficientes de troca (Sij) entre fonte/sumidouro, pela regra do
    trapézio estendida.
    Tradução literal de smatrix (module_canopy.F90, linhas 1947-2019).
    """
    Sij = np.zeros((num_layers + 3, num_layers + 3))
    integrand = np.zeros((num_layers + 3, num_layers + 3, max_thetar + 1))

    h = 1.0 / nthetar
    for i in range(1, num_layers + 3):
        for j in range(1, num_layers + 3):
            for r in range(1, nthetar + 2):
                integrand[i, j, r] = sinr[r] * cosr[r] * Fijr[i, j, r]    # ângulos de integração

    for i in range(1, num_layers + 3):
        for j in range(1, num_layers + 3):
            sumf = 0.0
            for r in range(2, nthetar + 1):
                sumf = sumf + integrand[i, j, r]                          # integração sem os pontos extremos
            Sij[i, j] = PI * h * (0.5 * (integrand[i, j, 1] + integrand[i, j, nthetar + 1]) + sumf)
            Sij[i, j] = _anint(Sij[i, j], 20)

    return Sij


# ---------------------------------------------------------------------------
# 6/22  scalc_2a -- orquestra o calculo da matriz de fatores de forma Sij e Wir
# ---------------------------------------------------------------------------
def scalc_2a(num_layers, max_thetar, n_thetar, n_thetak, foliage_type, lai, sai, clump):
    """
    Monta a matriz de fatores de forma/visão Sij e os fatores de visão do
    dossel Wir, a partir da distribuição de ângulos foliares e das
    probabilidades de intervalo (gap probabilities) por camada.

    Tradução literal de scalc_2a (module_canopy.F90, linhas 1448-1520).
    Chama kernal, leaf_slope, mean_proj, gap_prob, fmatrix e smatrix.

    doc ponto 7: Sij/Wir/Fijr são dimensionados aqui a partir do parâmetro
    num_layers recebido nesta chamada (na prática, ncsnow — que pode ser
    menor que nclayers quando neve/gelo cobrem camadas do dossel), não a
    partir de uma constante fixa do módulo. Replica o dimensionamento
    assumido (assumed-size) do Fortran.
    """
    max_thetak = 30

    Wir = np.zeros((num_layers + 3, max_thetar + 1))

    # *** núcleo p/ projeção média (por camada) da distrib. de ângulos foliares
    thetar, sinr, cosr, thetak, fk = kernal(max_thetar, max_thetak, n_thetar, n_thetak)

    # *** distribuições de inclinação foliar
    f = leaf_slope(num_layers, max_thetak, n_thetak, foliage_type, thetak)

    # *** projeção média
    g_bar = mean_proj(num_layers, max_thetak, max_thetar, n_thetar, n_thetak, fk, f)

    # *** probabilidades de intervalo (gap probabilities)
    po = gap_prob(num_layers, max_thetar, n_thetar, lai, sai, clump, g_bar, cosr)

    # *** Fijr
    Fijr = fmatrix(num_layers, max_thetar, n_thetar, po)

    # *** extrai os fatores de visão do dossel, Wir, da primeira linha de Fijr
    for j in range(1, num_layers + 3):
        for r in range(1, n_thetar + 2):
            Wir[j, r] = _anint(Fijr[1, j, r], 20)

    # *** Sij
    Sij = smatrix(num_layers, max_thetar, n_thetar, sinr, cosr, Fijr)

    return Sij, Wir


# ---------------------------------------------------------------------------
# 20/22  ludcmp -- decomposicao LU com pivoteamento parcial (Numerical Recipes)
# ---------------------------------------------------------------------------
def ludcmp(num_layers, a):
    """
    Decomposição LU (método de Crout com pivoteamento parcial).
    Referência: Numerical Recipes (Fortran), p.35, Press et al. (1986).
    Tradução literal de ludcmp (module_canopy.F90, linhas 2534-2644).

    Diferente do Fortran original (que destrói `a` no lugar), aqui
    trabalhamos sobre uma cópia — nada no restante do código reaproveita a
    matriz original depois de chamar ludcmp/matinv, então isso não muda o
    comportamento, só evita efeitos colaterais indesejados em Python.

    ATENÇÃO — doc ponto 1, bug herdado do Fortran original, preservado aqui de propósito:
    a checagem de pivô quase-nulo usa `dabs(a(j,j))` em toda a coluna
    (linha 2630 do .F90), EXCETO na checagem final do último elemento da
    diagonal (linha 2641), que compara `a(num_layers,num_layers) <= TINY1`
    SEM valor absoluto. Como qualquer número negativo satisfaz essa
    comparação, sempre que o último pivô calculado for negativo (o que é
    comum, não incomum), ele é substituído por TINY1 (positivo), corrompendo
    a decomposição. Confirmei isso numericamente: sem essa linha, a inversa
    bate com a referência; com ela, diverge sempre que o último pivô é
    negativo. Mantive o comportamento exatamente como está no Fortran —
    corrigir isso seria uma mudança de resultado numérico, não uma correção
    de tradução, e cabe a quem mantém o modelo decidir se quer corrigir.
    """
    a = a.copy()
    TINY1 = 1e-10
    n = num_layers
    indx = np.zeros(n + 1, dtype=int)
    vv = np.zeros(n + 1)
    imax = 0
    d = 1.0

    for i in range(1, n + 1):
        aamax = 0.0
        for j in range(1, n + 1):
            if abs(a[i, j]) <= TINY1:
                # equivalente ao "write(*,*) + stop" do Fortran (erro fatal)
                raise ValueError(f"ludcmp: elemento quase nulo em a[{i},{j}] = {a[i, j]}")
            if abs(a[i, j]) > aamax:
                aamax = abs(a[i, j])
        if aamax <= TINY1:
            aamax = TINY1                                                # matriz singular
        vv[i] = 1.0 / aamax

    for j in range(1, n + 1):
        if j > 1:
            for i in range(1, j):
                sumf = a[i, j]
                if i > 1:
                    for k in range(1, i):
                        sumf = sumf - a[i, k] * a[k, j]
                    a[i, j] = sumf

        aamax = 0.0
        for i in range(j, n + 1):
            sumf = a[i, j]
            if j > 1:
                for k in range(1, j):
                    sumf = sumf - a[i, k] * a[k, j]
                a[i, j] = sumf
            dum = vv[i] * abs(sumf)

            if (dum >= aamax) or (abs(dum - aamax) < 1e-10):
                imax = i
                aamax = dum

        if j != imax:
            for k in range(1, n + 1):
                dum = a[imax, k]
                a[imax, k] = a[j, k]
                a[j, k] = dum
            d = -d
            vv[imax] = vv[j]

        indx[j] = imax if imax != 0 else 1

        if abs(a[j, j]) <= TINY1:
            a[j, j] = TINY1
        if j != n:
            dum = 1.0 / a[j, j]
            for i in range(j + 1, n + 1):
                a[i, j] = a[i, j] * dum

    d = _anint(d, 20)
    if a[n, n] <= TINY1:            # ver nota de bug herdado no docstring (falta dabs aqui, igual ao Fortran)
        a[n, n] = TINY1

    return a, indx, d


# ---------------------------------------------------------------------------
# 21/22  lubksb -- substituicao direta/reversa sobre a decomposicao LU
# ---------------------------------------------------------------------------
def lubksb(num_layers, indx, a, b):
    """
    Substituição direta/reversa (forward/back substitution) para o sistema
    LU x = b, usando a decomposição de ludcmp.
    Referência: Numerical Recipes (Fortran), p.36, Press et al. (1986).
    Tradução literal de lubksb (module_canopy.F90, linhas 2647-2703).
    """
    b = b.copy()
    TINY1 = 1e-10
    n = num_layers
    ii = 0

    for i in range(1, n + 1):
        ll = indx[i]
        sumf = b[ll]
        b[ll] = b[i]
        if ii != 0:
            for j in range(ii, i):
                sumf = sumf - a[i, j] * b[j]
        elif sumf != 0.0:
            ii = i
        b[i] = sumf

    for i in range(n, 0, -1):
        sumf = b[i]
        if i < n:
            for j in range(i + 1, n + 1):
                sumf = sumf - a[i, j] * b[j]
        if abs(a[i, i]) > TINY1:
            b[i] = sumf / a[i, i]
        else:
            b[i] = sumf * 1e10

    return b


# ---------------------------------------------------------------------------
# 19/22  matinv -- inversao de matriz via decomposicao LU (Crout)
# ---------------------------------------------------------------------------
def matinv(num_layers, a):
    """
    Inversão de matriz via decomposição LU (Crout) + substituição.
    Referência: Numerical Recipes (Fortran), p.38, Press et al. (1986).
    Tradução literal de matinv (module_canopy.F90, linhas 2481-2531).
    """
    n = num_layers
    y = np.zeros((n + 1, n + 1))
    for i in range(1, n + 1):
        y[i, i] = 1.0

    a_lu, indx, d = ludcmp(n, a)

    for j in range(1, n + 1):
        col = y[:, j].copy()
        col = lubksb(n, indx, a_lu, col)
        y[:, j] = col

    return y


# ---------------------------------------------------------------------------
# 18/22  solve -- metodo de Newton generalizado: dx = (J'J)^-1 J'(-y)
# ---------------------------------------------------------------------------
def solve(num_layers, y, a):
    """
    Resolve o sistema não linear y(x) = 0 via método de Newton generalizado:
    dx = (J'J)^-1 J'(-y), onde J = a é a matriz jacobiana.
    Tradução literal de solve (module_canopy.F90, linhas 2394-2478).
    Chama matinv/ludcmp/lubksb (traduzidas acima, neste mesmo lote).
    """
    n = num_layers
    ata = np.zeros((n + 1, n + 1))
    aty = np.zeros(n + 1)
    dx = np.zeros(n + 1)

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            ata[i, j] = 0.0
            for k in range(1, n + 1):
                ata[i, j] = ata[i, j] + a[k, i] * a[k, j]

    ata_inv = matinv(n, ata)

    for i in range(1, n + 1):
        aty[i] = 0.0
        for j in range(1, n + 1):
            aty[i] = aty[i] + a[j, i] * (-y[j])

    for i in range(1, n + 1):
        dx[i] = 0.0
        for j in range(1, n + 1):
            dx[i] = dx[i] + ata_inv[i, j] * aty[j]
            dx[i] = _anint(dx[i], 15)

    return dx


# ---------------------------------------------------------------------------
# 13/22  radiosity -- resolve o balanco de radiosidade (onda longa) via inversao de matriz
# ---------------------------------------------------------------------------
def radiosity(num_layers, Sij, psi):
    """
    Resolve o balanço de radiosidade entre céu, camadas foliares e solo,
    via inversão de matriz, com espalhamento ψ = ρ + τ por camada foliar e
    albedo no solo (propriedades de onda curta; `canopy_met` a chama com
    esses ψ e guarda o resultado em `sabs`/`aground`). Calcula a energia
    absorvida por camada e a fração de onda curta que chega ao solo.
    O comentário original desta função falava em "onda longa"; os
    coeficientes usados e os nomes das saídas indicam onda curta.

    Referência: Smith & Goltz (1994), IEEE Trans. Geosci. Remote Sens.,
    32(5), p.1060-1066.
    Tradução literal de radiosity (module_canopy.F90, linhas 2022-2126).
    """
    ETOTAL = 1.0
    d1 = num_layers + 1

    abs_coeff = np.zeros(d1 + 1)
    emission = np.zeros(d1 + 1)
    scatt = np.zeros(d1 + 1)
    absorb = np.zeros(d1 + 1)
    matrix = np.zeros((d1 + 1, d1 + 1))
    b = np.zeros(d1 + 1)

    albedoc = 0.0
    total_absorption = 0.0

    for i in range(1, d1 + 1):
        emission[i] = psi[i] * ETOTAL * Sij[i + 1, 1]                    # eq. (8) da ref.

    for j in range(1, d1 + 1):
        for i in range(1, d1 + 1):
            matrix[i, j] = -psi[i] * Sij[i + 1, j + 1]
            if i == j:
                matrix[i, j] = 1.0 + matrix[i, j]

    matinvs = matinv(d1, matrix)                                        # resolve a radiosidade

    for i in range(1, d1 + 1):
        b[i] = 0.0
        for j in range(1, d1 + 1):
            b[i] = b[i] + matinvs[i, j] * emission[j]                    # radiosidade, eq. (7) da ref.

    # *** calcula "abs_coeff(i)", o coeficiente de absorção
    for i in range(1, num_layers + 1):
        abs_coeff[i] = 1.0 - 2.0 * psi[i]                                # camadas foliares
    abs_coeff[d1] = 1.0 - psi[num_layers + 1]                            # camada do solo

    # *** contribuição do espalhamento múltiplo
    for i in range(1, d1 + 1):
        scatt[i] = 0.0
        for j in range(1, d1 + 1):
            scatt[i] = scatt[i] + b[j] * Sij[i + 1, j + 1]

    # fluxo solar total que chega ao solo = fluxo no topo do dossel * Sij(nclayers+2,CÉU) + scatt(nclayers+1)
    aground = max(0.0, Sij[num_layers + 2, 1] + scatt[d1])
    aground = _anint(aground, 20)

    # verificação do balanço de energia
    for i in range(1, d1 + 1):
        absorb[i] = abs_coeff[i] * (ETOTAL * Sij[i + 1, 1] + scatt[i])   # eq. (9) da ref.
        absorb[i] = _anint(absorb[i], 20)

        total_absorption = total_absorption + absorb[i]
        albedoc = albedoc + b[i] * Sij[1, i + 1]

    # total_absorption + albedoc é calculado no Fortran original só para conferência,
    # sem uso posterior — mantido aqui por fidelidade, mas não retornado.
    total = total_absorption + albedoc  # noqa: F841

    return absorb, aground


# ---------------------------------------------------------------------------
# 3/22  canopy_moist -- interceptacao de chuva/neve, derretimento na folhagem, evapotranspiracao
# ---------------------------------------------------------------------------
def canopy_moist(state: FasstState, ncsnow, lai, x, ful, rh_prof, tai):
    """
    Calcula evaporação/condensação do dossel e interceptação de precipitação
    (chuva/neve) pela folhagem.

    Tradução literal de canopy_moist (module_canopy.F90, linhas 863-1254).
    `ncsnow, lai, x, ful, rh_prof, tai` são intent(in); todo o resto que a
    subrotina original devolve (rain, snow, tot_ep, rl, ra, ef, def, chf,
    sigfhi, pheat, dpheat, stcl, stcs) é retornado aqui como tupla, na mesma
    ordem da lista de argumentos original.

    Nota sobre sp_humid: a sub-rotina Fortran declara 10 saídas posicionais
    para cada chamada de sp_humid, mas este código só nomeia 4 delas (mixr,
    dqdt, vpress, wetbulb); as outras 6 nunca são lidas depois de
    calculadas. sp_humid devolve um MoistAirState com 23 campos nomeados;
    os 4 usados aqui são lidos por nome, e os 6 não utilizados não são
    capturados.
    """
    f1i = 0
    d1i = 0
    pdens1 = 0.0
    f1a = 0.0
    f1 = 0.0
    f2 = 0.0
    f3 = 0.0
    c1 = 0.0
    c2 = 0.0
    t1 = 0.0
    d1 = 0.0
    qaf = 0.0
    f2a = 0.0
    mixra = 0.0
    vpressa = 0.0
    wetbulba = 0.0
    mixrf = 0.0
    dqdtf = 0.0
    vpressf = 0.0
    ep = 0.0
    cf = 0.0
    rpp = 0.0
    rpf = 0.0
    rptr = 0.0
    pdens = 0.0
    sheatcw = 0.0
    sheatcs = 0.0
    ff = 0.0
    etr = 0.0
    taf = 0.0
    min_wat = 0.0
    max_wat = 0.0
    kths = 0.0
    rhoaf = 0.0
    tot_ep = 0.0

    rl = np.zeros(NC)
    ra = np.zeros(NC)
    ef = np.zeros(NC)
    defv = np.zeros(NC)          # `def` é palavra reservada em Python
    chf = np.zeros(NC)
    sigfhi = np.zeros(NC)
    pheat = np.zeros(NC)
    dpheat = np.zeros(NC)
    precip = np.zeros(NC)
    precip1 = np.zeros(NC)
    interc = np.zeros(NC)
    interc1 = np.zeros(NC)
    drip = np.zeros(NC)
    drip1 = np.zeros(NC)
    meltc = np.zeros(NC)
    max_wet = np.zeros(NC)
    max_wetl = np.zeros(NC)
    max_wets = np.zeros(NC)
    stcs = np.zeros(NC)
    stcl = np.zeros(NC)

    for i in range(1, NCLAYERS + 1):
        stcs[i] = state.storcs[i]
        stcl[i] = state.storcl[i]

    # interceptação de precipitação, armazenamento máximo -------------------
    d1 = state.storcs[1] + state.storcs[2] + state.storcs[3]
    if (state.hsaccum > EPS or state.hi > EPS) or d1 > EPS:
        f1i = 1

    c1 = 0.0
    c2 = state.dmet1[state.iw, 5] * 1e-2
    moist = sp_humid(f1i, state.met[state.iw, MetCol.AP], state.dmet1[state.iw, 4], c2, c1)
    mixra, vpressa, wetbulba = moist.mixr, moist.vpress, moist.wetbulb

    # densidade e calor específico da chuva/neve
    d1i = 1
    d1 = 0.0
    pdens = dense(wetbulba, d1, Phase.WATER)                      # kg/m^3, água
    sheatcw = spheats(wetbulba, Phase.WATER)                      # J/kg*K

    if _aint(state.met[state.iw, MetCol.PT]) == 3 or _aint(state.met[state.iw, MetCol.PT2]) == 3:
        pdens1 = dense(wetbulba, state.dmet1[state.iw, 3], Phase.SNOW)  # kg/m^3, neve
        sheatcs = spheats(wetbulba, Phase.ICE)                   # J/kg*K
    else:
        pdens1 = SDENSW
        sheatcs = 2050.0

    # armazenamento de precipitação na folhagem
    # Referências: Aston (1979), J. Hydrology 42; Hoyningen-Huene (1983)
    for i in range(1, ncsnow + 1):
        sigfhi[i] = _anint(state.sigfh * state.laif[i], 20)
        max_wet[i] = state.veg_prp[state.vegh_type, 7] * (lai[i] + state.veg_prp[state.vegh_type, 8])  # R&S (mm)

        max_wetl[i] = max_wet[i] * 1e-3                                     # m
        max_wets[i] = (max_wet[i] * 1e-3) * (pdens1 / pdens)                 # m
        max_wetl[i] = _anint(max_wetl[i], 20)
        max_wets[i] = _anint(max_wets[i], 20)

        if i == 1:
            if _aint(state.met[state.iw, MetCol.PT]) == 2 or _aint(state.met[state.iw, MetCol.PT]) == 4:   # chuva
                precip[i] = state.dmet1[state.iw, 6] * state.timstep                   # m/passo (hora parcial)
            elif _aint(state.met[state.iw, MetCol.PT]) == 3:                          # neve
                precip1[i] = state.dmet1[state.iw, 6] * state.timstep
            elif _aint(state.met[state.iw, MetCol.PT2]) == 3:                         # neve
                precip1[i] = state.dmet1[state.iw, 7] * state.timstep
        else:
            precip[i] = max(0.0, precip[i - 1] - stcl[i - 1])
            precip1[i] = max(0.0, precip1[i - 1] - stcs[i - 1])
        precip[i] = _anint(precip[i], 20)
        precip1[i] = _anint(precip1[i], 20)

        t1 = 0.5 * (lai[i] + state.veg_prp[state.vegh_type, 8])
        t1 = _anint(t1, 20)
        if t1 > 50.0:
            t1 = 50.0

        if precip[i] > 0.0:                                                  # chuva
            interc[i] = precip[i] * (1.0 - math.exp(-t1)) * sigfhi[i]            # m
            interc[i] = max(0.0, min(interc[i], max_wetl[i]))
        elif precip1[i] > 0.0:                                               # neve
            interc1[i] = precip1[i] * (1.0 - math.exp(-t1)) * sigfhi[i]          # m
            interc1[i] = max(0.0, min(interc1[i], max_wets[i]))

        interc[i] = _anint(interc[i], 20)
        interc1[i] = _anint(interc1[i], 20)

        stcs[i] = stcs[i] + interc1[i]
        if stcs[i] > max_wets[i]:
            drip1[i] = stcs[i] - max_wets[i]                                 # m
            stcs[i] = max_wets[i]                                           # m
        drip1[i] = _anint(drip1[i], 20)
        stcs[i] = _anint(stcs[i], 20)

        if stcs[i] > 0.0:
            if pdens1 < 0.156:                                              # W/m*K, cond. térmica (Sturm et al.)
                kths = 0.023 + 0.234 * pdens1
                kths = min(max(0.023, kths), 1.0)
            else:
                kths = 0.138 - 1.01 * pdens1 + 3.233 * pdens1 * pdens1
                kths = min(max(0.138, kths), 1.0)
            meltc[i] = _fdiv(kths, stcs[i]) * tai[i] * state.timstep * 3.6e3 / (pdens1 * LHFUS)
            meltc[i] = max(0.0, min(stcs[i] * pdens / pdens1, meltc[i]))
            meltc[i] = _anint(meltc[i], 20)
            stcs[i] = stcs[i] - meltc[i]
        stcs[i] = _anint(stcs[i], 20)

        stcl[i] = stcl[i] + interc[i] + meltc[i]
        if stcl[i] > max_wetl[i]:
            drip[i] = stcl[i] - max_wetl[i]
            stcl[i] = max_wetl[i]
        drip[i] = _anint(drip[i], 20)
        stcl[i] = _anint(stcl[i], 20)

    # raízes ------------------------------------------------------------------
    state.trmhm = 0.0                                                          # taxa máx. de transpiração (m/s)
    for i in range(1, state.nnodes + 1):
        if state.stt[i] > TREF:
            state.frh[i] = state.rk[state.vegh_type, i] * (
                1.0 - (state.soil_moist[i] - state.nsoilp[i, 16]) / (state.nsoilp[i, 17] - state.nsoilp[i, 16]))
            state.frh[i] = max(0.0, min(1.0, _anint(state.frh[i], 20)))
        else:
            state.frh[i] = 0.0

        state.trmhm += state.frh[i]
        f2a += state.rk[state.vegh_type, i] * state.soil_moist[i]
        min_wat += state.nsoilp[i, 16]
        max_wat += state.nsoilp[i, 17]

    ff = 1.0
    state.trmhm = 1.5e-7 * state.sigfh * ff * state.trmhm                            # m/s
    state.trmhm = _anint(state.trmhm, 20)

    # resistência estomática de cada camada
    # REF: http://www.ecmwf.int/research/ifsdocs/CY25r1/PHYSICS/
    f1a = (4e-3 * state.dmet1[state.iw, 1] + 5e-3) / (0.81 * (4e-3 * state.dmet1[state.iw, 1] + 1.0))
    f1 = 1.0 / min(1.0, f1a)                                              # adimensional

    # _fdiv: max_wat==min_wat é possível -- ver a mesma protecao em module_lowveg.py
    f2 = max(0.0, min(1.0, _fdiv(f2a - min_wat, max_wat - min_wat)))
    if f2 > EPS:
        f2 = 1.0 / f2

    for i in range(1, ncsnow + 1):
        c2 = rh_prof[i] * 1e-2
        moist = sp_humid(f1i, state.met[state.iw, MetCol.AP], tai[i] + TREF, c2, c1)
        mixra, vpressa, wetbulba = moist.mixr, moist.vpress, moist.wetbulb

        c2 = 1.0
        moist = sp_humid(f1i, state.met[state.iw, MetCol.AP], x[i] + TREF, c2, c1)
        mixrf, dqdtf, vpressf = moist.mixr, moist.dqdt, moist.vpress

        t1 = 3e-4 * (vpressf - vpressa)
        if vpressf != vpressa and abs(t1) < 50.0:
            f3 = math.exp(t1)                                                   # adimensional
        else:
            f3 = 1.0

        if lai[i] > EPS:
            rl[i] = max(0.0, (state.veg_prp[state.vegh_type, 1] / lai[i]) * f1 * f2 * f3)  # s/m
        else:
            rl[i] = 0.0
        rl[i] = _anint(rl[i], 20)

        c1 = vap_press(state, state.nnodes + 1, state.met[state.iw, MetCol.RH] * 1e-2, state.met[state.iw, MetCol.AP])
        taf = (1.0 - 0.65 * state.sigfh) * (tai[i] + TREF) + 0.65 * state.sigfh * (x[i] + TREF)
        rhoaf = (c1 / RV + (state.met[state.iw, MetCol.AP] * 1e2 - c1) / RD) / taf   # densidade do ar na folhagem
        rhoaf = _anint(rhoaf, 20)

        if ful[i] > EPS:
            cf = 1.0 + 0.3 / ful[i]                                          # coef. de transf. em bloco
            ra[i] = 1.0 / (cf * ful[i])                                      # resist. atmosférica (s/m)
        else:
            cf = 0.0
            ra[i] = 0.0
        cf = _anint(state.sigfh * cf, 20)
        ra[i] = _anint(ra[i], 20)

        if ra[i] + rl[i] > EPS:
            rpp = ra[i] / (ra[i] + rl[i])                                    # umectação da vegetação
        rpp = min(max(0.0, _anint(rpp, 20)), 1.0)

        d1 = 1.0 - 0.65 * sigfhi[i] * (1.0 - rpp)
        qaf = ((1.0 - 0.65 * sigfhi[i]) * mixra + 0.65 * sigfhi[i] * rpp * mixrf) / d1
        ep = lai[i] * cf * (rhoaf / dense(wetbulba, d1, Phase.WATER)) * ful[i] * (qaf - rpp * mixrf)  # evap. potencial (m/s)

        if rpp * mixrf >= qaf or state.trmhm <= EPS:
            if f2 > min_wat and ep <= EPS:
                ep = 0.0
        ep = _anint(ep, 20)

        if rl[i] + ra[i] > EPS:
            # doc ponto 4: ff/rpf/rptr NÃO são reiniciadas a cada i — se esta
            # condição for falsa, ficam com o valor da iteração anterior
            # (mesmo comportamento do Fortran original, preservado).
            if stcs[i] > 0.0:
                ff = _fdiv(stcs[i], max_wets[i]) ** (2.0 / 3.0)
            elif stcl[i] > 0.0:
                ff = _fdiv(stcl[i], max_wetl[i]) ** (2.0 / 3.0)
            rpf = 1.0 - (rl[i] / (rl[i] + ra[i])) * (1.0 - ff)
            rptr = (ra[i] / (rl[i] + ra[i])) * (1.0 - ff)

        rptr = min(max(0.0, rptr), 1.0)
        rptr = _anint(rptr, 20)
        rpf = min(max(0.0, rpf), 1.0)
        rpf = _anint(rpf, 20)

        ef[i] = rpf * ep                                                     # evap. potencial da folha (m/s)
        ef[i] = _anint(ef[i], 20)
        etr = min(state.trmhm, rptr * ep)                                     # transpiração da folha (m/s)
        etr = _anint(etr, 20)

        if stcs[i] > EPS:
            stcs[i] = stcs[i] - (ef[i] - etr) * (state.timstep * 3.6e3)         # m
            if stcs[i] > max_wets[i]:
                drip1[i] = drip1[i] + (stcs[i] - max_wets[i])                # m
                stcs[i] = max_wets[i]                                        # m
        elif stcl[i] > EPS:
            stcl[i] = stcl[i] - (ef[i] - etr) * (state.timstep * 3.6e3)         # m
            if stcl[i] > max_wetl[i]:
                drip[i] = drip[i] + (stcl[i] - max_wetl[i])                  # m
                stcl[i] = max_wetl[i]                                        # m

        drip[i] = max(0.0, _anint(drip[i], 20))
        stcl[i] = max(0.0, _anint(stcl[i], 20))
        drip1[i] = max(0.0, _anint(drip1[i], 20))
        stcs[i] = max(0.0, _anint(stcs[i], 20))

        ef[i] = rhoaf * lai[i] * cf * ful[i] * (qaf - rpp * mixrf)           # kg/m^2*s
        defv[i] = rhoaf * lai[i] * cf * ful[i] * rpp * dqdtf * (0.65 * sigfhi[i] / d1 - 1.0)

        ef[i] = _anint(ef[i], 20)
        defv[i] = _anint(defv[i], 20)

        tot_ep = tot_ep + ef[i] * rpf

        if abs(lai[i]) <= EPS:
            chf[i] = 0.0
        else:
            chf[i] = max(0.0, 1.1 * lai[i] * rhoaf * spheats(taf, Phase.DRY_AIR) * cf * ful[i])  # folhagem (W/m^2*K)
        chf[i] = _anint(chf[i], 20)

        pheat[i] = (1e-2 * sigfhi[i]
                    * (sheatcw * pdens * interc[i] + sheatcs * pdens1 * interc1[i])
                    * (wetbulba - (x[i] + TREF)))
        pheat[i] = _anint(pheat[i], 20)
        dpheat[i] = -1e-2 * sigfhi[i] * (sheatcw * pdens * interc[i] + sheatcs * pdens1 * interc1[i])
        dpheat[i] = _anint(dpheat[i], 20)

    rain = max(0.0, (precip[ncsnow] + drip[ncsnow]))                       # m
    rain = _anint(rain, 20)

    snow = max(0.0, (precip1[ncsnow] + drip1[ncsnow]))                     # m
    snow = _anint(snow, 20)

    tot_ep = _anint(tot_ep, 15)

    return rain, snow, tot_ep, rl, ra, ef, defv, chf, sigfhi, pheat, dpheat, stcl, stcs


# ---------------------------------------------------------------------------
# 15/22  bfunc -- fluxo emitido normalizado por sigma (epsilon*T^4) e sua derivada
# ---------------------------------------------------------------------------
def bfunc(epsi, xi):
    """
    Fluxo emitido, normalizado pela constante de Stefan-Boltzmann (sem o
    sigma). Tradução literal de bfunc (module_canopy.F90, linhas 2236-2267).
    """
    bxi = 0.0
    dbxi = 0.0

    tc = _anint((xi + TREF), 20)
    if tc > EPS:
        bxi = epsi * tc * tc * tc * tc                                   # emissividade * T^4
        bxi = _anint(bxi, 20)

        dbxi = 4.0 * epsi * tc * tc * tc                                 # derivada em relação a T
        dbxi = _anint(dbxi, 20)

    return bxi, dbxi


# ---------------------------------------------------------------------------
# 16/22  qfunc -- fluxo de calor sensivel e sua derivada
# ---------------------------------------------------------------------------
def qfunc(xi, chfi, sigfh, taii):
    """
    Fluxo de calor sensível.
    Referência: Tibbals et al. (1964), Am. J. Bot. 51:529-538.
    Tradução literal de qfunc (module_canopy.F90, linhas 2270-2320).

    doc ponto 8: o parâmetro `sigfh` aqui NÃO é o estado global state.sigfh
    (fração de cobertura do dossel como um todo) — no call site (dentro de
    feval), recebe sigfhi(i), a fração de cobertura ponderada por camada
    calculada em canopy_moist. Coincidência de nome herdada do Fortran
    (o parâmetro formal se chama sigfh mesmo recebendo sigfhi(i)).
    """
    tafi = (1.0 - 0.65 * sigfh) * taii + 0.65 * sigfh * xi
    qxi = chfi * (tafi - xi)
    qxi = _anint(qxi, 20)

    dqxi = chfi * (0.65 * sigfh - 1.0)
    dqxi = _anint(dqxi, 20)

    return qxi, dqxi


# ---------------------------------------------------------------------------
# 14/22  feval -- monta o residuo (fx) e a jacobiana (dfx) do balanco de energia por camada
# ---------------------------------------------------------------------------
def feval(num_layers, icall, Sij, x, epsc, alp, bta, btg, a, ef, defv, chf,
          sigfhi, tai, pheat, dpheat):
    """
    Monta o vetor de resíduos do balanço de energia (fx) e sua matriz
    jacobiana (dfx) para o método de Newton usado por `solve`.

    Convenção de sinais: Rn = LE + H (F = Rn - LE - H). B e A (radiativos)
    são positivos para fluxo entrando no dossel; LE, H são positivos para
    fluxo saindo do dossel.

    Tradução literal de feval (module_canopy.F90, linhas 2129-2233).
    """
    lecan = 0.0
    lwabs = 0.0
    lwemit = 0.0

    two_fac = np.zeros(num_layers + 1)
    drx = np.zeros(num_layers + 1)
    dbx = np.zeros(num_layers + 1)
    dqx = np.zeros(num_layers + 1)
    bx = np.zeros(num_layers + 1)
    qx = np.zeros(num_layers + 1)
    rx = np.zeros(num_layers + 1)
    fx = np.zeros(num_layers + 1)
    dfx = np.zeros((num_layers + 1, num_layers + 1))

    for i in range(1, num_layers + 1):
        # escala H e LE pelo LAI
        bx[i], dbx[i] = bfunc(epsc[i], x[i])
        qx[i], dqx[i] = qfunc(x[i], chf[i], sigfhi[i], tai[i])

        lhevap = 2500775.6 - 2369.729 * tai[i]                           # J/kg = (m/s)^2
        rx[i] = _anint(ef[i] * lhevap, 20)
        drx[i] = _anint(defv[i] * lhevap, 20)

    for il in range(1, num_layers + 1):
        lw_abs = 0.0
        lecan = lecan + rx[il]
        for j in range(1, num_layers + 1):
            lw_abs = lw_abs + bx[j] * Sij[il + 1, j + 1]
            two_fac[il] = two_fac[il] + Sij[j + 1, il + 1]
        two_fac[il] = two_fac[il] + Sij[1, il + 1] + Sij[num_layers + 2, il + 1]  # deveria ser exatamente 2.0

        fx[il] = (alp[il] * SIGMA * (bta * Sij[il + 1, 1] + lw_abs + btg * Sij[il + 1, num_layers + 2])
                  - two_fac[il] * SIGMA * bx[il]
                  + a[il] + qx[il] + rx[il] + pheat[il])
        fx[il] = _anint(fx[il], 20)

        if icall == 1:
            lwabs = lwabs + alp[il] * SIGMA * (bta * Sij[il + 1, 1] + lw_abs + btg * Sij[il + 1, num_layers + 2])
            lwemit = lwemit + two_fac[il] * SIGMA * bx[il]

    for il in range(1, num_layers + 1):
        for j in range(1, num_layers + 1):
            if j == il:
                dfx[il, j] = (alp[il] * SIGMA * dbx[il] * Sij[il + 1, j + 1]
                              - two_fac[il] * SIGMA * dbx[il] + dqx[il] + drx[il] + dpheat[il])
            else:
                dfx[il, j] = alp[il] * SIGMA * dbx[j] * Sij[il + 1, j + 1]
            dfx[il, j] = _anint(dfx[il, j], 20)

    return bx, qx, rx, fx, dfx


# ---------------------------------------------------------------------------
# 17/22  output -- fluxos de onda longa ceu/dossel/solo; energia total que chega ao solo (bgr)
# ---------------------------------------------------------------------------
def output(num_layers, Sij, bx, btg, bta):
    """
    Fluxos de onda longa trocados entre céu, camadas foliares e solo (para
    diagnóstico/saída) e a energia total de onda longa que chega ao solo (bgr).

    Tradução literal de output (module_canopy.F90, linhas 2323-2391).
    Nota: na matriz Sij, o índice num_layers+2 é o solo e o índice 1 é o céu
    (camada 2 é o topo do dossel); nos demais arrays, num_layers+1 é o solo.
    """
    bx_tognd = np.zeros(num_layers + 1)
    bx_fromgnd = np.zeros(num_layers + 1)
    bx_tosky = np.zeros(num_layers + 2)
    bx_fromsky = np.zeros(num_layers + 2)

    bgr = 0.0

    for i in range(1, num_layers + 1):
        bx_tosky[i] = SIGMA * bx[i] * Sij[1, i + 1]                   # da camada i para o céu
        bx_tognd[i] = SIGMA * bx[i] * Sij[num_layers + 2, i + 1]      # da camada i para o solo
        bgr = bgr + bx_tognd[i]
        bx_fromsky[i] = SIGMA * bta * Sij[i + 1, 1]                   # do céu para a camada i
        bx_fromgnd[i] = SIGMA * btg * Sij[i + 1, num_layers + 2]      # do solo para a camada i

    bx_tosky[num_layers + 1] = SIGMA * btg * Sij[1, num_layers + 2]           # do solo para o céu
    bx_fromsky[num_layers + 1] = SIGMA * bta * Sij[num_layers + 2, 1]         # do céu para o solo
    bgr = bgr + bx_fromsky[num_layers + 1]                               # total que chega ao solo
    bgr = _anint(bgr, 20)

    return bgr


# ---------------------------------------------------------------------------
# canopy_met  (module_canopy.F90, linhas 5-558)
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 1/22  canopy_met -- ponto de entrada publico -- orquestra um passo de tempo completo do dossel
# ---------------------------------------------------------------------------
def canopy_met(state: FasstState, canopy: CanopySavedState, alb: AlbedoEmisState, oldsd):
    """
    Calcula os perfis de temperatura e fluxos de energia do dossel para um
    passo de tempo. Tradução literal de canopy_met (Fortran), preservando
    indexação 1-based (arrays com folga no índice 0) e arredondamento
    Fortran (_anint/_aint) nos mesmos pontos do original.

    Nota sobre a chamada a albedo_emis: segue os três argumentos do
    comentário Fortran original (`call albedo_emis(oldsd,dmet1(iw,1),
    dmet1(iw,8))`), sem reatribuição — albedo_emis não devolve nada,
    atualiza state.albedo_fasst/state.emis como efeito colateral, lidos
    normalmente mais adiante (ex.: state.albedo_fasst em canopy_prop).
    `alb` carrega a variável `save` de albedo_emis (`albedoo`) entre
    chamadas sucessivas.
    """
    # ---- variáveis SAVE (persistem entre chamadas) ------------------------
    x, epsc, rho = canopy.x, canopy.epsc, canopy.rho
    clump, tau, alp, psi = canopy.clump, canopy.tau, canopy.alp, canopy.psi

    # ---- variáveis locais (zeradas a cada chamada, como no Fortran) -------
    iterations = 0
    ncsnow = 0
    icall = 1
    ll = 0
    dzflag = 0
    zos = 0.0
    fmu = 0.0
    fmuh = 0.0
    glb = 0.0
    bta = 0.0
    btg = 0.0
    ta = 0.0
    tg = 0.0
    rh = 0.0
    epstg = 0.0
    sdepth = 0.0
    htemp = 0.0
    lw = 0.0
    ftg = 0.0
    laim = 0.0
    wscorr = 0.0
    zu = 0.0
    rain = 0.0
    snow = 0.0
    temp = 0.0
    tot_ep = 0.0
    dtop = 0.0
    aground = 0.0
    zdc = 0.0
    zo = 0.0
    bgr = 0.0
    sumdx = 0.0
    tsign = 0.0
    ii = 0
    j = 0

    sabs = np.zeros(NC)
    a = np.zeros(NC)
    lai = np.zeros(NC)
    ful = np.zeros(NC)
    dzveg1 = np.zeros(NC)
    sai = np.zeros(NC)

    foliage_type = np.zeros(NC, dtype=int)
    dx = np.zeros(NC)
    fx = np.zeros(NC)
    sigfhi = np.zeros(NC)
    rl = np.zeros(NC)
    wind_prof = np.zeros(NC)
    ra = np.zeros(NC)
    ef = np.zeros(NC)
    defv = np.zeros(NC)          # `def` é palavra reservada em Python -> defv
    chf = np.zeros(NC)
    rh_prof = np.zeros(NC)
    tai = np.zeros(NC)
    pheat = np.zeros(NC)
    dpheat = np.zeros(NC)
    xco = np.zeros(NC)
    stcl = np.zeros(NC)
    stcs = np.zeros(NC)
    Sij = np.zeros((NC2, NC2))
    Wir = np.zeros((NC2, MAX_THETAR + 1))
    bx = np.zeros(NC)

    ttemp = state.toptemp
    state.toptemp = _anint(state.toptemp, 20)

    # call albedo_emis(oldsd,dmet1(iw,1),dmet1(iw,8)) -- ver nota no docstring
    albedo_emis(state, alb, oldsd, state.dmet1[state.iw, 1], state.dmet1[state.iw, 8])

    # call canopy_prop(season,foliage_type,epsc,alp,rho,tau,clump,psi)
    season, foliage_type, epsc, alp, rho, tau, clump, psi = canopy_prop(state, epsc, alp, rho, tau, clump, psi)
    canopy.epsc, canopy.alp = epsc, alp
    canopy.rho, canopy.tau = rho, tau
    canopy.clump, canopy.psi = clump, psi

    # Read met data (time, air tmp, grnd tmp, wind speed, rh, total solar, IR)
    fmu = max(0.1, state.dmet1[state.iw, 3])                                    # wind speed (m/s)
    ta = _anint((state.dmet1[state.iw, 4] - TREF), 20)               # air temp (C)
    tg = _anint((state.toptemp - TREF), 20)                       # surface temp (C)
    rh = state.dmet1[state.iw, 5]                                                 # relative humidity (decimal)
    glb = state.dmet1[state.iw, 1]                                                # total incoming solar (W/m^2)
    lw = state.dmet1[state.iw, 2]                                                 # incoming IR (W/m^2)
    bta = _anint((lw / SIGMA), 20)                             # normalized by St.Bolz. constant

    if _aint(abs(state.met[state.iw, MetCol.IRUP] - state.mflag) * 1e5) * 1e-5 <= EPS:
        if ttemp > EPS and state.ftemp > EPS:
            btg = ((1.0 - state.sigfl) * ((1.0 - state.emis) * state.dmet1[state.iw, 2]
                                        + state.emis * SIGMA * ttemp ** 4.0)
                   + state.sigfl * ((1.0 - state.epf) * state.dmet1[state.iw, 2]
                                 + state.epf * SIGMA * state.ftemp ** 4.0))     # W/m^2 (-)
            btg = _anint(btg, 20)
    else:
        btg = state.met[state.iw, MetCol.IRUP]
    btg = _anint((btg / SIGMA), 20)

    fmuh = fmu

    if state.iw == state.istart:
        for i in range(1, NCLAYERS + 1):
            x[i] = 0.9 * ta                                                 # inicializa camadas na temp. do ar (C)
            state.storcl[i] = 0.0                                    # água armazenada nas folhas (m)
            state.storcs[i] = 0.0                                    # neve armazenada nas folhas
            if state.infer_test == 0:
                state.avect[0, i] = x[i]

        if state.infer_test == 1:
            for i in range(1, NCLAYERS + 1):
                x[i] = state.canopy_temp[i, state.oldpos] - TREF               # C
                state.storcl[i] = state.stor[i]
                state.storcs[i] = state.stor[i + NCLAYERS]

    if state.iw >= 2:
        for i in range(1, NCLAYERS + 1):
            x[i] = state.canopy_temp[i, state.iw - 1] - TREF                   # C

    for i in range(1, NCLAYERS + 1):
        xco[i] = x[i]

    zu = state.zh + state.iheight                                                 # altura de medição do vento (m)
    sdepth = state.hsaccum + state.hi                                              # espessura combinada neve+gelo (m)

    htemp = state.hfol_tot if state.hfol_tot >= sdepth else sdepth

    if sdepth == 0.0:                                                       # sem neve no solo
        zos = state.rough if state.vegl_type == 0 else state.z0l
    else:                                                                   # com neve no solo
        if state.vegl_type == 0 or (htemp <= sdepth or abs(state.hfol_tot - sdepth) < 1e-10):
            zos = 7.775e-3                                                  # rugosidade da neve (m)
        else:
            zos = state.z0l

    ftg = 1.0 - 1.6e-3 * (298.0 - state.stt[state.refn]) * (298.0 - state.stt[state.refn])
    ftg = _anint(ftg, 20)

    laim = state.veg_prp[state.vegh_type, 6] + ftg * (state.veg_prp[state.vegh_type, 5] - state.veg_prp[state.vegh_type, 6])
    if laim > state.veg_prp[state.vegh_type, 5]:
        laim = state.veg_prp[state.vegh_type, 5]
    if laim < state.veg_prp[state.vegh_type, 6]:
        laim = state.veg_prp[state.vegh_type, 6]
    laim = _anint(laim, 20)

    # sigfh: Ramirez e Senarath (2000), J. Climate, 13(22), p.4050-
    sigfh = state.veg_prp[state.vegh_type, 2] - (1.0 - ftg) * (state.veg_prp[state.vegh_type, 2] - state.veg_prp[state.vegh_type, 3])
    if sigfh < state.veg_prp[state.vegh_type, 2]:
        sigfh = state.veg_prp[state.vegh_type, 2]
    if sigfh > state.veg_prp[state.vegh_type, 3]:
        sigfh = state.veg_prp[state.vegh_type, 3]
    sigfh = sigfh * 1e-2
    if season == state.iseason and abs(state.isigfh - SPFLAG) > EPS:
        sigfh = state.isigfh
    if abs(state.isigfh - SPFLAG) > EPS and (sigfh > state.isigfh and state.iseason >= season):
        sigfh = state.isigfh
    sigfh = _anint(sigfh, 20)
    state.sigfh = sigfh

    ncsnow = NCLAYERS
    for i in range(1, NCLAYERS + 1):
        lai[i] = laim * state.laif[i]                                         # divide entre as camadas
        if (state.iw == state.istart and state.infer_test == 0) and state.ilai[i] != SPFLAG:
            lai[i] = state.ilai[i]
        lai[i] = _anint(lai[i], 20)

        dzveg1[i] = state.dzveg[i]
        if i == NCLAYERS:
            dzveg1[i] = state.dzveg[i] - htemp
            if dzveg1[i] <= EPS:
                ncsnow -= 1
                dzveg1[i] = 0.0
                if i != 1:
                    dzveg1[i - 1] = dzveg1[i - 1] - (htemp - state.dzveg[i])
                    if dzveg1[i - 1] < 0.0:
                        ncsnow -= 1
                        dzveg1[i - 1] = 0.0
                        dzveg1[i - 2] = dzveg1[i - 2] - (htemp - (state.dzveg[i] + state.dzveg[i - 1]))
                        if dzveg1[i - 2] < 0.0:
                            ncsnow -= 1
                            dzveg1[i - 2] = 0.0
        dzveg1[i] = _anint(dzveg1[i], 20)
        if dzveg1[i] <= EPS:
            dzflag = 1
        sai[i] = state.veg_prp[state.vegh_type, 7] * (1.0 - state.laif[i])
        sai[i] = _anint(sai[i], 20)

    if ncsnow != 0:
        # determina o material da superfície (solo) e ajusta parâmetros
        psi[ncsnow + 1] = (1.0 - state.sigfl) * psi[NCLAYERS + 1] + state.sigfl * state.albf   # albedo da superfície
        epstg = state.emis                                                     # emissividade da superfície

        Sij, Wir = scalc_2a(ncsnow, MAX_THETAR, N_THETAR, N_THETAK, foliage_type, lai, sai, clump)

        # desloc. altura e comprimento de rugosidade do dossel
        sum_lai = sum(lai[1:ncsnow + 1])
        sum_sai = sum(sai[1:ncsnow + 1])

        zdc, zo = roughness(sum_lai, sum_sai, state.zh, htemp)

        # perfil de vento e umidade relativa no dossel
        wind_prof, rh_prof, tai = windprofile(ncsnow, state.zh, glb, sigfh, dzveg1)

        # fator de correção logarítmico do vento
        wscorr = math.log((state.zh - zdc) / zo) / math.log((zu - zdc) / zo)
        fmuh = fmu * wscorr                                                 # vento no topo do dossel
        for i in range(1, ncsnow + 1):
            ful[i] = _anint(fmuh * wind_prof[i], 20)

            rh_prof[i] = 1e2 + (rh - 1e2) * rh_prof[i]
            if rh_prof[i] > 1e2:
                rh_prof[i] = 1e2
            if rh_prof[i] < 1e-1:
                rh_prof[i] = 1e-1
            rh_prof[i] = _anint(rh_prof[i], 20)

            tai[i] = _anint(ta * tai[i], 20)                      # C

        sabs, aground = radiosity(ncsnow, Sij, psi)

        sumdx = 0.0
        for i in range(1, ncsnow + 1):
            a[i] = max(0.0, glb * sabs[i])                                # energia de onda curta absorvida
            a[i] = _anint(a[i], 20)
            sumdx = sumdx + 1.1 * TOL

        while abs(sumdx) > TOL:
            rain, snow, tot_ep, rl, ra, ef, defv, chf, sigfhi, pheat, dpheat, stcl, stcs = canopy_moist(
                state, ncsnow, lai, x, ful, rh_prof, tai)

            bx, qx, rx, fx, dfx = feval(
                ncsnow, icall, Sij, x, epsc, alp, bta, btg, a, ef, defv, chf,
                sigfhi, tai, pheat, dpheat)

            dx = solve(ncsnow, fx, dfx)

            sumdx = 0.0
            dtop = 2.5
            for i in range(1, ncsnow + 1):
                tsign = 0.0
                if abs(dx[i]) > EPS:
                    tsign = _anint(abs(dx[i]) / dx[i])
                x[i] = x[i] + dx[i]                                          # C
                x[i] = _anint(x[i], 15)

                sumdx = sumdx + abs(dx[i])

            iterations += 1
            if iterations >= MAX_ITER:
                break

        if iterations >= MAX_ITER:
            for i in range(1, ncsnow + 1):
                x[i] = tai[i]
            state.error_code = 1
            state.error_type = 2

        bx, qx, rx, fx, dfx = feval(
            ncsnow, icall, Sij, x, epsc, alp, bta, btg, a, ef, defv, chf,
            sigfhi, tai, pheat, dpheat)

        icall = 0
        bgr = output(ncsnow, Sij, bx, btg, bta)

    if ncsnow < NCLAYERS:
        if ncsnow == 0:
            aground = 1.0
            bgr = lw
        for i in range(ncsnow + 1, NCLAYERS + 1):
            state.storcl[i] = 0.0
            state.storcs[i] = 0.0

            if state.hfol_tot >= sdepth and state.hfol_tot > EPS:
                x[i] = state.ftemp - TREF                                   # C
            elif sdepth > EPS:
                x[i] = state.stt[state.nnodes + 1] - TREF
            else:
                x[i] = ta
    elif dzflag == 1:
        for i in range(1, ncsnow + 1):
            if dzveg1[i] <= EPS:
                state.storcl[i] = 0.0
                state.storcs[i] = 0.0
                x[i] = tai[i] - TREF

    ll = 2
    if state.timstep >= 2:
        ll = 1
    for i in range(1, NCLAYERS + 1):
        state.storcl[i] = _anint(stcl[i], 20)
        state.storcs[i] = _anint(stcs[i], 20)

        if state.iw > ll or (state.iw == ll and state.infer_test == 1):
            j = 0
            temp = state.avect[j + 1, i]
            while j < ll:
                state.avect[j, i] = temp
                j += 1
                if j <= ll - 1:
                    temp = state.avect[j + 1, i]
            state.avect[ll, i] = x[i]
            ii = ll
            j = ll
        else:
            j = state.iw
            state.avect[j, i] = x[i]
            ii = state.iw
            j = state.iw

        sumt = 0.0
        while j > -1:
            sumt = sumt + state.avect[j, i]
            j -= 1

        x[i] = sumt / float(ii + 1)
        x[i] = _anint(x[i], 10)

        state.canopy_temp[i, state.iw] = _anint((x[i] + TREF), 10)   # K
        state.stor[i] = state.storcl[i]
        state.stor[i + NCLAYERS] = state.storcs[i]

    # usa efeitos calculados do dossel p/ estimar temp. do ar, solar, IR, vento
    state.dmet1[state.iw, 1] = aground * glb                                     # solar total acima do dossel (W/m^2)
    if state.dmet1[state.iw, 1] < 1e-2:
        state.dmet1[state.iw, 1] = 0.0
    state.dmet1[state.iw, 2] = bgr                                               # IR acima do dossel (W/m^2)

    if ncsnow != 0:
        state.dmet1[state.iw, 3] = ful[ncsnow]                                   # vento acima do dossel (m/s)
        state.dmet1[state.iw, 4] = tai[ncsnow] + TREF                         # temp. do ar acima do dossel (K)
        state.dmet1[state.iw, 5] = rh_prof[ncsnow]                               # umidade relativa acima do dossel (%)

    state.dmet1[state.iw, 6] = rain
    state.dmet1[state.iw, 7] = snow
    if state.dmet1[state.iw, 6] <= EPS and _aint(state.met[state.iw, MetCol.PT]) != 1:
        state.met[state.iw, MetCol.PT] = float(1)
    if state.dmet1[state.iw, 6] > EPS and _aint(state.met[state.iw, MetCol.PT]) != 2:
        state.met[state.iw, MetCol.PT] = float(2)
    if state.dmet1[state.iw, 7] <= EPS and _aint(state.met[state.iw, MetCol.PT2]) != 1:
        state.met[state.iw, MetCol.PT2] = float(1)
    if state.dmet1[state.iw, 7] > EPS and _aint(state.met[state.iw, MetCol.PT2]) != 3:
        state.met[state.iw, MetCol.PT2] = float(3)

    if state.met[state.iw, MetCol.TSOL] > EPS and state.dmet1[state.iw, 1] > EPS:
        if _aint(abs(state.dmet1[state.iw, 8] - state.mflag) * 1e5) * 1e-5 > EPS:
            state.dmet1[state.iw, 8] = state.dmet1[state.iw, 8] * state.dmet1[state.iw, 1] / state.met[state.iw, MetCol.TSOL]  # solar refletido (W/m^2)
        state.dmet1[state.iw, 9] = state.dmet1[state.iw, 9] * state.dmet1[state.iw, 1] / state.met[state.iw, MetCol.TSOL]        # solar direto (W/m^2)
        state.dmet1[state.iw, 10] = state.dmet1[state.iw, 10] * state.dmet1[state.iw, 1] / state.met[state.iw, MetCol.TSOL]      # solar difuso (W/m^2)
    else:
        if _aint(abs(state.dmet1[state.iw, 8] - state.mflag) * 1e5) * 1e-5 > EPS:
            state.dmet1[state.iw, 8] = 0.0
        state.dmet1[state.iw, 9] = 0.0
        state.dmet1[state.iw, 10] = 0.0

    for i in range(1, 14):
        state.dmet1[state.iw, i] = _anint(state.dmet1[state.iw, i], 15)

    sumdx = 0.0
    for i in range(1, ncsnow + 1):
        sumdx = sumdx + epsc[i]
    if ncsnow != 0:
        state.surfemisc[state.iw] = sumdx / float(ncsnow)
    else:
        state.surfemisc[state.iw] = sumdx

    state.toptemp = ttemp
    state.cevap[state.iw] = tot_ep                                               # evaporação/condensação total do dossel (kg/m^2*s)

    # grava de volta o estado persistente (SAVE)
    canopy.x = x
    canopy.epsc = epsc


# ---------------------------------------------------------------------------
# veg_proph  (module_canopy.F90, linhas 2706-3010)
# ---------------------------------------------------------------------------
# cinfo(linha,coluna): igual ao df_alp do Lote 1, o DATA do Fortran preenche
# column-major — cada bloco de 7 valores com comentário de tipo de vegetação
# no original é uma COLUNA aqui, indexada por `vind`.
_CINFO = np.zeros((8, 6))
_CINFO[1, 1], _CINFO[2, 1], _CINFO[3, 1], _CINFO[4, 1], _CINFO[5, 1], _CINFO[6, 1], _CINFO[7, 1] = (
    25.5, 0.439215686, 0.462745098, 0.098039216, 0.4, 0.5, 0.1)      # agulha perene       (vegh_type 3, vind 1)
_CINFO[1, 2], _CINFO[2, 2], _CINFO[3, 2], _CINFO[4, 2], _CINFO[5, 2], _CINFO[6, 2], _CINFO[7, 2] = (
    36.0, 0.452777778, 0.477777778, 0.069444444, 0.4, 0.5, 0.1)      # folha larga perene  (vegh_type 6, vind 2)
_CINFO[1, 3], _CINFO[2, 3], _CINFO[3, 3], _CINFO[4, 3], _CINFO[5, 3], _CINFO[6, 3], _CINFO[7, 3] = (
    21.0, 0.428571429, 0.452380952, 0.119047619, 0.4, 0.5, 0.1)      # agulha decídua      (vegh_type 4, vind 3)
_CINFO[1, 4], _CINFO[2, 4], _CINFO[3, 4], _CINFO[4, 4], _CINFO[5, 4], _CINFO[6, 4], _CINFO[7, 4] = (
    31.5, 0.434920635, 0.46984127, 0.095238095, 0.4, 0.5, 0.1)       # folha larga decídua (vegh_type 5, vind 4)
_CINFO[1, 5], _CINFO[2, 5], _CINFO[3, 5], _CINFO[4, 5], _CINFO[5, 5], _CINFO[6, 5], _CINFO[7, 5] = (
    30.0, 0.45, 0.466666667, 0.083333333, 0.4, 0.5, 0.1)             # mista               (vegh_type 18, vind 5)


def _read_veg_table_row(file_unit, target_vid):
    """
    TODO: leitura de arquivo externo — no Fortran, `read(file_unit,*,iostat=io)
    vid,srmax,srmin,cmax,cmin,rl,laimax,laimin,ddmax,sai,ar,br,emissmin,
    emissmax,folamin,folamax,heigmin,heigmax`, repetida linha a linha até
    encontrar `vid == target_vid` ou chegar ao fim do arquivo.

    As unidades 31 (tabela padrão), 32 (Modis-NOAH) e 33 (UMD) são abertas
    fora de module_canopy.F90 (provavelmente em fasst_main.F90 ou
    fasst_driver.F90, ainda não traduzidos). Quando essa parte for
    implementada, esta função deve devolver uma tupla de 17 valores
    (srmax,srmin,cmax,cmin,rl,laimax,laimin,ddmax,sai,ar,br,emissmin,
    emissmax,folamin,folamax,heigmin,heigmax) — ou None se target_vid não
    for encontrado até o fim do arquivo (equivalente a io == -1 no Fortran).
    """
    raise NotImplementedError(
        "leitura das tabelas de propriedades de vegetação (units 31/32/33) "
        "ainda não foi traduzida — arquivos abertos fora de module_canopy.F90"
    )


# ---------------------------------------------------------------------------
# 22/22  veg_proph -- propriedades de vegetacao, fracao de raiz e geometria do dossel (2a porta de entrada)
# ---------------------------------------------------------------------------
def veg_proph(state: FasstState, biome_source, new_vt, veg_type):
    """
    Define as propriedades de vegetação (veg_prp) a partir da tabela padrão
    e, opcionalmente, de uma fonte de bioma alternativa (MODIS-NOAH ou UMD);
    calcula a fração de raiz por nó do perfil de solo (rk); e, quando for
    vegetação alta, define a geometria do dossel (zh, dzveg, laif).

    Tradução literal de veg_proph (module_canopy.F90, linhas 2706-3010).

    ATENÇÃO — doc ponto 3, dois comportamentos herdados do Fortran original, preservados
    aqui de propósito (não são erro de tradução):

    1) Quando `biome_source > 0`, o Fortran grava a coluna 1 (resistência
       estomática mínima) direto em `veg_prp(new_vt,1)` (linha 2875),
       enquanto todas as outras colunas (2 a 17) vão para o array local
       `newveg_prp`. Isso parece uma troca de nome por engano no original.
       O problema é que, logo depois (linha ~2904), o laço que decide o
       valor final de `veg_prp` para cada coluna 1..17 sempre LÊ de
       `newveg_prp(new_vt,i)` — inclusive para i=1, que nunca foi escrito
       nesse array (ficou no valor zerado da alocação). Ou seja, o valor de
       `srmin` lido da tabela de bioma para a coluna 1 é descartado e
       substituído por 0.0 logo em seguida, de forma silenciosa. Mantive
       exatamente essa cadeia de atribuições.
    2) doc ponto 5 — se `biome_source > 0` mas não é 1000 (Modis-NOAH) nem
       2000 (UMD), o Fortran tentaria alocar/ler de uma unidade de arquivo
       não definida (fid ficaria 0). Esse caso não está representado nos
       dados de teste e o comportamento original ali é essencialmente
       indefinido; aqui simplesmente não populamos `newveg_prp` nesse caso.

    Chama `_read_veg_table_row` (ainda stub — leitura de arquivo externo).
    """
    zdw = 0.0
    zup = 0.0
    sumf = 0.0

    defveg_prp = np.zeros((19, 18))    # (tipo de vegetação 1..18, coluna 1..17)
    newveg_prp = None

    if state.infer_test == 0:
        state.zh = 0.0
        for i in range(1, NCLAYERS + 1):
            state.dzveg[i] = 0.0
            state.laif[i] = 0.0

    # ---- tabela padrão (unit 31) -------------------------------------------
    row = _read_veg_table_row(31, veg_type)
    if row is not None:
        (srmax, srmin, cmax, cmin, rl, laimax, laimin, ddmax, sai, ar, br,
         emissmin, emissmax, folamin, folamax, heigmin, heigmax) = row
        defveg_prp[veg_type, 1] = srmin
        defveg_prp[veg_type, 2] = cmax
        defveg_prp[veg_type, 3] = cmin
        defveg_prp[veg_type, 4] = rl
        defveg_prp[veg_type, 5] = laimax
        defveg_prp[veg_type, 6] = laimin
        defveg_prp[veg_type, 7] = ddmax
        defveg_prp[veg_type, 8] = sai
        defveg_prp[veg_type, 9] = ar
        defveg_prp[veg_type, 10] = br
        defveg_prp[veg_type, 11] = srmax
        defveg_prp[veg_type, 12] = emissmin * 1e-2
        defveg_prp[veg_type, 13] = emissmax * 1e-2
        defveg_prp[veg_type, 14] = folamin * 1e-2
        defveg_prp[veg_type, 15] = folamax * 1e-2
        defveg_prp[veg_type, 16] = heigmin * 1e-2
        defveg_prp[veg_type, 17] = heigmax * 1e-2

    # ---- fonte de bioma alternativa (unit 32 ou 33), se houver -------------
    if biome_source > 0:
        fid = None
        ntypes = 0
        if biome_source == 1000:                                          # Modis-NOAH
            fid = 32
            ntypes = 20
        elif biome_source == 2000:                                        # UMD
            fid = 33
            ntypes = 14

        if fid is not None:
            newveg_prp = np.zeros((ntypes + 1, 18))
            row = _read_veg_table_row(fid, new_vt)
            if row is not None:
                (srmax, srmin, cmax, cmin, rl, laimax, laimin, ddmax, sai, ar, br,
                 emissmin, emissmax, folamin, folamax, heigmin, heigmax) = row
                # ver nota de bug herdado no docstring: coluna 1 vai para state.veg_prp,
                # não para newveg_prp — igual ao Fortran original.
                state.veg_prp[new_vt, 1] = srmin
                newveg_prp[new_vt, 2] = cmax
                newveg_prp[new_vt, 3] = cmin
                newveg_prp[new_vt, 4] = rl
                newveg_prp[new_vt, 5] = laimax
                newveg_prp[new_vt, 6] = laimin
                newveg_prp[new_vt, 7] = ddmax
                newveg_prp[new_vt, 8] = sai
                newveg_prp[new_vt, 9] = ar
                newveg_prp[new_vt, 10] = br
                newveg_prp[new_vt, 11] = srmax
                newveg_prp[new_vt, 12] = emissmin * 1e-2
                newveg_prp[new_vt, 13] = emissmax * 1e-2
                newveg_prp[new_vt, 14] = folamin * 1e-2
                newveg_prp[new_vt, 15] = folamax * 1e-2
                newveg_prp[new_vt, 16] = heigmin * 1e-2
                newveg_prp[new_vt, 17] = heigmax * 1e-2

    # ---- decide o valor final de veg_prp para cada coluna 1..17 -----------
    if biome_source == 0:
        for i in range(1, 18):
            state.veg_prp[veg_type, i] = defveg_prp[veg_type, i]
            state.veg_prp[veg_type, i] = _anint(state.veg_prp[veg_type, i], 20)
    else:
        for i in range(1, 18):
            nv = newveg_prp[new_vt, i] if newveg_prp is not None else SPFLAG
            if abs(nv - SPFLAG) <= EPS:
                state.veg_prp[veg_type, i] = defveg_prp[veg_type, i]
            else:
                state.veg_prp[veg_type, i] = nv
            state.veg_prp[veg_type, i] = _anint(state.veg_prp[veg_type, i], 20)

    # ---- fração de raiz na camada i, por tipo de vegetação -----------------
    for i in range(1, state.nnodes + 1):
        if i == 1:
            zdw = state.elev - state.nzi[i]
            zup = state.elev - (state.nzi[i] + state.nzi[i + 1]) * 0.5
        elif i == state.nnodes:
            zdw = state.elev - (state.nzi[i] + state.nzi[i - 1]) * 0.5
            zup = state.elev - state.nzi[i]
        else:
            zdw = state.elev - (state.nzi[i] + state.nzi[i - 1]) * 0.5
            zup = state.elev - (state.nzi[i] + state.nzi[i + 1]) * 0.5

        t1 = state.veg_prp[veg_type, 9] * zdw
        if t1 > 50.0:
            t1 = 50.0
        t2 = state.veg_prp[veg_type, 10] * zdw
        if t2 > 50.0:
            t2 = 50.0
        t3 = state.veg_prp[veg_type, 9] * zup
        if t3 > 50.0:
            t3 = 50.0
        t4 = state.veg_prp[veg_type, 10] * zup
        if t4 > 50.0:
            t4 = 50.0

        state.rk[veg_type, i] = -0.5 * (math.exp(-t1) + math.exp(-t2) - math.exp(-t3) - math.exp(-t4))
        state.rk[veg_type, i] = _anint(state.rk[veg_type, i], 20)
        if abs(state.rk[veg_type, i]) < EPS:
            state.rk[veg_type, i] = 0.0

    # ---- somente dossel (vegetação alta) -----------------------------------
    if state.infer_test == 0:
        vind = 0
        if state.vegh_type == 3 or state.vegh_type == 4:                          # agulha
            vind = 1 if state.vegh_type == 3 else 3
            state.zh = state.veg_prp[state.vegh_type, 17]
            for j in range(1, NCLAYERS + 1):
                state.dzveg[j] = _anint(_CINFO[j + 1, vind] * state.veg_prp[state.vegh_type, 17], 2)
                state.laif[j] = _CINFO[j + 4, vind]
        elif state.vegh_type == 5 or state.vegh_type == 6:                        # folha larga
            vind = 4 if state.vegh_type == 5 else 2
            state.zh = state.veg_prp[state.vegh_type, 17]
            for j in range(1, NCLAYERS + 1):
                state.dzveg[j] = _anint(_CINFO[j + 1, vind] * state.veg_prp[state.vegh_type, 17], 2)
                state.laif[j] = _CINFO[j + 4, vind]
        elif state.vegh_type == 18:                                           # mista
            vind = 5
            state.zh = state.veg_prp[state.vegh_type, 17]
            for j in range(1, NCLAYERS + 1):
                state.dzveg[j] = _anint(_CINFO[j + 1, vind] * state.veg_prp[state.vegh_type, 17], 2)
                state.laif[j] = _CINFO[j + 4, vind]

        if (3 <= state.vegh_type <= 6) or state.vegh_type == 18:
            for i in range(1, NCLAYERS + 1):
                if abs(state.izh - SPFLAG) > EPS and abs(state.idzveg[i] - SPFLAG) <= EPS:
                    state.dzveg[i] = state.izh * state.dzveg[i] / state.zh
                elif abs(state.izh - SPFLAG) > EPS and abs(state.idzveg[i] - SPFLAG) > EPS:
                    state.dzveg[i] = state.idzveg[i]
                state.dzveg[i] = _anint(state.dzveg[i], 5)
            if abs(state.izh - SPFLAG) > EPS:
                state.zh = state.izh

            # verifica o fechamento da espessura das camadas (dossel + solo = altura total)
            for j in range(1, NCLAYERS + 1):
                sumf = sumf + state.dzveg[j]

            if abs(sumf - state.zh) > 1e-3:
                if state.single_multi_flag == 0:
                    print(" canopy + ground thickness not equal to total height")
                # equivalente ao "stop" do Fortran (interrompe a execução)
                raise SystemExit("veg_proph: canopy + ground thickness not equal to total height")
