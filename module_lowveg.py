"""Balanço hídrico e propriedades sazonais da vegetação baixa (FASST).

Fonte Fortran: module_lowveg.F90 (798 linhas, 3 sub-rotinas: lowveg_met,
low_veg_prop, veg_propl — nenhuma auxiliar interna além dessas).

O que faz
---------
:func:`lowveg_met` resolve, a cada passo de tempo, a absorção radicular e o
balanço hídrico/energético da folhagem baixa: velocidade do vento e
resistência aerodinâmica dentro do dossel (dependentes de ``sqrt_chnf`` —
ver nota de dependência abaixo), densidade do ar e razão de mistura na
folhagem (via :func:`~fasst.functions.vap_press`/:func:`~fasst.sp_humid.sp_humid`),
resistência estomática por um modelo de três fatores de estresse
(temperatura, umidade do solo ponderada pela distribuição de raízes,
déficit de pressão de vapor), interceptação de chuva/neve limitada pela
capacidade de armazenamento do dossel, degelo da neve interceptada
(dependente da densidade da neve via ``kths``), e a partição final entre
evaporação direta e transpiração, com o excesso devolvido como gotejamento
ao fluxo de precipitação efetiva.
:func:`low_veg_prop` calcula as propriedades estruturais/ópticas sazonais da
vegetação baixa (LAI, fração de cobertura do solo, altura da folhagem,
rugosidade aerodinâmica, altura de deslocamento, emissividade e albedo) por
interpolação entre valores mínimo/máximo tabelados conforme um fator de
crescimento térmico (``ftg``) e a estação do ano (hemisfério + dia do ano).
:func:`veg_propl` carrega os 17 parâmetros físicos de um tipo de vegetação de
uma tabela padrão (e, opcionalmente, de uma fonte de bioma alternativa —
MODIS-NOAH ou UMD) e calcula a fração de distribuição de raízes por camada de
solo (``rk``) via uma função de decaimento duplo-exponencial.

Arquitetura
-------------
As três funções recebem `state: FasstState` como argumento explícito.
Constantes físicas (RD, RV, TREF, EPS, LHFUS, VK, SPFLAG) vêm de
`fasst.constants`; índices de coluna do meteorológico (`ip_doy`, `ip_ws`,
`ip_sd`, `ip_hi` no Fortran) são `MetCol.DOY`/`WS`/`SD`/`HI`, de
`fasst.indices`. `dense`/`spheats` (funções puras, sem estado) e
`vap_press` (recebe `state` explícito) vêm de `fasst.functions`;
`sp_humid` (pura) vem de `fasst.sp_humid`.

As 4 variáveis `save::` do Fortran (`stlli`, `stlsi`, `interc1o`,
`interco`) formam o dataclass :class:`LowVegState`, passado explicitamente
a `lowveg_met` e mutado in-place.

Indexação: índice de nó do perfil (`nsoilp`, `soil_moist`, `stt`, `frl`,
`nzi`, e a segunda dimensão de `rk`) é 1-based, `for i in range(1,
nnodes+1)`. Colunas de `veg_prp` (as 17 descritas no cabeçalho de
:func:`veg_propl`) e o **tipo de vegetação** (1..18, primeira dimensão de
`veg_prp`/`rk`) também são 1-based, direto como no Fortran — mesma
convenção usada por `initprofile.py`/`soil_moisture.py` para essas mesmas
arrays em `fasst/state.py` (`veg_prp`, `rk` e `nsoilp` são alocadas com
folga de índice 0 em todas as dimensões). `dmet` (lista local, não parte
de `FasstState`) é 0-based (`dmet[0]` == Fortran `dmet(1)`), por não ser
uma array compartilhada com o resto do projeto.

`_anint` (duplicada por arquivo — mesma convenção de
`fasst/functions.py`/`sp_humid.py`/`snow.py`, não um `fortran_compat`
compartilhado) arredonda para o inteiro mais próximo, metade para longe
de zero (via `math.copysign`), reproduzindo `ANINT` do Fortran inclusive
para valores negativos. `_fdiv` (mesma função de `sp_humid.py`/`snow.py`)
protege duas divisões cujo denominador pode chegar a 0.0: `f2 =
(f2a-min_wat)/(max_wat-min_wat)` (se `max_wat==min_wat`) e `ff =
(stls/max_wets)`/`(stll/max_wetl)` dentro do expoente 2/3 (se
`max_wets`/`max_wetl` ficarem em 0.0 por `pdens1<=eps` nunca ter
disparado a atribuição real) — ver função para o comentário pontual.

Estado persistente entre chamadas (equivalente ao ``save`` do Fortran):
``stlli``, ``stlsi``, ``interco``, ``interc1o`` são campos de
:class:`LowVegState`, passada explicitamente a :func:`lowveg_met` e mutada
in-place — substituindo a cláusula ``save:: stlli,stlsi,interc1o,interco``
de ``lowveg_met`` no original.

Dependência de ordem de chamada (comportamento herdado do Fortran):
``lowveg_met`` usa ``state.sqrt_chnf`` logo no cálculo do vento na folhagem
(``uaf``), mas essa variável só é recalculada dentro de ``low_veg_prop``.
``low_veg_prop`` precisa ser chamada antes de ``lowveg_met`` a cada passo
de tempo para que ``sqrt_chnf`` esteja atualizado.

Anomalia herdada do Fortran (comentário, não bug funcional): o cabeçalho de
:func:`veg_propl` no ``.F90`` rotula a coluna 12 de ``veg_prp`` como
"emissivity max" e a 13 como "emissivity min" — mas o código de leitura e de
uso (em ambas as sub-rotinas) sempre grava/lê ``emissmin`` na coluna 12 e
``emissmax`` na coluna 13, o oposto do rótulo. O comportamento numérico é
consistente (mín/máx corretos como limites de interpolação em
``low_veg_prop``); é só a legenda da coluna que está trocada no original.
Preservado aqui na mesma posição, com nota explícita no ponto do código.

Ponto em aberto: o quinto argumento passado a `sp_humid` em `lowveg_met`
(`ph`, fisicamente uma cabeça de pressão em metros) recebe `c1`, que é a
pressão de vapor em Pa calculada pela linha anterior (`vap_press`). A
correspondência com o argumento que o `module_lowveg.F90` original passa
nessa posição não está confirmada contra o arquivo fonte — ver nota no
ponto do código.
"""

import math
from dataclasses import dataclass

import numpy as np

from .constants import EPS, LHFUS, RD, RV, SPFLAG, TREF, VK
from .functions import Phase, dense, spheats, vap_press
from .indices import MetCol
from .sp_humid import sp_humid
from .state import FasstState

__all__ = ["LowVegState", "lowveg_met", "low_veg_prop", "veg_propl"]


@dataclass
class LowVegState:
    """Variáveis `save::` de `lowveg_met` (module_lowveg.F90, linhas 31-32).

    Persistem entre chamadas sucessivas para o mesmo passo/iteração (ii,
    iter_idx) — só são reatribuídas no primeiro passo (ii==0 and
    iter_idx==0) ou zeradas quando icaseo==4 (ver :func:`lowveg_met`).
    """

    stlli: float = 0.0
    stlsi: float = 0.0
    interc1o: float = 0.0
    interco: float = 0.0


def _anint(x: float, p: int = 0) -> float:
    """Fortran ``ANINT(x*10**p)*10**-p`` -- round half *away* from zero.

    Duplicada por arquivo (não importada de um `fortran_compat`
    compartilhado -- este projeto não tem um; ver `fasst/functions.py`,
    `sp_humid.py` e `snow.py`, que fazem o mesmo). NaN/inf passam
    inalterados, como o `ANINT` do Fortran.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def _fdiv(x: float, y: float) -> float:
    """Divisão IEEE-754 pura: ``x / 0.0`` -> +-inf / nan, sem levantar.

    Mesma função de `sp_humid.py`/`snow.py`, duplicada aqui pela mesma
    convenção. Usada nos dois pontos deste arquivo onde o Fortran tem uma
    divisão sem guarda que o Python, ao contrário do Fortran, faria
    estourar -- ver docstring do módulo.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


# ---------------------------------------------------------------------------
# lowveg_met  (module_lowveg.F90, linhas 5-334)
# ---------------------------------------------------------------------------
def lowveg_met(
    state: FasstState,
    lowveg: LowVegState,
    ii: int,
    iter_idx: int,
    pt1: int,
    pt2: int,
    mixrgr: float,
    wetness: float,
    taf: float,
    dmet: list,
    wetbulba: float,
):
    """
    Calcula a absorção radicular de vegetação baixa, interceptação de
    precipitação, evaporação e transpiração.

    Saídas devolvidas em tupla (rhoaf, cf, pheatf, rpp, mixra, mixrf,
    dqdtf) porque no Fortran são todas intent(out); dmet é intent(inout)
    no original e aqui é mutado in-place na própria lista recebida;
    wetbulba não é reatribuído dentro desta função.
    """
    f1i = 0
    c1 = 0.0
    c2 = 0.0
    vpressa = 0.0
    pdens = 0.0
    pdens1 = 0.0
    mixra = 0.0
    qaf = 0.0
    rhoaf = 0.0
    rpf = 0.0
    rptr = 0.0
    ra = 0.0
    ep = 0.0
    ef = 0.0
    etr = 0.0
    ff = 0.0
    d1 = 0.0
    rpp = 0.0
    t1 = 0.0
    cf = 0.0
    spheat = 0.0
    spheat1 = 0.0
    min_wat = 0.0
    max_wat = 0.0
    kths = 0.0
    f1 = 0.0
    f1a = 0.0
    f2a = 0.0
    f2 = 0.0
    f3 = 0.0
    max_wet = 0.0
    max_wetl = 0.0
    max_wets = 0.0
    interc = 0.0
    interc1 = 0.0
    drip = 0.0
    drip1 = 0.0
    meltlv = 0.0
    pheatf = 0.0

    state.stll = state.storll
    state.stls = state.storls

    if state.icaseo != 4 and (ii == 0 and iter_idx == 0):
        lowveg.stlli = state.stll
        lowveg.stlsi = state.stls
    elif state.icaseo == 4:
        state.stll = 0.0
        state.stls = 0.0

    # Velocidade do vento na folhagem (m/s). sqrt_chnf vem de low_veg_prop
    # (ver nota de dependência de ordem de chamada no docstring do módulo).
    state.uaf = 0.83 * state.sigfl * dmet[2] * state.sqrt_chnf + (1.0 - state.sigfl) * dmet[2]

    if state.uaf > 0.2:
        cf = (1.0 + 0.3 / state.uaf) * 0.01
        ra = 1.0 / (cf * state.uaf)
    else:
        cf = 0.01
        ra = 0.0

    cf = _anint(state.sigfl * cf, 20)

    c1 = vap_press(state, state.nnodes + 1, dmet[4] * 0.01, dmet[10])
    rhoaf = (c1 / RV + (dmet[10] * 1e2 - c1) / RD) / taf
    rhoaf = _anint(rhoaf, 20)

    if state.hsaccum > EPS or state.hi > EPS or state.stls > EPS:
        f1i = 1

    c2 = 1.0
    # ph=c1 (pressão de vapor, Pa) -- ver "Ponto em aberto" no docstring
    # do módulo: correspondência com o argumento do Fortran original não
    # confirmada.
    moist = sp_humid(f1i, dmet[10], state.ftemp, c2, c1)
    mixrf = moist.mixr
    dqdtf = moist.dqdt
    vpressf = moist.vpress
    wetbulbf = moist.wetbulb  # calculado, não utilizado adiante (ver Fortran)

    # Resistência estomática rs e transpiração máxima
    # (dmet(N) do Fortran == dmet[N-1] aqui, em todo o arquivo: dmet(1)->dmet[0]
    # acima, dmet(3)->dmet[2] no vento, dmet(5)/(11)->dmet[4]/[10] a seguir, etc.)
    f1a = (4e-3 * dmet[0] + 5e-3) / (0.81 * (4e-3 * dmet[0] + 1.0))
    f1 = 1.0 / min(1.0, max(EPS, f1a))

    # Laço sobre os nós do perfil de solo, 1-based (Fortran: do i=1,nnodes).
    # rk[vegl_type, i]: linha = tipo de vegetação, coluna = nó do solo,
    # ambos 1-based (ver docstring do módulo, "Indexação").
    state.trmlm = 0.0
    for i in range(1, state.nnodes + 1):
        if state.stt[i] > TREF:
            state.frl[i] = state.rk[state.vegl_type, i] * (
                1.0 - (state.soil_moist[i] - state.nsoilp[i, 16]) / (state.nsoilp[i, 17] - state.nsoilp[i, 16])
            )
            state.frl[i] = max(0.0, min(1.0, _anint(state.frl[i], 20)))
        else:
            state.frl[i] = 0.0

        state.trmlm += state.frl[i]
        f2a += state.rk[state.vegl_type, i] * state.soil_moist[i]
        min_wat += state.nsoilp[i, 16]
        max_wat += state.nsoilp[i, 17]

    # _fdiv: max_wat==min_wat é possível (todas as camadas com o mesmo
    # ponto de murcha/capacidade de campo) -- o Fortran deixaria a divisão
    # virar Infinity silenciosamente; ver docstring do módulo.
    f2 = max(0.0, min(1.0, _fdiv(f2a - min_wat, max_wat - min_wat)))
    if f2 > EPS:
        f2 = 1.0 / f2

    ff = 1.0
    state.trmlm = 1.5e-7 * ff * state.trmlm

    t1 = 3e-4 * (vpressf - vpressa)
    if abs(vpressf - vpressa) > EPS and t1 < 50.0:
        f3 = math.exp(t1)
    else:
        f3 = 1.0

    # veg_prp[vegl_type, 1] == veg_prp(vegl_type,1) no Fortran: resistência
    # estomática mínima (srmin), coluna 1 da tabela de veg_propl.
    if state.lail > EPS and state.ftemp > TREF:
        state.rsl = max(0.0, (state.veg_prp[state.vegl_type, 1] / state.lail) * f1 * f2 * f3)
    else:
        state.rsl = 0.0

    if (ra + state.rsl) > EPS:
        rpp = ra / (ra + state.rsl)
    rpp = min(max(0.0, _anint(rpp, 20)), 1.0)

    d1 = 1.0 - state.sigfl * (0.6 * (1.0 - rpp) + 0.1 * (1.0 - wetness))
    qaf = (
        (1.0 - 0.7 * state.sigfl) * mixra
        + 0.6 * state.sigfl * rpp * mixrf
        + 0.1 * state.sigfl * wetness * mixrgr
    ) / d1

    ep = state.lail * cf * (rhoaf / dense(wetbulba, 0.0, Phase.WATER)) * state.uaf * (qaf - rpp * mixrf)

    pdens = dense(wetbulba, 0.0, Phase.WATER)
    spheat = spheats(wetbulba, Phase.WATER)

    pdens1 = dense(wetbulba, dmet[2], Phase.SNOW)
    spheat1 = spheats(wetbulba, Phase.ICE)

    # veg_prp[...,7]/[...,8] == veg_prp(vegl_type,7)/(vegl_type,8) no Fortran:
    # capacidade de armazenamento por LAI (ddmax) e índice de área de caule
    # (sai) — não é um fator fixo, vem da tabela por tipo de vegetação.
    max_wet = state.veg_prp[state.vegl_type, 7] * (state.lail + state.veg_prp[state.vegl_type, 8])
    max_wetl = _anint(max_wet * 1e-3, 20)

    if pdens1 > EPS:
        max_wets = (max_wet * 1e-3) * dense(state.ftemp, 0.0, Phase.WATER) / pdens1
        max_wets = _anint(max_wets, 20)

    t1 = 0.5 * (state.lail + state.veg_prp[state.vegl_type, 8])
    if t1 > 50.0:
        t1 = 50.0

    if pt1 == 2 or pt1 == 4:  # Chuva
        interc = dmet[5] * (1.0 - math.exp(-t1))
        interc = max(0.0, min(interc, max_wetl))
    elif pt1 == 3:  # Neve
        interc1 = dmet[5] * (1.0 - math.exp(-t1))
        interc1 = max(0.0, min(interc1, max_wets))
    elif pt2 == 3:  # Neve
        interc1 = dmet[6] * (1.0 - math.exp(-t1))
        interc1 = max(0.0, min(interc1, max_wets))

    interc = _anint(interc, 20)
    interc1 = _anint(interc1, 20)
    lowveg.interco = interc
    lowveg.interc1o = interc1

    state.stls += interc1
    if state.stls > max_wets:
        drip1 = state.stls - max_wets
        state.stls = max_wets

    drip1 = _anint(drip1, 20)
    state.stls = _anint(state.stls, 20)
    if state.stls < 1e-4:
        state.stls = 0.0

    # Degelo da neve interceptada no dossel. kths: condutividade térmica
    # efetiva da neve em função da densidade (pdens1), duas correlações
    # empíricas conforme a faixa de densidade (Sturm et al., citado no
    # Fortran original); meltlv é limitado ao estoque de neve disponível
    # (convertido pela razão de densidades pdens/pdens1) logo abaixo.
    meltlv = 0.0
    if state.stls > 0.0:
        if pdens1 < 0.156:
            kths = 0.023 + 0.234 * pdens1
            kths = min(max(0.023, kths), 1.0)
        else:
            kths = 0.138 - 1.01 * pdens1 + 3.233 * pdens1 * pdens1
            kths = min(max(0.138, kths), 1.0)

        meltlv = max(
            0.0,
            (kths / state.stls) * (dmet[3] - TREF) * (state.deltat_fasst / (LHFUS * pdens1)),
        )
        if meltlv > state.stls * pdens / pdens1:
            meltlv = state.stls * pdens / pdens1

        meltlv = _anint(meltlv, 20)
        state.stls -= meltlv
        state.stls = _anint(state.stls, 20)
        if state.stls < 1e-4:
            state.stls = 0.0

    state.stll += interc + meltlv
    drip = 0.0
    if state.stll > max_wetl:
        drip = state.stll - max_wetl
        state.stll = max_wetl

    drip = _anint(drip, 20)
    state.stll = _anint(state.stll, 20)
    if state.stll < 1e-6:
        state.stll = 0.0

    if rpp * mixrf >= qaf or state.trmlm <= EPS:
        if f2 > min_wat and ep <= EPS:
            ep = 0.0

    if (state.rsl + ra) > EPS:
        if state.stls > EPS:
            # _fdiv: max_wets pode ficar em 0.0 (nunca atribuído de fato,
            # se pdens1<=eps) enquanto stls>eps -- ver docstring do módulo.
            ff = _fdiv(state.stls, max_wets) ** (2.0 / 3.0)
        elif state.stll > EPS:
            ff = _fdiv(state.stll, max_wetl) ** (2.0 / 3.0)

        rpf = 1.0 - (state.rsl / (state.rsl + ra)) * (1.0 - ff)
        rptr = (ra / (state.rsl + ra)) * (1.0 - ff)

    rptr = min(max(0.0, _anint(rptr, 20)), 1.0)
    rpf = min(rpf, 1.0)
    rpf = _anint(rpf, 20)

    ef = rpf * ep
    etr = min(state.trmlm, rptr * ep)

    if state.stls > EPS:
        state.stls -= (ef - etr) * state.deltat_fasst
        if state.stls < 1e-4:
            state.stls = 0.0
        if state.stls > max_wets:
            drip1 += state.stls - max_wets
            state.stls = max_wets
    elif state.stll > EPS:
        state.stll -= (ef - etr) * state.deltat_fasst
        if state.stll < 1e-6:
            state.stll = 0.0
        if state.stll > max_wetl:
            drip += state.stll - max_wetl
            state.stll = max_wetl

    drip = max(0.0, _anint(drip, 20))
    state.stll = max(0.0, _anint(state.stll, 20))
    drip1 = max(0.0, _anint(drip1, 20))
    state.stls = max(0.0, _anint(state.stls, 20))

    # Gotejamento redistribuído de volta ao fluxo de precipitação efetiva
    # que chega ao solo (dmet(6)/dmet(7) no Fortran -> dmet[5]/dmet[6] aqui).
    # Nota de precisão: aqui o arredondamento usa k=15, não k=20 como no
    # restante da função — reproduz o Fortran original (anint(...*1d15)),
    # não é uma inconsistência da tradução.
    dmet[5] += state.sigfl * drip / (state.timstep * 3600.0)
    dmet[6] += state.sigfl * drip1 / (state.timstep * 3600.0)

    dmet[5] = _anint(dmet[5], 15)
    dmet[6] = _anint(dmet[6], 15)

    pheatf = state.sigfl * (spheat * (pdens * interc) + spheat1 * (pdens1 * interc1))
    pheatf = _anint(pheatf, 20)

    return rhoaf, cf, pheatf, rpp, mixra, mixrf, dqdtf


# ---------------------------------------------------------------------------
# low_veg_prop  (module_lowveg.F90, linhas 337-567)
# ---------------------------------------------------------------------------
def low_veg_prop(state: FasstState, fatemp: float, c3_hgt: float):
    """
    Calcula as propriedades estruturais e ópticas da vegetação baixa
    conforme a estação.

    c3_hgt é recebido mas não utilizado: no Fortran original a única
    referência a esse parâmetro está comentada (linhas 470-471, ajuste de
    hfol contra a altura da camada basal do dossel alto), então ele
    permanece um argumento morto também aqui — preservado por fidelidade
    de assinatura, não por engano.
    """
    season = 1
    ftg = 0.0
    hgt1 = 0.0

    # Classificação de estação por hemisfério (sinal de lat) e dia do ano;
    # a numeração 1=inverno,2=primavera,3=verão,4=outono e os limiares de
    # doy batem exatamente com o Fortran (aint(met(iw,ip_doy)) lá, int()
    # aqui — mesmo efeito de truncamento em direção a zero para doy>=0).
    doy = int(state.met[state.iw, MetCol.DOY])
    if state.lat >= 0.0:  # Hemisfério Norte
        if doy >= 335 or doy <= 80:
            season = 1  # Inverno
        elif 81 <= doy <= 151:
            season = 2  # Primavera
        elif 152 <= doy <= 243:
            season = 3  # Verão
        elif 244 <= doy <= 334:
            season = 4  # Outono
    else:  # Hemisfério Sul
        if doy >= 335 or doy <= 80:
            season = 3  # Verão
        elif 81 <= doy <= 151:
            season = 4  # Outono
        elif 152 <= doy <= 243:
            season = 1  # Inverno
        elif 244 <= doy <= 334:
            season = 2  # Primavera

    if state.iw == state.istart and state.infer_test == 0:
        state.iseason = season

    ftg = 1.0 - 1.6e-3 * (298.0 - state.stt[state.refn]) * (298.0 - state.stt[state.refn])
    ftg = _anint(ftg, 20)

    sigfl_old = state.sigfl
    lail_old = state.lail
    hfol_old = state.hfol
    hfoltot_old = state.hfol_tot
    epf_old = state.epf
    fola_old = state.fola

    # Split por tipo de vegetação, igual ao Fortran: para vtype==2 (short
    # grass) o LAI é calculado já aqui e sigfl deriva dele por uma relação
    # exponencial; para os demais tipos, sigfl vem direto da tabela
    # (colunas cmax/cmin) e o LAI correspondente só é calculado mais abaixo,
    # dentro do bloco "if sigfl > eps" (branch "if vtype != 2"). Não é
    # código duplicado nem morto — são os dois ramos mutuamente exclusivos
    # do mesmo cálculo, na mesma ordem em que aparecem no .F90.
    vtype = state.vegl_type
    if vtype == 2:
        state.lail = state.veg_prp[vtype, 6] + ftg * (state.veg_prp[vtype, 5] - state.veg_prp[vtype, 6])
        if state.lail > state.veg_prp[vtype, 5]:
            state.lail = state.veg_prp[vtype, 5]
        if state.lail < state.veg_prp[vtype, 6] and state.lail < state.ilail:
            state.lail = state.veg_prp[vtype, 6]
        if state.iw == state.istart and state.infer_test == 0:
            state.ilail = state.lail
        if season < state.iseason and state.lail > state.ilail:
            state.lail = state.ilail
        if abs(-0.75 * state.lail) < 50.0:
            state.sigfl = 1.0 - math.exp(-0.75 * state.lail)
    else:
        state.sigfl = state.veg_prp[vtype, 2] - (1.0 - ftg) * (state.veg_prp[vtype, 2] - state.veg_prp[vtype, 3])
        if state.sigfl < state.veg_prp[vtype, 2]:
            state.sigfl = state.veg_prp[vtype, 2]
        if state.sigfl > state.veg_prp[vtype, 3]:
            state.sigfl = state.veg_prp[vtype, 3]
        state.sigfl *= 1e-2

    if season == state.iseason and abs(state.isigfl - SPFLAG) > EPS:
        state.sigfl = state.isigfl
    if (season < state.iseason and state.sigfl > state.isigfl) and abs(state.isigfl - SPFLAG) > EPS:
        state.sigfl = state.isigfl

    if state.sigfl > EPS:
        # Altura da folhagem por grupo de tipo (mesma partição do Fortran):
        # vtype 2/8 (grama curta/deserto) e 16/17 (arbustos) têm regras
        # próprias; os demais tipos (culturas, grama alta, tundra etc.)
        # caem no "else", com a mesma fórmula usada para 2/8 mas indexada
        # de novo por conveniência de leitura, não por acidente.
        if vtype in (2, 8):
            state.hfol = state.veg_prp[vtype, 17] - (1.0 - ftg) * (state.veg_prp[vtype, 17] - state.veg_prp[vtype, 16])
            if state.hfol < state.veg_prp[vtype, 16]:
                state.hfol = state.veg_prp[vtype, 16]
            if state.hfol > state.veg_prp[vtype, 17]:
                state.hfol = state.veg_prp[vtype, 17]
        elif vtype in (16, 17):
            state.hfol = state.veg_prp[vtype, 16]
        else:
            state.hfol = state.veg_prp[vtype, 17] - (1.0 - ftg) * (state.veg_prp[vtype, 17] - state.veg_prp[vtype, 16])
            if state.hfol < state.veg_prp[vtype, 16]:
                state.hfol = state.veg_prp[vtype, 16]
            if state.hfol > state.veg_prp[vtype, 17]:
                state.hfol = state.veg_prp[vtype, 17]

        if season == state.iseason and abs(state.ihfol - SPFLAG) > EPS:
            state.hfol = state.ihfol
        if season < state.iseason and state.hfol > state.ihfol:
            state.hfol = state.ihfol

        state.hfol_tot = state.hfol * 0.01
        if state.iw != 0:
            state.hfol_tot = 0.5 * (state.hfol_tot + hfoltot_old)
        state.hfol_tot = _anint(state.hfol_tot, 20)

        state.hfol = max(0.0, state.hfol * 0.01 - state.hm)

        if state.iheightn <= state.hfol:
            state.hfol = state.iheightn * 0.5

        if state.iw >= 1:
            if abs(state.met[state.iw - 1, MetCol.SD] + state.met[state.iw - 1, MetCol.HI]) < EPS:
                state.hfol = 0.5 * (state.hfol + hfol_old)

        state.hfol = _anint(state.hfol, 20)

        hgt1 = state.iheightn - state.hm
        if hgt1 <= state.hfol:
            hgt1 = state.hfol

        # Rugosidade aerodinâmica (z0l), altura de deslocamento (zd) e
        # sqrt_chnf: usados por lowveg_met no cálculo do vento na folhagem
        # (uaf) — ver nota de dependência de ordem de chamada no docstring
        # do módulo. z0l = 7.775e-3 m é o valor fixo de neve/solo nu.
        if state.hfol <= 0.0:  # Neve no solo
            state.z0l = 7.775e-3
            state.zd = state.z0l
            state.sqrt_chnf = 0.0
        else:
            state.z0l = 0.131 * (state.hfol**0.997)
            state.z0l = min(state.z0l, state.veg_prp[vtype, 4])
            state.zd = 0.701 * (state.hfol**0.979)
            state.sqrt_chnf = VK / math.log((hgt1 - state.zd) / state.z0l)

        state.z0l = _anint(state.z0l, 20)
        state.zd = _anint(state.zd, 20)
        state.sqrt_chnf = _anint(state.sqrt_chnf, 20)

        if vtype != 2:
            state.lail = state.veg_prp[vtype, 6] + ftg * (state.veg_prp[vtype, 5] - state.veg_prp[vtype, 6])
            if state.lail > state.veg_prp[vtype, 5]:
                state.lail = state.veg_prp[vtype, 5]
            if state.lail < state.veg_prp[vtype, 6] and state.lail < state.ilail:
                state.lail = state.veg_prp[vtype, 6]
            if state.iw == state.istart and state.infer_test == 0:
                state.ilail = state.lail
            if season < state.iseason and state.lail > state.ilail:
                state.lail = state.ilail

        # veg_prp[...,12]/[...,13] == colunas 12/13 do Fortran. O cabeçalho
        # de veg_propl rotula a col.12 como "emissivity max" e a 13 como
        # "emissivity min", mas o código sempre grava emissmin na 12 e
        # emissmax na 13 (ver docstring do módulo) — por isso aqui col.12
        # funciona como piso (mín) e col.13 como teto (máx), como abaixo.
        state.epf = state.veg_prp[vtype, 12] + ftg * (state.veg_prp[vtype, 13] - state.veg_prp[vtype, 12])
        if state.epf < state.veg_prp[vtype, 12]:
            state.epf = state.veg_prp[vtype, 12]
        if state.epf > state.veg_prp[vtype, 13]:
            state.epf = state.veg_prp[vtype, 13]
        if season == state.iseason and abs(state.iepf - SPFLAG) > EPS:
            state.epf = state.iepf
        if season < state.iseason and state.epf > state.iepf:
            state.epf = state.iepf

        # fola: componente de reflectância da folhagem (colunas 14/15 do
        # Fortran = albedo min/max, sem a inversão de rótulo do bloco acima
        # — aqui a legenda bate com o uso). albf = 1 - fola é o albedo
        # final da folhagem, calculado logo após o "else" deste if/else.
        state.fola = state.veg_prp[vtype, 15] - (1.0 - ftg) * (state.veg_prp[vtype, 15] - state.veg_prp[vtype, 14])
        if state.fola < state.veg_prp[vtype, 14]:
            state.fola = state.veg_prp[vtype, 14]
        if state.fola > state.veg_prp[vtype, 15]:
            state.fola = state.veg_prp[vtype, 15]

        if state.iw == state.istart and state.infer_test == 0:
            if abs(state.ifola - SPFLAG) > EPS:
                state.fola = state.ifola
            state.ftemp = fatemp
            state.uaf = state.met[state.iw, MetCol.WS]
        else:
            state.lail = 0.5 * (lail_old + state.lail)
            state.sigfl = 0.5 * (sigfl_old + state.sigfl)
            state.fola = 0.5 * (fola_old + state.fola)
            state.epf = 0.5 * (epf_old + state.epf)

        if state.hfol <= 0.0 and state.hm <= EPS:
            state.sigfl = 0.0
            state.lail = 0.0
            state.fola = 0.0
            state.epf = 0.0

        state.albf = 1.0 - state.fola
    else:
        state.epf = 0.0
        state.albf = 0.0
        state.lail = 0.0
        state.hfol = 0.0
        state.hfol_tot = 0.0
        hgt1 = state.iheightn - state.hm
        state.z0l = 0.0
        state.zd = 0.0
        state.sqrt_chnf = 0.0

    state.sigfl = _anint(state.sigfl, 20)
    state.lail = _anint(state.lail, 20)
    state.fola = _anint(state.fola, 20)
    state.albf = _anint(state.albf, 20)
    state.epf = _anint(state.epf, 20)


# ---------------------------------------------------------------------------
# veg_propl  (module_lowveg.F90, linhas 570-795)
# ---------------------------------------------------------------------------
def veg_propl(
    state: FasstState,
    biome_source: int,
    new_vt: int,
    veg_type: int,
    file_unit_31_data: list,
    file_unit_biome_data: list = None,
):
    """
    Carrega parâmetros de vegetação de arquivos de dados/tabelas e calcula
    a fração de distribuição de raízes por camada de solo (rk).

    file_unit_31_data / file_unit_biome_data substituem os arquivos por
    unidade lógica (31, e 32 ou 33 conforme biome_source) do Fortran: aqui
    chegam como listas de linhas de texto já lidas, em vez de leitura
    direta por unit number.
    """
    # Tabelas locais 1-based (índice 0 sem uso, tipo de vegetação 1..18 ou
    # 1..ntypes, coluna de propriedade 1..17) — mesma folga de índice 0 em
    # toda dimensão que fasst/state.py usa para veg_prp/rk (ver docstring
    # do módulo, "Indexação").
    defveg_prp = np.zeros((19, 18))

    # Lê os cabeçalhos de entrada (linhas 1 a 36 no Fortran)
    lines = file_unit_31_data[36:]
    for line in lines:
        parts = line.split()
        if not parts:
            continue
        vid = int(parts[0])
        srmax, srmin, cmax, cmin, rl = map(float, parts[1:6])
        laimax, laimin, ddmax, sai, ar, br = map(float, parts[6:12])
        emissmin, emissmax, folamin, folamax, heigmin, heigmax = map(float, parts[12:18])

        if veg_type == vid:
            # Colunas 1..17, identicas ao Fortran (ver nota sobre a coluna
            # 12/13 abaixo e o docstring do módulo).
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
            # Note col. 12/13: gravam emissmin/emissmax nessa ordem, embora
            # o cabeçalho de colunas do .F90 rotule a col.12 como "max" e a
            # 13 como "min" — rótulo trocado no original, comportamento
            # preservado como está (ver docstring do módulo).
            defveg_prp[veg_type, 12] = emissmin * 1e-2
            defveg_prp[veg_type, 13] = emissmax * 1e-2
            defveg_prp[veg_type, 14] = folamin * 1e-2
            defveg_prp[veg_type, 15] = folamax * 1e-2
            defveg_prp[veg_type, 16] = heigmin
            defveg_prp[veg_type, 17] = heigmax
            break

    newveg_prp = None
    if biome_source > 0 and file_unit_biome_data is not None:
        if biome_source == 1000:  # Modis_NOAH
            hlines, ntypes = 45, 20
        elif biome_source == 2000:  # UMD
            hlines, ntypes = 31, 14
        else:
            hlines, ntypes = 0, 0

        newveg_prp = np.zeros((ntypes + 1, 18))
        lines_biome = file_unit_biome_data[hlines:]

        for line in lines_biome:
            parts = line.split()
            if not parts:
                continue
            vid = int(parts[0])
            srmax, srmin, cmax, cmin, rl = map(float, parts[1:6])
            laimax, laimin, ddmax, sai, ar, br = map(float, parts[6:12])
            emissmin, emissmax, folamin, folamax, heigmin, heigmax = map(float, parts[12:18])

            if new_vt == vid:
                newveg_prp[new_vt, 1] = srmin
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
                newveg_prp[new_vt, 16] = heigmin
                newveg_prp[new_vt, 17] = heigmax
                break

    if biome_source == 0:
        for i in range(1, 18):
            state.veg_prp[veg_type, i] = _anint(defveg_prp[veg_type, i], 10)
    else:
        for i in range(1, 18):
            if abs(newveg_prp[new_vt, i] - SPFLAG) <= EPS:
                state.veg_prp[veg_type, i] = defveg_prp[veg_type, i]
            else:
                state.veg_prp[veg_type, i] = newveg_prp[new_vt, i]
            state.veg_prp[veg_type, i] = _anint(state.veg_prp[veg_type, i], 10)

    # Cálculo da fração de raiz por camada (BATS: ar/br são os coeficientes
    # de decaimento por tipo de vegetação, colunas 9/10 do Fortran).
    # Nota de precisão: arredondamento em k=10 neste método inteiro
    # (tabela de propriedades e rk), diferente do k=20 usado no resto do
    # arquivo — reproduz o Fortran (anint(...*1d10)), não é inconsistência.
    for i in range(1, state.nnodes + 1):
        if i == 1:
            zdw = state.elev - state.nzi[1]
            zup = state.elev - (state.nzi[1] + state.nzi[2]) * 0.5
        elif i == state.nnodes:
            zdw = state.elev - (state.nzi[i] + state.nzi[i - 1]) * 0.5
            zup = state.elev - state.nzi[i]
        else:
            zdw = state.elev - (state.nzi[i] + state.nzi[i - 1]) * 0.5
            zup = state.elev - (state.nzi[i] + state.nzi[i + 1]) * 0.5

        t1 = min(50.0, state.veg_prp[veg_type, 9] * zdw)
        t2 = min(50.0, state.veg_prp[veg_type, 10] * zdw)
        t3 = min(50.0, state.veg_prp[veg_type, 9] * zup)
        t4 = min(50.0, state.veg_prp[veg_type, 10] * zup)

        state.rk[veg_type, i] = -0.5 * (math.exp(-t1) + math.exp(-t2) - math.exp(-t3) - math.exp(-t4))
        state.rk[veg_type, i] = _anint(state.rk[veg_type, i], 10)

        if abs(state.rk[veg_type, i]) < EPS:
            state.rk[veg_type, i] = 0.0
