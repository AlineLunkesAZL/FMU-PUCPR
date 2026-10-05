"""Radiação de onda longa (IR) e onda curta (solar) incidentes sobre a superfície (FASST).

Fonte Fortran: module_radiation.F90 (724 linhas, 6 sub-rotinas, todas
internas ao módulo — nenhuma chamada externa). Contém: fluxo IV
descendente e altura de base de nuvem (:func:`dnirflx` +
:func:`cloudbase`), emissividade atmosférica (:func:`emisatm`), posição
solar zênite/azimute (:func:`sol_zen`) e radiação solar direta/difusa
(:func:`Solflx` + :func:`insol`, método de Shapiro).

Arquitetura
-------------
Cada uma das 6 sub-rotinas é uma função que recebe `state: FasstState`
como primeiro argumento. `lat`, `mlong`, `timeoffset`, `albedo_fasst` e
`mflag` são campos de `state` (mflag é lido do cabeçalho do arquivo
meteorológico em tempo de execução); `EPS`, `TREF`, `SIGMA` e `PI` são
constantes de `fasst.constants`. Este é o único arquivo do projeto sem
nenhuma dependência de outro arquivo do FASST.

Parâmetros ``intent(inout)`` do Fortran (``mcld``/``hcld`` em
``dnirflx``; ``icld``/``cover``/``hgt``/``tsol`` em ``Solflx``;
``icld``/``cover``/``hgt`` em ``insol``) são tratados por cópia local +
retorno em tupla, e não por mutação in-place — quem chama precisa
reatribuir o retorno.

As matrizes de coeficientes do antigo ``blockdata dinsol`` (``r1``,
``r2``, ``t1``, ``t2``, ``wt``) foram transcritas linha a linha na
mesma ordem em que aparecem no ``DATA`` do Fortran; como o Fortran
preenche essas matrizes em ordem column-major, o array ``numpy``
resultante é a transposta da matriz Fortran — o padrão de acesso em
:func:`insol` (``r1[coef_idx, layer_idx]`` em vez de
``r1[layer_idx, coef_idx]``) compensa essa transposição.

Dois bugs do Fortran original, preservados de propósito
----------------------------------------------------------
* :func:`cloudbase` — no cálculo de ``ilat``, o Fortran testa
  ``plat >= dabs(25d0)``, isto é, ``dabs`` de uma constante (25.0), não
  de ``plat``. Isso equivale a ``plat >= 25`` puro e simples, ignorando
  o sinal da latitude. Para latitudes negativas (hemisfério sul) abaixo
  de -25°, o Fortran nunca seleciona o conjunto de coeficientes de
  ``ilat`` de latitude alta — mesmo estando geograficamente lá.
* :func:`cloudbase` — no cálculo de ``isean``, a condição para inverno
  no hemisfério sul (``int(doy) > 150 .or. int(doy) < 250``) usa
  ``.or.`` onde aparentemente deveria haver um ``.and.``: como 150 <
  250, a condição é sempre verdadeira para qualquer dia do ano, então
  o hemisfério sul é sempre classificado como "inverno" nessa
  parametrização climatológica, independente da estação real.

Decisão de tradução (não é bug do Fortran): :func:`sol_zen`
-------------------------------------------------------------
O Fortran calcula ``phi = acos(cosphi)`` sem checar o domínio de
``cosphi``/``costheta``; se o argumento escapar de [-1, 1] por erro de
ponto flutuante, o comportamento típico é retornar ``NaN`` e seguir a
execução. O equivalente direto em Python, ``math.acos``, lança
``ValueError`` nesse caso em vez de propagar ``NaN`` — por isso aqui
usamos ``numpy.arccos`` (sem clamping), que replica o comportamento
"NaN silencioso" do Fortran em vez de interromper a execução.

`_fdiv` (mesma função em `fasst/sp_humid.py`/`snow.py`/
`module_lowveg.py`/`module_canopy.py`) substitui a divisão comum em
dois pontos onde o denominador pode ser exatamente 0.0: a divisão
final de `sdowns` por ela mesma em :func:`insol` (o acumulado de
reflectâncias/transmitâncias pode chegar a 0.0 em condições de
noite/pôr-do-sol, com `sdown0` também 0.0) e `costheta` em
:func:`sol_zen` (denominador `cos(zlat)*sin(phi)`, que pode ser 0.0 no
polo ou com o sol exatamente no zênite/nadir). O Fortran propagaria
+-inf/nan silenciosamente nesses casos; a divisão comum do Python
levantaria `ZeroDivisionError`.
"""

import math
import numpy as np

from .constants import EPS, PI, SIGMA, TREF
from .state import FasstState


def _aint(x: float, p: int = 0) -> float:
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direção a zero.

    Duplicada por arquivo -- este projeto não tem um `fortran_compat`
    compartilhado (ver `fasst/functions.py`, `sp_humid.py`, `snow.py`,
    `module_lowveg.py`, `module_canopy.py`).
    """
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _fdiv(x: float, y: float) -> float:
    """Divisão IEEE-754 pura: ``x / 0.0`` -> +-inf / nan, sem levantar.

    Mesma função de `sp_humid.py`/`snow.py`/`module_lowveg.py`/
    `module_canopy.py`.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


# ---------------------------------------------------------------------------
# dnirflx  (module_radiation.F90, linhas 5-69)
# ---------------------------------------------------------------------------
def dnirflx(state: FasstState, ematm: float, ta: float, doy: float, lcldbse: float,
            mcldbse: float, hcldbse: float, lcld: float,
            mcld: float, hcld: float):
    """
    Calcula o fluxo de radiação infravermelha (onda longa) descendente (flxir).

    `mcld`/`hcld` são `intent(inout)` no Fortran: a sub-rotina
    primeiro usa os valores originais para achar a base das nuvens
    (via `cloudbase`) e só depois os sobrescreve com a cobertura
    efetiva (overlap aleatório), na mesma ordem do Fortran. Aqui
    isso vira retorno em tupla (`mcld_eff`, `hcld_eff`) em vez de
    mutação in-place.
    """
    # Substitui a modificação inout do mcld e hcld sem alterar o escopo externo
    mcld_local = float(mcld)
    hcld_local = float(hcld)

    # 1. Calcula as alturas das bases das nuvens
    zlcld, zmcld, zhcld, lcld_out, mcld_local, hcld_local = cloudbase(
        state, doy, state.lat, lcldbse, mcldbse, hcldbse, lcld, mcld_local, hcld_local
    )

    # 2. Coberturas efetivas de nuvens médias e altas (overlap aleatório)
    # (usa os valores de mcld/hcld anteriores à sobreposição, na ordem do Fortran)
    hcld_eff = hcld_local * (1.0 - mcld_local) * (1.0 - lcld)
    mcld_eff = mcld_local * (1.0 - lcld)

    # 3. Fluxo de céu limpo (flxclr) e fluxo com nuvens (flxcld)
    flxclr = 0.0
    if (ta + TREF) > EPS:
        flxclr = ematm * SIGMA * ((ta + TREF) ** 4.0)

    flxcld = (lcld * (94.0 - 5.8 * zlcld) +
              mcld_eff * (94.0 - 5.8 * zmcld) +
              hcld_eff * (94.0 - 5.8 * zhcld))

    # 4. Fluxo descendente total
    flxir = flxclr + flxcld

    return flxir, mcld_eff, hcld_eff


# ---------------------------------------------------------------------------
# cloudbase  (module_radiation.F90, linhas 73-155)
# ---------------------------------------------------------------------------
def cloudbase(state: FasstState, doy: float, plat: float, lcldbse: float, mcldbse: float,
              hcldbse: float, lcld: float, mcld: float, hcld: float):
    """
    Calcula a altura das bases das nuvens (zlcld, zmcld, zhcld) com base em
    parametrização climatológica (latitude/estação).

    Preserva de propósito dois bugs do Fortran original (ver
    docstring do módulo, "Dois bugs do Fortran original"): o
    cálculo de `ilat` ignora o sinal de `plat`, e a condição de
    `isean` para o hemisfério sul é sempre verdadeira.
    """
    # Coeficientes da matriz (coef(2, 2, 3, 4) no Fortran)
    coef = np.array([
        [
            [[1.05, 0.6, 5.0, 25.0], [4.1, 0.3, 4.0, 25.0], [7.0, 1.5, 3.0, 30.0]],
            [[1.05, 0.6, 1.5, 25.0], [4.1, 2.0, 1.7, 25.0], [7.0, 1.5, 3.0, 30.0]]
        ],
        [
            [[1.15, 0.45, 5.0, 25.0], [4.4, 0.3, 4.0, 25.0], [7.0, 1.5, 3.0, 30.0]],
            [[1.15, 0.6, 1.5, 25.0], [4.4, 1.2, 3.0, 25.0], [7.0, 1.5, 3.0, 30.0]]
        ]
    ])

    # BUG PRESERVADO (doc "cloudbase — ilat"): o Fortran testa
    # `plat >= dabs(25d0)`, ou seja, compara `plat` com a constante
    # 25.0 — não com `dabs(plat)`. Latitudes negativas abaixo de
    # -25° nunca caem no `ilat` de latitude alta. NÃO trocar por
    # `abs(plat) >= 25.0` sem decisão explícita do projeto.
    ilat = 0
    if plat >= abs(25.0):
        ilat = 1

    isean = 1  # não inverno (índice 1 no Python, representa 2 no Fortran)

    doy_int = int(doy)
    if plat > 0.0 and (doy_int > 330 or doy_int < 65):
        isean = 0  # inverno, hemisfério norte

    # BUG PRESERVADO (doc "cloudbase — isean"): a condição abaixo
    # (`doy > 150 or doy < 250`) é sempre verdadeira para qualquer
    # dia do ano (150 < 250), então o hemisfério sul é classificado
    # como "inverno" o ano inteiro nesta parametrização.
    if plat < 0.0 and (doy_int > 150 or doy_int < 250):
        isean = 0  # inverno, hemisfério sul

    # Base da Nuvem Baixa
    if (_aint(abs(lcldbse - state.mflag) * 1e5) * 1e-5 <= EPS) and (lcld > EPS):
        a, b, c, d = coef[isean, ilat, 0, :]
        zlcld = a - b * (1.0 - abs(math.cos(c * (state.lat - d))))
    else:
        zlcld = lcldbse
    if zlcld <= EPS:
        zlcld = 0.0

    # Base da Nuvem Média
    if (_aint(abs(mcldbse - state.mflag) * 1e5) * 1e-5 <= EPS) and (mcld > EPS):
        a, b, c, d = coef[isean, ilat, 1, :]
        zmcld = a - b * (1.0 - abs(math.cos(c * (state.lat - d))))
    else:
        zmcld = mcldbse
    if zmcld <= EPS:
        zmcld = 0.0

    # Base da Nuvem Alta
    if (_aint(abs(hcldbse - state.mflag) * 1e5) * 1e-5 <= EPS) and (hcld > EPS):
        a, b, c, d = coef[isean, ilat, 2, :]
        zhcld = a - b * (1.0 - abs(math.cos(c * (state.lat - d))))
    else:
        zhcld = hcldbse
    if zhcld <= EPS:
        zhcld = 0.0

    return zlcld, zmcld, zhcld, lcld, mcld, hcld


# ---------------------------------------------------------------------------
# emisatm  (module_radiation.F90, linhas 158-220)
# ---------------------------------------------------------------------------
def emisatm(state: FasstState, tmp: float, rh: float):
    """
    Calcula a emissividade atmosférica de ondas longas (ematm).

    Formulação de Crawford & Duchon (1999): ematm = 1.24*(ea/Ta)^(1/7),
    com a pressão de vapor `ea` obtida da umidade relativa via
    Clausius-Clapeyron.
    """
    rv = 461.0
    eso = 6.13
    ta = tmp + TREF
    ea = 0.0
    ematm = 0.0

    l_val = (-2.43e-3 * ta + 3.166659) * 1e6
    if abs(l_val / rv * (1.0 / TREF - 1.0 / ta)) < 50.0:
        ea = (eso * math.exp(l_val / rv * (1.0 / TREF - 1.0 / ta))) * rh * 1e-2

    if ea > 0.0:
        ematm = 1.24 * ((ea / ta) ** (1.0 / 7.0))

    return ematm


# ---------------------------------------------------------------------------
# sol_zen  (module_radiation.F90, linhas 223-315)
# ---------------------------------------------------------------------------
def sol_zen(state: FasstState, year: float, doy: float, hour: float, minute: float):
    """
    Calcula os ângulos zênite (szen) e azimute solar (saz) em graus.

    Baseado no algoritmo da NOAA (Meeus, "Astronomical Algorithms").

    Decisão de tradução (ver docstring do módulo): o Fortran chama
    `acos(cosphi)`/`acos(costheta)` sem checar o domínio [-1, 1].
    Usamos `numpy.arccos` (sem clamping) em vez de `math.acos` para
    preservar o comportamento de "NaN silencioso" do Fortran, já
    que `math.acos` lançaria `ValueError` fora do domínio. A divisão
    de `costheta` usa `_fdiv` pelo mesmo motivo (denominador pode ser
    0.0 no polo ou com o sol no zênite/nadir).
    """
    f1 = 1.0 / 180.0
    f2 = 1.0 / PI
    zlat = PI * state.lat * f1
    zlong = -PI * state.mlong * f1  # radianos; calculada e nao usada adiante (igual ao Fortran)

    thour = 0.0 if abs(hour - 24.0) <= EPS else hour

    # Ano fracionário
    mody = year - int(year * 0.25) * 4.0
    if abs(mody) > EPS:
        fyear = (2.0 * PI / 365.0) * (doy - 1.0 + (thour - 12.0) / 24.0)
    else:
        fyear = (2.0 * PI / 366.0) * (doy - 1.0 + (thour - 12.0) / 24.0)

    # Equação do tempo (minutos)
    eqtime = 229.18 * (7.5e-5 + 1.868e-3 * math.cos(fyear) -
                       3.2077e-2 * math.sin(fyear) - 1.4615e-2 * math.cos(2.0 * fyear) -
                       4.0849e-2 * math.sin(2.0 * fyear))

    # Declinação solar (radianos)
    decl = (6.918e-3 - 3.99912e-1 * math.cos(fyear) + 7.0257e-2 * math.sin(fyear) -
            6.758e-3 * math.cos(2.0 * fyear) + 9.07e-4 * math.sin(2.0 * fyear) -
            2.697e-3 * math.cos(3.0 * fyear) + 1.48e-3 * math.sin(3.0 * fyear))

    # Tempo Solar Verdadeiro
    time_offset = eqtime + 4.0 * state.mlong - 60.0 * state.timeoffset
    tst = thour * 60.0 + minute + time_offset

    # Ângulo de Hora Solar (radianos)
    ha = (tst * 0.25 - 180.0) * PI * f1

    # Ângulos de Zênite e Azimute
    cosphi = math.sin(zlat) * math.sin(decl) + math.cos(zlat) * math.cos(decl) * math.cos(ha)
    # sem clamping — fiel ao Fortran; np.arccos devolve nan em vez de lançar exceção
    phi = np.arccos(cosphi)
    szen = phi * 180.0 * f2

    costheta = _fdiv(math.sin(zlat) * cosphi - math.sin(decl), math.cos(zlat) * math.sin(phi))
    # sem clamping aqui também — mesmo racional acima

    if zlat >= 0.0:
        if ha < 0.0:
            saz = 180.0 - np.arccos(costheta) * 180.0 * f2
        else:
            saz = 180.0 + np.arccos(costheta) * 180.0 * f2
    else:
        if ha < 0.0:
            saz = np.arccos(costheta) * 180.0 * f2
        else:
            saz = 360.0 - np.arccos(costheta) * 180.0 * f2

    return szen, saz


# ---------------------------------------------------------------------------
# Solflx  (module_radiation.F90, linhas 318-402)
# ---------------------------------------------------------------------------
def Solflx(state: FasstState, icld: list, szen: float, doy: float, prcp: float,
           tsol: float, cover: list, hgt: list):
    """
    Calcula os componentes de radiação solar direta (sdir) e difusa (sdif).

    `icld`/`cover`/`hgt`/`tsol` são `intent(inout)` no Fortran;
    aqui viram cópias locais + retorno em tupla (ver docstring do
    módulo).
    """
    icld_local = list(icld)
    cover_local = list(cover)
    hgt_local = list(hgt)

    ctest = 1 if (cover_local[0] > 0.5 or cover_local[1] > 0.5 or cover_local[2] > 0.5) else 0
    prcp1 = prcp * 1e-3
    frad = 0.017453292
    cosz = math.cos(szen * frad)

    jday = int(doy)

    # Chama o algoritmo de insolação (Shapiro)
    sdowns, direct, diffuse, icld_local, cover_local, hgt_local = insol(
        state, jday, icld_local, cosz, prcp1, hgt_local, cover_local
    )

    tsol_out = tsol
    if _aint(abs(tsol_out - state.mflag) * 1e5) * 1e-5 <= EPS:
        tsol_out = sdowns

    test = direct + diffuse
    f1 = tsol_out / sdowns if sdowns > EPS else 0.0

    if abs(test - tsol_out) > 1e-3:
        if ctest == 1:
            direct = direct * f1
            diffuse = max(0.0, tsol_out - direct)
            if diffuse < EPS:
                direct = tsol_out
        else:
            diffuse = diffuse * f1
            direct = max(0.0, tsol_out - diffuse)
            if direct < 0.0:
                diffuse = tsol_out

    sdir = direct
    sdif = diffuse

    return sdir, sdif, tsol_out, icld_local, cover_local, hgt_local


# ---------------------------------------------------------------------------
# insol  (module_radiation.F90, linhas 405-721)
# ---------------------------------------------------------------------------
# Constantes de dados (antigo `blockdata dinsol`) ----------------------------
# ATENÇÃO: o Fortran preenche essas matrizes em ordem column-major a
# partir de uma lista `DATA` plana. As linhas abaixo replicam a lista
# plana na mesma ordem de leitura do Fortran; o resultado é a
# TRANSPOSTA da matriz Fortran. O padrão de acesso em `insol`
# (`r1[coef_idx, layer_idx]`, e não `r1[layer_idx, coef_idx]`)
# compensa essa transposição — conferido numericamente contra o
# Fortran, não trocar a ordem dos índices sem revalidar.
def insol(state: FasstState, jday: int, icld: list, cosz: float, prcp1: float, hgt: list, cover: list):
    """
    Calcula os componentes de radiação direta e difusa método de R. Shapiro.
    """
    r1 = np.array([
        [0.12395,  0.15325,  0.15946,  0.27436],
        [-0.34765, -0.3962,  -0.42185, -0.43132],
        [0.39478,  0.42095,  0.488,    0.2692],
        [-0.14627, -0.142,   -0.18492, -0.00447]
    ])

    r2 = np.array([
        [0.25674,  0.42111,  0.61394,  0.69143],
        [-0.18077, -0.04002, -0.01469, -0.14419],
        [-0.21961, -0.51833, -0.174,   -0.051],
        [0.25272,  0.4054,   0.14215,  0.06682]
    ])

    t1 = np.array([
        [0.76977,  0.69318,  0.68679,  0.55336],
        [0.49407,  0.68227,  0.71012,  0.61511],
        [-0.44647, -0.64289, -0.71463, -0.29816],
        [0.11558,  0.1791,   0.22339,  -0.06663]
    ])

    t2 = np.array([
        [0.63547,  0.43562,  0.23865,  0.15785],
        [0.35229,  0.26094,  0.20143,  0.3241],
        [0.08709,  0.36428, -0.01183, -0.14458],
        [-0.22902, -0.38556, -0.07892,  0.01457]
    ])

    wt = np.array([
        [0.675,   1.552,   1.429,   1.512],
        [-3.432, -1.957,  -1.207,  -1.176],
        [1.929,  -1.762,  -2.008,  -2.16],
        [0.842,   2.067,   0.853,   1.42],
        [2.693,   0.448,   0.324,  -0.032],
        [-1.354,  0.932,   1.582,   1.422]
    ])

    icld_l = list(icld)
    cover_l = list(cover)
    hgt_l = list(hgt)

    sdowne = 2.0 * PI * float(jday - 2) / 365.242
    sdowne = (1.0001399 + 1.67261e-2 * math.cos(sdowne)) ** 2

    coszsq = cosz * cosz
    coszcube = coszsq * cosz
    sdown0 = 1369.2 * sdowne * cosz

    icla = [icld_l[2], icld_l[1], icld_l[0]]
    covera = [cover_l[2], cover_l[1], cover_l[0]]

    if icla[2] != 0:
        icla[2] = 4
    if icla[1] != 0:
        icla[1] = 3

    if icla[0] != 0:
        if cover_l[0] > 0.0 or icla[0] == 5:
            icla[0] = 1
        if cover_l[0] > 0.4 or icla[0] == 7:
            icla[0] = 2

    if prcp1 > EPS:
        covera = [1.0, 1.0, 1.0]
        icla = [2, 3, 4]

        if (abs(cover_l[2]) <= EPS) and (abs(cover_l[1]) <= EPS) and (abs(cover_l[0]) <= EPS):
            cover_l[2] = covera[0]
            hgt_l[2] = 8.0
            icld_l[2] = 5

            cover_l[1] = covera[1]
            hgt_l[1] = 4.0
            icld_l[1] = 3

            cover_l[0] = covera[2]
            hgt_l[0] = 1.5
            icld_l[0] = 6

    rks = [0.0, 0.0, 0.0]
    tk = [0.0, 0.0, 0.0]
    tdk = [0.0, 0.0, 0.0]

    for i in range(3):
        j = icla[i]
        if j > 4:
            j = 4
        fr = covera[i]
        l = i + 1

        wgt = 0.0
        rcld = 0.0
        tcld = 0.0
        rclr = 0.0
        tclr = 0.0

        j_idx = j - 1
        l_idx = l - 1

        if j != 0:
            wgt = (wt[0, j_idx] + wt[1, j_idx] * cosz + wt[2, j_idx] * fr +
                   wt[3, j_idx] * cosz * fr + wt[4, j_idx] * coszsq + wt[5, j_idx] * fr * fr)
            wgt = wgt * fr

            if wgt > 0.0:
                rcld = (r2[0, j_idx] + r2[1, j_idx] * cosz +
                        r2[2, j_idx] * coszsq + r2[3, j_idx] * coszcube)
                tcld = (t2[0, j_idx] + t2[1, j_idx] * cosz +
                        t2[2, j_idx] * coszsq + t2[3, j_idx] * coszcube)

        if wgt < 1.0:
            rclr = (r1[0, l_idx] + r1[1, l_idx] * cosz +
                    r1[2, l_idx] * coszsq + r1[3, l_idx] * coszcube)
            tclr = (t1[0, l_idx] + t1[1, l_idx] * cosz +
                    t1[2, l_idx] * coszsq + t1[3, l_idx] * coszcube)

        rks[i] = wgt * rcld + (1.0 - wgt) * rclr
        tk[i] = wgt * tcld + (1.0 - wgt) * tclr
        tdk[i] = max(0.0, tk[i] - rks[i])

    rg = state.albedo_fasst
    d1 = 1.0 - (rks[0] * rks[1])
    d2 = 1.0 - (rks[1] * rks[2])
    d3 = 1.0 - (rks[2] * rg)

    sdowns = d1 * d2 - (rks[0] * rks[2] * tk[1] * tk[1])
    sdowns = d3 * sdowns - (d1 * rks[1] * rg * tk[2] * tk[2])
    sdowns = sdowns - (rks[0] * rg * tk[1] * tk[1] * tk[2] * tk[2])
    # _fdiv: sdowns (denominador) pode ser exatamente 0.0 em condicoes de
    # noite/por-do-sol (sdown0=0 tambem) -- ver docstring do modulo.
    sdowns = _fdiv(tk[0] * tk[1] * tk[2] * sdown0, sdowns)

    if sdowns <= 0.0:
        sdowns = 0.0
        direct = 0.0
        diffuse = 0.0
    else:
        direct = max(0.0, tdk[0] * tdk[1] * tdk[2] * sdown0)
        diffuse = max(0.0, sdowns - direct)

    return sdowns, direct, diffuse, icld_l, cover_l, hgt_l

