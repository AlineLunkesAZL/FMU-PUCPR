"""
flow_param.py -- traducao de flow_param.F90

Fonte Fortran: flow_param.F90 (472 linhas). Uma unica sub-rotina publica
(`flow_param`). Calcula os parametros de transporte de agua liquida e de
vapor em cada no do perfil, e a partir deles os fluxos advectivos que
`soil_tmp`/`soil_moisture` consomem:

  * `klh` -- condutividade hidraulica por gradiente de pressao,
    van Genuchten (1980), com reducao por gelo segundo
    Niu & Yang (2006) JHM 7(5) 937-952 [array LOCAL]
  * `klt` -- condutividade hidraulica por gradiente de temperatura,
    Hansson et al. (2004) [array LOCAL]
  * `kvh` -- condutividade de vapor por gradiente de pressao,
    UNSAT-H / Fayer (2000) [array LOCAL]
  * `kvt` -- condutividade de vapor por gradiente de temperatura,
    UNSAT-H / Fayer (2000) [array LOCAL]
  * `dsmdh` -- d(umidade do solo)/d(cabeca de pressao) [campo de estado]

Assinatura original:
    subroutine flow_param(isn,qtop,qtopv,simeltsm,zt,delz,qbot,
                          klhtop,kvhtop,kvttop)
    integer(ip),intent(in):: isn
    real(dp),intent(in):: qtop,qtopv,simeltsm
    real(dp),intent(in):: zt(ntot),delz(ntot)
    real(dp),intent(out):: qbot,klhtop,kvhtop,kvttop

Arquitetura
-------------
`flow_param(state, isn, qtop, qtopv, simeltsm, zt, delz)` devolve
`(qbot, klhtop, kvhtop, kvttop)` -- os quatro `intent(out)`, nessa
ordem. E DROP-IN para o contrato `FlowParamFn` de
`fasst.new_profile` (`new_profile.py:92-101, 575-578`), sem adaptador,
ao contrario de `th_param` (ver PENDENCIAS.md §2.1). Os 4 call sites
sao new_profile.F90:937 (`isn = nnodes`), 1041, 1240 e 1562
(`isn = 1`) -- em new_profile.py, linhas 1479, 1583, 1763 e 1969.

DOIS MODOS, selecionados por `isn`
------------------------------------
`isn` nao e so um limite de laco, e um seletor de modo:

  * `isn == nnodes` (1o call site): o laco de condutividades roda para
    UM unico no (o de topo). O bloco de medias (linhas 259-470) e
    inteiramente PULADO (`if(isn == 1)`), e `qbot` sai SEMPRE 0 porque
    `klh(1)`/`kvh(1)` continuam com o zero da zeragem inicial. Ou seja,
    nesse modo o unico resultado util sao `klhtop`/`kvhtop`/`kvttop`.
  * `isn == 1` (os outros 3): perfil completo, laco de condutividades
    em 1..nnodes e depois as medias de Cherry & Freeze em 1..ntemp,
    escrevendo `khu`, `khl`, `vin`, `fv1`, `flowu`, `flowl`.

Note a assimetria de faixas: as condutividades sao calculadas em
`isn..nnodes`, mas as medias percorrem `1..ntemp`, e `ntemp` pode
exceder `nnodes` (nos de neve/vegetacao/ar acima do perfil de solo).
Para `nnodes < i <= ntemp`, `klh(i)`/`kvh(i)` valem 0, logo
`kdenom == 0` e os caminhos khll/khlv caem nos seus `else` -- e por
isso que os guards `dabs(kdenom) /= 0d0` existem (ver ponto 6).

Estado mutado (nenhum aparece na lista de argumentos formais)
--------------------------------------------------------------
`dsmdh`, `vin`, `fv1`, `flowu`, `flowl`, `khu`, `khl` -- todos zerados
em 1..ntot no inicio e preenchidos depois. E TAMBEM `stt`, via o
"safety catch" da linha 112 (ver ponto 1) -- o unico write em `stt`
aqui, facil de nao ver ao perguntar "o que flow_param muta?".

`klh`, `klt`, `kvh`, `kvt` e `impflag` sao arrays LOCAIS (declarados
`(ntot)` no Fortran, dimensionados `MAXN + EXTRAN` aqui por
consistencia com os demais modulos do projeto), nao campos de estado --
apesar de `klh`/`kvh` parecerem irmaos de `khl`/`khu`, que sao estado.

Comparacao de `node_type`
--------------------------
`node_type` e array de strings de largura fixa; o Fortran compara com
padding (`'AI' == 'AI '`). Aqui as comparacoes usam
`str(state.node_type[i]).strip()` uma vez por iteracao, inclusive para
o teste de BRANCO da linha 262 (`node_type(i) /= '  '`), que em NumPy
pode aparecer como `''` ou como espacos -- ver o docstring de
`fasst.state`, ponto 10. (`soil_tmp.py` compara direto, sem `.strip()`,
porque so testa rotulos nao-brancos.)
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

from .constants import EPS, EXTRAN, GRAV, MAXN, RV, TREF
from .functions import Phase, dense, head, soilhumid, spheats
from .state import FasstState


_NODESZ = MAXN + EXTRAN

# Quarta/quinta copia manual do ajuste WES de calor latente de vaporizacao
# (linha 264). Ja sinalizada em sflux.py (_LHVAP_LITERAL), surfenergy.py e
# soil_tmp.py (_LHEVAP_INTERCEPT/_LHEVAP_SLOPE) e module_canopy.py:1553.
# Mantida local aqui tambem, sem importar de nenhum irmao -- ver ponto 8.
_LHEVAP_INTERCEPT = 2500775.6  # J/kg
_LHEVAP_SLOPE = 2369.729       # J/(kg*K)

# Linha 203: `610.78d0` escrito a mao; duplica constants.VPSAT0 (mesmo valor,
# duas fontes). Mantido literal, mesma convencao de sp_humid.py:331 e
# functions.py:491.
_VPSAT0_LITERAL = 610.78  # Pa


def _reals(n: int) -> np.ndarray:
    return np.zeros(n + 1, dtype=np.float64)


def _ints(n: int) -> np.ndarray:
    return np.zeros(n + 1, dtype=np.int64)


def _anint(x: float, p: int) -> float:
    """Fortran ``ANINT(x * 10**p) * 10**-p`` -- arredonda para longe de zero.

    NAO usa o idioma ``floor(|xs| + 0.5)`` dos demais modulos: para
    2**52 <= |xs| < 2**53 o espacamento entre doubles e 1, e ``|xs| + 0.5``
    arredonda (half-even) para o PROXIMO inteiro par -- um quantum acima do
    ANINT do Fortran (C ``round()``), que devolve o inteiro intacto. Medido
    contra gfortran na validacao (tests/test_flow_param_ref.py): o idioma
    antigo diverge em 39% dos valores dessa faixa; esta versao, em 0%. A parte
    fracionaria ``a - floor(a)`` e exata em qualquer faixa.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    a = abs(xs)
    r = math.floor(a)
    if a - r >= 0.5:
        r += 1.0
    return math.copysign(r, xs) * (10.0**-p)


def _round20(x: float) -> float:
    """O unico idioma de arredondamento deste arquivo: ``anint(x*1d20)*1d-20``."""
    return _anint(x, 20)


def flow_param(
    state: FasstState,
    isn: int,
    qtop: float,
    qtopv: float,
    simeltsm: float,
    zt: np.ndarray,
    delz: np.ndarray,
) -> Tuple[float, float, float, float]:
    """Parametros de fluxo de agua/vapor do perfil.

    Devolve `(qbot, klhtop, kvhtop, kvttop)`. Muta `state.dsmdh`,
    `state.vin`, `state.fv1`, `state.flowu`, `state.flowl`, `state.khu`,
    `state.khl` e -- inesperadamente -- `state.stt` (ponto 1).
    """
    ntot = state.ntot
    ntemp = state.ntemp
    nnodes = state.nnodes
    ntype = state.ntype
    nsoilp = state.nsoilp
    stt = state.stt

    # ---------------- initialize variables (linhas 41-89) ----------------
    # Transcrito literalmente (mesma disciplina de icethick.py/soil_moisture.py).
    # Linha 41-42: `d1i = 0` e imediatamente `d1i = 1` -- a primeira atribuicao
    # e morta (ponto 9).
    d1i = Phase.WATER          # 1
    d2i = Phase.WATER_VAPOR    # 0
    w = 0.0
    hc1 = 0.0
    c1 = 0.0
    c2 = 0.0
    h0 = 0.0                   # declarado, zerado, NUNCA lido (ponto 10)
    fvl = 0.0
    fvu = 0.0
    fll = 0.0
    flu = 0.0
    cw = 0.0
    cv = 0.0
    kv = 0.0
    kl = 0.0
    dcoef = 0.0                # Fortran `D` (renomeado: D nao e nome util em Python)
    eta = 0.0
    hr = 0.0                   # Fortran `Hr`
    a = 0.0
    b = 0.0
    t1 = 0.0
    vpsat = 0.0
    desdt = 0.0
    rhovs = 0.0
    drvsdt = 0.0
    kdenom = 0.0
    vinlw = 0.0
    vinuw = 0.0
    vinlv = 0.0
    vinuv = 0.0
    dh = 0.0
    dtg = 0.0                  # Fortran `dT` (gradiente de temperatura, K/m)
    lhe = 0.0
    qbot = 0.0
    fvl1 = 0.0
    fvu1 = 0.0
    khll = 0.0
    khlv = 0.0
    khul = 0.0
    khuv = 0.0
    phtemp = 0.0               # declarado, zerado, NUNCA lido (ponto 10)
    ptemp1 = 0.0
    slpr = 0.0
    sph1 = 0.0
    sph2 = 0.0
    rd1 = 0.0
    tsign = 0.0                # declarado, zerado, NUNCA lido (ponto 10)
    frz = 0.0                  # Fortran `Frz`

    impflag = _ints(_NODESZ)
    klh = _reals(_NODESZ)
    klt = _reals(_NODESZ)
    kvh = _reals(_NODESZ)
    kvt = _reals(_NODESZ)

    for i in range(1, ntot + 1):
        impflag[i] = 0
        state.dsmdh[i] = 0.0
        klh[i] = 0.0
        klt[i] = 0.0
        kvh[i] = 0.0
        kvt[i] = 0.0
        state.vin[i] = 0.0                                               # m/s
        state.fv1[i] = 0.0                                               # W/m^2
        state.flowu[i] = 0.0                                             # m/s
        state.flowl[i] = 0.0                                             # m/s
        state.khu[i] = 0.0                                               # m/s
        state.khl[i] = 0.0                                               # m/s

    # Linha 106: o fator de projecao de rampa esta FIXO em 1 -- a formula real
    # (`dmax1(0d0,dmin1(1d0,dcos(slope_fasst*pi/1.8d2)))`) esta comentada no
    # original. Logo as 4 multiplicacoes por slpr das linhas 234-237 sao
    # no-ops. Preservado, ver ponto 5.
    slpr = 1.0

    # ---- solve for dsmdh, kvh, kvt, klh, klt (linhas 108-238) ----
    # Linha 109: `d2i = 2` (ICE) -- ATRIBUICAO MORTA. Nenhuma leitura de d2i
    # ocorre antes da linha 260, que o devolve a 0 (WATER_VAPOR). Ver ponto 9:
    # quem ler so esta linha conclui, errado, que `cv`/`sph2` usam GELO.
    d2i = Phase.ICE

    for i in range(isn, nnodes + 1):
        # SAFETY CATCH QUE MUTA ESTADO -- ver ponto 1. `dabs(stt) <= eps`
        # significa ~0 K (impossivel fisicamente), ou seja no nao inicializado;
        # o 1.0 K escrito aqui persiste para todos os consumidores seguintes.
        if abs(stt[i]) <= EPS:
            stt[i] = 1.0                                                 # safety catch to prevent /0

        if ((ntype[i] == 26 or ntype[i] == 27) or ntype[i] == 0) or (
            nsoilp[i][11] <= EPS or nsoilp[i][12] <= EPS
        ):
            impflag[i] = 1

        if impflag[i] == 0:
            rhow = dense(stt[i], 0.0, d1i)                               # kg/m^3 (water density)

            # Sem guard contra nsoilp[i][9] == nsoilp[i][8] -- ver ponto 11.
            w = (state.soil_moist[i] - nsoilp[i][8]) / (
                nsoilp[i][9] - nsoilp[i][8]
            )                                                            # unitless

            if w > EPS:
                if abs(state.phead[i]) <= EPS:
                    # trailing `!1d0 !0` no original: alternativas vestigiais
                    ptemp1 = abs(head(state, i, 0.99 * nsoilp[i][9])) * 1e2
                else:
                    ptemp1 = abs(state.phead[i]) * 1e2                   # cm

                c1 = (1.0 + (nsoilp[i][10] * ptemp1) ** nsoilp[i][11]) ** (
                    -nsoilp[i][12] - 1.0
                )                                                        # unitless
                c2 = nsoilp[i][12] * nsoilp[i][11] * (
                    nsoilp[i][10] ** nsoilp[i][11]
                )                                                        # unitless

                state.dsmdh[i] = (
                    (nsoilp[i][9] - nsoilp[i][8])
                    * c2
                    * c1
                    * (ptemp1 ** (nsoilp[i][11] - 1.0))
                    * 1e2
                )                                                        # 1/m

            if w <= 0.99 and w > 0.0:
                hc1 = 1.0 - (1.0 - (w ** (1.0 / nsoilp[i][12]))) ** nsoilp[i][12]  # unitless
                klh[i] = nsoilp[i][7] * math.sqrt(w) * hc1 * hc1 * 1e-2   # m/s (hydraulic conductivity due to waterflow)
                klh[i] = min(klh[i], nsoilp[i][7] * 1e-2)
                if state.ice[i] > EPS:
                    frz = math.exp(-3e-1 * (1.0 - state.ice[i] / nsoilp[i][9])) - math.exp(-3e-1)
                    klh[i] = klh[i] * (1.0 - frz)
                    # alternativa comentada no original:
                    #   klh(i) = klh(i)*(1d0 + 8d0*ice(i))**2d0   !m/s
                    # (note que ela AUMENTARIA klh com gelo -- ver ponto 4)
                elif state.ice[i] <= EPS and stt[i] > TREF:
                    klh[i] = klh[i] * (5.3888e-1 + 2.096e-2 * (stt[i] - TREF))
            elif w > 0.99:
                klh[i] = nsoilp[i][7] * 1e-2
                if state.ice[i] > EPS:
                    frz = math.exp(-3e-1 * (1.0 - state.ice[i] / nsoilp[i][9])) - math.exp(-3e-1)
                    klh[i] = klh[i] * (1.0 - frz)
                    #   klh(i) = klh(i)*(1d0 + 8d0*ice(i))**2d0   !m/s  (comentado)
                elif state.ice[i] <= EPS and stt[i] > TREF:
                    klh[i] = klh[i] * (5.3888e-1 + 2.096e-2 * (stt[i] - TREF))
            # w <= 0 ou NaN: klh[i] fica com o 0 da zeragem

            klh[i] = max(0.0, klh[i])

            # Bloco comentado no original (arvore de 5 ramos por ntype):
            #   if(ntype(i) <= 4)      G = 5d0   !2d0
            #   else if(5..8)          G = 7d0
            #   else if(ntype(i)==15)  G = 10d0
            #   else                   G = 9d0
            # O codigo ATIVO abaixo resulta em G == 10 em todos os casos
            # alcancaveis -- ver ponto 2.
            if ntype[i] != 15:
                g = 2.5e1 * (nsoilp[i][20] + 1.0)
                if g < 5.0:
                    g = 5.0
                if g > 1e1:
                    g = 1e1
            else:
                g = 1e1

            t1 = stt[i] - TREF
            if t1 <= EPS:
                klt[i] = 0.0
            else:
                gamma = 1e-3 * (7.56e1 - 1.425e-1 * t1 - 2.38e-4 * t1 * t1)   # surface tension (N/m)
                dgdt = 1e-3 * (-1.425e-1 - 2.0 * 2.38e-4 * t1)                # d(gamma)/dT (N/m*K)
                if abs(gamma) > EPS:
                    # gamma ~ 0: klt[i] fica com o 0 da zeragem (nao e reatribuido)
                    klt[i] = max(0.0, klh[i] * (g * state.phead[i] * dgdt / gamma))  # m^2/s*K

            if state.wvc[i] >= 1e-15:
                hr = soilhumid(state, i, state.phead[i], state.soil_moist[i], stt[i])

                if abs(state.ice[i]) <= EPS:                             # over water
                    a = 17.269
                    b = 35.86
                else:                                                    # over ice/snow
                    a = 21.8745
                    b = 7.66

                t1 = stt[i] - b                                          # K
                if abs(a * (stt[i] - TREF) / t1) > 5e1:
                    vpsat = 0.0
                else:
                    vpsat = _VPSAT0_LITERAL * math.exp(a * (stt[i] - TREF) / t1)  # Pa (saturation vapor pressure)

                rhovs = vpsat / (RV * stt[i])                             # kg/m^3 (saturated water vapor density)
                desdt = vpsat * a * (1.0 - (stt[i] - TREF) / t1) / t1     # Pa/K (saturated)
                drvsdt = (desdt - vpsat / stt[i]) / (RV * stt[i])         # kg/m^3*K (drhovs/dtemp1)

                t1 = max(0.0, nsoilp[i][2] - (state.soil_moist[i] + state.ice[i]))
                dcoef = 0.0
                if t1 > 0.0:
                    dcoef = (t1 ** (5.0 / 3.0)) * (2.12e-5 * (stt[i] / TREF) ** 2.0)  # m^2/s (diff. coeff. water vapor in air)
                    dcoef = max(0.0, dcoef)

                kvh[i] = dcoef * rhovs * GRAV * hr / (rhow * RV * stt[i])  # m/s
                kvh[i] = max(0.0, kvh[i])

                if nsoilp[i][20] > EPS:
                    # nsoilp[:,20] e FRACAO de argila (initprofile.py:267 faz
                    # *1e-2); o *1d2 aqui volta para porcento. Ver ponto 3.
                    t1 = (
                        (1.0 + 2.6 / math.sqrt(nsoilp[i][20] * 1e2)) * state.soil_moist[i]
                    ) ** 4.0
                    if t1 > 5e1:
                        t1 = 5e1
                    eta = 9.5 + 6.0 * state.soil_moist[i] - 8.5 * math.exp(-t1)  # unitless
                else:
                    eta = 0.0

                kvt[i] = max(0.0, dcoef * eta * hr * drvsdt / rhow)      # m^2/s*K
            # fim de wvc(i) >= 1d-15
        # fim de impflag(i) == 0

        state.dsmdh[i] = _round20(state.dsmdh[i])                        # 1/m
        klh[i] = _round20(klh[i] * slpr)                                 # m/s   (slpr == 1, ponto 5)
        kvh[i] = _round20(kvh[i] * slpr)                                 # m/s
        klt[i] = _round20(klt[i] * slpr)                                 # m^2/s*K
        kvt[i] = _round20(kvt[i] * slpr)                                 # m^2/s*K

    # ---- condicao de contorno inferior e saidas de topo (linhas 240-256) ----
    if abs(zt[1] - state.gwl) > EPS:
        dh = (state.phead[1] - 0.0) / (zt[1] - state.gwl)                # head at gwl = 0
        dtg = 0.0  # formula real comentada: dexp(-2.08d0*(zt(1)-gwl))/(zt(1)-gwl)
    else:
        dh = 0.0
        dtg = 0.0

    klhtop = klh[nnodes]
    kvhtop = kvh[nnodes]
    kvttop = kvt[nnodes]

    klbot = -(klh[1] * (dh + 1.0) + klt[1] * dtg)                        # m/s
    kvbot = -(kvh[1] * dh + kvt[1] * dtg)                                # m/s
    # ORDEM IMPORTA (ponto 7): qbot soma os valores NAO quantizados; so depois
    # klbot/kvbot sao quantizados, e e a versao QUANTIZADA que o no de base le
    # mais abaixo (vinlw = klbot, vinlv = kvbot). Nao unificar num passo.
    qbot = klbot + kvbot

    klbot = _round20(klbot)
    kvbot = _round20(kvbot)

    # ---- medias de klh e kvh (Cherry & Freeze, fluxo perpendicular) ----
    if isn == 1:
        d2i = Phase.WATER_VAPOR   # linha 260 -- e ESTE o valor que sph2/cv usam
        for i in range(1, ntemp + 1):  # `do i=1,ntemp  !ntot` no original
            nt = str(state.node_type[i]).strip()
            if (nt != "AI" and nt != "") and impflag[i] == 0:             # not air and not impervious
                lhe = _LHEVAP_INTERCEPT - _LHEVAP_SLOPE * (stt[i] - TREF)  # J/kg = (m/s)^2

                # ---- top node, snow, veg ----
                if i >= nnodes:
                    if nt == "HM":                                        # snow on ground
                        fvl1 = 0.0
                        fvu1 = 0.0
                        fvl = 0.0
                        fvu = 0.0
                        vinlw = 0.0
                        vinuw = 0.0
                        vinlv = 0.0
                        vinuv = 0.0

                        fll = -state.grspheat[i] * simeltsm / state.deltat_fasst          # W/m^2*K
                        flu = -state.grspheat[i] * state.atopf / state.deltat_fasst / float(state.step)  # W/m^2*K
                    elif nt == "VG":
                        fvl1 = 0.0
                        fvu1 = 0.0
                        fvl = 0.0
                        fvu = 0.0
                        vinlw = 0.0
                        vinuw = 0.0
                        vinlv = 0.0
                        vinuv = 0.0
                        fll = 0.0
                        flu = 0.0
                    else:
                        rd1 = dense(stt[i], 0.0, d1i)
                        sph1 = spheats(stt[i], d1i) * rd1
                        sph2 = spheats(stt[i], d2i) * dense(stt[i], 0.0, d2i)

                        knum = delz[i - 1] + delz[i]                              # m
                        kdenom = klh[i - 1] * delz[i] + klh[i] * delz[i - 1]      # m^2/s
                        if abs(kdenom) != 0.0:                                    # ponto 6: teste de zero exato
                            khll = klh[i] * klh[i - 1] * knum / kdenom            # m/s
                        else:
                            khll = 0.0

                        kdenom = kvh[i - 1] * delz[i] + kvh[i] * delz[i - 1]      # m^s/s [sic]
                        if abs(kdenom) != 0.0:
                            khlv = kvh[i] * kvh[i - 1] * knum / kdenom            # m/s
                        else:
                            khlv = 0.0
                        state.khl[i] = (khll + khlv) / (zt[i] - zt[i - 1])         # 1/s

                        kv = 5e-1 * (kvt[i - 1] + kvt[i])                          # m^2/s*K
                        kl = 5e-1 * (klt[i - 1] + klt[i])                          # m^2/s*K
                        cw = 5e-1 * (
                            sph1 + spheats(stt[i - 1], d1i) * dense(stt[i - 1], 0.0, d1i)
                        )                                                          # J/m^3*K
                        cv = 5e-1 * (
                            sph2 + spheats(stt[i - 1], d2i) * dense(stt[i - 1], 0.0, d2i)
                        )                                                          # J/m^3*K

                        dh = (state.phead[i] - state.phead[i - 1]) / (zt[i] - zt[i - 1])  # m/m
                        dtg = (stt[i] - stt[i - 1]) / (zt[i] - zt[i - 1])                 # K/m

                        vinlw = -(khll * (dh + 1.0) + kl * dtg)                    # m/s
                        fll = cw * vinlw                                           # W/m^2*K
                        vinlv = -(khlv * dh + kv * dtg)                            # m/s
                        fvl = cv * vinlv                                           # W/m^2*K
                        fvl1 = lhe * rd1 * vinlv                                   # W/m^2

                        khul = 0.0                                                 # m/s
                        khuv = 0.0                                                 # m/s
                        state.khu[i] = 0.0                                         # 1/s

                        vinuw = -qtop - qtopv                                      # m/s
                        vinuv = 0.0  # `0d0 !-qtopv` no original

                        if state.hm > EPS:
                            vinuv = 0.0   # ja e 0 na linha acima -- no-op preservado
                        flu = sph1 * vinuw                                         # W/m^2*K
                        fvu = sph2 * vinuv                                         # W/m^2*K
                        fvu1 = lhe * rd1 * vinuv                                   # W/m^2

                # ---- bottom node ----
                elif i == 1:
                    rd1 = dense(stt[i], 0.0, d1i)
                    sph1 = spheats(stt[i], d1i) * rd1
                    sph2 = spheats(stt[i], d2i) * dense(stt[i], 0.0, d2i)

                    khll = 0.0                                                     # m/s
                    khlv = 0.0                                                     # m/s
                    state.khl[i] = 0.0                                             # 1/s

                    vinlw = klbot                                                  # m/s (JA quantizado -- ponto 7)
                    fll = sph1 * vinlw                                             # W/m^2*K
                    vinlv = kvbot                                                  # m/s
                    fvl = sph2 * vinlv                                             # W/m^2*K
                    fvl1 = lhe * rd1 * vinlv                                       # W/m^2

                    knum = delz[i] + delz[i + 1]                                   # m
                    kdenom = klh[i] * delz[i + 1] + klh[i + 1] * delz[i]           # m^2/s
                    if abs(kdenom) != 0.0:
                        khul = klh[i] * klh[i + 1] * knum / kdenom                 # m/s
                    else:
                        khul = 0.0
                    kdenom = kvh[i] * delz[i + 1] + kvh[i + 1] * delz[i]           # m^2/s
                    if abs(kdenom) != 0.0:
                        khuv = kvh[i] * kvh[i + 1] * knum / kdenom                 # m/s
                    else:
                        khuv = 0.0
                    state.khu[i] = (khul + khuv) / (zt[i + 1] - zt[i])             # 1/s

                    kv = 5e-1 * (kvt[i] + kvt[i + 1])                              # m^2/s*K
                    kl = 5e-1 * (klt[i] + klt[i + 1])                              # m^2/s*K
                    # literais 1 e 0 no original em vez de d1i/d2i -- mesmos
                    # valores (d1i==1, d2i==0 aqui), ver ponto 9
                    cw = 5e-1 * (
                        sph1 + spheats(stt[i + 1], Phase.WATER) * dense(stt[i + 1], 0.0, d1i)
                    )                                                              # J/m^3*K
                    cv = 5e-1 * (
                        sph2 + spheats(stt[i + 1], Phase.WATER_VAPOR) * dense(stt[i + 1], 0.0, d2i)
                    )                                                              # J/m^3*K

                    dh = (state.phead[i + 1] - state.phead[i]) / (zt[i + 1] - zt[i])  # m/m
                    dtg = (stt[i + 1] - stt[i]) / (zt[i + 1] - zt[i])                 # K/m

                    vinuw = -(khul * (dh + 1.0) + kl * dtg)                        # m/s
                    flu = cw * vinuw                                               # W/m^2*K
                    vinuv = -(khuv * dh + kv * dtg)                                # m/s
                    fvu = cv * vinuv                                               # W/m^2*K
                    fvu1 = lhe * rd1 * vinuv                                       # W/m^2

                # ---- interior nodes ----
                elif i > 1 and i < nnodes:
                    rd1 = dense(stt[i], 0.0, d1i)
                    sph1 = spheats(stt[i], d1i) * rd1
                    sph2 = spheats(stt[i], d2i) * dense(stt[i], 0.0, d2i)

                    knum = delz[i - 1] + delz[i]                                   # m
                    kdenom = klh[i - 1] * delz[i] + klh[i] * delz[i - 1]           # m^2/s
                    if abs(kdenom) != 0.0:
                        khll = klh[i] * klh[i - 1] * knum / kdenom                 # m/s
                    else:
                        khll = 0.0
                    kdenom = kvh[i - 1] * delz[i] + kvh[i] * delz[i - 1]           # m^s/s [sic]
                    if abs(kdenom) != 0.0:
                        khlv = kvh[i] * kvh[i - 1] * knum / kdenom                 # m/s
                    else:
                        khlv = 0.0
                    state.khl[i] = (khll + khlv) / (zt[i] - zt[i - 1])             # 1/s

                    kv = 5e-1 * (kvt[i - 1] + kvt[i])                              # m^2/s*K
                    kl = 5e-1 * (klt[i - 1] + klt[i])                              # m^2/s*K
                    cw = 5e-1 * (
                        sph1 + spheats(stt[i - 1], d1i) * dense(stt[i - 1], 0.0, d1i)
                    )                                                              # J/m^3*K
                    cv = 5e-1 * (
                        sph2 + spheats(stt[i - 1], d2i) * dense(stt[i - 1], 0.0, d2i)
                    )                                                              # J/m^3*K

                    dh = (state.phead[i] - state.phead[i - 1]) / (zt[i] - zt[i - 1])  # m/m
                    dtg = (stt[i] - stt[i - 1]) / (zt[i] - zt[i - 1])                 # K/m

                    vinlw = -(khll * (dh + 1.0) + kl * dtg)                        # m/s
                    fll = cw * vinlw                                               # W/m^2*K
                    vinlv = -(khlv * dh + kv * dtg)                                # m/s
                    fvl = cv * vinlv                                               # W/m^2*K
                    fvl1 = lhe * rd1 * vinlv                                       # W/m^2

                    knum = delz[i] + delz[i + 1]                                   # m
                    kdenom = klh[i] * delz[i + 1] + klh[i + 1] * delz[i]           # m^2/s
                    if abs(kdenom) != 0.0:
                        khul = klh[i] * klh[i + 1] * knum / kdenom                 # m/s
                    else:
                        khul = 0.0
                    kdenom = kvh[i] * delz[i + 1] + kvh[i + 1] * delz[i]           # m^2/s
                    if abs(kdenom) != 0.0:
                        khuv = kvh[i] * kvh[i + 1] * knum / kdenom                 # m/s
                    else:
                        khuv = 0.0
                    state.khu[i] = (khul + khuv) / (zt[i + 1] - zt[i])             # 1/s

                    kv = 5e-1 * (kvt[i] + kvt[i + 1])                              # m^2/s*K
                    kl = 5e-1 * (klt[i] + klt[i + 1])                              # m^2/s*K
                    cw = 5e-1 * (
                        sph1 + spheats(stt[i + 1], Phase.WATER) * dense(stt[i + 1], 0.0, d1i)
                    )                                                              # J/m^3*K
                    cv = 5e-1 * (
                        sph2 + spheats(stt[i + 1], Phase.WATER_VAPOR) * dense(stt[i + 1], 0.0, d2i)
                    )                                                              # J/m^3*K
                    dh = (state.phead[i + 1] - state.phead[i]) / (zt[i + 1] - zt[i])  # m/m
                    dtg = (stt[i + 1] - stt[i]) / (zt[i + 1] - zt[i])                 # K/m

                    vinuw = -(khul * (dh + 1.0) + kl * dtg)                        # m/s
                    flu = cw * vinuw                                               # W/m^2*K
                    vinuv = -(khuv * dh + kv * dtg)                                # m/s
                    fvu = cv * vinuv                                               # W/m^2*K
                    fvu1 = lhe * rd1 * vinuv                                       # W/m^2

                state.vin[i] = (vinuw + vinuv) - (vinlw + vinlv)                   # m/s

                state.fv1[i] = fvu1 - fvl1                                         # W/m^2
                state.flowu[i] = flu + fvu
                state.flowl[i] = fll + fvl                                         # W/m^2*K

                state.khu[i] = _round20(state.khu[i])
                state.khl[i] = _round20(state.khl[i])
                state.vin[i] = _round20(state.vin[i])
                state.fv1[i] = _round20(state.fv1[i])
                state.flowu[i] = _round20(state.flowu[i])
                state.flowl[i] = _round20(state.flowl[i])

    return qbot, klhtop, kvhtop, kvttop
