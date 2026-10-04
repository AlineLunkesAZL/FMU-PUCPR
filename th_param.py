"""
th_param.py -- traducao de th_param.F90

Fonte Fortran: th_param.F90 (379 linhas). Uma unica sub-rotina publica
(`th_param`), sem nenhuma sub-rotina chamada ("no subroutines called" no
cabecalho do .F90). Calcula as propriedades termicas de cada no do
perfil -- condutividade termica (`grthcond`, W/m*K) e capacidade
calorifica volumetrica (`grspheat`, J/m^3*K) -- para solo, materiais de
superficie (concreto, asfalto, rocha, agua, ar, gelo glacial) e camadas
de neve/gelo/vegetacao acima do solo.

Assinatura original:
    subroutine th_param(n,pdens,pdensnew,rhov,rhoda,rhotot)
    integer(ip),intent(in):: n
    real(dp),intent(in):: pdens,pdensnew
    real(dp),intent(in):: rhov(n),rhoda(n)
    real(dp),intent(out):: rhotot

Arquitetura
-------------
`th_param` recebe `state: FasstState` como argumento explicito.
`EPS`, `SDENSW`, `SPFLAG`, `TREF`, `KVEG` e `SPHVEG` vem de
`fasst.constants`; `ntype`, `nsoilp`, `node_type`, `iw`, `nnodes`,
`ntemp`, `ice`, `soil_moist`, `hi`, `hsaccum`, `newsd`, `sigfl`, `stt`,
`wvc`, `sdens`, `grspheat`, `grthcond`, `sphm` e `km` sao campos de
`FasstState`. `dense` (pura), `spheats` e `thconds` sao chamadas reais,
importadas de `fasst.functions`, com `Phase` no lugar dos inteiros
locais `d0i`/`d1i`/`d2i`/`d3i` do Fortran (`WATER_VAPOR`/`WATER`/`ICE`/
`DRY_AIR`).

As saidas reais desta sub-rotina (`grthcond`, `grspheat`, `km`, `sphm`)
sao campos de `FasstState`, atualizados como efeito colateral; a unica
saida `intent(out)` (`rhotot`, densidade media da camada de neve/gelo)
e devolvida diretamente. `n` (dimensao de `rhov`/`rhoda` no Fortran) e
mantido na assinatura por fidelidade, mas nao e lido no corpo.
`rhov`/`rhoda` sao listas 1-based (indice 0 sem uso), uma entrada por no.

Referencias citadas no codigo Fortran: Farouki, CRREL Monograph 81-1,
p.112-116 (condutividade termica de solos); Lu et al. (2007) com
coeficientes de Tarnawski et al. (2009) (numero de Kersten); Sturm et
al. (1997), J. Glaciology 43(143), p.26-41 (condutividade termica da
neve).

Pontos de atencao
-------------------
* Provavel bug do Fortran original, preservado sem correcao: nos tres
  blocos que calculam a temperatura media da interface neve/gelo
  (`tice`, linhas 302-304, 314-316 e 329-330 do .F90), a fracao esta
  invertida. A expressao e
      tice = (a + b) / (b*T_base + a*T_topo)
  com a = condutancia da neve (k/espessura) e b = condutancia do gelo.
  Uma media ponderada de temperatura por condutancia seria
  (b*T_base + a*T_topo) / (a + b) -- o inverso. Como escrito, `tice`
  vale aproximadamente 1/T (~0.004), nao T (~270 K). `tice` e depois
  passado a `thconds`/`spheats`/`dense` para a fase gelo; as clampagens
  internas dessas funcoes (fasst.functions) limitam o resultado a
  valores finitos, mas diferentes dos que uma temperatura de ~270 K
  produziria. Exemplo numerico (neve de 0.2 m com k=0.1 W/m*K sobre gelo
  de 0.1 m, T_base=273 K, T_topo=270 K): `tice` = 0.0037 (media
  ponderada correta: 272.9 K); condutividade do gelo = 3.41 W/m*K
  (contra 2.22, +54%); calor especifico do gelo = 1389 J/kg*K (contra
  2050, -32%); densidade do gelo praticamente inalterada (916.20 contra
  916.21 kg/m^3). Traduzido literalmente.
* Provavel bug do Fortran original, preservado sem correcao: no
  calculo da condutividade dos solidos do solo, `ksq = 8.4d0*pq`
  (multiplicacao) enquanto os outros dois fatores do mesmo produto
  usam exponenciacao (`kso = 0.25d0**nsoilp(i,14)`,
  `ksn = 2.9d0**(...)`), como em um modelo de media geometrica. A
  versao anterior, deixada comentada na linha 112, usava `8.4d0**`.
* `ka` e `kv` (condutividade do ar e do vapor d'agua) sao calculadas
  por `thconds` para todo no de solo mas nunca lidas depois; `cs` e
  zerada no topo e nunca usada. Codigo morto, preservado por
  fidelidade.
* Tipos de material sem ramo proprio (`ntype` 19, 22, 23, 24, 28, 29)
  nao recebem nenhum valor de `grthcond` nesta chamada -- mantem o
  valor anterior de `state.grthcond` (a capacidade calorifica
  `grspheat`, ao contrario, e sempre calculada).
* `tdens` (densidade da neve acumulada) e atribuida so quando
  `hsaccum > eps` e persiste entre iteracoes do laco `do i`; nos
  demais casos permanece com o valor da iteracao anterior (0.0 se
  nenhuma iteracao anterior a atribuiu). Quando `iw >= 2`, le
  `state.sdens[iw-1]`, que precisa estar alocado (`FasstState.sdens`
  e `None` por padrao).
* `rhotot` so e atribuida em nos de neve/gelo (`node_type` `HM`/`MX`);
  cada no sobrescreve o anterior, e o valor devolvido e o do ultimo.
  Se nao houver nenhum, devolve 0.0.
* No ramo de neve sem gelo, o caso `else` (que cobre "neve antiga e nova
  juntas", mas tambem "nem uma nem outra") divide por
  `newsd + hsaccum` e por `hsaccum/kms + newsd/kmsn`, que valem 0/0
  quando `hsaccum` e `newsd` sao zero. O Fortran propaga NaN
  silenciosamente; aqui `_fdiv` (mesma funcao de `sp_humid.py`/
  `snow.py`/`module_lowveg.py`/`module_canopy.py`/`module_radiation.py`/
  `soil_strength.py`) reproduz esse comportamento em vez de levantar
  `ZeroDivisionError`.
"""

import math

import numpy as np

from .constants import EPS, KVEG, SDENSW, SPFLAG, SPHVEG, TREF
from .functions import Phase, dense, spheats, thconds
from .state import FasstState


def _anint(x, p=0):
    """Fortran ``ANINT(x*10**p)*10**-p`` -- arredonda para longe de zero."""
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def _fdiv(x, y):
    """Divisao IEEE-754 pura: ``x / 0.0`` -> +-inf / nan, sem levantar."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


def _sturm_conductivity(densgc):
    """Condutividade termica da neve (W/m*K), Sturm et al. (1997), em
    funcao da densidade em g/cm^3, limitada a [limite inferior, 1.0]."""
    if densgc < 1.56e-1:
        k = 2.3e-2 + 2.34e-1 * densgc
        return min(max(2.3e-2, k), 1.0)
    k = 1.38e-1 - 1.01 * densgc + 3.233 * densgc * densgc
    return min(max(1.38e-1, k), 1.0)


def th_param(state: FasstState, n, pdens, pdensnew, rhov, rhoda):
    # ---------------- zero-out variables (linhas 28-66) ----------------
    sr = 0.0
    ke = 0.0
    ksat = 0.0
    ks = 0.0
    kw = 0.0
    ki = 0.0
    ka = 0.0  # calculada para solos, nunca lida -- ver docstring
    kv = 0.0  # calculada para solos, nunca lida -- ver docstring
    gs = 0.0
    gw = 0.0
    gi = 0.0
    ga = 0.0
    gv = 0.0
    ns = 0.0
    na = 0.0
    sphmi = 0.0
    sphms = 0.0
    sphmsn = 0.0
    kmi = 0.0
    kms = 0.0
    kmsn = 0.0
    tice = 0.0
    hms = 0.0
    tdens = 0.0
    densgc = 0.0
    ttemp = 0.0
    sphsnow = 0.0
    rhotot = 0.0
    cs = 0.0  # nunca usada -- ver docstring
    ksq = 0.0
    kso = 0.0
    ksn = 0.0
    pq = 0.0
    np_ = 0.0

    # ke = Kersten number; sr = saturation ratio (vol. water/vol. voids)
    # nsoilp(i,2) = porosity; nsoilp(i,5) = quartz content (fraction);
    # nsoilp(i,13) = specific heat of dry soil (J/kg*K);
    # nsoilp(i,14) = organic fraction of the soil
    # thermal diffusivity (m^2/s) = grthcond(i)/grspheat(i)
    # thermal inertia (J/m^2Ks^0.5) = sqrt(grspheat(i)*grthcond(i))

    for i in range(1, state.ntemp + 1):
        if i <= state.nnodes:  # soils, etc....
            # equations for soil thcond are based on Farouki, p.112-116,
            # CRREL Monograph 81-1
            if state.ice[i] + state.soil_moist[i] > state.nsoilp[i][2]:
                np_ = state.ice[i] + state.soil_moist[i]
            else:
                np_ = state.nsoilp[i][2]
            ns = 1.0 - np_

            if state.ntype[i] != 27 and state.ntype[i] != 26:  # not water or air
                na = max(0.0, np_ - (state.soil_moist[i] + state.ice[i] + state.wvc[i]))
            elif state.ntype[i] == 27:
                na = max(0.0, 1.0 - state.wvc[i])
            else:
                na = 0.0

            if state.ntype[i] <= 18:
                kw = thconds(state.stt[i], Phase.WATER)  # W/m*K, water
                ka = thconds(state.stt[i], Phase.DRY_AIR)  # W/m*K, air
                kv = thconds(state.stt[i], Phase.WATER_VAPOR)  # W/m*K, water vapor
                sr = max(EPS, min(1.0, state.soil_moist[i] / state.nsoilp[i][9]))  # unitless

                # thermal conductivity of soil solids (W/m*K)
                ksq = 1.0  # quartz
                # if(nsoilp(i,5) > 0d0) ksq = 8.4d0**nsoilp(i,5)  -- comentado no Fortran
                pq = 3.39e-1 + 4.17e-1 * (state.nsoilp[i][18] + state.nsoilp[i][25])  # Tarnawski et al. (2009)
                if pq > 0.0:
                    ksq = 8.4 * pq  # multiplicacao, nao exponenciacao -- ver docstring
                kso = 1.0  # organics
                if state.nsoilp[i][14] > 0.0:
                    kso = 0.25 ** state.nsoilp[i][14]
                ksn = 1.0  # all other solids
                # versao anterior (comentada no Fortran): usava nsoilp(i,5) em vez de pq
                if abs(1.0 - (pq + state.nsoilp[i][14])) > EPS:
                    ksn = 2.9 ** (1.0 - (pq + state.nsoilp[i][14]))  # 5d0
                ks = ksq * kso * ksn  # total solids

                # saturated thermal conductivity (ksat); wetness factor (ke)
                ke = sr
                if state.ntype[i] == 15:
                    if sr > 5e-2:
                        ke = 0.7 * math.log10(sr) + 1.0
                else:
                    if sr > 0.0:
                        if state.ntype[i] <= 8 or state.ntype[i] == 18:
                            # Lu et al (2007) with Tarnawski (2009) coeff.
                            ke = math.exp(7.28e-1 * (1.0 - sr ** (7.28e-1 - 1.165)))
                        else:
                            ke = math.exp(3.7e-1 * (1.0 - sr ** (3.7e-1 - 1.29)))
                ke = min(max(0.0, ke), 1.0)

                if state.ice[i] > EPS:  # frozen
                    ki = thconds(state.stt[i], Phase.ICE)
                    if state.soil_moist[i] > EPS:
                        ksat = (ks ** ns) * (ki ** state.ice[i]) * (kw ** (np_ - state.ice[i]))  # W/m*K
                    elif state.soil_moist[i] <= EPS:
                        ksat = (ks ** ns) * (ki ** state.ice[i])
                else:  # unfrozen
                    ksat = (kw ** np_) * (ks ** ns)  # W/m*K

                # effective thermal conductivity (W/m*K)
                state.grthcond[i] = state.nsoilp[i][6] + (ksat - state.nsoilp[i][6]) * ke

            elif state.ntype[i] == 20:  # concrete
                if abs(state.nsoilp[i][6] - SPFLAG) <= EPS:
                    state.grthcond[i] = 0.875  # range is 0.05 - 1.7 W/m*K
                else:
                    state.grthcond[i] = state.nsoilp[i][6]
            elif state.ntype[i] == 21:  # asphalt
                if abs(state.nsoilp[i][6] - SPFLAG) <= EPS:
                    state.grthcond[i] = 0.7  # range is 0.15 - 1.4 W/m*K
                else:
                    state.grthcond[i] = state.nsoilp[i][6]
            elif state.ntype[i] == 25:  # bed rock - assumed to be granite
                if abs(state.nsoilp[i][6] - SPFLAG) <= EPS:
                    state.grthcond[i] = 2.0  # range is 2.0 - 3.5 W/m*K
                else:
                    state.grthcond[i] = state.nsoilp[i][6]
            elif state.ntype[i] == 26:  # water
                if abs(state.nsoilp[i][6] - SPFLAG) <= EPS:
                    state.grthcond[i] = thconds(state.stt[i], Phase.WATER)  # W/m*K
                else:
                    state.grthcond[i] = state.nsoilp[i][6]
            elif state.ntype[i] == 27:  # air
                if abs(state.nsoilp[i][6] - SPFLAG) <= EPS:
                    state.grthcond[i] = thconds(state.stt[i], Phase.DRY_AIR)  # W/m*K
                else:
                    state.grthcond[i] = state.nsoilp[i][6]
            elif state.ntype[i] == 30:  # glaciers +/- permanent snow
                if abs(state.nsoilp[i][6] - SPFLAG) <= EPS:
                    state.grthcond[i] = thconds(state.stt[i], Phase.ICE)  # W/m*K
                else:
                    state.grthcond[i] = state.nsoilp[i][6]

            # specific heat (J/m^3*K)
            if state.node_type[i] != 'WA' and abs(state.nsoilp[i][13] - SPFLAG) > EPS:
                gs = ns * (state.nsoilp[i][1] * 1e3) * state.nsoilp[i][13]  # solids
            else:
                gs = 0.0

            if state.ntype[i] != 27 and state.ntype[i] != 26:
                gw = (state.soil_moist[i] * dense(state.stt[i], 0.0, Phase.WATER)
                      * spheats(state.stt[i], Phase.WATER))  # water
                gi = (state.ice[i] * dense(state.stt[i], 0.0, Phase.ICE)
                      * spheats(state.stt[i], Phase.ICE))  # ice
                ga = na * rhoda[i] * spheats(state.stt[i], Phase.DRY_AIR)  # air
                gv = state.wvc[i] * rhov[i] * spheats(state.stt[i], Phase.WATER_VAPOR)  # water vapor
            elif state.ntype[i] == 27:
                gw = 0.0
                gi = 0.0
                ga = na * rhoda[i] * spheats(state.stt[i], Phase.DRY_AIR)
                gv = state.wvc[i] * rhov[i] * spheats(state.stt[i], Phase.WATER_VAPOR)
            elif state.ntype[i] == 26:
                gw = (state.soil_moist[i] * dense(state.stt[i], 0.0, Phase.WATER)
                      * spheats(state.stt[i], Phase.WATER))
                gi = 0.0
                ga = 0.0
                gv = 0.0
            state.grspheat[i] = gs + gw + gi + ga + gv  # J/m^3*K

        elif i > state.nnodes:
            if state.node_type[i] == 'HM' or state.node_type[i] == 'MX':  # snow, ice
                j = i - state.nnodes
                # calculate snow/ice thermal conductivity
                # Sturm et al. (1997), J. Glaciology 43(143), pp.26-41
                kmi = 0.0
                kms = 0.0
                kmsn = 0.0
                sphmi = 0.0
                sphms = 0.0
                sphmsn = 0.0
                ttemp = state.stt[i]
                sphsnow = 2.09e3  # J/kg*K

                if state.stt[i] > TREF:
                    ttemp = TREF

                if state.hsaccum > EPS:
                    if state.iw >= 2:
                        tdens = state.sdens[state.iw - 1]  # kg/m^3
                        if abs(state.sdens[state.iw - 1]) <= EPS:
                            tdens = pdensnew
                    else:
                        tdens = SDENSW
                    densgc = tdens * 1e-3  # g/cm^3
                    kms = _sturm_conductivity(densgc)  # W/m*K
                    sphms = sphsnow * tdens  # J/m^3*K, specific heat of snow

                if state.newsd > EPS:
                    densgc = pdens * 1e-3
                    kmsn = _sturm_conductivity(densgc)
                    sphmsn = sphsnow * pdens  # J/m^3*K, specific heat of snow

                if state.hi > EPS:
                    if abs(state.hsaccum + state.newsd) <= EPS:
                        state.km[j] = thconds(ttemp, Phase.ICE)  # W/m*K, thermal conductivity ice
                        state.sphm[j] = spheats(ttemp, Phase.ICE) * dense(ttemp, 0.0, Phase.ICE)  # J/m^3*K
                    else:
                        if state.hsaccum > EPS and abs(state.newsd) <= EPS:
                            # `tice`: fracao invertida -- ver "Pontos de atencao"
                            tice = (
                                (kms / state.hsaccum + 2.29 / state.hi)
                                / ((2.29 / state.hi) * state.stt[state.nnodes]
                                   + (kms / state.hsaccum) * ttemp)
                            )
                            kmi = thconds(tice, Phase.ICE)
                            sphmi = spheats(tice, Phase.ICE) * dense(tice, 0.0, Phase.ICE)

                            rhotot = ((tdens * state.hsaccum + dense(tice, 0.0, Phase.ICE) * state.hi)
                                      / (state.hsaccum + state.hi))
                            state.km[j] = (state.hsaccum + state.hi) / (state.hi / kmi + state.hsaccum / kms)
                            state.sphm[j] = (
                                (sphmi * dense(tice, 0.0, Phase.ICE) + sphms * tdens) / rhotot
                            )  # J/m^3*K, snow/ice layer
                        elif abs(state.hsaccum) <= EPS and state.newsd > EPS:
                            tice = (
                                (kmsn / state.newsd + 2.29 / state.hi)
                                / ((2.29 / state.hi) * state.stt[state.nnodes]
                                   + (kmsn / state.newsd) * ttemp)
                            )
                            kmi = thconds(tice, Phase.ICE)
                            sphmi = spheats(tice, Phase.ICE) * dense(tice, 0.0, Phase.ICE)

                            rhotot = ((pdens * state.newsd + dense(tice, 0.0, Phase.ICE) * state.hi)
                                      / (state.newsd + state.hi))
                            state.km[j] = (state.hi + state.newsd) / (state.hi / kmi + state.newsd / kmsn)
                            state.sphm[j] = (
                                (sphmi * dense(tice, 0.0, Phase.ICE) + sphmsn * pdens) / rhotot
                            )  # J/m^3*K, snow/ice layer
                        else:  # both old and new snow
                            hms = state.hsaccum + state.newsd
                            kms = hms / (state.hsaccum / kms + state.newsd / kmsn)

                            tice = (
                                (kms / hms + 2.29 / state.hi)
                                / ((2.29 / state.hi) * state.stt[state.nnodes] + (kms / hms) * ttemp)
                            )
                            kmi = thconds(tice, Phase.ICE)
                            sphmi = spheats(tice, Phase.ICE) * dense(tice, 0.0, Phase.ICE)

                            rhotot = (
                                (tdens * state.hsaccum + dense(tice, 0.0, Phase.ICE) * state.hi
                                 + pdens * state.newsd)
                                / (state.hsaccum + state.hi + state.newsd)
                            )
                            state.km[j] = (state.hi + hms) / (
                                state.hi / kmi + state.hsaccum / kms + state.newsd / kmsn
                            )
                            state.sphm[j] = (
                                (sphmi * dense(tice, 0.0, Phase.ICE) + sphms * tdens + sphmsn * pdens)
                                / rhotot
                            )  # J/m^3*K, snow/ice layer
                else:  # no ice layer
                    if state.hsaccum > EPS and abs(state.newsd) <= EPS:
                        rhotot = tdens
                        state.km[j] = kms
                        state.sphm[j] = sphms
                    elif abs(state.hsaccum) <= EPS and state.newsd > EPS:
                        rhotot = pdens
                        state.km[j] = kmsn
                        state.sphm[j] = sphmsn
                    else:
                        # `_fdiv`: 0/0 quando hsaccum e newsd sao zero -- ver docstring
                        rhotot = _fdiv(pdens * state.newsd + tdens * state.hsaccum,
                                       state.newsd + state.hsaccum)
                        state.km[j] = _fdiv(
                            state.hsaccum + state.newsd,
                            _fdiv(state.hsaccum, kms) + _fdiv(state.newsd, kmsn),
                        )
                        state.sphm[j] = _fdiv(sphms * tdens + sphmsn * pdens, rhotot)  # J/m^3*K

                state.km[j] = _anint(state.km[j], 20)
                state.sphm[j] = _anint(state.sphm[j], 20)

                if state.node_type[i] == 'MX':
                    state.grthcond[i] = (1.0 - state.sigfl) * state.km[j] + state.sigfl * KVEG
                    state.grspheat[i] = (
                        (1.0 - state.sigfl) * state.sphm[j]
                        + state.sigfl * SPHVEG * dense(state.stt[i], 0.0, Phase.WATER)
                    )
                else:
                    state.grthcond[i] = state.km[j]
                    state.grspheat[i] = state.sphm[j]
            elif state.node_type[i] == 'VG':
                state.grthcond[i] = state.sigfl * KVEG
                state.grspheat[i] = state.sigfl * SPHVEG * dense(state.stt[i], 0.0, Phase.WATER)

        state.grthcond[i] = _anint(state.grthcond[i], 20)
        state.grspheat[i] = _anint(state.grspheat[i], 20)

    return rhotot
