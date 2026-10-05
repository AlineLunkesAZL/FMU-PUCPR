"""
soil_moisture.py -- traducao de soil_moisture.F90

Fonte Fortran: soil_moisture.F90 (402 linhas). Uma unica sub-rotina
publica (`soil_moisture`). Calcula o perfil de umidade do solo
(soil_moist) e a cabeca de pressao (phead) em cada no, resolvendo a
equacao de Richards linearizada (variavel dependente: cabeca de
pressao) pelo metodo de Newton-Raphson, com o incremento
delta-phead(i) obtido por eliminacao tridiagonal (algoritmo de Thomas).
Tambem calcula o escoamento superficial (runoff) por excesso de agua e
tres medidas de erro/diagnostico (rhs_errorm, errorm, smerror) usadas
pelo chamador para decidir se mais iteracoes de Newton-Raphson sao
necessarias.

Assinatura original:
    subroutine soil_moisture(sm_old,wvco,iceo,sourceo,sinkro,vino,
                             told,thvc,dthvdh,runoff,rhs_errorm,
                             errorm,smerror,sert1)
    real(dp),intent(in):: sm_old(ntot),wvco(ntot),iceo(ntot)
    real(dp),intent(in):: sourceo(ntot),sinkro(ntot),vino(ntot)
    real(dp),intent(in):: told(ntot),thvc(ntot),dthvdh(ntot)
    real(dp),intent(inout):: runoff(ntot)
    real(dp),intent(out):: rhs_errorm,errorm,smerror,sert1

Arquitetura
-------------
`soil_moisture` recebe `state: FasstState` como argumento explicito.
`EPS` e `TREF` vem de `fasst.constants`; `ip_pt`/`ip_rh` (indices de
coluna do meteorologico, no Fortran) sao `MetCol.PT`/`MetCol.RH`.
`dense` (pura) e `head(state,i,smt)` sao chamadas reais, importadas de
`fasst.functions`. O cabecalho do .F90 tambem lista `soilhumid` entre
as funcoes usadas ("uses the function: dense,head,soilhumid") -- mas
`soilhumid` nunca e de fato chamada em lugar nenhum do corpo desta
sub-rotina.

`runoff` e o unico parametro `intent(inout)` (um array, nao um
escalar) -- por consistencia com o resto do projeto (que sempre isola
copias locais de arrays `intent(inout)`/`out` em vez de depender de
mutacao in-place por referencia), a funcao recebe `runoff` e
imediatamente faz uma copia local (`runoff = list(runoff)`) antes de
mutar -- quem chamar deve reatribuir a partir do retorno em tupla,
exatamente como os quatro escalares `intent(out)` (rhs_errorm, errorm,
smerror, sert1):

    runoff, rhs_errorm, errorm, smerror, sert1 = soil_moisture(
        state, sm_old, wvco, iceo, sourceo, sinkro, vino, told, thvc,
        dthvdh, runoff,
    )

Os demais parametros (sm_old, wvco, iceo, sourceo, sinkro, vino, told,
thvc, dthvdh) sao `intent(in)` puros -- nunca escritos no Fortran,
lidos diretamente aqui sem copia.

Indexacao 1-based identica ao Fortran em todos os vetores/matrizes
(posicao 0 sem uso), para os arrays LOCAIS desta sub-rotina (A, B, C,
D, gam, delphead -- todos dimensionados por nnodes, nao ntot: o
sistema tridiagonal so cobre os nos de solo/agua, nunca os nos de ar
acima do perfil). As dezenas de campos de `FasstState` usados
(soil_moist, phead, ice, wvc, vin, source, sink, sinkr, khl, khu,
dsmdh, pheadmin, nsoilp, node_type, ntype, stt, met, delzs, frl, frh)
tem essa mesma convencao.

Pontos de atencao
-------------------
* `sd` nunca deixa de ser 0.0: o Fortran zera `sd` no fim de cada
  iteracao do laco principal (`sd = 0d0`) e so o alteraria de volta
  para um valor diferente de zero atraves de uma linha comentada
  (`!sd = excess`, dentro de um `if(i /= 1)` tambem comentado) -- como
  essa linha esta desativada, `sd` e sempre 0.0 quando lido no inicio
  da proxima iteracao (em phead(i), sms, smv). Na pratica, todo o
  excesso de umidade vira runoff incondicionalmente (linha ativa logo
  abaixo do bloco comentado), nunca e "empurrado" para o no de baixo
  via `sd`. Preservado fielmente -- inclusive o uso de `sd` nas
  formulas, mesmo sabendo que seu valor e sempre zero.
* `pho` e `d0i` sao calculados/zerados mas nunca lidos de novo em
  nenhum lugar do arquivo -- codigo morto inofensivo, preservado por
  fidelidade.
"""

import math

from .constants import EPS, TREF
from .indices import MetCol
from .state import FasstState
from .functions import Phase, dense, head


def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _anint(x, p=0):
    """Fortran ``ANINT(x*10**p)*10**-p`` -- arredonda para longe de zero."""
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def soil_moisture(state: FasstState, sm_old, wvco, iceo, sourceo, sinkro, vino, told,
                   thvc, dthvdh, runoff):
    nnodes = state.nnodes

    # copia local do unico parametro intent(inout) (ver docstring)
    runoff = list(runoff)

    # ---------------- zero-out variables (linhas 37-76) ----------------
    d0i = 0
    ni = 0
    nb = 0
    ni1 = 0
    nb1 = 0
    sn2 = 0
    sn3 = 0
    f1 = 0.0
    f2 = 0.0
    f3 = 0.0
    rsinkl = 0.0
    rsinkh = 0.0
    ff = 0.0
    w1 = 0.0
    c1 = 0.0
    bet = 0.0
    fts = 0.0
    smt = 0.0
    totmoisture = 0.0
    oldmoist = 0.0
    totsink = 0.0
    pho = 0.0
    smh = 0.0
    excess = 0.0
    t1 = 0.0
    t3 = 0.0
    sm_old1 = 0.0
    sd = 0.0
    rhow = 0.0
    rhoi = 0.0
    tsign = 0.0
    dwi = 0.0
    didh = 0.0
    dvdh = 0.0
    smv = 0.0
    sumvin = 0.0
    sms = 0.0
    sert1 = 0.0

    # ---------------- zero out arrays (linhas 78-88) ----------------
    delphead = _new_array(nnodes)
    A = _new_array(nnodes)
    B = _new_array(nnodes)
    C = _new_array(nnodes)
    D = _new_array(nnodes)
    gam = _new_array(nnodes)

    for i in range(1, nnodes + 1):
        delphead[i] = 0.0
        A[i] = 0.0
        B[i] = 0.0
        C[i] = 0.0
        D[i] = 0.0
        gam[i] = 0.0
        runoff[i] = 0.0
        state.sinkr[i] = 0.0

    # ---------------- initialize errors; constants (linhas 90-97) ----------------
    rhs_errorm = -abs(state.mflag)
    errorm = -abs(state.mflag)
    smerror = -abs(state.mflag)

    fts = 1.0 / state.deltat_fasst

    # ---------------- determine sources and sinks (linhas 99-132) ----------------
    for i in range(1, nnodes + 1):
        if state.ntype[i] != 27:
            if state.veg_flagl == 0 or abs(state.sigfl) <= EPS:  # low vegetation
                rsinkl = 0.0  # root uptake (unitless)
                state.frl[i] = 0.0
            else:
                if state.trmlm > EPS:
                    ff = 1.0
                    if state.stt[i] <= TREF:
                        ff = 0.0
                    rsinkl = 1.5e-7 * state.sigfl * ff * state.frl[i]  # m/s
                else:
                    rsinkl = 0.0

            if state.veg_flagh == 0:  # high vegetation (trees)
                rsinkh = 0.0
                state.frh[i] = 0.0
            else:
                if state.trmhm > EPS:
                    ff = 1.0
                    if state.stt[i] <= TREF:
                        ff = 0.0
                    rsinkh = 1.5e-7 * state.sigfh * ff * state.frh[i]  # m/s
                else:
                    rsinkh = 0.0

            state.sinkr[i] = (rsinkl + rsinkh) * state.deltat_fasst / state.delzs[i]  # unitless
        else:
            state.sinkr[i] = 0.0

    # ---------------- monta a matriz tridiagonal (linhas 134-183) ----------------
    # A(i) = coeff. para delphead(i-1); B(i) = coeff. para delphead(i);
    # C(i) = coeff. para delphead(i+1); D(i) = rhs.
    # z e positivo para cima a partir do fundo do perfil.
    for i in range(1, nnodes + 1):
        if state.ntype[i] != 27:
            dwi = 0.0
            diw = 0.0
            didh = 0.0
            dvdh = 0.0
            if state.ice[i] + iceo[i] > EPS:
                dwi = 0.5 * (
                    dense(state.stt[i], 0.0, Phase.WATER) / dense(state.stt[i], 0.0, Phase.ICE)
                    + dense(told[i], 0.0, Phase.WATER) / dense(told[i], 0.0, Phase.ICE)
                )  # unitless (water/ice)
                diw = 0.5 * (
                    dense(state.stt[i], 0.0, Phase.ICE) / dense(state.stt[i], 0.0, Phase.WATER)
                    + dense(told[i], 0.0, Phase.ICE) / dense(told[i], 0.0, Phase.WATER)
                )  # unitless (ice/water)
                didh = -dwi * state.dsmdh[i]

            if state.wvc[i] + wvco[i] > EPS:
                dvdh = (
                    dthvdh[i] * (state.nsoilp[i][2] - (state.soil_moist[i] + state.ice[i]))
                    - thvc[i] * (1 - dwi) * state.dsmdh[i]
                )  # unitless (vapor/water)

            f1 = state.deltat_fasst / (2.0 * state.delzs[i])  # s/m
            f2 = (
                (state.soil_moist[i] - sm_old[i])
                + (state.wvc[i] - wvco[i])
                + diw * (state.ice[i] - iceo[i])
            )  # unitless
            f3 = f1 * (state.source[i] + sourceo[i]) - 0.5 * (state.sinkr[i] + sinkro[i])  # unitless

            A[i] = -f1 * state.khl[i]  # 1/m
            B[i] = f1 * (state.khl[i] + state.khu[i]) + state.dsmdh[i] + didh + dvdh  # 1/m
            C[i] = -f1 * state.khu[i]  # 1/m
            D[i] = -f2 - f1 * (state.vin[i] + vino[i]) + f3  # unitless

            A[i] = _anint(A[i], 15)
            B[i] = _anint(B[i], 15)
            C[i] = _anint(C[i], 15)
            D[i] = _anint(D[i], 15)

            if abs(D[i]) > rhs_errorm:
                rhs_errorm = abs(D[i])
                sert1 = D[i]
        elif state.ntype[i] == 27:
            D[i] = 0.0
            if sn2 == 0:
                sn2 = i
            sn3 = i

    # ---------------- solve the matrix eqn. for delta phead(i) (linhas 185-203) ----------------
    # determina os intervalos [ni,nb] (e [ni1,nb1], se houver uma camada
    # de ar entre duas camadas de solo/agua) onde a matriz tridiagonal e
    # de fato resolvida.
    if sn2 == 0 and sn3 == 0:  # no air nodes
        ni = 1
        nb = nnodes
    elif sn2 == 1 and sn3 == nnodes:  # all air nodes
        ni = 0
        nb = 0
    elif sn2 > 1 and sn3 == nnodes:  # air on top
        ni = 1
        nb = sn2 - 1
    elif sn2 == 1 and sn3 < nnodes:  # air on bottom
        ni = sn3 + 1
        nb = nnodes
    elif sn2 > 1 and sn3 < nnodes:  # air layer
        ni = 1
        nb = sn2 - 1
        ni1 = sn3 + 1
        nb1 = nnodes

    if ni > 0 and nb > 0:
        bet = B[ni]  # 1/m
        if abs(bet) <= 1e-10:
            tsign = 1.0
            if bet < 0.0:
                tsign = -1.0
            bet = tsign * 1e-10
            delphead[ni] = tsign * abs(D[ni])
        else:
            delphead[ni] = D[ni] / bet  # m

        for j in range(ni + 1, nb + 1):
            gam[j] = C[j - 1] / bet  # unitless
            bet = B[j] - A[j] * gam[j]  # 1/m
            if abs(bet) <= 1e-10:
                tsign = 1.0
                if bet < 0.0:
                    tsign = -1.0
                bet = 1e-10 * tsign  # tridiag fails, not enough moisture flow
                delphead[j] = tsign * abs(D[j])
            else:
                delphead[j] = (D[j] - A[j] * delphead[j - 1]) / bet  # m

        for j in range(nb - 1, ni - 1, -1):
            delphead[j] = delphead[j] - gam[j + 1] * delphead[j + 1]  # m

    if ni1 > 0 and nb1 > 0:
        bet = B[ni1]  # 1/m
        if abs(bet) <= 1e-10:
            tsign = 1.0
            if bet < 0.0:
                tsign = -1.0
            bet = tsign * 1e-10
            delphead[ni1] = tsign * abs(D[ni1])
        else:
            delphead[ni1] = D[ni1] / bet  # m

        for j in range(ni1 + 1, nb1 + 1):
            gam[j] = C[j - 1] / bet  # unitless
            bet = B[j] - A[j] * gam[j]  # 1/m
            if abs(bet) <= 1e-10:
                tsign = 1.0
                if bet < 0.0:
                    tsign = -1.0
                bet = 1e-10 * tsign  # tridiag fails, not enough moisture flow
                delphead[j] = tsign * abs(D[j])
            else:
                delphead[j] = (D[j] - A[j] * delphead[j - 1]) / bet  # m

        for j in range(nb1 - 1, ni1 - 1, -1):
            delphead[j] = delphead[j] - gam[j + 1] * delphead[j + 1]  # m

    # ---------------- determine the maximum error, update phead(i); ----------------
    # ---------------- calculate the corresponding soil moisture (linhas 263-399) ----------------
    for i in range(nnodes, 0, -1):
        delphead[i] = _anint(delphead[i], 15)

        if state.ntype[i] != 27:
            rhow = dense(state.stt[i], 0.0, Phase.WATER)
            rhoi = dense(state.stt[i], 0.0, Phase.ICE)
            dwi = rhoi / rhow
            f1 = state.soil_moist[i] + state.ice[i] * dwi
            f3 = _anint(delphead[i], 15)

            tsign = 0.0
            if abs(state.vin[i]) > EPS:
                tsign = -abs(state.vin[i]) / state.vin[i]
            elif abs(state.vin[i]) <= EPS and abs(delphead[i]) > EPS:
                tsign = abs(delphead[i]) / delphead[i]

            if abs(delphead[i]) > abs(state.pheadmin[i]) * 2.5e-2:
                delphead[i] = abs(state.pheadmin[i]) * 2.5e-2 * tsign

            delphead[i] = _anint(delphead[i], 15)

            if abs(delphead[i]) > errorm:
                errorm = abs(delphead[i])  # m

            sm_old1 = state.soil_moist[i]
            pho = state.phead[i]  # nunca lida de novo (ver docstring)
            excess = 0.0

            if i < ni:
                state.phead[i] = state.phead[i]  # m -- no-op fiel ao Fortran
                smt = state.soil_moist[i]  # unitless
            else:
                state.phead[i] = state.phead[i] + delphead[i] + sd * state.delzs[i]  # m

                if state.phead[i] < 0.0 and state.phead[i] > state.pheadmin[i]:
                    iflag = 0
                    c1 = 1.0 + (state.nsoilp[i][10] * abs(state.phead[i] * 1e2)) ** state.nsoilp[i][11]
                    w1 = max(0.0, 1.0 / (c1 ** state.nsoilp[i][12]))
                    smh = state.nsoilp[i][8] + w1 * (state.nsoilp[i][9] - state.nsoilp[i][8])
                    smh = max(0.0, smh)
                else:
                    if state.phead[i] >= 0.0:
                        iflag = 2
                    else:
                        iflag = 1
                    smh = state.nsoilp[i][24] - state.ice[i] * rhoi / rhow

                t3 = state.dsmdh[i] * delphead[i]
                sms = max(0.0, state.soil_moist[i] + t3 + sd)  # unitless

                smv = sm_old1 - state.vin[i] * (state.deltat_fasst / state.delzs[i]) + sd
                smv = max(0.0, smv)

                if iflag == 0 and state.ice[i] + iceo[i] < EPS:
                    smt = smh
                    if smt > state.nsoilp[i][24] or abs(f3 - delphead[i]) > EPS:
                        smt = sms
                else:
                    if iflag == 2 and _anint(state.met[state.iw][MetCol.PT]) == 2:
                        smt = smh
                    else:
                        smt = sms
                        # comentado no Fortran: if(ice(i) > eps) smt = 0.5*(sms + smv)

                if (state.node_type[i] == 'CO' or state.node_type[i] == 'AS'
                        or state.node_type[i] == 'RO'):
                    smt = smv

                smt = _anint(smt, 10)
                f2 = smt + state.ice[i] * dwi

                # check for continuity in time
                t1 = abs(state.vin[i] * state.deltat_fasst / state.delzs[i]) + sd
                if t1 < EPS:
                    t1 = 1e-2

                if abs(f2 - f1) > t1:
                    if f2 > f1:
                        f2 = f1 + t1
                    else:
                        f2 = f1 - t1
                    smt = f2 - state.ice[i] * dwi
                smt = _anint(smt, 10)

                # check for over/under allowed water content
                excess = 0.0
                if state.ice[i] <= EPS:
                    if smt - state.nsoilp[i][15] < EPS:
                        smt = state.nsoilp[i][15]
                    elif smt - state.nsoilp[i][24] > EPS:
                        excess = smt - state.nsoilp[i][24]  # unitless
                        smt = state.nsoilp[i][24]
                elif state.ice[i] > EPS:
                    if smt < EPS:
                        smt = 0.0
                    elif smt > state.nsoilp[i][24]:
                        excess = smt - state.nsoilp[i][24]
                        smt = state.nsoilp[i][24]
                smt = _anint(smt, 10)
                excess = _anint(excess, 10)

                # sd sempre fica 0.0 daqui em diante (ver docstring do
                # modulo): o unico jeito de virar != 0 esta comentado no
                # Fortran original.
                sd = 0.0
                if excess > 0.0:
                    # bloco comentado no Fortran, preservado como comentario:
                    #   if(i /= 1): sd = excess
                    #   else if(i == 1): runoff(i) = runoff(i) + excess*delzs(i)*fts
                    # como o if/else esta comentado, a linha do runoff roda
                    # incondicionalmente, para qualquer i:
                    runoff[i] = runoff[i] + excess * state.delzs[i] * fts  # m/s

            state.sink[i] = state.sinkr[i] * state.delzs[i] * fts  # m/s
            state.sink[i] = _anint(state.sink[i], 10)
            runoff[i] = _anint(runoff[i], 10)

            state.soil_moist[i] = smt
            state.phead[i] = head(state, i, state.soil_moist[i])

            totmoisture = totmoisture + state.soil_moist[i] * state.delzs[i]  # m
            totsink = totsink + state.sinkr[i] * state.delzs[i] + runoff[i] * state.deltat_fasst  # m
            oldmoist = oldmoist + sm_old1 * state.delzs[i]  # m
            sumvin = sumvin + state.vin[i] * state.deltat_fasst  # m
            smerror = abs((totmoisture - oldmoist) - (totsink + sumvin))  # m
        else:
            if abs(delphead[i]) > errorm:
                errorm = abs(delphead[i])  # m
            state.soil_moist[i] = 1e-2 * state.met[state.iw][MetCol.RH]

    return runoff, rhs_errorm, errorm, smerror, sert1
