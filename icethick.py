"""
icethick.py -- traducao de icethick.F90

Fonte Fortran: icethick.F90 (119 linhas). Uma unica sub-rotina publica
(`icethick`). Calcula a acrecao/depleicao de gelo de superficie (gelo
de estrada/pista, nao gelo de lago) em um unico timestep: espessura
final `hfinal` a partir da espessura inicial `hinit`, mais dois
diagnosticos escritos em `FasstState` (`vimelt`, a agua de degelo
liberada, e `refreezei`, a fracao de recongelamento). Se ha neve
acumulada, nada acontece -- o modelo de neve (module_snow) e quem cuida
do congelamento/degelo naquele caso.

Assinatura original:
    subroutine icethick(hinit,hfinal)
    real(dp),intent(in):: hinit
    real(dp),intent(out):: hfinal

Arquitetura
-------------
`icethick` recebe `state: FasstState` como argumento explicito e
devolve `hfinal` (o unico `intent(out)`). As duas outras saidas da
sub-rotina sao variaveis de modulo no Fortran e continuam sendo campos
de estado aqui, mutados no lugar: `state.vimelt` e `state.refreezei`.
Mesma convencao de `soil_moisture` (escalares `intent(out)` retornados,
estado de modulo mutado no lugar).

Fan-out ZERO: nao chama nenhuma sub-rotina e nenhuma funcao do projeto
-- so intrinsecas do Fortran (`aint`, `int`, `dabs`, `dmax1`, `dmin1`,
`anint`). `eps`, `Tref`, `idens`, `lhfus` sao `parameter` de
fasst_global e vieram para `fasst.constants` (EPS, TREF, IDENS, LHFUS);
`ip_pt`/`ip_prec` sao `MetCol.PT`/`MetCol.PREC`. Isso torna `icethick`
validavel isoladamente HOJE, sem esperar nenhum outro modulo (ver o doc
de transcricao, secao VALIDACAO).

Fan-in 1: o unico chamador e fasst_main.F90:343, dentro de
`if(node_type(nnodes) /= 'WA')` e de
`if(hi > eps .or. pt == 2 .or. pt == 4)`:

    if(dabs(hsaccum) <= eps) call icethick(oldhi,hi)

ou seja, `hinit` = `oldhi` e `hfinal` -> `hi`. Como `fasst_main.F90`
ainda nao foi transcrito, o call site Python ainda nao existe.

A ZERAGEM INICIAL E OBRIGATORIA (linhas 25-34), NAO E CERIMONIA
------------------------------------------------------------------
Varios caminhos nao entram em NENHUM dos dois ramos externos e
dependem literalmente desses zeros:

  * `pt == 3` (neve) -- nem o ramo de chuva (2/4) nem o de "sem
    chuva" (1). Alcancavel pelo call site quando `hi > eps`.
  * `pt` 2 ou 4 com `dabs(hsaccum) > eps` -- excluido pelo guard
    interno (hoje inalcancavel, ver ponto 1 do doc).
  * `pt` 2 ou 4, `hsaccum` zero, mas `hinit` NEGATIVO -- passa o teste
    externo e falha os dois internos (`dabs(hinit) <= eps` e
    `hinit > eps`).

Em Fortran esses caminhos caem nos zeros das linhas 25-34; em Python,
sem a zeragem explicita, seriam `UnboundLocalError`. Por isso o bloco
de zeragem abaixo e transcrito literalmente, um a um, e NAO deve ser
"limpo".

Medido por mutacao na validacao (tests/test_icethick_ref.py): dos 10
zeros, exatamente QUATRO sao load-bearing -- `hnew1`, `hnew2`,
`vimelt1`, `vimelt2`, os lidos depois dos ramos (linhas 103, 107, 111,
113); remover qualquer um deles da `UnboundLocalError` nos caminhos de
vao. Os outros seis (`qtop`, `freeze_frac`, `qb_f1`, `qbot`, `hfinal`,
`vimelt`) sao redundantes: so sao lidos onde ja foram atribuidos, ou
(caso de `hfinal`/`vimelt`) sao sobrescritos incondicionalmente nas
linhas 103/107. Mantidos mesmo assim, por fidelidade 1-para-1.

Convencao de sinal de `melt(iw)` / `lheat(iw)`
-----------------------------------------------
`lheat1 = melt(iw)` (fluxo de energia de degelo, W/m^2). Positivo =
energia disponivel para DERRETER; negativo = energia sendo retirada,
disponivel para CONGELAR (`lheat1 < 0d0` e exatamente o teste de
congelamento nas linhas 41 e 53). Confirmado num segundo call site
independente: module_snow.F90:716 condiciona a `melt(iw) > 0d0` e as
linhas 727/732 fazem `atop` -- uma profundidade POSITIVA de
degelo/sublimacao -- proporcional a `melt(iw)`.

`lheat(iw)` (usado em `freeze_frac`) e um array DIFERENTE de
`melt(iw)`: `freeze_frac = lheat(iw)/melt(iw)` le como uma razao de uma
grandeza consigo mesma, e nao e -- ver ponto 3 do doc.
"""

import math

from .constants import EPS, IDENS, LHFUS, TREF
from .indices import MetCol
from .state import FasstState


# icethick.F90:21 -- `real(dp),parameter:: pdens = 1d3` (densidade da
# precipitacao, kg/m^3). Local da sub-rotina no Fortran; nao existe
# constante de densidade de agua em fasst.constants, mantida local aqui.
_PDENS = 1e3


def _anint(x, p=0):
    """Fortran ``ANINT(x*10**p)*10**-p`` -- arredonda para longe de zero.

    Identica a `soil_moisture._anint` / `snow._anint`; `fasst.functions`
    so expoe `anint10`/`anint20` e este arquivo precisa tambem de 1d15.
    """
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


def icethick(state: FasstState, hinit: float) -> float:
    """Espessura de gelo de superficie apos um timestep.

    Devolve `hfinal` (m). Muta `state.vimelt` (m de agua de degelo) e
    `state.refreezei` (fracao de recongelamento, adimensional).
    """
    iw = state.iw
    nnodes = state.nnodes

    # ---------------- initialize variables (linhas 25-34) ----------------
    # NAO REMOVER -- ver "A ZERAGEM INICIAL E OBRIGATORIA" no docstring do
    # modulo: ha caminhos que nao entram em nenhum ramo.
    hnew1 = 0.0
    qtop = 0.0
    vimelt1 = 0.0
    freeze_frac = 0.0
    qb_f1 = 0.0
    qbot = 0.0
    hnew2 = 0.0
    vimelt2 = 0.0
    hfinal = 0.0
    state.vimelt = 0.0

    lheat1 = state.melt[iw]  # W/m^2; > 0 derrete, < 0 congela (ver docstring)

    # `aint` na linha 38 e `int` na linha 75 no original -- ambos truncam
    # para zero, igual a int() aqui; a inconsistencia e cosmetica.
    pt = int(state.met[iw][MetCol.PT])

    if (pt == 2 or pt == 4) and abs(state.hsaccum) <= EPS:   # rain or freezing rain
        if abs(hinit) <= EPS:                                # no ice to start
            # NOTA: aqui o teste de congelamento usa stt(nnodes) (temperatura
            # do no de topo do perfil); no ramo `hinit > eps` abaixo usa
            # toptemp. Assimetria do original, preservada -- ponto 2 do doc.
            if state.stt[nnodes] <= TREF and lheat1 < 0.0:
                freeze_frac = max(0.0, min(1.0, state.lheat[iw] / lheat1))  # based on crrel rep 96-2
                hnew1 = (freeze_frac * state.met[iw][MetCol.PREC] * 1e-3
                         * state.timstep * (_PDENS / IDENS))                # m
                hnew2 = 0.0
            else:
                hnew1 = 0.0
                hnew2 = 0.0
            vimelt1 = 0.0                                                   # m
            vimelt2 = 0.0                                                   # m
        elif hinit > EPS:                                                   # ice to start
            if state.toptemp <= TREF and lheat1 < 0.0:                      # freezing
                freeze_frac = max(0.0, min(1.0, state.lheat[iw] / lheat1))  # based on crrel rep 96-2
                hnew1 = (freeze_frac * state.met[iw][MetCol.PREC] * 1e-3
                         * state.timstep * (_PDENS / IDENS))                # m
                vimelt1 = 0.0
            else:
                qtop = -lheat1                                              # W/m^2
                # SEM clamp em 0 aqui (ao contrario da linha 79) e sem limite
                # contra hinit: hnew1 pode ficar mais negativo que -hinit, e e
                # dai que vem o hfinal NEGATIVO que a linha 104 apara.
                hnew1 = (state.timstep * 3.6e3) * qtop / (IDENS * LHFUS)     # m
                vimelt1 = abs(hnew1) * IDENS * 1e-3                          # m (from the top)

            if state.stt[nnodes] > TREF:                                    # bottom temp > 0
                qb_f1 = (((state.grthcond[nnodes] + state.grthcond[nnodes - 1]) * 5e-1)
                         / ((state.nz[nnodes] - state.nz[nnodes - 1]) * 5e-1))  # W/m^2 K
                qbot = qb_f1 * (state.stt[nnodes] - TREF)                    # W/m^2
                hnew2 = max(0.0, min(hinit - abs(hnew1),
                                     (state.timstep * 3.6e3) * qbot / (IDENS * LHFUS)))  # m, bottom melt depth
            else:
                hnew2 = 0.0
            vimelt2 = hnew2 * IDENS * 1e-3
    elif pt == 1 and abs(state.hsaccum) <= EPS:                             # no rain
        if hinit > EPS:                                                     # ice to start
            if state.toptemp > TREF:
                qtop = -lheat1                                              # W/m^2
                # dmax1(0d0,...): com toptemp > Tref o esperado e derreter
                # (lheat1 > 0 => qtop < 0 => hnew1 < 0), e o clamp zera isso.
                # Preservado literalmente -- ver ponto 2 do doc.
                hnew1 = max(0.0, (state.timstep * 3.6e3) * qtop / (IDENS * LHFUS))
            else:
                hnew1 = 0.0
            vimelt1 = abs(hnew1) * IDENS * 1e-3                             # m (from the top)

            if state.stt[nnodes] > TREF:                                    # bottom temp > 0
                qb_f1 = (((state.grthcond[nnodes] + state.grthcond[nnodes - 1]) * 5e-1)
                         / ((state.nz[nnodes] - state.nz[nnodes - 1]) * 5e-1))  # W/m^2 K
                qbot = qb_f1 * (state.stt[nnodes] - TREF)                    # W/m^2
                # `3600d0` no original aqui, `3.6d3` na linha 69 -- mesmo
                # valor, literal diferente; preservado como esta.
                hnew2 = max(0.0, min(hinit - abs(hnew1),
                                     (state.timstep * 3600.0) * qbot / (IDENS * LHFUS)))  # melt depth bottom temp >0
            else:
                hnew2 = 0.0
            vimelt2 = hnew2 * IDENS * 1e-3
        else:                                                               # no ice, no melting
            hnew1 = 0.0
            hnew2 = 0.0
            vimelt1 = 0.0
            vimelt2 = 0.0

    hfinal = hinit + hnew1 - hnew2
    # `<`, nao `dabs(...) <`: apara tanto valores minusculos positivos quanto
    # o hfinal GENUINAMENTE NEGATIVO do degelo de topo sem limite (ver acima).
    if hfinal < 1e-5:
        hfinal = 0.0
    hfinal = _anint(hfinal, 10)

    state.vimelt = vimelt1 + vimelt2
    state.vimelt = _anint(state.vimelt, 15)

    # ORDEM IMPORTA: refreezei divide pelo hfinal JA aparado e JA quantizado.
    if hinit > EPS:
        state.refreezei = (hnew1 - hnew2) / hinit
    elif abs(hinit) <= EPS and hfinal > EPS:
        state.refreezei = (hnew1 - hnew2) / hfinal
    else:
        state.refreezei = 0.0
    state.refreezei = _anint(state.refreezei, 15)

    return hfinal
