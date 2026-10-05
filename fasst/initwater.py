"""
initwater.py -- traducao de initwater.F90

Fonte Fortran: initwater.F90 (203 linhas). Uma unica sub-rotina publica
(`initwater`). Inicializa o estado de neve/gelo acumulados no inicio de
cada execucao (ou a cada timestep, dependendo de `sprint`/`infer_test`):
zera variaveis de controle de vegetacao/estacao, decide se os dados vem
de uma execucao anterior (`infer_test==1`, cenario multi-run) ou sao
lidos direto do arquivo meteorologico (`infer_test==0`), inicializa a
profundidade e densidade de neve acumulada, escreve o cabecalho do
arquivo de saida de neve (unidade 55) quando solicitado, e varre o
arquivo meteorologico em busca de passos de tempo faltantes.

Assinatura original:
    subroutine initwater(sprint,phie,pdens,oldsd,oldhi,sdensi)
    integer(ip),intent(in):: sprint
    real(dp),intent(inout):: phie,pdens,oldsd,oldhi,sdensi

Arquitetura
-------------
`initwater` recebe `state: FasstState` como argumento explicito.
Constantes fisicas (`EPS`, `TREF`, `SPFLAG`, `SDENSD`, `SDENSW`) e
tamanhos fixos (`MAXCOL`, `MOVERLAP`, `NCLAYERS`) vem de
`fasst.constants` -- `sdensd`/`sdensw` nao sao campos de `FasstState`.
`moverlap` (usado aqui para `wstart` e para dimensionar `timeo`) tambem
e constante -- diferente de `initsurface.py`, que usa `moverlapr`
(campo real de estado) para `wstart` (ver "Nao confundir com
initsurface.py" abaixo). Indices de coluna do meteorologico (`ip_tmp`,
`ip_sd` no Fortran) sao `MetCol.TMP`/`MetCol.SD`. `mflag` e campo de
`FasstState`, lido do cabecalho do arquivo meteorologico em tempo de
execucao.

`met_date` e real e pura, importada de `fasst.functions`. `missing_met`
tambem e real, importada de `fasst.missing_met` (as sub-rotinas que
chama -- sol_zen, Solflx, emisatm, dnirflx -- ja estao resolvidas em
`fasst.module_radiation`). `read_old_data` continua sem implementacao
conhecida -- e um stub de modulo (levanta `NotImplementedError`),
duplicado neste arquivo e nao importado de `initsurface.py` (mesma
convencao de nao compartilhar infraestrutura entre arquivos deste
projeto, ver `_anint`/`_fdiv`).

Indexacao 1-based identica ao Fortran em todos os vetores/matrizes
(posicao 0 sem uso), com os mesmos indices (i, j, ...) do original.
Inclui `met1`, que no Fortran e declarada `met1(1,maxcol)` (matriz
1×maxcol) mas aqui vira um vetor 1D 1-based, ja que a primeira
dimensao e sempre 1 em todo o arquivo.

Nao confundir com initsurface.py
-----------------------------------
Este arquivo e initsurface.py se parecem muito, mas tem diferencas
numericas/de nomenclatura reais -- ver o docstring de initsurface.py
para a lista completa (a condicao extra `i>4`, `moverlap` vs
`moverlapr`, o fator de arredondamento de hsaccum/hm, o valor padrao de
pdens, veg_propl/veg_proph comentadas aqui mas reais la, o FORMAT 92
do cabecalho de neve).

Pontos de atencao
-------------------
* A linha de cabecalho do arquivo de saida de neve (unidade 55) e
  reconstruida a partir do FORMAT 92 do Fortran (contagem exata de
  espacos entre rotulos de coluna: 150 caracteres, com espacos antes
  de add/meta/dwind/atop/abot/sd) -- conferida byte a byte contra o
  FORMAT original.
"""

import math

from .constants import EPS, MAXCOL, MOVERLAP, NCLAYERS, SDENSD, SDENSW, SPFLAG, TREF
from .indices import MetCol
from .state import FasstState
from .functions import met_date
from .missing_met import missing_met


def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _aint(x, p=0):
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direção a zero.

    Duplicada por arquivo -- ver `module_radiation.py`/`module_lowveg.py`/
    `module_canopy.py` para a mesma convenção.
    """
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


# ---------------------------------------------------------------------------
# Dependencias externas ainda sem implementacao conhecida (ver docstring)
# ---------------------------------------------------------------------------
def read_old_data(d1i, mpos, wstart, sdensi, phie, timeo):
    """TODO: funcao externa -- recupera sdensi/phie/timeo de uma execucao anterior (multi-run)."""
    raise NotImplementedError("read_old_data ainda nao foi traduzida (funcao externa)")


def initwater(state: FasstState, sprint, phie, pdens, oldsd, oldhi, sdensi):
    # ---------------- STEP 1 (initwater.F90, linhas 31-51) ----------------
    # initialize certain check/test variables
    state.iseason = 0

    wcheck = 0
    wstart = 0
    wend = 0
    d1i = 0

    oldsd = 0.0
    oldhi = 0.0
    sdensi = 0.0
    state.zh = 0.0
    pdens = 0.0

    for i in range(1, NCLAYERS + 1):
        state.laif[i] = 0.0
        state.dzveg[i] = 0.0

    state.iheightn = state.iheight

    # ---------------- STEP 2 (initwater.F90, linhas 54-166) ----------------
    # initialize surface conditions and soil profile, max/min values, etc
    met1 = _new_array(MAXCOL)  # met1(1,maxcol) -- linha unica, 1 x maxcol

    if state.infer_test == 1:
        # ramo infer_test==1: dados de uma execucao anterior (multi-run)
        # (initwater.F90, linhas 56-105)
        for j in range(1, MAXCOL + 1):
            met1[j] = state.mflag  # oldpos = 1

        d1i = 1
        wstart = max(MOVERLAP, state.istart)
        timeo = [[0.0] * 5 for _ in range(MOVERLAP + 1)]  # timeo(moverlap,4)
        sdensi, phie = read_old_data(d1i, state.mpos, wstart, sdensi, phie, timeo)

        for i in range(1, MAXCOL + 1):
            met1[i] = state.met[state.istart - 1][i]

        d1i = 1
        missing_met(state, d1i, state.oldpos, state.oldpos, met1)

        # double check met parameters for missing values
        for j in range(state.istart, state.iend + 1):
            for i in range(1, MAXCOL - 5 + 1):
                if i > 4 and (i != 26 and i != 28):
                    if (_aint(abs(state.met[j][i] - state.mflag) * 1e5) * 1e-5 <= EPS
                            or _aint(abs(state.met[j][i] - _aint(state.mflag)) * 1e5) * 1e-5 <= EPS):
                        wcheck = 1
                        if wstart == 0:
                            wstart = j
                        wend = j

        if wcheck == 1:
            for j in range(1, MAXCOL + 1):
                met1[j] = state.mflag

            d1i = 2
            missing_met(state, d1i, wstart, wend, met1)

        state.sdens[state.oldpos] = sdensi

        state.ft[state.oldpos] = state.ftemp
        state.tt[state.oldpos] = state.toptemp
        state.airt[1][state.oldpos] = state.met[state.oldpos][MetCol.TMP] + TREF
        state.airt[2][state.oldpos] = state.dmet1[state.oldpos][4] + TREF

        # (bloco original comentado no Fortran -- mantido comentado)
        # if vegl_type != 0: veg_propl(...)
        # if vegh_type != 0: veg_proph(...)

    elif state.infer_test == 0:
        # ramo infer_test==0: dados lidos direto do arquivo meteorologico
        # (initwater.F90, linhas 106-166)
        # double check met parameters for missing values
        for j in range(state.istart, state.iend + 1):
            for i in range(1, MAXCOL - 5 + 1):
                if i != 26 and i != 28:
                    if (_aint(abs(state.met[j][i] - state.mflag) * 1e5) * 1e-5 <= EPS
                            or _aint(abs(state.met[j][i] - _aint(state.mflag)) * 1e5) * 1e-5 <= EPS):
                        wcheck = 1
                        if wstart == 0:
                            wstart = j
                        wend = j

        if wcheck == 1:
            for j in range(1, MAXCOL + 1):
                met1[j] = state.mflag

            d1i = 2
            missing_met(state, d1i, wstart, wend, met1)

        # initialize snow depth
        if (state.met[state.istart][MetCol.SD] > EPS
                and _aint(abs(state.met[state.istart][MetCol.SD] - state.mflag) * 1e5) * 1e-5 > EPS):
            if (abs(state.hsaccum) <= EPS
                    or _aint(abs(state.hsaccum - SPFLAG) * 1e5) * 1e-5 <= EPS):
                state.hsaccum = state.met[state.istart][MetCol.SD]
            elif (abs(state.hsaccum) >= EPS
                    and _aint(abs(state.hsaccum - SPFLAG) * 1e5) * 1e-5 > EPS):
                state.hsaccum = max(state.hsaccum, state.met[state.istart][MetCol.SD])

        state.hsaccum = _anint(state.hsaccum, 5)
        state.hm = _anint((state.hsaccum + state.newsd + state.hi), 5)  # m

        # initialize soil profile
        # (bloco original comentado no Fortran -- mantido comentado)
        # if veg_flagl == 1: veg_propl(...)
        # if veg_flagh == 1: veg_proph(...)

        if state.hsaccum > EPS or state.hi > EPS:
            state.toptemp = TREF
        if state.hsaccum > EPS:
            if (abs(state.iswe) >= EPS
                    and _aint(abs(state.iswe - SPFLAG) * 1e5) * 1e-5 > EPS):
                sdensi = (state.iswe / state.hsaccum) * 1e3
            else:
                sdensi = 0.5 * (SDENSD + SDENSW)
        else:
            sdensi = 0.0
        state.sdens[state.istart] = sdensi
        state.storll = 0.0
        state.storls = 0.0
    # end if(infer_test == 1)

    # more snow depth initialization (initwater.F90, linhas 168-172)
    oldsd = state.hsaccum
    oldhi = state.hi
    if abs(pdens) <= EPS:
        pdens = 500.0
    if oldsd + oldhi > EPS:
        phie = (1.0 - pdens * 1e-3) * 0.95

    # write the snow output file header line (initwater.F90, linhas 174-184)
    if sprint == 1:
        state.unit55.write(
            "{:10d} {:6d} {:3d} {:10.6f} {:11.6f} {:11.6f} {:8d} {:8d} "
            "{:6.2f} {:5.2f} {:8.2f}\n".format(
                state.freq_id, state.iend, 9, state.lat, state.mlong, state.elev, state.vitd_index,
                state.met_count, state.timeoffset, state.timstep, state.mflag,
            )
        )
        state.unit55.write(
            # Reconstruida diretamente do FORMAT 92 do Fortran:
            # 6x,'doy',8x,'sdold',10x,'add',12x,'meta',11x,'dwind',10x,'atop',
            # 11x,'abot',11x,'sd',13x,'delta sd',1x,'mode',7x,'diameter'
            # Total: 150 caracteres.
            "      doy        sdold          add            meta           dwind"
            "          atop           abot           sd             delta sd mode"
            "       diameter\n"
        )

    # check met file for missing time steps (initwater.F90, linhas 186-201)
    for i in range(state.istart, state.iend + 1):
        if i != 1:
            time1 = met_date(state.met[i - 1][1], state.met[i - 1][2],
                              state.met[i - 1][3], state.met[i - 1][4])
            time2 = met_date(state.met[i][1], state.met[i][2],
                              state.met[i][3], state.met[i][4])
            if abs(state.timstep - (time2 - time1) * 365.0 * 24.0) > 1e-1:
                for j in range(1, MAXCOL + 1):
                    met1[j] = state.met[i - 1][j]  # oldpos = 1

                d1i = 3
                j = i - 1
                missing_met(state, d1i, j, i, met1)

    return phie, pdens, oldsd, oldhi, sdensi
