"""
initsurface.py -- traducao de initsurface.F90

Fonte Fortran: initsurface.F90 (326 linhas). Uma unica sub-rotina
publica (`initsurface`). Orquestra a inicializacao de uma execucao do
FASST: zera variaveis de controle de vegetacao/estacao, decide se os
dados de neve/solo vem de uma execucao anterior (infer_test==1, cenario
multi-run) ou sao lidos/calculados do zero (infer_test==0 -- inclusive
chamando initprofile para montar o perfil de solo), inicializa a
rugosidade da superficie, a profundidade/densidade de neve acumulada, e
escreve os cabecalhos dos arquivos de saida de nos (unidade 3) e de
neve (unidade 55) quando solicitado.

E a irma mais completa de initwater.py: cobre praticamente o mesmo
terreno (mesma estrutura STEP 1 / STEP 2, mesmas checagens de dado
meteorologico faltante, mesma inicializacao de neve), mas soma a isso a
leitura de propriedades de solo por camada e a chamada a initprofile.
Varios detalhes numericos diferem dos dois arquivos mesmo em trechos
aparentemente identicos -- ver "Nao confundir com initwater.py" abaixo.

Assinatura original:
    subroutine initsurface(nt0,nm0,nprint,sprint,sname,ttest,mtest,
                            nstr_flag,code1,phie,pdens,oldsd,oldhi,sdensi)

Arquitetura
-------------
`initsurface` recebe `state: FasstState` como argumento explicito.
Constantes fisicas (`EPS`, `TREF`, `SPFLAG`, `SDENSW`) e tamanhos fixos
(`MAXL`, `MAXN`, `MAXP`, `MAXCOL`, `NCLAYERS`) vem de `fasst.constants`
-- nao sao estado. `MOVERLAP` (sem "r", usado so para dimensionar
`timeo`) tambem e constante; `moverlapr` (usado no calculo de `wstart`)
E campo real de `FasstState` -- os dois continuam distintos aqui,
exatamente como no Fortran (ver "Nao confundir com initwater.py").
`ip_tmp`/`ip_sd` (indices de coluna do meteorologico, no Fortran) sao
`MetCol.TMP`/`MetCol.SD`. `mflag` e campo real de `FasstState`, lido do
cabecalho do arquivo meteorologico em tempo de execucao.

`veg_propl`/`veg_proph` sao as funcoes reais, importadas de
`fasst.module_lowveg`/`fasst.module_canopy`. `veg_propl` espera as
linhas do arquivo de propriedades ja lidas (`file_unit_31_data`),
enquanto `veg_proph` usa um stub interno de leitura linha a linha
(`_read_veg_table_row`, em module_canopy.py) -- as duas dependencias de
leitura de arquivo (units 31/32/33) nao seguem a mesma convencao entre
si. Aqui, isso fica exposto via um stub local (`_read_table_lines`),
usado na chamada de `veg_propl` -- ver "Dependencias externas".

`read_old_data` nao tem implementacao conhecida -- e um stub de modulo
(levanta `NotImplementedError`), no mesmo padrao de `_read_veg_table_row`
em `module_canopy.py`. `met_date`, `missing_met`, `upr_case` e
`get_user_soil_params` sao reais, importadas de `fasst.functions`/
`fasst.missing_met`/`fasst.us_soil_tools`. `get_user_soil_params` nao usa
`state`: recebe o arquivo da unidade 90 (`state.unit90`) e o nome do solo,
e devolve os 23 valores de saida, atribuidos aqui a `sname[i]`,
`state.sclass[i]`, `state.soiltype[i]` e `state.isoilp[i][...]`.

`state.unit3`, `state.unit10`, `state.unit55`, `state.unit90`,
`state.unit31`, `state.unit32`, `state.unit33` sao usados aqui como
identificadores de arquivo ja abertos -- essa convencao nao esta
confirmada contra `fasst/state.py`, que nao declara campos de
arquivo/unidade logica. Ver Pendencias.

Estrutura desta transcricao (parametros)
------------------------------------------
Parametros Fortran e seu destino em Python:
  nt0, nm0             intent(in)     -> ver nota abaixo: initprofile os
                                         trata como inout (interface
                                         implicita, sem checagem de
                                         intent pelo compilador), entao
                                         voltam alterados via aliasing
                                         de memoria -- devolvidos na
                                         tupla de retorno por fidelidade
  nprint, sprint       intent(in)     -> parametros normais de entrada
  sname(maxl)          intent(in)     -> lista de strings (1-indexada)
  ttest                intent(in)     -> parametro normal de entrada
                                         (nunca modificado por initprofile)
  mtest                intent(in)     -> mesma situacao de nt0/nm0: alterado
                                         dentro de initprofile via aliasing,
                                         volta na tupla de retorno
  nstr_flag(maxn)       intent(out)    -> devolvido na tupla de retorno
  code1                 intent(inout)  -> devolvido na tupla de retorno
  phie,pdens,oldsd,
  oldhi,sdensi          intent(inout)  -> devolvido na tupla de retorno

Nota sobre nt0/nm0/mtest: no Fortran, initsurface declara os tres como
intent(in) -- mas os passa para initprofile (um arquivo externo, sem
interface explicita), cuja propria assinatura os trata como
intent(inout). Fortran passa argumentos por referencia independente do
intent, e sem uma INTERFACE explicita o compilador nao valida essa
inconsistencia entre as duas assinaturas -- entao os valores QUE
initprofile escreve realmente propagam de volta para initsurface (e
para quem chamou initsurface), apesar do intent(in) declarado aqui.
Python nao tem esse aliasing por acidente, entao a unica forma fiel de
reproduzir o comportamento e devolver os tres na tupla de retorno de
initsurface, mesmo eles sendo "so entrada" na assinatura Fortran local.

Dependencias externas
-----------------------
Nenhuma implementada por suposicao: `read_old_data` e um stub local
(levanta `NotImplementedError`); `met_date`, `missing_met`, `upr_case` e
`get_user_soil_params` sao reais (`fasst.functions`/`fasst.missing_met`/
`fasst.us_soil_tools`); `veg_propl`/`veg_proph` sao as funcoes reais (em
`module_lowveg.py`/`module_canopy.py`).

Nao confundir com initwater.py
--------------------------------
Varios trechos deste arquivo se parecem MUITO com initwater.py, mas tem
diferencas numericas ou de condicao que a traducao preserva
corretamente -- nao usar um arquivo como referencia visual para "corrigir"
o outro sem checar o .F90 de cada um:
* A condicao extra `i>4` na checagem de met faltante esta no ramo
  OPOSTO em cada arquivo: aqui (initsurface) e o ramo infer_test==0 que
  tem `i>4`; em initwater e o ramo infer_test==1.
* wstart usa `moverlapr` aqui (initwater usa `moverlap`, sem "r"); o
  array timeo, no mesmo trecho, usa `moverlap` (sem "r") nos dois arquivos.
* hsaccum/hm sao arredondados com fator 1e15 aqui (initwater usa 1e5).
* O valor padrao de pdens e SDENSW aqui (initwater usa a constante 500).
* veg_propl/veg_proph sao chamadas reais aqui, nos dois ramos; em
  initwater.py as mesmas chamadas aparecem comentadas (codigo morto no
  proprio Fortran original).
* O FORMAT 92 (cabecalho do arquivo de neve, unidade 55) tem campos
  extras aqui (year/swe/%wet) e espacamento proprio -- ver "Pontos de
  atencao".

Pontos de atencao
-------------------
* `code1 = 1`, `error_code = 1` e a mensagem "Used default soil type" (so
  quando `single_multi_flag == 0`) sao executados para TODA camada que nao
  e solo definido pelo usuario (`soiltype != -1`), inclusive quando o tipo
  de solo informado e valido: no Fortran (linhas 179-183), essas tres
  instrucoes ficam depois da cadeia `if/else if` que substitui os tipos
  0, -2, -3 e -4, e nao dentro dela. A mensagem sugere que a intencao era
  registrar apenas a substituicao por um tipo padrao. Preservado como no
  original.
* A linha de cabecalho do arquivo de saida de neve (unidade 55) e
  reconstruida a partir do FORMAT 92 DESTE arquivo (nao do FORMAT 92 de
  initwater.F90, que e diferente): 192 caracteres, contagem exata de
  espacos entre rotulos de coluna -- conferida byte a byte contra a
  reconstrucao do FORMAT original.
"""

from .constants import EPS, MAXCOL, MAXL, MAXN, MAXP, MOVERLAP, NCLAYERS, SDENSD, SDENSW, SPFLAG, TREF
from .indices import MetCol
from .state import FasstState
from .functions import met_date
from .missing_met import missing_met
from .us_soil_tools import get_user_soil_params, upr_case
from .module_lowveg import veg_propl
from .module_canopy import veg_proph
from .initprofile import initprofile


def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _aint(x, p=0):
    import math
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _anint(x, p=0):
    import math
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


def _read_table_lines(unit):
    """
    TODO: leitura do arquivo de propriedades de vegetacao (units 31/32/33),
    ja aberto em algum lugar nao identificado (fasst_main.F90?). Deveria
    devolver a lista de linhas do arquivo (para `veg_propl`) -- ver nota
    de inconsistencia com `veg_proph`/`_read_veg_table_row`, em
    "Reconciliação com a arquitetura real" no docstring do módulo.
    """
    raise NotImplementedError(
        f"leitura do arquivo de propriedades de vegetacao (unit {unit}) "
        "ainda nao foi traduzida"
    )


def initsurface(state: FasstState, nt0, nm0, nprint, sprint, sname, ttest, mtest,
                 code1, phie, pdens, oldsd, oldhi, sdensi):

    # ---------------- STEP 1 (initsurface.F90, linhas 45-81) ----------------
    state.iseason = 0
    code1 = 0
    tstps = 0

    wcheck = 0
    wstart = 0
    wend = 0
    d1i = 0

    str_flag = _new_array(MAXL, 0)
    nstr_flag = _new_array(MAXN, 0)

    oldsd = 0.0
    oldhi = 0.0
    sdensi = 0.0
    state.zh = 0.0
    pdens = 0.0

    for i in range(1, NCLAYERS + 1):
        state.laif[i] = 0.0
        state.dzveg[i] = 0.0

    state.iheightn = state.iheight
    rtstps = 24.0 / state.timstep
    tstps = int(rtstps)  # INT trunca em direcao a zero == int() do Python

    # ---------------- STEP 2 (initsurface.F90, linhas 84-247) ----------------
    met1 = _new_array(MAXCOL)

    if state.infer_test == 1:
        # ramo infer_test==1: dados de uma execucao anterior (multi-run)
        # (initsurface.F90, linhas 86-135)
        for j in range(1, MAXCOL + 1):
            met1[j] = state.mflag  # oldpos = 1

        d1i = 1
        wstart = max(state.moverlapr, state.istart)
        timeo = [[0.0] * 5 for _ in range(MOVERLAP + 1)]
        sdensi, phie = read_old_data(d1i, state.mpos, wstart, sdensi, phie, timeo)

        for i in range(1, MAXCOL + 1):
            met1[i] = state.met[state.istart - 1][i]

        d1i = 1
        missing_met(state, d1i, state.oldpos, state.oldpos, met1)

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

        state.sdens[state.oldpos] = sdensi
        state.ft[state.oldpos] = state.ftemp
        state.tt[state.oldpos] = state.toptemp
        state.airt[1][state.oldpos] = state.met[state.oldpos][MetCol.TMP] + TREF
        state.airt[2][state.oldpos] = state.dmet1[state.oldpos][4] + TREF

        if state.vegl_type != 0:
            veg_propl(state, state.biome_source, state.new_vtl, state.vegl_type,
                      _read_table_lines(31))
        if state.vegh_type != 0:
            veg_proph(state, state.biome_source, state.new_vth, state.vegh_type)

    elif state.infer_test == 0:
        # ramo infer_test==0: propriedades de solo lidas/calculadas do
        # zero, perfil montado via initprofile (initsurface.F90, linhas
        # 136-247)
        for i in range(1, state.nlayers + 1):
            for j in range(1, MAXP + 1):
                state.isoilp[i][j] = SPFLAG
                state.soilp[i][j] = SPFLAG

        # atribui/valida tipo de solo por camada, ou aceita parametros
        # informados pelo usuario (initsurface.F90, linhas 144-185)
        for i in range(1, state.nlayers + 1):
            str_flag[i] = 0
            if state.soiltype[i] == -1:  # user supplied soil parameters
                str_flag[i] = 3
                d1i = len(sname[i].rstrip())
                sname[i] = upr_case(d1i, sname[i])

                (sname[i], state.sclass[i], state.soiltype[i],
                 state.isoilp[i][1], state.isoilp[i][2], state.isoilp[i][3],
                 state.isoilp[i][4], state.isoilp[i][5], state.isoilp[i][6],
                 state.isoilp[i][7], state.isoilp[i][8], state.isoilp[i][9],
                 state.isoilp[i][10], state.isoilp[i][11], state.isoilp[i][13],
                 state.isoilp[i][12], state.isoilp[i][14], state.isoilp[i][18],
                 state.isoilp[i][19], state.isoilp[i][20], state.isoilp[i][21],
                 state.isoilp[i][22], state.isoilp[i][23]) = get_user_soil_params(
                    state.unit90, sname[i])
                state.unit90.seek(0)
            else:
                if state.soiltype[i] == 0:  # unknown soil type
                    state.soiltype[i] = 7  # default: dirty sand
                    str_flag[i] = 1
                elif state.soiltype[i] == -2:  # unknown fine grained
                    state.soiltype[i] = 12
                    str_flag[i] = 2
                elif state.soiltype[i] == -3:  # unknown coarse grained
                    state.soiltype[i] = 6
                    str_flag[i] = 1
                elif state.soiltype[i] == -4:  # unknown disturbed
                    state.soiltype[i] = 7
                    str_flag[i] = 0
                    state.rho_fac[i] = 0.8

                code1 = 1
                state.error_code = 1
                if state.single_multi_flag == 0:
                    state.unit10.write(" Used default soil type. \n\n")

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

        # initialize snow depth
        if (state.met[state.istart][MetCol.SD] > EPS
                and _aint(abs(state.met[state.istart][MetCol.SD] - state.mflag) * 1e5) * 1e-5 > EPS):
            if (abs(state.hsaccum) <= EPS
                    or _aint(abs(state.hsaccum - SPFLAG) * 1e5) * 1e-5 <= EPS):
                state.hsaccum = state.met[state.istart][MetCol.SD]
            elif (abs(state.hsaccum) >= EPS
                    and _aint(abs(state.hsaccum - SPFLAG) * 1e5) * 1e-5 > EPS):
                state.hsaccum = max(state.hsaccum, state.met[state.istart][MetCol.SD])

        state.hsaccum = _anint(state.hsaccum, 15)
        state.hm = _anint((state.hsaccum + state.newsd + state.hi), 15)  # m

        # initialize soil profile
        (nt0, nm0, str_flag, mtest, nstr_flag) = initprofile(
            state, ttest, nt0, mtest, nm0, str_flag
        )

        if state.veg_flagl == 1:
            veg_propl(state, state.biome_source, state.new_vtl, state.vegl_type,
                      _read_table_lines(31))
        if state.veg_flagh == 1:
            veg_proph(state, state.biome_source, state.new_vth, state.vegh_type)

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

    # initialize soil roughness length (initsurface.F90, linhas 249-261)
    if state.rough <= EPS:
        if state.ntype[state.nnodes] >= 19:  # concrete, asphalt, bedrock, glaciers
            state.rough = 1e-3
        elif state.ntype[state.nnodes] <= 4 or state.ntype[state.nnodes] == 15:  # gravels, peat
            state.rough = 5e-2
        elif 4 < state.ntype[state.nnodes] <= 6:  # rough sands
            state.rough = 1e-2
        else:  # all else
            state.rough = 1e-3
    state.rough = _anint(state.rough, 15)

    # initialize snow/ice (initsurface.F90, linhas 263-267)
    oldsd = state.hsaccum
    oldhi = state.hi
    if abs(pdens) <= EPS:
        pdens = SDENSW
    if oldsd + oldhi > EPS:
        phie = (1.0 - pdens * 1e-3) * 0.95

    # write the header lines for the node input file (initsurface.F90,
    # linhas 269-276)
    if state.single_multi_flag == 0 and nprint == 1:
        state.unit3.write(
            " Total Number of Nodes: {:4d},  freq_id:{:10d}\n".format(
                state.nnodes, state.freq_id
            )
        )
        state.unit3.write(
            " Year   JD   Hr    M    N  USCS    Depth   Grtemp   "
            "Water         Ice       W + I       Vapor     F/T\n"
        )
        state.unit3.write("                                     m       K\n")

    # write the snow output file header line (initsurface.F90, linhas
    # 278-288 -- FORMAT 92 proprio deste arquivo, ver docstring do modulo)
    if sprint == 1:
        state.unit55.write(
            "{:10d} {:6d} {:3d} {:10.6f} {:11.6f} {:11.6f} {:8d} {:8d} "
            "{:6.2f} {:5.2f} {:8.2f}\n".format(
                state.freq_id, state.iend, 9, state.lat, state.mlong, state.elev, state.vitd_index,
                state.met_count, state.timeoffset, state.timstep, state.mflag,
            )
        )
        state.unit55.write(
            # Reconstruida diretamente do FORMAT 92 do Fortran (initsurface.F90):
            # 6x,'year',8x,'doy',8x,'sdold',10x,'add',12x,'meta',11x,'dwind',
            # 10x,'atop',11x,'abot',11x,'sd',11x,'delta sd',3x,'mode',
            # 7x,'diameter',12x,'swe',11x,'%wet'
            # Total: 192 caracteres. Nao e o mesmo cabecalho de initwater.py
            # (outro FORMAT, sem os campos 'year'/'swe'/'%wet', 150 caracteres).
            "      year        doy        sdold          add            meta"
            "           dwind          atop           abot           sd     "
            "      delta sd   mode       diameter            swe           "
            "%wet\n"
        )

    # check met file for missing time steps (initsurface.F90, linhas
    # 290-324)
    for i in range(state.istart, state.iend + 1):
        if i != 1:
            time1 = met_date(state.met[i - 1][1], state.met[i - 1][2],
                              state.met[i - 1][3], state.met[i - 1][4])
            time2 = met_date(state.met[i][1], state.met[i][2],
                              state.met[i][3], state.met[i][4])

            mody = state.met[i][1] - _aint(state.met[i][1] * 0.25) * 4.0
            # Fortran tambem calcula `modyo` aqui (mesma formula, para
            # met(i-1,1)), mas ela nunca e lida depois -- a unica linha
            # que a usaria (`if(mody /= modyo) daylim = 1.07d0`) esta
            # comentada no proprio .F90. Omitida: codigo morto, sem efeito
            # observavel.
            daylim = 1.042

            if state.met[i - 1][1] != state.met[i][1]:
                daylim = 25.07

            if abs(mody) <= EPS:
                if (time2 - time1) * 367.0 * tstps >= daylim:
                    for j in range(1, MAXCOL + 1):
                        met1[j] = state.met[i - 1][j]
                    d1i = 3
                    j = i - 1
                    missing_met(state, d1i, j, i, met1)
            else:
                if (time2 - time1) * 366.0 * tstps >= daylim:
                    for j in range(1, MAXCOL + 1):
                        met1[j] = state.met[i - 1][j]
                    d1i = 3
                    j = i - 1
                    missing_met(state, d1i, j, i, met1)

    return nstr_flag, code1, phie, pdens, oldsd, oldhi, sdensi, mtest, nt0, nm0
