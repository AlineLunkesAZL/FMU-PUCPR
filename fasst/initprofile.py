"""
initprofile.py -- traducao de initprofile.F90 do FASST

Fonte Fortran: initprofile.F90 (1156 linhas). Uma unica sub-rotina publica
(`initprofile`) mais duas sub-rotinas internas apendicadas ao mesmo
arquivo (`sort`, `locate`, ao final deste arquivo -- ver nota abaixo).

O que faz
---------
Posiciona os nos do perfil vertical (solo + camadas de neve/vegetacao,
se houver), a partir da espessura de cada camada e, opcionalmente, de
profundidades de medicao de temperatura/umidade (`zti`/`zm`); le e ajusta
as propriedades de cada camada de solo a partir de uma tabela em disco
(unidade 30); calcula albedo/emissividade de neve; inicializa
temperatura, cabeca de pressao (`phead`), umidade (`soil_moist`), teor de
vapor d'agua/gelo (`wvc`/`ice`) em cada no; e por fim calcula o passo de
tempo estavel (`deltat_fasst`/`step`) usado pelo resto do FASST na
integracao temporal.

Arquitetura
-------------
`initprofile` recebe `state: FasstState` como argumento explicito.
Constantes fisicas/tamanhos fixos (`EPS`, `TREF`, `SPFLAG`, `PI`, `GRAV`,
`RV`, `RD`, `MAXL`, `MAXN`, `MAXP`, `EXTRAN`, `HT_MIN`, `HT_MINM`) vem de
`fasst.constants`; indices de coluna do meteorologico (`ap`, `prec`,
`prec2`, `pt`, `pt2`, `rh`, `tmp`, `tsoil` no Fortran) sao
`MetCol.AP/PREC/PREC2/PT/PT2/RH/TMP/TSOIL`. `mflag` e campo de
`FasstState`, lido do cabecalho do arquivo meteorologico em tempo de
execucao.

`dense`, `head`, `soilhumid` e `vap_press` sao chamadas reais, importadas
de `fasst.functions` (`dense(temp1,wind1,phase)` -- pura --,
`head(state,i,smt)`, `soilhumid(state,i,ph,sms,st)`,
`vap_press(state,i,rh,ap)`). `th_param` tambem e real, importada de
`fasst.th_param` (nao e a mesma coisa que `thconds`, de `fasst.functions`,
que tem assinatura de 2 argumentos, contra os 6 de `th_param`). Como
`th_param` le `nnodes`/`ntemp` de `state` (variaveis globais no Fortran),
`initprofile` os grava em `state` imediatamente antes dessa chamada, alem
de novamente no final da funcao.

`nt0`, `nm0`, `str_flag`, `mtest` e `nstr_flag` nao pertencem a
`fasst_global`/`FasstState` -- continuam sendo passados e devolvidos por
valor, via argumentos e retorno em tupla. `sid` e `totthick` sao
variaveis puramente locais desta sub-rotina, tambem sem campo
correspondente em `FasstState`.

Indexacao 1-based com folga no indice 0, igual ao restante do projeto:
todo array cujo indice Fortran e usado literalmente no codigo (nz(i),
soilp(i,j), zti(j)...) e alocado aqui com uma posicao a mais em cada
dimensao 1-based, indice 0 sempre sem uso, indices 1..N correspondem
exatamente aos do Fortran.

`sort` e `locate`, ao final deste arquivo, sao duas sub-rotinas internas
do mesmo arquivo Fortran (apendice de initprofile.F90, logo depois de
`end subroutine initprofile`) -- mantidas aqui, sem modulo separado,
porque e assim que existem no original.

Pontos de atencao herdados do Fortran original
-----------------------------------------------
* `sort`: o teste `if(i==1.and.arr(i)>a) i=0`, dentro do insertion sort,
  nunca dispara de fato -- se i==1 e arr(i)>a, o `while` anterior ainda
  estaria rodando. Preservado por fidelidade (ver docstring de `sort`).
* Bloco de estabilidade do passo de tempo (calculo de `mindeltat`): a
  variavel local `deltatt` NAO e reiniciada a cada iteracao do laço
  `do i=1,nnodes` -- se `dabs(grthcond(i))<=eps`, o valor calculado na
  iteracao anterior e reaproveitado silenciosamente. Igual ao Fortran
  (variavel local ao escopo da sub-rotina, nao do laço); preservado.
* `mindeltatm` (criterio de estabilidade da umidade, 0.45*dz^2*theta_max/
  (K_sat*|h(theta_max)|)) e calculado em todos os nos de solo, mas nunca
  entra em `deltat_fasst` nem em outra expressao: `deltat_fasst` vem so
  do criterio termico (`mindeltat`). Igual ao Fortran (linhas 1055-1061
  do .F90), onde `mindeltatm` tambem so e atribuida.
* `locate`, no ramo `ttest=='y'`: a interpolacao de `stt(i)` a partir de
  `zti`/`tm` e feita sem nenhuma checagem de que `jj` esteja estritamente
  dentro de `[1, nt0-1]` -- se `locate` devolver `jj==nt0` (`z` fora do
  intervalo coberto por `zti`), o Fortran usa `zti(nt0+1)`/`tm(nt0+1)`
  mesmo assim, sem checagem de limite. Preservado por fidelidade.
* Leitura de propriedades de solo: quando `soilp(i,8)` (umidade minima) e
  `soilp(i,9)` (umidade maxima) da camada ficam iguais dentro de `eps`, o
  Fortran interrompe a execucao (`write` + `stop`, incondicional -- nao
  depende de `single_multi_flag`) -- levanta `SystemExit` aqui, igual ao
  segundo `stop` desta mesma sub-rotina (umidade inicial fora da faixa).
"""

import math

from .constants import (
    EPS, EXTRAN, GRAV, HT_MIN, HT_MINM, IALBEDO, IEMIS, MAXL, MAXN, MAXP, PI, RD, RV,
    SEMIS, SNALBEDO, SOALBEDO, SPFLAG, TREF,
)
from .indices import MetCol
from .state import FasstState
from .functions import Phase, dense, head, soilhumid, vap_press
from .th_param import th_param

# ------------------------------------------------------------------
# DATA ics / ... /   (classificacao USCS -> tipo grosso(1)/fino(2)/outro(0))
# ------------------------------------------------------------------
ICS = [0,  # posicao 0 nao usada (placeholder de indexacao 1-based)
       2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

MINDEPTH = 2.0          # minimum allowed soil depth
WLTPOINT = -1.5e4       # wilting point (cm)
FLDCAP = -3.4e2         # field capacity (cm)
MINWAT = -3.1e4         # minimum allowed water (cm)
AIRDRY = -1.0e6         # minimum head due to air drying (cm)


def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _aint(x, p=0):
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direção a zero.

    Duplicada por arquivo -- ver `module_radiation.py`/`module_lowveg.py`/
    `module_canopy.py`/`initsurface.py`/`initwater.py` para a mesma
    convenção (este projeto não tem um `fortran_compat` compartilhado).
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


def _read_soil_record(state):
    """Le um registro da tabela de propriedades de solo (unidade 30)."""
    line = state.unit30.readline()
    if not line:
        return None
    parts = line.split()
    if len(parts) < 16:
        return None
    sid = int(float(parts[0]))
    ssname = parts[1]
    vals = [float(x) for x in parts[2:16]]
    (dens, pors, ssemis, ssalb, shc, smin, smax, salpha, svgn,
     sspheat, sorgan, spsand, spsilt, spclay) = vals
    return (sid, ssname, dens, pors, ssemis, ssalb, shc, smin, smax,
            salpha, svgn, sspheat, sorgan, spsand, spsilt, spclay)


def initprofile(
    state: FasstState,
    ttest='n',
    nt0=0,
    mtest='n',
    nm0=0,
    str_flag=None,
    *args,
    **kwargs
):
    """
    Traducao fiel de `subroutine initprofile(ttest,nt0,mtest,nm0,str_flag,
    nstr_flag)` (initprofile.F90, linhas 1-1076).

    Posiciona os nos do perfil, le/ajusta as propriedades de solo por
    camada, inicializa temperatura/umidade/vapor-dagua-e-gelo em cada no,
    e calcula o passo de tempo estavel para integracao (deltat_fasst,
    step) -- ver docstring do modulo para a visao geral e os pontos de
    atencao herdados do Fortran original.
    """
    if str_flag is None:
        str_flag = [0] * (MAXL + 1)

    # ---- zero out local variables (arrays) (initprofile.F90, linhas 117-131) ----
    nz1 = _new_array(MAXN)
    nz3 = _new_array(MAXN)
    nz4 = _new_array(MAXN)
    rhoda = _new_array(MAXN)
    rhov = _new_array(MAXN)
    nzt = _new_array(MAXN)
    nz2 = _new_array(MAXL)
    sph = _new_array(MAXL)
    nstr_flag = _new_array(MAXN, 0)

    nmtest = nm0 + 1
    nttest = nt0 + 1
    state.toptemp = 0.0

    # =================================================================
    # initialize the layer properties (le tabela de propriedades, unit 30)
    # (initprofile.F90, linhas 140-320)
    # =================================================================
    for i in range(1, state.nlayers + 1):
        io_eof = False

        # pula as 20 linhas de cabecalho
        for _k in range(20):
            state.unit30.readline()

        jj = 0

        while (not io_eof) and jj == 0:
            rec = _read_soil_record(state)
            if rec is None:
                io_eof = True
                continue
            (sid, ssname, dens, pors, ssemis, ssalb, shc, smin, smax,
             salpha, svgn, sspheat, sorgan, spsand, spsilt, spclay) = rec

            if state.soiltype[i] == sid:
                jj = 1
                state.soiltype[i] = sid
                state.stype[i] = ssname
                state.soilp[i][1] = dens
                state.soilp[i][2] = pors
                state.soilp[i][3] = ssalb
                state.soilp[i][4] = ssemis
                state.soilp[i][7] = shc
                state.soilp[i][8] = smin
                state.soilp[i][9] = smax
                state.soilp[i][11] = svgn
                state.soilp[i][13] = sspheat
                state.soilp[i][14] = sorgan
                state.soilp[i][18] = spsand
                state.soilp[i][19] = spsilt
                state.soilp[i][20] = spclay

                if state.soilp[i][2] <= EPS:
                    state.soilp[i][2] = 1e-3

                if sid <= 2 or sid in (5, 6):
                    state.soilp[i][23] = 2.5
                elif sid in (3, 4):
                    state.soilp[i][23] = 30.0
                else:
                    state.soilp[i][23] = 100.0 - state.soilp[i][18]

                if abs(state.rho_fac[i] - 1.0) > EPS:
                    state.soilp[i][1] = state.soilp[i][1] * (
                        (1.0 - state.soilp[i][2] / state.rho_fac[i])
                        / (1.0 - state.soilp[i][2])
                    )
                    state.soilp[i][2] = state.soilp[i][2] / state.rho_fac[i]
                    state.soilp[i][7] = state.soilp[i][7] / state.rho_fac[i]
                    state.soilp[i][9] = state.soilp[i][2]
                    state.soilp[i][11] = state.soilp[i][11] / state.rho_fac[i]
                    salpha = salpha * state.rho_fac[i] * 0.5

                if abs(state.soilp[i][5] - SPFLAG) <= EPS:
                    state.soilp[i][5] = 1e-2 * state.soilp[i][18]
                if abs(salpha - SPFLAG) > EPS:
                    state.soilp[i][10] = 1.0 / salpha
                if abs(state.soilp[i][11] - SPFLAG) > EPS:
                    state.soilp[i][12] = 1.0 - 1.0 / state.soilp[i][11]
                if abs(state.soilp[i][18] - SPFLAG) > EPS:
                    state.soilp[i][18] = 1e-2 * state.soilp[i][18]
                if abs(state.soilp[i][19] - SPFLAG) > EPS:
                    state.soilp[i][19] = 1e-2 * state.soilp[i][19]
                if abs(state.soilp[i][20] - SPFLAG) > EPS:
                    state.soilp[i][20] = 1e-2 * state.soilp[i][20]
                if abs(state.soilp[i][21] - SPFLAG) > EPS:
                    state.soilp[i][21] = 0.5 * state.soilp[i][14]
                if abs(state.soilp[i][23] - SPFLAG) > EPS:
                    state.soilp[i][23] = 1e-2 * state.soilp[i][23]

        state.unit30.seek(0)  # rewind(30)

        if state.soiltype[i] >= 19:
            str_flag[i] = 4

        if state.soiltype[i] != 26:
            if state.soilp[i][9] > state.soilp[i][2]:
                state.soilp[i][9] = state.soilp[i][2]
        else:
            if state.soilp[i][9] > state.soilp[i][1]:
                state.soilp[i][9] = state.soilp[i][1]

        for j in range(1, MAXP + 1):
            f1 = _anint(state.isoilp[i][j], 10)
            if abs(f1 - SPFLAG) > EPS:
                state.soilp[i][j] = state.isoilp[i][j]

        if state.soilp[i][2] <= EPS:
            state.soilp[i][2] = 1e-3

        if (abs(state.isoilp[i][1] - SPFLAG) <= EPS
                and abs(state.isoilp[i][2] - SPFLAG) > EPS):
            state.soilp[i][1] = state.soilp[i][2] * 2.7

        state.soilp[i][25] = max(0.0, (1.0 - state.soilp[i][23]) - state.soilp[i][18])

        if ((state.soiltype[i] != 26 and state.soiltype[i] != 27)
                and (state.soilp[i][10] > EPS
                     and abs(state.soilp[i][10] - SPFLAG) > EPS)):
            c1 = 1.0 + (state.soilp[i][10] * abs(AIRDRY)) ** state.soilp[i][11]
            w1 = 1.0 / (c1 ** state.soilp[i][12])
            temp = state.soilp[i][8] + w1 * (state.soilp[i][9] - state.soilp[i][8])
            state.soilp[i][8] = max(1e-3, min(temp, state.soilp[i][8]))
            if abs(state.soilp[i][8] - state.soilp[i][9]) <= EPS:
                state.soilp[i][8] = 0.1 * state.soilp[i][9]

            c1 = 1.0 + (state.soilp[i][10] * abs(MINWAT)) ** state.soilp[i][11]
            w1 = 1.0 / (c1 ** state.soilp[i][12])
            state.soilp[i][15] = state.soilp[i][8] + w1 * (state.soilp[i][9] - state.soilp[i][8])
            state.soilp[i][15] = max(state.soilp[i][15], state.soilp[i][8] * 1.001)

            w = (state.soilp[i][15] - state.soilp[i][8]) / (state.soilp[i][9] - state.soilp[i][8])
            sph[i] = -(1.0 / state.soilp[i][10]) * (
                ((w ** (-1.0 / state.soilp[i][12])) - 1.0) ** (1.0 / state.soilp[i][11])
            )
            sph[i] = 1e-2 * sph[i]
            sph[i] = _anint(sph[i], 10)

            c1 = 1.0 + (state.soilp[i][10] * abs(WLTPOINT)) ** state.soilp[i][11]
            w1 = 1.0 / (c1 ** state.soilp[i][12])
            state.soilp[i][16] = state.soilp[i][8] + w1 * (state.soilp[i][9] - state.soilp[i][8])
            state.soilp[i][16] = max(state.soilp[i][16], state.soilp[i][15] * 1.01)

            c1 = -0.6 * (2.0 + math.log10(state.soilp[i][7] * 8.64e4))
            w1 = state.soilp[i][11] ** c1
            state.soilp[i][17] = state.soilp[i][8] + w1 * (state.soilp[i][9] - state.soilp[i][8])

            c1 = 1.0 + (state.soilp[i][10] * abs(FLDCAP)) ** state.soilp[i][11]
            w1 = 1.0 / (c1 ** state.soilp[i][12])
            temp = state.soilp[i][8] + w1 * (state.soilp[i][9] - state.soilp[i][8])
            state.soilp[i][17] = max(0.75 * state.soilp[i][9], state.soilp[i][17], temp)
            state.soilp[i][24] = 0.999 * state.soilp[i][9]
        else:
            sph[i] = -1e-5
            state.soilp[i][15] = state.soilp[i][8]
            state.soilp[i][16] = state.soilp[i][15]
            state.soilp[i][17] = state.soilp[i][9]
            state.soilp[i][24] = state.soilp[i][9]

        if abs(state.soilp[i][8] - state.soilp[i][9]) <= EPS:
            print(" minimum and maximum water contents are equal")
            # equivalente ao "stop" do Fortran (interrompe a execucao,
            # incondicional -- nao depende de single_multi_flag, ao contrario
            # do stop de umidade inicial mais adiante nesta mesma sub-rotina).
            raise SystemExit(
                "initprofile: minimum and maximum water contents are equal "
                f"(camada {i}, soiltype={state.soiltype[i]})"
            )

        if state.soiltype[i] != 26 and state.soiltype[i] != 27:
            if abs(state.soilp[i][6] - SPFLAG) <= EPS:
                if state.soiltype[i] <= 4:
                    state.soilp[i][6] = 0.039 * (state.soilp[i][2] ** (-2.2))
                    state.soilp[i][6] = -0.56 * state.soilp[i][2] + 0.51
                elif 4 < state.soiltype[i] <= 18:
                    f1 = 1.0 / (1.0 - state.soilp[i][2])
                    state.soilp[i][6] = (0.135 * state.soilp[i][1] * 1e3 + 64.7) / (
                        state.soilp[i][1] * 1e3 * (f1 - 0.947)
                    )
                    if state.soiltype[i] != 15:
                        state.soilp[i][6] = -0.56 * state.soilp[i][2] + 0.51
                    else:
                        state.soilp[i][6] = 2.0 * state.soilp[i][6]

        for j in range(1, MAXP + 1):
            state.soilp[i][j] = _anint(state.soilp[i][j], 10)

    # =================================================================
    # snow/ice on the ground: albedo and emissivity
    # (initprofile.F90, linhas 322-348)
    # =================================================================
    state.sgralbedo = state.soilp[1][3]
    state.sgremis = state.soilp[1][4]
    state.albedo_fasst = 1.0
    newsnow = 0.0

    if ((abs(state.hi) <= EPS and abs(state.hsaccum) <= EPS)
            and (_aint(state.met[1][MetCol.PT]) != 3 or _aint(state.met[1][MetCol.PT2]) != 3)):
        state.albedo_fasst = min(state.albedo_fasst, state.sgralbedo)
        state.emis = state.sgremis
        state.hsaccum = 0.0
        state.hi = 0.0
        newsnow = 0.0
    elif state.hsaccum > EPS:
        state.emis = SEMIS
        if _aint(state.met[1][MetCol.PT]) != 3 or _aint(state.met[1][MetCol.PT2]) != 3:
            state.albedo_fasst = min(state.albedo_fasst, SOALBEDO)
            newsnow = 0.0
        else:
            state.albedo_fasst = min(state.albedo_fasst, SNALBEDO)
            newsnow = (state.met[1][MetCol.PREC] + state.met[1][MetCol.PREC2]) * 1e-3
    elif state.hi > EPS and abs(state.hsaccum) <= EPS:
        state.albedo_fasst = min(state.albedo_fasst, IALBEDO)
        state.emis = IEMIS
        state.hsaccum = 0.0
        newsnow = 0.0

    htot = state.hsaccum + state.hi + newsnow

    # =================================================================
    # node placement based on # of layers, yes/no to temp/moist meas.
    # (initprofile.F90, linhas 350-462, inclui a chamada a `sort`)
    # =================================================================
    totthick = 0.0
    for i in range(1, state.nlayers + 1):
        totthick += state.lthick[i]
        nz2[i] = totthick
    sumthick = totthick
    if totthick < MINDEPTH:
        totthick = MINDEPTH

    ttest_y = ttest in ('y', 'Y')
    mtest_y = mtest in ('y', 'Y')

    nnodesf = 0

    if ttest_y and mtest_y:
        nnodesi = nt0
        if abs(state.zti[1]) > EPS:
            for j in range(2, nnodesi + 2):
                nz1[j] = state.zti[j - 1]
            nnodesi = nt0 + 1
        else:
            for j in range(1, nnodesi + 1):
                nz1[j] = state.zti[j]
            nnodesi = nt0

        nnodesi2 = nm0
        if abs(state.zm[1]) > EPS:
            for j in range(2, nnodesi2 + 2):
                nz4[j] = state.zm[j - 1]
            nnodesi2 = nm0 + 1
        else:
            for j in range(1, nnodesi2 + 1):
                nz4[j] = state.zm[j]
            nnodesi2 = nm0

        nnodesf = nnodesi + nnodesi2 + state.nlayers
        for i in range(1, nnodesf + 1):
            if i <= nnodesi:
                nz3[i] = nz1[i]
            elif nnodesi < i <= nnodesi2 + nnodesi:
                nz3[i] = nz4[i - nnodesi]
            else:
                nz3[i] = nz2[i - (nnodesi + nnodesi2)]

    elif ttest_y and (mtest in ('n', 'N')):
        nnodesi = nt0
        if abs(state.zti[1]) > EPS:
            for j in range(2, nnodesi + 2):
                nz1[j] = state.zti[j - 1]
            nnodesi = nt0 + 1
        else:
            for j in range(1, nnodesi + 1):
                nz1[j] = state.zti[j]
            nnodesi = nt0

        nnodesf = nnodesi + state.nlayers
        for i in range(1, nnodesf + 1):
            if i <= nnodesi:
                nz3[i] = nz1[i]
            else:
                nz3[i] = nz2[i - nnodesi]

    elif (ttest in ('n', 'N')) and mtest_y:
        nnodesi = nm0
        if abs(state.zm[1]) > EPS:
            for j in range(2, nnodesi + 2):
                nz1[j] = state.zm[j - 1]
            nnodesi = nm0 + 1
        else:
            for j in range(1, nnodesi + 1):
                nz1[j] = state.zm[j]
            nnodesi = nm0

        nnodesf = nnodesi + state.nlayers
        for i in range(1, nnodesf + 1):
            if i <= nnodesi:
                nz3[i] = nz1[i]
            else:
                nz3[i] = nz2[i - nnodesi]

    elif (ttest in ('n', 'N')) and (mtest in ('n', 'N')):
        nz1[2] = totthick
        nnodesi = 2

        nnodesf = nnodesi + state.nlayers
        for i in range(1, nnodesf + 1):
            if i <= nnodesi:
                nz3[i] = nz1[i]
            else:
                nz3[i] = nz2[i - nnodesi]

    nnodesi, state.nz = sort(nnodesf, nz3)

    # =================================================================
    # add more nodes if spacing is too large
    # (initprofile.F90, linhas 464-533)
    # =================================================================
    nz1 = _new_array(MAXN)
    nz1[1] = state.nz[1]
    icount = 1
    if state.water_flag == 0:
        while nz1[icount] < state.nz[nnodesi]:
            icount += 1
            if nz1[icount - 1] <= 1.0:
                nz1[icount] = nz1[icount - 1] + 1.5e-2 * (2.0 * float(icount) - 3.0)
            else:
                nz1[icount] = nz1[icount - 1] + 3e-2 * (2.0 * float(icount) - 3.0)
        if nz1[icount] > state.nz[nnodesi]:
            icount -= 1

        for i in range(2, icount + 1):
            for k in range(2, nnodesi + 1):
                if abs(nz1[i] - state.nz[k]) <= 0.5 * (1.5e-2 * (2.0 * float(i) - 3.0)):
                    nz1[i] = state.nz[k]

        for i in range(2, nnodesi + 1):
            icount += 1
            nz1[icount] = state.nz[i]

        nnodesi, state.nz = sort(icount, nz1)

    elif state.water_flag in (1, 2):
        while nz1[icount] < state.nz[nnodesi]:
            icount += 1
            if nz1[icount - 1] <= 1.0:
                nz1[icount] = nz1[icount - 1] + 1.5e-2 * (2.0 * float(icount) - 3.0)
            else:
                nz1[icount] = nz1[icount - 1] + 3e-2 * (2.0 * float(icount) - 3.0)
        if nz1[icount] > state.nz[nnodesi]:
            icount -= 1

        icount += 1
        if state.vegh_type > 0:
            nz1[icount] = 3.0
        elif state.vegl_type > 0:
            nz1[icount] = 1.5
        else:
            nz1[icount] = 10.0
        icount += 1
        nz1[icount] = nz1[icount - 1] + 1.0

        for i in range(2, icount + 1):
            for k in range(2, nnodesi + 1):
                if abs(nz1[i] - state.nz[k]) <= 1.5e-2 * (2.0 * float(i) - 3.0):
                    nz1[i] = state.nz[k]

        for i in range(2, nnodesi + 1):
            icount += 1
            nz1[icount] = state.nz[i]

        nnodesi, state.nz = sort(icount, nz1)

    # add more nodes if total thickness < 1m (initprofile.F90, linhas 535-570)
    if state.nz[nnodesi] - totthick < 0.0:
        delzz = min(0.3, 2.0 * (state.nz[nnodesi] - state.nz[nnodesi - 1]))
        rnnodes = (totthick - state.nz[nnodesi]) / delzz
        nnodes = nnodesi + int(rnnodes) + 1

        nz1 = _new_array(MAXN)
        for j in range(1, nnodesi + 1):
            nz1[j] = state.nz[j]

        for j in range(nnodesi + 1, nnodes + 1):
            nz1[j] = nz1[j - 1] + delzz

        while nz1[nnodes] > totthick:
            nnodes -= 1

        while nz1[nnodes] < totthick:
            nnodes += 1
            nz1[nnodes] = nz1[nnodes - 1] + delzz

        nnodesf, state.nz = sort(nnodes, nz1)
        nnodes = nnodesf
    else:
        nnodes = nnodesi

    # =================================================================
    # adjust number of nodes for presence of snow and/or low veg
    # (initprofile.F90, linhas 572-640; inclui ntot/totthick/test_thick)
    # =================================================================
    ntemp = nnodes
    icase = 0
    node_type = state.node_type

    if state.veg_flagl == 0 or state.sigfl <= EPS:
        if state.hm > EPS:
            icase = 1
            if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                ntemp = nnodes + 2
            else:
                ntemp = nnodes + 1
            for i in range(nnodes + 1, ntemp + 1):
                node_type[i] = 'HM'
    else:
        if state.hm <= EPS:
            icase = 2
            ntemp = nnodes + 1
            node_type[ntemp] = 'VG'
        else:
            if state.hfol_tot - state.hm <= EPS:
                icase = 4
                if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                    if (abs(state.hfol_tot - state.hm) <= EPS
                            or abs(state.hfol_tot - HT_MIN) <= EPS):
                        ntemp = nnodes + 2
                        if abs(state.hfol_tot - state.hm) <= EPS:
                            node_type[ntemp] = 'MX'
                        else:
                            node_type[ntemp] = 'HM'
                    else:
                        ntemp = nnodes + 3
                        node_type[ntemp] = 'HM'
                        if state.hfol_tot > HT_MIN:
                            node_type[nnodes + 2] = 'MX'
                        else:
                            node_type[nnodes + 2] = 'HM'
                    node_type[nnodes + 1] = 'MX'
                elif state.hm <= HT_MIN:
                    if abs(state.hfol_tot - state.hm) <= EPS:
                        ntemp = nnodes + 1
                    else:
                        ntemp = nnodes + 2
                    node_type[nnodes + 1] = 'MX'
                    if ntemp == nnodes + 2:
                        node_type[nnodes + 2] = 'MX'
            elif state.hfol_tot - state.hm > EPS:
                icase = 3
                if state.hm > HT_MIN and state.hm - HT_MIN > HT_MINM:
                    ntemp = nnodes + 3
                else:
                    ntemp = nnodes + 2
                node_type[ntemp] = 'VG'
                node_type[nnodes + 1] = 'MX'
                if ntemp == nnodes + 3:
                    node_type[nnodes + 2] = 'MX'

    ntot = nnodes + EXTRAN
    totthick = state.nz[nnodes]
    test_thick = totthick

    # reverse nodes so that 1 is at the bottom and nnodes at the surface
    # (initprofile.F90, linhas 642-656)
    for i in range(1, nnodes + 1):
        j = nnodes - i + 1
        nzt[i] = state.nz[j]
    if state.nz[nnodes] < EPS:
        state.nz[nnodes] = 0.0

    refn = 0
    for i in range(1, nnodes + 1):
        state.nz[i] = _anint(nzt[i], 5)
        if 0.1 < state.nz[i] < 0.3:
            refn = i

    if sumthick < totthick:
        state.lthick[state.nlayers] = state.lthick[state.nlayers] + (totthick - sumthick)

    # =================================================================
    # initialize soil type at a node
    # (initprofile.F90, linhas 658-707)
    # =================================================================
    delzi = totthick
    rotest = 'n'
    for i in range(state.nlayers, 0, -1):
        botz = _anint((delzi - state.lthick[i]), 5)
        for j in range(1, nnodes + 1):
            if j == 1 and i == state.nlayers:
                node_type[j] = state.stype[state.nlayers]
                state.ntype[j] = state.soiltype[state.nlayers]
                for k in range(1, MAXP + 1):
                    state.nsoilp[j][k] = state.soilp[state.nlayers][k]
                state.pheadmin[j] = sph[state.nlayers]
                nstr_flag[j] = str_flag[state.nlayers]
                if state.ntype[j] >= 19:
                    rotest = 'y'
                state.icourse[j] = ICS[state.ntype[j]]
            elif (state.nz[j] - botz) > EPS and (state.nz[j] - delzi) <= EPS:
                node_type[j] = state.stype[i]
                state.ntype[j] = state.soiltype[i]
                for k in range(1, MAXP + 1):
                    state.nsoilp[j][k] = state.soilp[i][k]
                state.pheadmin[j] = sph[i]
                nstr_flag[j] = str_flag[i]
                if state.ntype[j] >= 19:
                    rotest = 'y'
                state.icourse[j] = ICS[state.ntype[j]]
            elif j == nnodes or state.nz[j] < botz:
                node_type[j] = state.stype[1]
                state.ntype[j] = state.soiltype[1]
                for k in range(1, MAXP + 1):
                    state.nsoilp[j][k] = state.soilp[1][k]
                state.pheadmin[j] = sph[1]
                nstr_flag[j] = str_flag[1]
                if state.ntype[j] >= 19:
                    rotest = 'y'
                state.icourse[j] = ICS[state.soiltype[1]]

        for l in range(1, nm0 + 1):
            if botz <= state.zm[l] <= delzi:
                if str_flag[i] <= 2:
                    state.sm[l] = state.sm[l] * state.soilp[i][9]
                if state.sm[l] < 1.1 * state.soilp[i][8]:
                    state.sm[l] = 1.1 * state.soilp[i][8]
                if state.sm[l] > 0.99 * state.soilp[i][9]:
                    state.sm[l] = 0.99 * state.soilp[i][9]

        delzi = botz

    # =================================================================
    # initialize temperature and head above the "soil"
    # (initprofile.F90, linhas 709-756)
    # =================================================================
    sid = nnodes
    for i in range(1, nnodes + 1):
        if 0.9 - state.nz[i] <= EPS:
            sid = i

    pres = state.met[1][MetCol.AP]
    state.stt[nnodes + 1] = state.met[1][MetCol.TMP] + TREF
    vp = vap_press(state, nnodes + 1, 1e-2 * state.met[1][MetCol.RH], pres)
    rhoa = (pres * 1e2 - vp) / (RD * (state.met[1][MetCol.TMP] + TREF))
    hatm = -vp / (rhoa * GRAV)
    hatm = _anint(hatm, 10)

    for i in range(ntemp + 1, ntot + 1):
        state.phead[i] = hatm
        state.stt[i] = state.met[1][MetCol.TMP] + TREF
    state.ftemp = state.stt[ntemp + 1]

    if icase == 0:
        if (_aint(abs(state.met[1][MetCol.TSOIL] - state.mflag) * 1e5) * 1e-5 > EPS
                and state.mstflag == 1):
            for i in range(ntemp + 1, ntot + 1):
                state.stt[i] = state.met[1][MetCol.TSOIL]
            state.ftemp = state.stt[ntemp]
    elif icase in (1, 4):
        for i in range(nnodes + 1, ntemp + 1):
            state.stt[i] = min(TREF, state.met[1][MetCol.TMP] + TREF)
            state.phead[i] = 0.0
        state.tmelt[state.istart] = state.stt[ntemp]
        state.ftemp = state.stt[ntemp]
    elif icase == 2:
        for i in range(nnodes + 1, ntemp + 1):
            state.stt[i] = state.met[1][MetCol.TMP] + TREF
        state.ftemp = state.stt[ntemp]
    elif icase == 3:
        state.stt[nnodes + 1] = min(TREF, state.met[1][MetCol.TMP] + TREF)
        state.stt[nnodes + 2] = state.met[1][MetCol.TMP] + TREF
        state.tmelt[state.istart] = state.stt[nnodes + 1]
        state.ftemp = state.stt[nnodes + 2]
        state.phead[nnodes + 1] = hatm
        state.phead[nnodes + 2] = 0.0

    # =================================================================
    # initialize soil temperature
    # (initprofile.F90, linhas 758-810; ver docstring do modulo sobre o
    # caso jj==nt0 em locate(), ramo ttest=='y')
    # =================================================================
    fr_test = 0

    if ttest_y:
        if state.zti[nt0] < test_thick:
            nt0 = nttest
            state.zti[nt0] = test_thick
            state.tm[nt0] = state.tm[nt0 - 1]

        if abs(state.zti[1]) > EPS:
            state.stt[nnodes] = (
                (state.met[1][MetCol.TMP] + TREF) * state.iheight
                + state.nz[nnodes - 1] * state.tm[1]
            ) / (state.nz[nnodes - 1] + state.iheight)
            if state.nz[nnodes - 1] < 1e-2:
                state.stt[nnodes] = state.tm[1]
        else:
            state.stt[nnodes] = state.tm[1]

        if state.stt[nnodes] <= TREF:
            fr_test = 1

        for i in range(1, nnodes + 1):
            jj = locate(state.zti, nt0, state.nz[i])
            # Sem guarda: o Fortran interpola incondicionalmente, mesmo se
            # locate devolver jj no limite (jj==nt0) -- ver docstring do
            # modulo. state.zti tem folga suficiente (tamanho MAXN) para
            # que state.zti[jj+1] sempre exista.
            state.stt[i] = (
                ((state.nz[i] - state.zti[jj]) / (state.zti[jj] - state.zti[jj + 1]))
                * (state.tm[jj] - state.tm[jj + 1]) + state.tm[jj]
            )

            if state.stt[i] <= TREF and state.ntype[i] != 27:
                fr_test = 1
            if node_type[i] == 'AI':
                state.stt[i] = state.met[1][MetCol.TMP] + TREF
    else:
        if _aint(abs(state.met[1][MetCol.TSOIL] - state.mflag) * 1e5) * 1e-5 > EPS:
            state.stt[nnodes] = state.met[1][MetCol.TSOIL]
        elif abs(htot) <= EPS and node_type[nnodes] != 'SN':
            state.stt[nnodes] = state.met[1][MetCol.TMP] + TREF
        elif abs(htot) > EPS or node_type[nnodes] == 'SN':
            state.stt[nnodes] = min(TREF, state.met[1][MetCol.TMP] + TREF)

        if state.stt[nnodes] <= TREF and state.ntype[nnodes] != 27:
            fr_test = 1

        state.stt[sid] = state.stt[nnodes] * (1.0 + 2.5e-2 * state.nz[sid])
        if htot > EPS:
            state.stt[sid] = state.stt[nnodes] * (1.0 + 1.5e-2 * state.nz[sid])

        for i in range(1, nnodes):
            state.stt[i] = state.stt[nnodes] * (1.0 + 2.5e-2 * state.nz[i])
            if htot > EPS:
                state.stt[i] = state.stt[nnodes] * (1.0 + 1.5e-2 * state.nz[i])
            if state.nz[i] >= state.nz[sid]:
                state.stt[i] = state.stt[sid]

            if node_type[i] == 'SN':
                state.stt[i] = min(TREF, state.met[1][MetCol.TMP] + TREF)
            if node_type[i] == 'AI':
                state.stt[i] = state.met[1][MetCol.TMP] + TREF
            if state.stt[i] <= TREF and state.ntype[i] != 27:
                fr_test = 1

    if abs(htot) <= EPS and abs(state.zti[1]) <= EPS:
        state.toptemp = state.stt[nnodes]
    else:
        state.toptemp = min(TREF, state.met[1][MetCol.TMP] + TREF)

    # =================================================================
    # initialize soil moisture
    #     determine if a buried rock layer is present
    # (initprofile.F90, linhas 812-949)
    # =================================================================
    tcount = 0
    lcount = 0
    i = nnodes
    while i >= 1:
        if state.ntype[i] >= 19:
            tcount = i                     # node of top of rock or air layer
            lcount = max(lcount, tcount)
        i -= 1
    tcount = lcount

    lcount = 0
    i = 1
    while i < nnodes:
        if state.ntype[i] < 19 and state.ntype[i + 1] >= 19:
            lcount = i                     # node of bottom of rock layer
        i += 1
    if lcount != 0:
        lcount = lcount + 1
    else:
        if tcount != 0:
            lcount = 1
        else:
            lcount = tcount

    # see if more buried layers; extent of rock layer
    rocount = 0                            # number of impermeable layers
    if rotest == 'y':
        for i in range(1, nnodes + 1):
            if state.ntype[i] >= 19:
                rocount += 1
    if rocount == nnodes or state.water_flag in (1, 2):
        mtest = '3'

    if tcount == 0:
        if rocount == nnodes:
            tcount = nnodes
        else:
            tcount = 1

    # determine position of groundwater table if not given
    # based on Beven, K., R. Lamb, P. Quinn, R. Romanowicz, J. Freer (1995)
    # in Singh, V.P. (Ed.), Computer Models of Watershed Hydrology. WRR
    # Publications, PP. 627-668.
    if state.gwl < 0.0:
        i = max(1, lcount - 1)
        if state.slope_fasst > EPS and state.slope_fasst < 90.0:
            aslope = max(0.0, min(1.57, state.slope_fasst * PI / 1.8e2))
            state.gwl = 10.0 * (1e-1 / (state.nsoilp[i][9] - state.nsoilp[i][8])) * math.log(
                abs((-state.pheadmin[i] / totthick + math.cos(aslope))
                    / (totthick * math.tan(aslope)))
            )
            if state.gwl < EPS:
                state.gwl = 10.0
        else:
            state.gwl = 10.0 * (1e-1 / (state.nsoilp[i][9] - state.nsoilp[i][8])) * math.log(
                abs(-state.pheadmin[i] / totthick)
            )
            if state.gwl < EPS:
                state.gwl = 10.0
        if state.vegl_type == 8 or state.vegl_type == 11:
            state.gwl = 50.0

    if mtest in ('y', 'Y'):                # some measured moistures
        if tcount != nnodes:
            if state.zm[nm0] < state.nz[tcount]:
                nm0 = nmtest
                state.zm[nm0] = test_thick
                state.sm[nm0] = (state.nsoilp[tcount][8] + state.nsoilp[tcount][9]) * 0.5
                if state.zm[nm0] >= state.gwl:
                    state.sm[nm0] = state.nsoilp[tcount][9]

        stemp = 0.0
        for i in range(1, nnodes + 1):
            smmax = state.nsoilp[i][9] - state.nsoilp[i][8]

            jj = locate(state.zm, nm0, state.nz[i])

            if abs(state.zm[jj] - state.zm[jj + 1]) > EPS and state.zm[jj + 1] > state.zm[jj]:
                state.soil_moist[i] = (
                    ((state.nz[i] - state.zm[jj]) / (state.zm[jj] - state.zm[jj + 1]))
                    * (state.sm[jj] - state.sm[jj + 1]) + state.sm[jj]
                )

                if state.soil_moist[i] < state.nsoilp[i][8]:
                    state.soil_moist[i] = 1.1 * state.nsoilp[i][8]
                if state.nz[i] >= state.gwl:
                    state.soil_moist[i] = state.nsoilp[i][9]
                if state.soil_moist[i] > state.nsoilp[i][9]:
                    state.soil_moist[i] = 0.9 * state.nsoilp[i][9]
                stemp = state.soil_moist[i]
                if state.ntype[i] > 18 and (state.soil_moist[i] < state.nsoilp[i][8]
                                             or state.soil_moist[i] > state.nsoilp[i][9]):
                    state.soil_moist[i] = state.nsoilp[i][8] + 0.3 * smmax
            elif abs(state.zm[jj] - state.nz[i]) <= EPS:
                state.soil_moist[i] = state.sm[jj]
                if state.soil_moist[i] < state.nsoilp[i][8]:
                    state.soil_moist[i] = 1.1 * state.nsoilp[i][8]
                if state.nz[i] >= state.gwl:
                    state.soil_moist[i] = state.nsoilp[i][9]
                if state.soil_moist[i] > state.nsoilp[i][9]:
                    state.soil_moist[i] = 0.9 * state.nsoilp[i][9]
            else:
                if state.ntype[i] <= 18 and abs(stemp) > EPS:
                    state.soil_moist[i] = stemp + state.sm[nm0] * smmax
                else:
                    state.soil_moist[i] = state.nsoilp[i][8] + 0.3 * smmax

            if state.ntype[i] == 26:
                state.soil_moist[i] = 1.0
            if state.nz[i] >= state.gwl or (state.water_flag in (1, 2) and state.ntype[i] != 26):
                state.soil_moist[i] = state.nsoilp[i][9]

            if state.ntype[i] == 27:
                state.soil_moist[i] = 1e-2 * state.met[1][MetCol.RH]

            if state.soil_moist[i] < state.nsoilp[i][8] or state.soil_moist[i] > state.nsoilp[i][9]:
                if state.single_multi_flag == 0:
                    print(" Initial soil moistures out of range for user "
                          "soil type. Fix input file.")
                # equivalente ao "stop" do Fortran (interrompe a execucao)
                raise SystemExit(
                    "initprofile: initial soil moisture out of range for "
                    f"node {i} (soil type {state.ntype[i]})"
                )
    else:                                   # no measured moistures
        for i in range(1, nnodes + 1):
            smmax = state.nsoilp[i][9] - state.nsoilp[i][8]

            if node_type[i] != 'SN':
                state.soil_moist[i] = state.nsoilp[i][8] + 0.3 * smmax
                if state.ntype[i] == 27:
                    state.soil_moist[i] = 1e-2 * state.met[1][MetCol.RH]
            else:
                state.soil_moist[i] = state.nsoilp[i][8] + 0.98 * smmax
            if state.ntype[i] == 26:
                state.soil_moist[i] = 1.0
            if state.nz[i] >= state.gwl or (state.water_flag in (1, 2) and state.ntype[i] != 26):
                state.soil_moist[i] = state.nsoilp[i][9]

    # =================================================================
    # determine the head, adjust ksat for depth
    # initialize total moisture and ice content (fractions), node state
    # (initprofile.F90, linhas 951-986)
    # =================================================================
    for i in range(nnodes + 1, ntot + 1):
        state.ice[i] = 0.0

    for i in range(1, nnodes + 1):
        state.ice[i] = 0.0
        if i == nnodes and (state.ntype[i] == 26 and state.hm > EPS):
            state.ice[i] = 1.0
            state.soil_moist[i] = 0.0

        state.soil_moist[i] = _anint(state.soil_moist[i], 10)

        if state.ntype[i] < 19:
            deep = 1.0
            if state.elev - state.nz[i] < 1.0:
                if abs(2e-2 * (state.nz[i] - state.elev)) < 50.0:
                    deep = max(0.5, math.exp(2e-2 * (state.nz[i] - state.elev)))
                state.nsoilp[i][7] = deep * state.nsoilp[i][7]

        state.soil_moist[i] = _anint(state.soil_moist[i], 10)
        state.stt[i] = _anint(state.stt[i], 10)
        state.ice[i] = _anint(state.ice[i], 10)

        if state.ntype[i] != 26 and state.ntype[i] != 27:
            state.phead[i] = head(state, i, state.soil_moist[i])
        else:
            state.phead[i] = 0.0

    # =================================================================
    # determine the water vapor content; adjust nz(i) for elevation
    # (initprofile.F90, linhas 988-1031)
    # =================================================================
    state.gwl = state.elev - state.gwl
    for i in range(nnodes + 1, ntot + 1):
        state.wvc[i] = 0.0

    for i in range(1, nnodes + 1):
        if state.ntype[i] != 26 and state.ntype[i] != 27:
            rh = soilhumid(state, i, state.phead[i], state.soil_moist[i], state.stt[i])
            p = state.met[1][MetCol.AP] + 1e-2 * (state.nz[i] + abs(state.phead[i])) \
                * dense(state.stt[i], 0.0, Phase.WATER) * GRAV                    # mbar
        elif state.ntype[i] == 26:
            rh = 1.0
            p = state.met[1][MetCol.AP]
        elif state.ntype[i] == 27:
            rh = state.met[state.istart][MetCol.RH] * 1e-2
            p = state.met[1][MetCol.AP]

        rhow = dense(state.stt[i], 0.0, Phase.WATER)
        vpress = vap_press(state, i, rh, p)                                    # Pa
        mixr = 0.622 * vpress / (p * 1e2 - vpress)                        # kg/kg
        rhov[i] = 0.622 * vpress / (RV * state.stt[i])                      # kg/m^3
        rhoda[i] = max(0.95, min(2.8, (p * 1e2 - vpress) / (RD * state.stt[i])))

        t1 = min(1.0, max(mixr * rhoda[i] / (rhow + mixr * rhoda[i]), 0.0))
        state.wvc[i] = t1 * (state.nsoilp[i][2] - (state.soil_moist[i] + state.ice[i]))
        if state.ntype[i] == 27:
            state.wvc[i] = t1

        state.wvc[i] = max(0.0, min(state.wvc[i], state.nsoilp[i][2] - (state.soil_moist[i] + state.ice[i])))
        if state.ntype[i] == 26:
            state.wvc[i] = 0.0
        state.wvc[i] = _anint(state.wvc[i], 15)

        state.nz[i] = state.elev - state.nz[i]           # place z-axis so m.s.l.=0, z+ upward
        state.nz[i] = _anint(state.nz[i], 5)
        state.nzi[i] = state.nz[i]

    # =================================================================
    # diffusion*dt/dz^2 < 0.5 for stability in finite diff calculations
    # (initprofile.F90, linhas 1033-1076; ver docstring do modulo sobre
    # `deltatt` nao ser reiniciado a cada iteracao, e sobre `th_param`)
    # =================================================================
    state.newsd = 0.0
    # th_param nao e mencionada no comentario de cabecalho do .F90 original
    # ("calls the following subroutines: sort, locate, thcond") -- discrepancia
    # do proprio Fortran. Atualiza state.grthcond/grspheat/km/sphm; `rhotot`
    # (densidade media da camada de neve/gelo) nao e lido de novo neste arquivo.
    # nnodes/ntemp sao locais aqui, mas variaveis globais no Fortran; th_param
    # le os dois de `state`, entao precisam estar gravados antes da chamada.
    state.nnodes = nnodes
    state.ntemp = ntemp
    th_param(state, MAXN, 0.0, 0.0, rhov, rhoda)

    mindeltat = 99999.0                                                   # seconds
    mindeltatm = 99999.0                                                  # seconds
    deltatt = 0.0
    deltatm = 0.0
    for i in range(1, nnodes + 1):
        if i == 1:
            state.delzs[i] = (state.nz[i + 1] - state.nz[i])                          # m
        elif i == nnodes:
            state.delzs[i] = (state.nz[i] - state.nz[i - 1])                          # m
        else:
            state.delzs[i] = 0.5 * (state.nz[i + 1] - state.nz[i - 1])                # m
        state.delzs[i] = _anint(state.delzs[i], 5)
        state.delzsi[i] = state.delzs[i]

        # deltatt NAO e reiniciado a cada iteracao (igual ao Fortran): se
        # dabs(grthcond[i])<=eps, o valor da iteracao anterior e reaproveitado.
        if abs(state.grthcond[i]) > EPS:
            deltatt = 0.45 * state.delzs[i] * state.delzs[i] * state.grspheat[i] / state.grthcond[i]
        mindeltat = min(mindeltat, deltatt)                             # s

        if abs(state.nsoilp[i][7]) > EPS and state.ntype[i] < 20:
            deltatm = 0.45 * state.delzs[i] * state.delzs[i] * state.nsoilp[i][9] / (
                state.nsoilp[i][7] * 1e-2 * abs(head(state, i, state.nsoilp[i][17]))
            )                                                              # s
            mindeltatm = min(mindeltatm, deltatm)

    state.deltat_fasst = mindeltat
    if state.timstep >= 1.0:
        state.deltat_fasst = min(max(state.deltat_fasst, 3e2), state.timstep * 3.6e3)
    else:
        state.deltat_fasst = min(min(state.deltat_fasst, 1e1), state.timstep * 3.6e3)

    rstep = state.timstep * 3.6e3 / state.deltat_fasst
    step = int(rstep)
    step = max(1, step)

    state.deltat_fasst = state.timstep * 3.6e3 / float(step)
    state.deltati = state.deltat_fasst
    state.stepi = step

    # Atualizacoes finais no objeto de estado do modelo. nt0/nm0/sid/totthick
    # nao sao campos de FasstState -- nt0/nm0 voltam na tupla de retorno;
    # sid/totthick sao puramente locais desta sub-rotina.
    state.nnodes = nnodes
    state.ntemp = ntemp
    state.ntot = ntot
    state.icase = icase
    state.refn = refn

    return nt0, nm0, str_flag, mtest, nstr_flag


# ******************************************************************************
# sort  (initprofile.F90, linhas 1079-1121 -- "subroutine sort(n,arr,ncount,b)")
# ******************************************************************************
def sort(n, arr):
    """
    Traducao fiel de `subroutine sort(n,arr,ncount,b)`.

    Ordena arr[1..n] em ordem crescente, IN-PLACE (insertion sort, igual
    ao algoritmo Fortran, inclusive o teste `if(i==1.and.arr(i)>a) i=0`,
    que nunca dispara de fato -- se i==1 e arr(i)>a, o `while` anterior
    ainda estaria rodando -- mas foi preservado por fidelidade). Depois,
    deduplica valores adjacentes cuja diferenca (escalada por 1000) fique
    abaixo de 1e-4; a comparacao e sempre contra o elemento anterior do
    array JA ORDENADO (arr[j-1]), nunca contra o ultimo valor mantido em
    `b` -- exatamente como no original.

    Retorna (ncount, b), com b[1..ncount] os valores deduplicados.
    """
    for j in range(2, n + 1):
        a = arr[j]
        i = j - 1
        while i >= 1 and arr[i] > a:
            arr[i + 1] = arr[i]
            i -= 1
        if i == 1 and arr[i] > a:  # nunca dispara (ver docstring); preservado por fidelidade
            i = 0
        arr[i + 1] = a

    b = _new_array(n)
    ncount = 1
    b[ncount] = arr[1]
    for j in range(2, n + 1):
        a = arr[j] * 1000.0
        if abs(a - arr[j - 1] * 1000.0) > 1e-4:
            ncount += 1
            b[ncount] = arr[j]

    return ncount, b


# ******************************************************************************
# locate  (initprofile.F90, linhas 1124-1156 -- "subroutine locate(zmi,n0,z,jj)")
# ******************************************************************************
def locate(zmi, n0, z):
    """
    Traducao fiel de `subroutine locate(zmi,n0,z,jj)`.

    Busca por bissecao (Numerical Recipes `locate`): funciona tanto para
    zmi[1..n0] crescente quanto decrescente, decidindo a direcao pelo
    sinal de (zmi[n0] - zmi[1]). Devolve jj tal que z esta (ou deveria
    estar) entre zmi[jj] e zmi[jj+1] -- incluindo jj==0 ou jj==n0 quando
    z esta fora do intervalo coberto por zmi, sem nenhuma checagem
    adicional, igual ao original. Quem chama locate() e responsavel por
    garantir que zmi tenha um elemento valido em jj+1 (o Fortran tampouco
    checa isso).
    """
    jl = 1
    ju = n0 + 1
    while ju - jl > 1:
        jm = int((ju + jl) * 0.5)  # aint trunca em direcao a zero; sempre positivo aqui
        if (zmi[n0] > zmi[1]) == (z > zmi[jm]):
            jl = jm
        else:
            ju = jm
    return jl


