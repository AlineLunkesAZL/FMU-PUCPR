"""
missing_met.py -- traducao de missing_met.F90

Fonte Fortran: missing_met.F90 (658 linhas). Uma unica sub-rotina
publica (`missing_met`), com tres modos de operacao (`wmode`):
  wmode==1: infer_test==1, continuacao de uma execucao anterior --
            concatena o registro inicial de uma nova execucao aos
            dados ja carregados, preenchendo passos de tempo
            faltantes entre as duas fontes, se houver.
  wmode==2: preenche parametros meteorologicos faltantes (valor ==
            mflag) dentro de um intervalo de registros ja carregados,
            por interpolacao linear entre o registro anterior e o
            seguinte (ou por persistencia/valor padrao nas bordas).
  wmode==3: preenche uma lacuna de passos de tempo inteiros no meio de
            uma serie ja carregada, entre `ipos` e `ie`, usando o
            mesmo esquema de preenchimento de wmode==1 -- mas com uma
            estrategia adicional de persistencia por hora-do-dia (ver
            a variavel `ib`/`ib0` na Secao 1 da documentacao).

Em todos os casos que recalculam radiacao solar/infravermelha para um
registro preenchido, chama sol_zen, Solflx, emisatm e dnirflx (todas
de module_radiation.py) e a funcao met_date (fasst.functions).

Assinatura original:
    subroutine missing_met(wmode,ipos,ie,met1)
    integer(ip),intent(in):: ipos,ie
    integer(ip),intent(in):: wmode
    real(dp),intent(in):: met1(1,maxcol)

Arquitetura
-------------
`missing_met` recebe `state: FasstState` como argumento explicito.
`EPS` vem de `fasst.constants`; `mflag`, `istart`, `iend`, `timstep`,
`met`, `maxlines`, `mgap`, `ncols`, `unit10` sao campos de `FasstState`.
Os indices de coluna do meteorologico (`ip_year`, `ip_doy`, ...) sao os
membros correspondentes de `MetCol`. `sol_zen`/`Solflx`/`emisatm`/
`dnirflx` sao chamadas reais, importadas de `fasst.module_radiation`;
`met_date`, de `fasst.functions`.

`met1`, declarada `met1(1,maxcol)` no Fortran (uma matriz com a
primeira dimensao forcada a 1), chega aqui como uma lista 1D 1-based
(`met1[coluna]`) -- mesma simplificacao que `initsurface.py`/
`initwater.py` ja aplicam ao montar esse array antes de chamar
`missing_met`. `met2`, variavel de trabalho local (nao um parametro),
e uma matriz 2D 1-based (linha, coluna), dimensionada por
`state.maxlines` x `MAXCOL`.

`missing_met` nao devolve nada: `state.met` e mutado em memoria
(sobrescrito linha a linha, dentro do numero de linhas efetivamente
usado) e `state.iend` e atualizado (wmode 1 e 3; wmode 2 preenche
lacunas sem mudar o tamanho da serie).

Pontos de atencao
-------------------
* Provavel bug do Fortran original, preservado sem correcao: no laco
  de preenchimento de registro novo (usado em wmode 1 e 3), a condicao
  `if(j /= 29 .or. j /= 30)` -- com o comentario `!ip_zen, ip_az`,
  sugerindo a intencao de PULAR as colunas zenite/azimute solar (que
  ja foram calculadas por `sol_zen` poucas linhas antes) -- e sempre
  verdadeira, para qualquer valor inteiro de `j`: se `j==29`, entao
  `j/=30` e verdadeiro; se `j==30`, entao `j/=29` e verdadeiro; para
  qualquer outro `j`, os dois lados sao verdadeiros. E uma tautologia
  booleana -- `.or.` onde muito provavelmente deveria haver um
  `.and.`, o mesmo tipo de troca ja encontrada em cloudbase
  (module_radiation.py). Na pratica, isso significa que as colunas
  zenite e azimute, calculadas explicitamente por `sol_zen`, sao
  imediatamente sobrescritas pela formula de interpolacao linear do
  `else` mais abaixo (j nao cai em nenhum dos outros casos especiais).
  Preservado literalmente (`jcol != 29 or jcol != 30`, sempre
  verdadeiro), sem filtrar nada -- corrigir mudaria o resultado
  numerico do modelo.
* A mesma condicao (`daylim = 1.07d0` se `mody != modyo`) esta
  comentada em tres pontos do arquivo (wmode 1, e as duas ocorrencias
  em wmode 3) -- `daylim` fica sempre em 1.042, exceto quando o ano
  muda entre os dois registros (25.07). Codigo morto preservado por
  fidelidade.
* wmode==1 interrompe a execucao (`stop`) se os arquivos meteorologicos
  estiverem afastados demais (`msteps > mgap*tstps`); wmode==3, na
  mesma situacao, so imprime um aviso (para a tela e para
  `state.unit10`) e continua, assumindo persistencia. Essa diferenca e
  do proprio Fortran original, nao desta traducao -- os dois modos
  tratam o mesmo problema (lacuna grande demais) de formas diferentes,
  coerente com wmode==1 ser sobre continuar uma execucao (dado novo
  pode simplesmente estar incompleto) e wmode==3 ser sobre uma lacuna
  no meio de dados que ja existem (nao ha novo arquivo para rejeitar).
"""

import math

from .constants import EPS, MAXCOL
from .indices import MetCol
from .state import FasstState
from .module_radiation import sol_zen, Solflx, emisatm, dnirflx
from .functions import met_date


def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _new_matrix(nrows, ncols, fill=0.0):
    """1-based com folga no indice 0, nas duas dimensoes."""
    return [[fill] * (ncols + 1) for _ in range(nrows + 1)]


def _aint(x, p=0):
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direção a zero."""
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _fill_solar_ir(state, met2, c3):
    """
    Recalcula radiacao solar (TSOL/DIR/DIF/UPSOL, via Solflx) e
    infravermelha (IR, via emisatm+dnirflx) para o registro `c3` de
    `met2`, in-place -- bloco repetido em todos os tres modos sempre
    que um registro novo e criado. Corresponde as linhas 185-231
    (wmode==1), 386-442 (wmode==2) e 587-633 (wmode==3) do .F90 -- as
    tres copias sao identicas no Fortran original; fundidas aqui num
    unico helper para reduzir risco de erro de transcricao, nao por
    escolha do Fortran.
    """
    row = met2[c3]

    if row[MetCol.ZEN] >= 90.0:
        row[MetCol.TSOL] = 0.0
        row[MetCol.DIR] = 0.0
        row[MetCol.DIF] = 0.0
        row[MetCol.UPSOL] = 0.0
    else:
        row[MetCol.TSOL] = state.mflag
        row[MetCol.DIR] = state.mflag
        row[MetCol.DIF] = state.mflag
        row[MetCol.UPSOL] = state.mflag

        cover = _new_array(3)
        hgt = _new_array(3)
        icld = [0, 0, 0, 0]
        cover[1] = row[MetCol.LCD]
        cover[2] = row[MetCol.MCD]
        cover[3] = row[MetCol.HCD]
        hgt[1] = row[MetCol.LHGT]
        hgt[2] = row[MetCol.MHGT]
        hgt[3] = row[MetCol.HHGT]
        icld[1] = int(row[MetCol.LCT])
        icld[2] = int(row[MetCol.MCT])
        icld[3] = int(row[MetCol.HCT])

        prcp = row[MetCol.PREC] + row[MetCol.PREC2]

        sdir, sdif, tsol_new, icld, cover, hgt = Solflx(
            state, icld, row[MetCol.ZEN], row[MetCol.DOY], prcp,
            row[MetCol.TSOL], cover, hgt,
        )
        row[MetCol.DIR] = sdir
        row[MetCol.DIF] = sdif
        row[MetCol.TSOL] = tsol_new

        if abs(row[MetCol.TSOL]) <= EPS:
            row[MetCol.UPSOL] = 0.0
        row[MetCol.LCD] = cover[1]
        row[MetCol.MCD] = cover[2]
        row[MetCol.HCD] = cover[3]
        row[MetCol.LHGT] = hgt[1]
        row[MetCol.MHGT] = hgt[2]
        row[MetCol.HHGT] = hgt[3]
        row[MetCol.LCT] = float(icld[1])
        row[MetCol.MCT] = float(icld[2])
        row[MetCol.HCT] = float(icld[3])

    mcd = row[MetCol.MCD]
    hcd = row[MetCol.HCD]
    ematm = emisatm(state, row[MetCol.TMP], row[MetCol.RH])
    flxir, _mcld_eff, _hcld_eff = dnirflx(
        state, ematm, row[MetCol.TMP], row[MetCol.DOY],
        row[MetCol.LHGT], row[MetCol.MHGT], row[MetCol.HHGT],
        row[MetCol.LCD], mcd, hcd,
    )
    row[MetCol.IR] = flxir


def missing_met(state: FasstState, wmode, ipos, ie, met1):
    maxlines = state.maxlines
    ncols = state.ncols

    c3 = 0
    msteps = 0
    tstps = 0
    ib = 0
    ib0 = 0

    mody = 0.0
    modyo = 0.0
    daylim = 0.0
    t1 = 0.0
    t2o = 0.0
    f1 = 0.0

    rtstps = 24.0 / state.timstep
    tstps = int(rtstps)

    # =================================================================
    # wmode == 1: infer_test == 1, continuacao de execucao anterior
    # (linhas 55-257)
    # =================================================================
    if wmode == 1:
        c3 = 0
        met2 = _new_matrix(maxlines, MAXCOL, state.mflag)

        # determine initial times of new and old data
        t1 = met_date(met1[MetCol.YEAR], met1[MetCol.DOY], met1[MetCol.HR], met1[MetCol.MIN])
        t2o = met_date(
            state.met[state.istart][MetCol.YEAR], state.met[state.istart][MetCol.DOY],
            state.met[state.istart][MetCol.HR], state.met[state.istart][MetCol.MIN],
        )

        if t2o <= t1:  # files overlap
            for jj in range(state.istart, state.iend + 1):
                c3 += 1
                for k in range(1, MAXCOL + 1):
                    met2[c3][k] = state.met[jj][k]
        else:
            c3 += 1
            for k in range(1, MAXCOL + 1):
                met2[c3][k] = met1[k]

            mody = state.met[state.istart][MetCol.YEAR] - _aint(state.met[state.istart][MetCol.YEAR] * 0.25) * 4.0
            modyo = met1[MetCol.YEAR] - _aint(met1[MetCol.YEAR] * 0.25) * 4.0  # noqa: F841 (nunca lida de novo -- ver docstring)
            daylim = 1.042  # 1.003
            # if(mody /= modyo) daylim = 1.07d0  -- comentado no Fortran, ver docstring
            if state.met[state.istart][MetCol.YEAR] != state.met[1][MetCol.YEAR]:
                daylim = 25.07

            if abs(mody) <= EPS:
                msteps = int((t2o - t1) * 367.0 * tstps)
                if (t2o - t1) * 367.0 * tstps < daylim:
                    msteps = 0
            else:
                msteps = int((t2o - t1) * 366.0 * tstps)
                if (t2o - t1) * 366.0 * tstps < daylim:
                    msteps = 0

            if msteps > int(state.mgap * tstps):
                print("Met files too far apart, stopping!")
                print(
                    f"New; year: {_aint(state.met[state.istart][MetCol.YEAR])} "
                    f"day: {_aint(state.met[state.istart][MetCol.DOY])} "
                    f"hour: {_aint(state.met[state.istart][MetCol.HR])} "
                    f"minute: {_aint(state.met[state.istart][MetCol.MIN])}"
                )
                print(
                    f"Old; year: {_aint(state.met[1][MetCol.YEAR])} "
                    f"day: {_aint(state.met[1][MetCol.DOY])} "
                    f"hour: {_aint(state.met[1][MetCol.HR])} "
                    f"minute: {_aint(state.met[1][MetCol.MIN])}"
                )
                raise SystemExit("missing_met: met files too far apart")

            elif msteps != 0 and msteps <= int(2.0 * tstps):
                for i in range(1, msteps + 1):
                    c3 += 1

                    met2[c3][1] = met1[1]
                    met2[c3][2] = met1[2]

                    mody = state.timstep * 10.0 - _aint(state.timstep) * 10.0
                    if abs(mody) <= EPS:
                        met2[c3][3] = met1[3] + i * state.timstep
                        met2[c3][4] = met1[4]
                    else:
                        met2[c3][3] = met1[3] + i * _aint(state.timstep)
                        met2[c3][4] = met1[4] + i * (mody * 0.1)

                    if met2[c3][4] >= 60.0:
                        met2[c3][3] = met1[3] + _aint(met2[c3][4] / 60.0)
                        met2[c3][4] = met2[c3][4] - _aint(met2[c3][4] / 60.0) * 60.0

                    if met2[c3][3] >= 24.0:
                        met2[c3][2] = met1[2] + _aint(met2[c3][3] / 24.0)
                        met2[c3][3] = met2[c3][3] - _aint(met2[c3][3] / 24.0) * 24.0

                    mody = met1[1] - _aint(met1[1] * 0.25) * 4.0
                    if met2[c3][2] > 366.0 and abs(mody) <= EPS:
                        met2[c3][2] = 1.0
                        met2[c3][1] = met2[c3][1] + 1.0
                    elif met2[c3][2] > 365.0 and abs(mody) > EPS:
                        met2[c3][2] = 1.0
                        met2[c3][1] = met2[c3][1] + 1.0

                    met2[c3][MetCol.ZEN], met2[c3][MetCol.AZ] = sol_zen(
                        state, met2[c3][MetCol.YEAR], met2[c3][MetCol.DOY],
                        met2[c3][MetCol.HR], met2[c3][MetCol.MIN],
                    )

                    f1 = 1.0 / (msteps + 1.0)
                    for jcol in range(5, ncols + 1):
                        # ver "Pontos de atencao" no docstring do modulo:
                        # esta condicao e sempre verdadeira.
                        if jcol != 29 or jcol != 30:  # ip_zen, ip_az
                            if jcol == MetCol.PT or jcol == MetCol.PT2:
                                met2[c3][jcol] = met1[jcol]
                            elif jcol == MetCol.LCT or jcol == MetCol.MCT:
                                met2[c3][jcol] = met1[jcol]
                            elif jcol == MetCol.HCT or jcol == 35:
                                met2[c3][jcol] = met1[jcol]
                            else:
                                met2[c3][jcol] = met1[jcol] + i * (state.met[1][jcol] - met1[jcol]) * f1

                    met2[c3][MetCol.SD] = state.mflag

                    if abs(met2[c3][MetCol.PREC]) <= EPS:
                        met2[c3][MetCol.PT] = 1
                    elif met2[c3][MetCol.PREC] > EPS and met2[c3][MetCol.TMP] > 0.0:
                        met2[c3][MetCol.PT] = 2
                    elif met2[c3][MetCol.PREC] > EPS and met2[c3][MetCol.TMP] <= 0.0:
                        met2[c3][MetCol.PT] = 3

                    if abs(met2[c3][MetCol.PREC2]) <= EPS:
                        met2[c3][MetCol.PT2] = 1
                    else:
                        met2[c3][MetCol.PT2] = 3

                    _fill_solar_ir(state, met2, c3)

                for jj in range(state.istart, state.iend + 1):
                    c3 += 1
                    for k in range(1, MAXCOL + 1):
                        met2[c3][k] = state.met[jj][k]

            elif msteps == 0:
                for jj in range(state.istart, state.iend + 1):
                    c3 += 1
                    for k in range(1, MAXCOL + 1):
                        met2[c3][k] = state.met[jj][k]
                # c3 = c3 - 1  -- comentado no Fortran

        for i in range(1, c3 + 1):
            for j in range(1, MAXCOL + 1):
                state.met[i][j] = met2[i][j]

        state.iend = c3

    # =================================================================
    # wmode == 2: preenche parametros faltantes num intervalo ja carregado
    # (linhas 260-443)
    # =================================================================
    elif wmode == 2:
        for i in range(ipos, ie + 1):
            # time
            if _aint(abs(state.met[i][3] - state.mflag) * 1e5) * 1e-5 <= EPS:  # hour
                if i != state.istart:
                    state.met[i][3] = state.met[ipos][3] + (i - ipos) * state.timstep
                else:
                    if abs(state.met[ipos + 1][3]) <= EPS:
                        state.met[ipos + 1][3] = 24.0
                    state.met[i][3] = state.met[ipos + 1][3] - (i - ipos + 1) * state.timstep

            if _aint(abs(state.met[i][4] - state.mflag) * 1e5) * 1e-5 <= EPS:  # minute
                if i != state.istart:
                    state.met[i][4] = state.met[ipos][4] + (i - ipos) * (state.timstep - _aint(state.timstep))
                else:
                    state.met[i][4] = state.met[ipos + 1][4] - (i - ipos + 1) * (state.timstep - _aint(state.timstep))

            if state.met[i][4] >= 60.0:
                state.met[i][3] = state.met[ipos][3] + _aint(state.met[i][4] / 60.0)
                state.met[i][4] = state.met[i][4] - _aint(state.met[i][4] / 60.0) * 60.0

            if state.met[i][3] >= 24.0:
                state.met[i][2] = state.met[ipos][2] + _aint(state.met[i][3] / 24.0)
                state.met[i][3] = state.met[i][3] - _aint(state.met[i][3] / 24.0) * 24.0

            mody = state.met[i][1] - _aint(state.met[i][1] * 0.25) * 4.0
            if _aint(state.met[i][2]) > _aint(366.0) and abs(mody) <= EPS:
                state.met[i][2] = 1.0
                state.met[i][1] = state.met[i][1] + 1.0
            elif _aint(state.met[i][2]) > _aint(365.0) and abs(mody) > EPS:
                state.met[i][2] = 1.0
                state.met[i][1] = state.met[i][1] + 1.0

            # wind, pressure, temp, precipitation, rh
            for j in range(5, 13):
                if _aint(abs(state.met[i][j] - state.mflag) * 1e5) * 1e-5 <= EPS:
                    if i == state.istart:
                        state.met[i][j] = state.met[ipos + 1][j]
                        if _aint(abs(state.met[i][j] - state.mflag) * 1e5) * 1e-5 <= EPS:
                            if j == 5:
                                state.met[i][j] = 1e3
                            if j == 7:
                                state.met[i][j] = 60.0
                            if j == 8:
                                state.met[i][j] = 1.2
                            if j == 9:
                                state.met[i][j] = 0.0
                            if j == 10 or j == 12:
                                state.met[i][j] = 0.0
                            if j == 11:  # precip type
                                if abs(state.met[i][MetCol.PREC]) <= EPS:
                                    state.met[i][MetCol.PT] = 1
                                elif state.met[i][MetCol.PREC] > EPS and state.met[i][MetCol.TMP] > 0.0:
                                    state.met[i][MetCol.PT] = 2
                                elif state.met[i][MetCol.PREC] > EPS and state.met[i][MetCol.TMP] <= 0.0:
                                    state.met[i][MetCol.PT] = 3
                            if j == 13:  # precip2 type
                                if state.met[i][12] > EPS:
                                    state.met[i][j] = 3.0
                                else:
                                    state.met[i][j] = 1.0
                    elif i >= ie:
                        state.met[i][j] = state.met[ie - 1][j]
                    else:
                        if _aint(abs(state.met[i + 1][j] - state.mflag) * 1e5) * 1e-5 > EPS:
                            state.met[i][j] = 0.5 * (state.met[i - 1][j] + state.met[i + 1][j])
                        else:
                            state.met[i][j] = state.met[i - 1][j]

            # clouds
            for j in range(13, 23):
                if _aint(abs(state.met[i][j] - state.mflag) * 1e5) * 1e-5 <= EPS:
                    if i == state.istart:
                        state.met[i][j] = state.met[ipos + 1][j]
                        if _aint(abs(state.met[i][j] - state.mflag) * 1e5) * 1e-5 <= EPS:
                            if state.met[i][10] + state.met[i][12] > EPS:
                                state.met[i][14] = 1.0
                                state.met[i][16] = 6.0
                                state.met[i][17] = 1.0
                                state.met[i][19] = 3.0
                                state.met[i][20] = 1.0
                                state.met[i][22] = 5.0
                            else:
                                state.met[i][14] = 0.5
                                state.met[i][16] = 6.0
                                state.met[i][17] = 0.0
                                state.met[i][19] = 0.0
                                state.met[i][20] = 0.0
                                state.met[i][22] = 0.0
                    elif i >= ie:
                        state.met[i][j] = state.met[ie - 1][j]
                    else:
                        if _aint(abs(state.met[i + 1][j] - state.mflag) * 1e5) * 1e-5 > EPS:
                            state.met[i][j] = 0.5 * (state.met[i - 1][j] + state.met[i + 1][j])
                        else:
                            state.met[i][j] = state.met[i - 1][j]

            for j in (29, 30):
                if _aint(abs(state.met[i][j] - state.mflag) * 1e5) * 1e-5 <= EPS:
                    state.met[i][MetCol.ZEN], state.met[i][MetCol.AZ] = sol_zen(
                        state, state.met[i][MetCol.YEAR], state.met[i][MetCol.DOY],
                        state.met[i][MetCol.HR], state.met[i][MetCol.MIN],
                    )

            # solar
            for j in range(23, 27):
                if (_aint(abs(state.met[i][j] - state.mflag) * 1e5) * 1e-5 <= EPS
                        and j != 26):
                    if state.met[i][MetCol.ZEN] >= 90.0:
                        state.met[i][MetCol.TSOL] = 0.0
                        state.met[i][MetCol.DIR] = 0.0
                        state.met[i][MetCol.DIF] = 0.0
                        state.met[i][MetCol.UPSOL] = 0.0
                    else:
                        state.met[i][MetCol.TSOL] = state.mflag
                        state.met[i][MetCol.DIR] = state.mflag
                        state.met[i][MetCol.DIF] = state.mflag
                        state.met[i][MetCol.UPSOL] = state.mflag

                        cover = _new_array(3)
                        hgt = _new_array(3)
                        icld = [0, 0, 0, 0]
                        cover[1] = state.met[i][MetCol.LCD]
                        cover[2] = state.met[i][MetCol.MCD]
                        cover[3] = state.met[i][MetCol.HCD]
                        hgt[1] = state.met[i][MetCol.LHGT]
                        hgt[2] = state.met[i][MetCol.MHGT]
                        hgt[3] = state.met[i][MetCol.HHGT]
                        icld[1] = int(state.met[i][MetCol.LCT])
                        icld[2] = int(state.met[i][MetCol.MCT])
                        icld[3] = int(state.met[i][MetCol.HCT])

                        prcp = state.met[i][MetCol.PREC] + state.met[i][MetCol.PREC2]

                        sdir, sdif, _tsol_new, icld, cover, hgt = Solflx(
                            state, icld, state.met[i][MetCol.ZEN], state.met[i][MetCol.DOY],
                            prcp, state.met[i][MetCol.TSOL], cover, hgt,
                        )
                        state.met[i][MetCol.DIR] = sdir
                        state.met[i][MetCol.DIF] = sdif

                        state.met[i][MetCol.LCD] = cover[1]
                        state.met[i][MetCol.MCD] = cover[2]
                        state.met[i][MetCol.HCD] = cover[3]
                        state.met[i][MetCol.LHGT] = hgt[1]
                        state.met[i][MetCol.MHGT] = hgt[2]
                        state.met[i][MetCol.HHGT] = hgt[3]
                        state.met[i][MetCol.LCT] = float(icld[1])
                        state.met[i][MetCol.MCT] = float(icld[2])
                        state.met[i][MetCol.HCT] = float(icld[3])

            # ir
            for j in (27, 28):
                if _aint(abs(state.met[i][j] - state.mflag) * 1e5) * 1e-5 <= EPS:
                    if j == 27:
                        mcd = state.met[i][MetCol.MCD]
                        hcd = state.met[i][MetCol.HCD]
                        ematm = emisatm(state, state.met[i][MetCol.TMP], state.met[i][MetCol.RH])
                        flxir, _mcld_eff, _hcld_eff = dnirflx(
                            state, ematm, state.met[i][MetCol.TMP], state.met[i][MetCol.DOY],
                            state.met[i][MetCol.LHGT], state.met[i][MetCol.MHGT],
                            state.met[i][MetCol.HHGT], state.met[i][MetCol.LCD], mcd, hcd,
                        )
                        state.met[i][MetCol.IR] = flxir
                    # j == 28 (ip_irup): ramo vazio no Fortran original -- nunca faz nada.

    # =================================================================
    # wmode == 3: preenche uma lacuna de passos de tempo no meio da serie
    # (linhas 447-656)
    # =================================================================
    elif wmode == 3:
        c3 = 0
        met2 = _new_matrix(maxlines, MAXCOL, state.mflag)

        # determine gap in data
        t1 = met_date(met1[MetCol.YEAR], met1[MetCol.DOY], met1[MetCol.HR], met1[MetCol.MIN])  # beginning
        t2o = met_date(
            state.met[ie][MetCol.YEAR], state.met[ie][MetCol.DOY],
            state.met[ie][MetCol.HR], state.met[ie][MetCol.MIN],
        )  # end

        mody = state.met[ie][MetCol.YEAR] - _aint(state.met[ie][MetCol.YEAR] * 0.25) * 4.0
        modyo = met1[MetCol.YEAR] - _aint(met1[MetCol.YEAR] * 0.25) * 4.0  # noqa: F841 (nunca lida de novo -- ver docstring)
        daylim = 1.042  # 1.003
        # if(mody /= modyo) daylim = 1.07d0  -- comentado no Fortran, ver docstring
        if met1[MetCol.YEAR] != state.met[ie][MetCol.YEAR]:
            daylim = 25.07

        if abs(mody) <= EPS:
            msteps = int((t2o - t1) * 367.0 * tstps)
            if (t2o - t1) * 367.0 * tstps < daylim:
                msteps = 0
        else:
            msteps = int((t2o - t1) * 366.0 * tstps)
            if (t2o - t1) * 366.0 * tstps < daylim:
                msteps = 0

        if msteps > int(state.mgap * tstps):  # >= 2 (user defined) days missing
            msg = (
                f"Missing met data over {state.mgap:6.2f} days apart. Filled in "
                f"parameters except solar and IR are dubious. Persistence assumed.\n"
                f" Year {int(state.met[ie][MetCol.YEAR]):5d}, DOY {int(state.met[ie][MetCol.DOY]):4d}, "
                f"Hour {int(state.met[ie][MetCol.HR]):3d}, Minute {int(state.met[ie][MetCol.MIN]):3d}\n"
            )
            print(msg)
            if hasattr(state, "unit10"):
                state.unit10.write(msg)

        ib = ipos
        if ib != 1:
            ib = ib - 1
        while abs(
            (state.met[ib][MetCol.HR] + state.met[ib][MetCol.MIN] / 60.0)
            - (state.met[ipos][MetCol.HR] + state.met[ipos][MetCol.MIN] / 60.0)
        ) > EPS:
            ib = ib - 1
            if ib == 0:
                break
        if ib != 0:
            ib = ib + 1
        ib0 = ib

        for i in range(1, ipos + 1):
            c3 += 1
            for j in range(1, MAXCOL + 1):
                met2[c3][j] = state.met[i][j]

        for i in range(1, msteps + 1):
            c3 += 1

            met2[c3][1] = met1[1]
            met2[c3][2] = met1[2]

            mody = state.timstep * 10.0 - _aint(state.timstep) * 10.0
            if abs(mody) <= EPS:
                met2[c3][3] = met1[3] + i * state.timstep
                met2[c3][4] = met1[4]
            else:
                met2[c3][3] = met1[3] + i * _aint(state.timstep)
                met2[c3][4] = met1[4] + i * (mody / 10.0)

            if met2[c3][4] >= 60.0:
                met2[c3][3] = met1[3] + _aint(met2[c3][4] / 60.0)
                met2[c3][4] = met2[c3][4] - _aint(met2[c3][4] / 60.0) * 60.0

            if met2[c3][3] >= 24.0:
                met2[c3][2] = met1[2] + _aint(met2[c3][3] / 24.0)
                met2[c3][3] = met2[c3][3] - _aint(met2[c3][3] / 24.0) * 24.0

            mody = met1[1] - _aint(met1[1] * 0.25) * 4.0
            if met2[c3][2] > 366.0 and abs(mody) <= EPS:
                met2[c3][2] = 1.0
                met2[c3][1] = met2[c3][1] + 1.0
            elif met2[c3][2] > 365.0 and abs(mody) > EPS:
                met2[c3][2] = 1.0
                met2[c3][1] = met2[c3][1] + 1.0

            met2[c3][MetCol.ZEN], met2[c3][MetCol.AZ] = sol_zen(
                state, met2[c3][MetCol.YEAR], met2[c3][MetCol.DOY],
                met2[c3][MetCol.HR], met2[c3][MetCol.MIN],
            )

            f1 = 1.0 / (msteps + 1.0)
            for jcol in range(5, ncols + 1):
                # ver "Pontos de atencao" no docstring do modulo: esta
                # condicao e sempre verdadeira.
                if jcol != 29 or jcol != 30:  # ip_zen, ip_az
                    if ib == 0:
                        if jcol == MetCol.PT or jcol == MetCol.PT2:
                            met2[c3][jcol] = met1[jcol]
                        elif jcol == MetCol.LCT or jcol == MetCol.MCT:
                            met2[c3][jcol] = met1[jcol]
                        elif jcol == MetCol.HCT or jcol == 35:
                            met2[c3][jcol] = met1[jcol]
                        else:
                            met2[c3][jcol] = met1[jcol] + i * (state.met[ie][jcol] - met1[jcol]) * f1
                    else:
                        met2[c3][jcol] = state.met[ib][jcol]

            met2[c3][MetCol.SD] = state.mflag

            if abs(met2[c3][MetCol.PREC]) <= EPS:
                met2[c3][MetCol.PT] = 1
            elif met2[c3][MetCol.PREC] > EPS and met2[c3][MetCol.TMP] > 0.0:
                met2[c3][MetCol.PT] = 2
            elif met2[c3][MetCol.PREC] > EPS and met2[c3][MetCol.TMP] <= 0.0:
                met2[c3][MetCol.PT] = 3

            if abs(met2[c3][MetCol.PREC2]) <= EPS:
                met2[c3][MetCol.PT2] = 1
            else:
                met2[c3][MetCol.PT2] = 3

            _fill_solar_ir(state, met2, c3)

            if ib != 0:
                ib = ib + 1
                if ib > ipos:
                    ib = ib0

        for i in range(ie, state.iend + 1):
            c3 += 1
            for j in range(1, MAXCOL + 1):
                met2[c3][j] = state.met[i][j]

        for i in range(1, c3 + 1):
            for j in range(1, MAXCOL + 1):
                state.met[i][j] = met2[i][j]

        state.iend = c3
