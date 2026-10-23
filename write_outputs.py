"""
write_outputs.py -- traducao de write_outputs.F90

Fonte Fortran: write_outputs.F90 (375 linhas). Uma unica sub-rotina
publica (`write_outputs`). Escreve os arquivos de saida do FASST ao
final da simulacao (ou de um trecho dela):

* unidade 2  -- arquivo principal: um registro por passo de tempo, com
  48 colunas (meteorologia, radiacao, estado da superficie, resistencia
  do solo, indicadores de escorregamento);
* unidade 28 -- vegetacao e temperatura de superficie (so se
  `vegint > eps`): 21 colunas, com interpolacao linear entre passos
  quando o intervalo de saida da vegetacao e menor que o passo do
  modelo;
* unidade 24 -- fluxos de energia no solo (so se `water_flag <= 1` e
  `flprint == 1`): 17 colunas.

Assinatura original:
    subroutine write_outputs(flprint,keep,vegint)
    integer(ip),intent(in):: flprint,keep
    real(dp),intent(in):: vegint

Nenhuma sub-rotina e chamada (o proprio .F90 diz "no subroutines
called").

Arquitetura
-------------
`write_outputs` recebe `state: FasstState` como primeiro argumento. As
constantes `EPS`, `TREF` e `NCLAYERS` vem de `fasst.constants` e os
indices de coluna (`ip_*`) de `fasst.indices.MetCol`. Nao ha retorno:
a sub-rotina Fortran so tem argumentos `intent(in)` e seu efeito e (a)
escrever nas unidades 2, 28 e 24 e (b) alterar, in-place, campos do
estado -- tudo preservado:

    met[i][PT], met[i][PT2]  (troca de tipo de precipitacao)
    met[i][SD]               (dado ausente -> 0)
    sdens[i]                 (zerada quando nao ha neve)
    frthick[i], twthick[i], sstate[i]   (so se water_flag > 1)

Arquivos de saida: seguindo a convencao dos demais modulos
(`initsurface.py` usa `state.unit3`, `state.unit55` etc.), as unidades
sao objetos tipo-arquivo em `state.unit2`, `state.unit28` e
`state.unit24`, que recebem `.write(str)`. A convencao de nomes ainda
nao foi confirmada contra o codigo do driver (pendente).

Escrita formatada
-------------------
O Fortran usa FORMAT (1000, 1001, 1002, 1004) com larguras fixas. Para
a saida ser identica byte a byte, este modulo tem emuladores locais dos
descritores usados -- `_fmt_i` (Iw), `_fmt_f` (Fw.d), `_fmt_e` (Ew.d) e
`_fmt_a1` (A1); `nX` vira espacos. Comportamentos do gfortran
reproduzidos: campo estourado vira `*` repetido; zero a esquerda cai
quando falta largura (`0.5` -> `.5`); `-0.00` para negativos que
arredondam a zero; `NaN`/`Infinity` alinhados a direita; expoente de
tres digitos em E perde a letra (`0.1234+100`). Os cabecalhos de texto
(`_HDR_*`) foram extraidos do proprio fonte, trecho a trecho.

Helpers locais (`_aint`, `_fdiv`, `_fint`) nao sao compartilhados com
outros arquivos, igual aos demais modulos do projeto.

Pontos de atencao (preservados do Fortran, nao corrigidos)
------------------------------------------------------------
* `ncol = 43` e o numero escrito no cabecalho da unidade 2, mas cada
  registro tem 48 campos (44 sem as quatro colunas de data/hora). O
  numero `16` da unidade 28 tambem esta defasado (21 campos; 17 sem a
  data/hora); so a unidade 24 (13) bate com as colunas de dados. Quem
  ler esses arquivos pelo cabecalho deve saber disso.
* Unidade 28: o numero de registros declarado no cabecalho e
  `intval*(iend - istartp + 1)`, mas o laco escreve 1 registro para o
  primeiro passo e `intval` registros para cada um dos demais, ou seja,
  `1 + intval*(n - 1)` (igual ao declarado so quando `intval == 1`).
* `fw_flag` e `slippery` sao variaveis locais que NAO sao zeradas a
  cada passo: no ramo `sd > eps`, `hi <= eps`, `i /= 1`, se a neve nao
  diminuiu, `fw_flag` mantem o valor do passo anterior (0 no primeiro).
  Os ramos `slushy(i) == 1` e `fw_flag = 5` so sao alcancaveis com
  `i == 1`, porque `else if (i /= 1)` os antecede.
* `met(i-1,...)` com `i = 1` e `iw /= 1` leria fora dos limites no
  Fortran. Aqui le a linha 0 "fantasma" (zeros), sem erro; so ocorre se
  `keep` fizer `istartp` valer 1 com `iw /= 1`.
* A deteccao de ano bissexto na interpolacao da unidade 28 usa o ano de
  `met(i, ip_year)` (passo atual) e nao o de `t1(j,1)` (passo anterior):
  na virada de ano o teste pode usar o calendario errado.
* Hora (`t1(j,3)`) so e normalizada uma vez (`>= 24` -> `-24`); minutos
  usam laco `while`.
* `sstate` vazio (valor padrao do estado) sai como espaco em A1.
* `int(x)` de NaN/Inf/valor fora de 32 bits: o gfortran/x86 devolve
  -2147483648; `_fint` imita isso (comportamento de plataforma, nao do
  padrao Fortran).
* O bloco final comentado ("temporary for stdmod readin file", formatos
  3000/3001) nao foi traduzido: e codigo morto no original. As variaveis
  `ice_code`, `ttemp` (usada), `max_rain` etc. seguem o original.
"""

import math

import numpy as np

from .constants import EPS, NCLAYERS, TREF
from .indices import MetCol
from .state import FasstState


def _aint(x: float) -> float:
    """Fortran ``AINT(x)``: trunca em direcao a zero (NaN/Inf passam)."""
    if not math.isfinite(x):
        return x
    return float(math.trunc(x))


def _fdiv(x: float, y: float) -> float:
    """Divisao IEEE-754 pura: ``x / 0.0`` -> +-inf / nan, sem levantar."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


def _fint(x: float) -> int:
    """Fortran ``INT(x)`` (inteiro de 32 bits): trunca para zero.

    NaN, Inf ou fora de 32 bits -> -2147483648 (resultado do gfortran em
    x86-64, `cvttsd2si`).
    """
    if x != x or abs(x) >= 2147483648.0:
        return -2147483648
    return int(x)


# ---------------------------------------------------------------------------
# Emuladores dos descritores de edicao FORMAT do Fortran
# ---------------------------------------------------------------------------
def _special(x: float, w: int):
    """NaN / +-Inf em um campo de largura `w` (None se `x` e finito)."""
    if x != x:
        s = "NaN"
    elif x == math.inf:
        s = "Infinity" if w >= 8 else "Inf"
    elif x == -math.inf:
        s = "-Infinity" if w >= 9 else "-Inf"
    else:
        return None
    return s.rjust(w) if len(s) <= w else "*" * w


def _fmt_i(n: int, w: int) -> str:
    """Descritor ``Iw``."""
    s = str(int(n))
    return s.rjust(w) if len(s) <= w else "*" * w


def _fmt_f(x: float, w: int, d: int) -> str:
    """Descritor ``Fw.d``."""
    sp = _special(x, w)
    if sp is not None:
        return sp
    s = f"{x:.{d}f}"
    if len(s) > w:  # o zero a esquerda do ponto e opcional
        if s.startswith("0."):
            s = s[1:]
        elif s.startswith("-0."):
            s = "-" + s[2:]
    return s.rjust(w) if len(s) <= w else "*" * w


def _fmt_e(x: float, w: int, d: int) -> str:
    """Descritor ``Ew.d`` (forma ``0.dddd E+ee``)."""
    sp = _special(x, w)
    if sp is not None:
        return sp
    neg = math.copysign(1.0, x) < 0.0
    if x == 0.0:
        digits = "0" * d
        ex = 0
    else:
        m, e = f"{abs(x):.{d - 1}e}".split("e")
        digits = m.replace(".", "")
        ex = int(e) + 1
    if abs(ex) <= 99:
        es = f"E{ex:+03d}"
    elif abs(ex) <= 999:
        es = f"{ex:+04d}"
    else:
        return "*" * w
    s = ("-" if neg else "") + "0." + digits + es
    return s.rjust(w) if len(s) <= w else "*" * w


def _fmt_a1(ch) -> str:
    """Descritor ``A1`` aplicado a ``character(len=1)``."""
    ch = str(ch)
    return ch[0] if ch else " "


_HDR_MAIN_NOMES = (
    'Year  JD  Hr   M   APres   ATemp     RH     WnSp   '
    ' WDir    Prec   PT   Prec2   PT2  LCAmt    LCHt   LCT  '
    'MCAmt    MCHt   MCT   HCAmt   HCHt   HCT   STot    SDir    '
    'SDif   UpSol       IR     IRUp     SZen     SAz    Snow  '
    '   Ice  Grtemp  Grmoistv  Grmoistm    RCI6   RCI12    CBR '
    ' Frz Dp. T/U Dp. F/T/U    Sdens   Slpry FWType   SOdens     '
    '   Emis       Vis        TotD'
)

_HDR_MAIN_UNID = (
    '                    mbar     K               m/s   '
    '        mm/stp       mm/stp                 km           '
    '        km                   km         W/m^2   W/m^2   '
    'W/m^2   W/m^2    W/m^2    W/m^2                     m      '
    ' m       K        % vol    % mass                          '
    '    m      m            kg/m^3                 kg/m^3      '
    '             km           m'
)

_HDR_VEG_NOMES = (
    'Year  JD  Hr   M   LVegTyp   LVegTmp     LVegEms'
    '      LVegEvp  CTyp    CTopTmp    CMTmp    CLTmp   CAvEms  '
    '    CTotEvp   SurfTmp    Grtemp  Grmoistv         GrEvp  '
    'Snow+Ice     ATTmp    ASTmp'
)

_HDR_VEG_UNID = (
    '                                 K              '
    '    kg/m^2*s              K        K        K              '
    '  kg/m^2*s       K         K       % vol      kg/m^2*s     '
    ' m          K        K'
)

_HDR_FLUXO_NOMES = (
    'Year  JD  Hr   M     SDown       SUp    IRdown  '
    '    IRup     LHeat     SHeat     PHeat     CHeat    CHeat1'
    '       Sum    SnDpth     MlTmp     SfTmp'
)

_HDR_FLUXO_UNID = (
    '                     W/m^2     W/m^2     W/m^2  '
    '   W/m^2     W/m^2     W/m^2     W/m^2     W/m^2     W/m^2  '
    '   W/m^2       m         K         K'
)


def _fmt_1000(freq_id, nrec, ncol, lat, mlong, elev, vitd_index, met_count,
              timeoffset, timstep, mflag):
    """FORMAT 1000: i10,1x,i6,1x,i3,1x,f10.6,1x,f11.6,1x,f11.6,1x,i8,1x,i8,
    1x,f6.2,1x,f5.2,1x,f8.2"""
    return (_fmt_i(freq_id, 10) + " " + _fmt_i(nrec, 6) + " "
            + _fmt_i(ncol, 3) + " " + _fmt_f(lat, 10, 6) + " "
            + _fmt_f(mlong, 11, 6) + " " + _fmt_f(elev, 11, 6) + " "
            + _fmt_i(vitd_index, 8) + " " + _fmt_i(met_count, 8) + " "
            + _fmt_f(timeoffset, 6, 2) + " " + _fmt_f(timstep, 5, 2) + " "
            + _fmt_f(mflag, 8, 2))


def _fmt_1004(t, vegl_type, vegh_type, v):
    """FORMAT 1004 (unidade 28). `t` = (ano, dia, hora, min); `v` = os 15
    valores reais na ordem do WRITE:
    ft, surfemisf, f2, ct1, ct2, ct3, surfemisc, cevap, tt, tsoil,
    swater, f3, t8, airt1, airt2 (15 valores).

    4(i4),3x,i4,5x,f8.2,3x,f9.3,2x,e11.4,1x,i4,3x,4(f9.2),2x,e11.4,1x,
    f9.2,2x,f8.2,2x,f8.5,3x,e11.4,1x,f8.4,2x,2(f9.2)
    """
    (ft, semf, f2, c1, c2, c3, semc, cev, tt, tsoil, swat, f3, t8,
     a1, a2) = v
    return (_fmt_i(t[0], 4) + _fmt_i(t[1], 4) + _fmt_i(t[2], 4)
            + _fmt_i(t[3], 4) + "   " + _fmt_i(vegl_type, 4) + "     "
            + _fmt_f(ft, 8, 2) + "   " + _fmt_f(semf, 9, 3) + "  "
            + _fmt_e(f2, 11, 4) + " " + _fmt_i(vegh_type, 4) + "   "
            + _fmt_f(c1, 9, 2) + _fmt_f(c2, 9, 2) + _fmt_f(c3, 9, 2)
            + _fmt_f(semc, 9, 2) + "  " + _fmt_e(cev, 11, 4) + " "
            + _fmt_f(tt, 9, 2) + "  " + _fmt_f(tsoil, 8, 2) + "  "
            + _fmt_f(swat, 8, 5) + "   " + _fmt_e(f3, 11, 4) + " "
            + _fmt_f(t8, 8, 4) + "  " + _fmt_f(a1, 9, 2) + _fmt_f(a2, 9, 2))


def write_outputs(state: FasstState, flprint, keep, vegint):
    with np.errstate(all="ignore"):  # NaN/Inf silenciosos, como no Fortran
        _write_outputs(state, flprint, keep, vegint)


def _write_outputs(state: FasstState, flprint, keep, vegint):
    """Escreve os arquivos de saida (unidades 2, 28 e 24). Sem retorno.

    flprint: 1 -> escreve tambem o arquivo de fluxos (unidade 24).
    keep:    < 0 -> escreve desde `istart`; senao, so os ultimos `keep`+1
             passos (`istartp = iend - keep`, minimo 1).
    vegint:  intervalo de saida da vegetacao em minutos (<= eps -> sem
             unidade 28).
    """
    maxlines = state.maxlines
    istart = state.istart
    iend = state.iend
    freq_id = state.freq_id
    iw = state.iw
    lat = state.lat
    mlong = state.mlong
    elev = state.elev
    vitd_index = state.vitd_index
    met_count = state.met_count
    timeoffset = state.timeoffset
    timstep = state.timstep
    mflag = state.mflag
    water_flag = state.water_flag
    vegh_type = state.vegh_type
    vegl_type = state.vegl_type

    met = state.met
    sdens = state.sdens
    sstate = state.sstate
    slushy = state.slushy
    surfice = state.surfice
    surfmoist = state.surfmoist
    surficep = state.surficep
    surfmoistp = state.surfmoistp
    frthick = state.frthick
    twthick = state.twthick
    surfci = state.surfci
    surfrci = state.surfrci
    surfcbr = state.surfcbr
    surfd = state.surfd
    surfemis = state.surfemis
    tot_thick = state.tot_thick
    lhes = state.lhes
    lheatf = state.lheatf
    lheat = state.lheat
    ft = state.ft
    surfemisf = state.surfemisf
    canopy_temp = state.canopy_temp
    surfemisc = state.surfemisc
    cevap = state.cevap
    tt = state.tt
    airt = state.airt
    tmelt = state.tmelt
    sheat = state.sheat
    cheat = state.cheat
    cheat1 = state.cheat1
    pheat1 = state.pheat1
    sdown = state.sdown
    sup = state.sup
    irdown = state.irdown
    irup = state.irup
    evap_heat = state.evap_heat

    unit2 = state.unit2

    # zero-out parameters (persistem entre passos de i: ver docstring)
    slippery = 0
    fw_flag = 0
    intval = 0
    ptime = 0.0
    swater = 0.0
    swaterp = 0.0
    f2 = 0.0
    f3 = 0.0

    ncol = 43  # fixo no original (os registros tem 48 campos)

    if keep < 0:
        istartp = istart
    else:
        istartp = iend - keep
        if istartp < 1:
            istartp = 1

    # t1(maxlines,4) inteiro e t2(maxlines,15) real, zerados
    t1 = np.zeros((maxlines + 1, 5), dtype=np.int64)
    t2 = np.zeros((maxlines + 1, 16), dtype=np.float64)

    # ******************************************************************
    # Main FASST output file
    # ******************************************************************
    unit2.write(_fmt_1000(freq_id, iend - istartp + 1, ncol, lat, mlong, elev,
                          vitd_index, met_count, timeoffset, timstep, mflag)
                + "\n")
    unit2.write(_HDR_MAIN_NOMES + "\n")
    unit2.write(_HDR_MAIN_UNID + "\n")

    # print the met data
    max_rain = 0.25 * 25.4 * timstep * 2.0  # (mm) needed for slippery

    for i in range(istartp, iend + 1):
        ttime = met[i][MetCol.DOY] + met[i][MetCol.HR] / 24.0
        if met[i][MetCol.PREC] >= max_rain:
            ptime = ttime

        if met[i][MetCol.PREC] > EPS and _aint(met[i][MetCol.PT]) == 1:
            met[i][MetCol.PT] = met[i][MetCol.PT2]
            met[i][MetCol.PT2] = 1.0

        if _aint(abs(met[i][MetCol.SD] - mflag) * 1e5) * 1e-5 <= EPS:
            met[i][MetCol.SD] = 0.0

        if water_flag == 0:
            if met[i][MetCol.SD] > EPS:
                slippery = 3
                if met[i][MetCol.HI] > EPS:
                    fw_flag = 6
                    if iw != 1:
                        if met[i][MetCol.HI] < met[i - 1][MetCol.HI]:
                            fw_flag = 1
                elif i != 1:
                    if met[i][MetCol.SD] < met[i - 1][MetCol.SD]:
                        fw_flag = 1
                elif slushy[i] == 1:
                    fw_flag = 4
                else:
                    fw_flag = 5
            elif met[i][MetCol.HI] > EPS and met[i][MetCol.SD] <= EPS:
                slippery = 4  # 2
                fw_flag = 0
            elif ttime - ptime <= 0.25:
                slippery = 2  # 1
                fw_flag = 3
            else:
                slippery = 1  # 0
                fw_flag = 3
        else:
            slippery = 1
            fw_flag = 3
            if met[i][MetCol.HI] > EPS:
                slippery = 4
                fw_flag = 0
                if iw != 1:
                    if met[i][MetCol.HI] < met[i - 1][MetCol.HI]:
                        fw_flag = 1
            elif met[i][MetCol.SD] > EPS:
                slippery = 3
                fw_flag = 6
                if iw != 1:
                    if met[i][MetCol.SD] < met[i - 1][MetCol.SD]:
                        fw_flag = 1

        if abs(met[i][MetCol.SD]) <= EPS:
            sdens[i] = 0.0

        swater = surfmoist[i] + surfice[i]
        swaterp = surfmoistp[i] + surficep[i]

        # determine surface freeze or thaw depth and state
        if water_flag > 1:
            frthick[i] = met[i][MetCol.HI]
            twthick[i] = 0.0
            sstate[i] = "u"
            if frthick[i] > EPS or met[i][MetCol.TSOIL] <= TREF:
                sstate[i] = "f"

        m = met[i]
        # FORMAT 1001: 4(i4),6(f8.2),i5,f8.2,i5,3(2(f8.2),i5),4(f8.2),1x,
        # f8.2,1x,f8.2,1x,2(f8.2),2(f8.5),f8.2,2x,f8.5,2x,f8.5,3(f8.2),
        # f8.4,f8.4,3x,a1,3x,f8.2,7x,i1,3x,i1,4x,f8.3,3x,f9.3,3x,f8.2,
        # 3x,f8.4
        line = (
            _fmt_i(_fint(m[MetCol.YEAR]), 4) + _fmt_i(_fint(m[MetCol.DOY]), 4)
            + _fmt_i(_fint(m[MetCol.HR]), 4) + _fmt_i(_fint(m[MetCol.MIN]), 4)
            + _fmt_f(m[MetCol.AP], 8, 2)
            + _fmt_f(m[MetCol.TMP] + TREF, 8, 2)
            + _fmt_f(m[MetCol.RH], 8, 2) + _fmt_f(m[MetCol.WS], 8, 2)
            + _fmt_f(m[MetCol.WDIR], 8, 2) + _fmt_f(m[MetCol.PREC], 8, 2)
            + _fmt_i(_fint(m[MetCol.PT]), 5)
            + _fmt_f(m[MetCol.PREC2], 8, 2)
            + _fmt_i(_fint(m[MetCol.PT2]), 5)
            + _fmt_f(m[MetCol.LCD], 8, 2) + _fmt_f(m[MetCol.LHGT], 8, 2)
            + _fmt_i(_fint(m[MetCol.LCT]), 5)
            + _fmt_f(m[MetCol.MCD], 8, 2) + _fmt_f(m[MetCol.MHGT], 8, 2)
            + _fmt_i(_fint(m[MetCol.MCT]), 5)
            + _fmt_f(m[MetCol.HCD], 8, 2) + _fmt_f(m[MetCol.HHGT], 8, 2)
            + _fmt_i(_fint(m[MetCol.HCT]), 5)
            + _fmt_f(m[MetCol.TSOL], 8, 2) + _fmt_f(m[MetCol.DIR], 8, 2)
            + _fmt_f(m[MetCol.DIF], 8, 2) + _fmt_f(m[MetCol.UPSOL], 8, 2)
            + " " + _fmt_f(m[MetCol.IR], 8, 2)
            + " " + _fmt_f(m[MetCol.IRUP], 8, 2)
            + " " + _fmt_f(m[MetCol.ZEN], 8, 2) + _fmt_f(m[MetCol.AZ], 8, 2)
            + _fmt_f(m[MetCol.SD], 8, 5) + _fmt_f(m[MetCol.HI], 8, 5)
            + _fmt_f(m[MetCol.TSOIL], 8, 2)
            + "  " + _fmt_f(swater, 8, 5)
            + "  " + _fmt_f(swaterp, 8, 5)
            + _fmt_f(surfci[i], 8, 2) + _fmt_f(surfrci[i], 8, 2)
            + _fmt_f(surfcbr[i], 8, 2)
            + _fmt_f(frthick[i], 8, 4) + _fmt_f(twthick[i], 8, 4)
            + "   " + _fmt_a1(sstate[i])
            + "   " + _fmt_f(sdens[i], 8, 2)
            + "       " + _fmt_i(slippery, 1)
            + "   " + _fmt_i(fw_flag, 1)
            + "    " + _fmt_f(surfd[i], 8, 3)
            + "   " + _fmt_f(surfemis[i], 9, 3)
            + "   " + _fmt_f(m[MetCol.VIS], 8, 2)
            + "   " + _fmt_f(tot_thick[i], 8, 4)
        )
        unit2.write(line + "\n")

    # ******************************************************************
    # vegetation and surface temperature output file
    # ******************************************************************
    if vegint > EPS:
        unit28 = state.unit28
        rintval = timstep * 60.0 / vegint
        intval = _fint(rintval)
        if intval < 1:
            intval = 1
        if intval > maxlines:
            intval = maxlines

        unit28.write(_fmt_1000(freq_id, intval * (iend - istartp + 1), 16, lat,
                               mlong, elev, vitd_index, met_count, timeoffset,
                               timstep, mflag) + "\n")
        unit28.write(_HDR_VEG_NOMES + "\n")
        unit28.write(_HDR_VEG_UNID + "\n")

        for i in range(istartp, iend + 1):
            swater = surfmoist[i] + surfice[i]

            if i == istartp or intval == 1:
                if (_aint(abs(met[i][MetCol.SD] - mflag) * 1e5) * 1e-5 > EPS
                        and _aint(abs(met[i][MetCol.HI] - mflag) * 1e5) * 1e-5 > EPS):
                    t2[i][8] = met[i][MetCol.SD] + met[i][MetCol.HI]
                else:
                    t2[i][8] = 0.0

                f2 = _fdiv(lheatf[i], lhes[i])
                f3 = _fdiv(lheat[i], lhes[i])
                unit28.write(_fmt_1004(
                    (_fint(met[i][MetCol.YEAR]), _fint(met[i][MetCol.DOY]),
                     _fint(met[i][MetCol.HR]), _fint(met[i][MetCol.MIN])),
                    vegl_type, vegh_type,
                    (ft[i], surfemisf[i], f2, canopy_temp[1, i],
                     canopy_temp[2, i], canopy_temp[3, i], surfemisc[i],
                     cevap[i], tt[i], met[i][MetCol.TSOIL], swater, f3,
                     t2[i][8], airt[1][i], airt[2][i])) + "\n")
            else:
                swater1 = surfmoist[i - 1] + surfice[i - 1]

                for j in range(1, intval + 1):
                    if j == intval:
                        if (_aint(abs(met[i][MetCol.SD] - mflag) * 1e5) * 1e-5 > EPS
                                and _aint(abs(met[i][MetCol.HI] - mflag) * 1e5) * 1e-5 > EPS):
                            t2[i][8] = met[i][MetCol.SD] + met[i][MetCol.HI]
                        else:
                            t2[i][8] = 0.0

                        f2 = _fdiv(lheatf[i], lhes[i])
                        f3 = _fdiv(lheat[i], lhes[i])
                        unit28.write(_fmt_1004(
                            (_fint(met[i][MetCol.YEAR]), _fint(met[i][MetCol.DOY]),
                             _fint(met[i][MetCol.HR]), _fint(met[i][MetCol.MIN])),
                            vegl_type, vegh_type,
                            (ft[i], surfemisf[i], f2, canopy_temp[1, i],
                             canopy_temp[2, i], canopy_temp[3, i],
                             surfemisc[i], cevap[i], tt[i],
                             met[i][MetCol.TSOIL], swater, f3, t2[i][8],
                             airt[1][i], airt[2][i])) + "\n")
                    else:
                        t1[j][1] = _fint(met[i - 1][MetCol.YEAR])
                        t1[j][2] = _fint(met[i - 1][MetCol.DOY])
                        t1[j][3] = _fint(met[i - 1][MetCol.HR])
                        t1[j][4] = _fint(met[i - 1][MetCol.MIN] + j * vegint)

                        # corrections may be necessary due to inconsistant
                        # data recording
                        k = 0
                        while int(t1[j][4]) >= 60:
                            t1[j][4] = t1[j][4] - 60
                            k = k + 1
                        t1[j][3] = t1[j][3] + k
                        if t1[j][3] >= 24:
                            t1[j][3] = t1[j][3] - 24
                            t1[j][2] = t1[j][2] + 1

                        mody = met[i][MetCol.YEAR] - _aint(met[i][MetCol.YEAR] * 2.5e-1) * 4.0
                        if abs(float(t1[j][2]) - 366.0) <= EPS and abs(mody) > EPS:
                            t1[j][2] = 1
                            t1[j][1] = t1[j][1] + 1
                        elif abs(float(t1[j][2]) - 367.0) <= EPS and abs(mody) <= EPS:
                            t1[j][2] = 1
                            t1[j][1] = t1[j][1] + 1

                        f1 = (j * 1e2 / intval) * 1e-2
                        t2[j][1] = ft[i - 1] + (ft[i] - ft[i - 1]) * f1  # low veg temp

                        for k in range(1, NCLAYERS + 1):
                            t2[j][k + 1] = (canopy_temp[k, i - 1]
                                            + (canopy_temp[k, i] - canopy_temp[k, i - 1])
                                            * f1)  # canopy temps

                        t2[j][5] = tt[i - 1] + (tt[i] - tt[i - 1]) * f1  # surface temp
                        t2[j][6] = (met[i - 1][MetCol.TSOIL]
                                    + (met[i][MetCol.TSOIL] - met[i - 1][MetCol.TSOIL])
                                    * f1)  # ground temp
                        t2[j][7] = swater1 + (swater - swater1) * f1  # ground ice + water

                        if (_aint(abs(met[i][MetCol.SD] - mflag) * 1e5) * 1e-5 > EPS
                                and _aint(abs(met[i][MetCol.HI] - mflag) * 1e5) * 1e-5 > EPS):
                            t2[j][8] = ((met[i - 1][MetCol.SD] + met[i - 1][MetCol.HI])
                                        + ((met[i][MetCol.SD] + met[i][MetCol.HI])
                                           - (met[i - 1][MetCol.SD] + met[i - 1][MetCol.HI]))
                                        * f1)  # surface icing + snow depth
                        else:
                            t2[j][8] = 0.0

                        t2[j][9] = airt[1][i - 1] + (airt[1][i] - airt[1][i - 1]) * f1  # air temp above canopy
                        t2[j][10] = airt[2][i - 1] + (airt[2][i] - airt[2][i - 1]) * f1  # air temp below canopy
                        t2[j][11] = surfemisf[i - 1] + (surfemisf[i] - surfemisf[i - 1]) * f1  # low veg emissivity
                        t2[j][12] = surfemisc[i - 1] + (surfemisc[i] - surfemisc[i - 1]) * f1  # canopy emissivity

                        f2 = _fdiv(lheatf[i - 1], lhes[i - 1])
                        t2[j][13] = f2 + (_fdiv(lheatf[i], lhes[i]) - f2) * f1  # low veg cond(+)/evap(-)
                        t2[j][14] = cevap[i - 1] + (cevap[i] - cevap[i - 1]) * f1  # total canopy cond/evap
                        f3 = _fdiv(lheat[i - 1], lhes[i - 1])
                        t2[j][15] = f3 + (_fdiv(lheat[i], lhes[i]) - f3) * f1  # ground cond/evap

                        unit28.write(_fmt_1004(
                            (int(t1[j][1]), int(t1[j][2]), int(t1[j][3]), int(t1[j][4])),
                            vegl_type, vegh_type,
                            (t2[j][1], t2[j][11], t2[j][13], t2[j][2],
                             t2[j][3], t2[j][4], t2[j][12], t2[j][14],
                             t2[j][5], t2[j][6], t2[j][7], t2[j][15],
                             t2[j][8], t2[j][9], t2[j][10])) + "\n")

    # ******************************************************************
    # print the ground flux information -> Does not include vegetation
    # component
    # ******************************************************************
    if water_flag <= 1 and flprint == 1:
        unit24 = state.unit24
        unit24.write(_fmt_1000(freq_id, iend - istartp + 1, 13, lat, mlong,
                               elev, vitd_index, met_count, timeoffset,
                               timstep, mflag) + "\n")
        unit24.write(_HDR_FLUXO_NOMES + "\n")
        unit24.write(_HDR_FLUXO_UNID + "\n")

        for i in range(istartp, iend + 1):
            ttemp = tt[i]
            sd = met[i][MetCol.SD]
            if i != istartp:
                sd = met[i - 1][MetCol.SD]

            # FORMAT 1002: 4(i4),10(f10.2),f10.5,2(f10.2)
            unit24.write(
                _fmt_i(_fint(met[i][MetCol.YEAR]), 4)
                + _fmt_i(_fint(met[i][MetCol.DOY]), 4)
                + _fmt_i(_fint(met[i][MetCol.HR]), 4)
                + _fmt_i(_fint(met[i][MetCol.MIN]), 4)
                + _fmt_f(sdown[i], 10, 2) + _fmt_f(sup[i], 10, 2)
                + _fmt_f(irdown[i], 10, 2) + _fmt_f(irup[i], 10, 2)
                + _fmt_f(lheat[i], 10, 2) + _fmt_f(sheat[i], 10, 2)
                + _fmt_f(pheat1[i], 10, 2) + _fmt_f(cheat[i], 10, 2)
                + _fmt_f(cheat1[i], 10, 2) + _fmt_f(evap_heat[i], 10, 2)
                + _fmt_f(sd, 10, 5) + _fmt_f(tmelt[i], 10, 2)
                + _fmt_f(ttemp, 10, 2) + "\n")
