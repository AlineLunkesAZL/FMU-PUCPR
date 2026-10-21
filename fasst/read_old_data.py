"""
read_old_data.py -- traducao de read_old_data.F90

Fonte Fortran: read_old_data.F90 (287 linhas). Uma unica sub-rotina
publica (`read_old_data`), que le e grava os arquivos binarios usados
para continuar uma simulacao a partir de uma execucao anterior
("hot-start", cenario multi-run). Tem tres modos (`iomode`):

* 1 -- le o estado salvo de uma execucao anterior (unidade 11): procura
  o bloco cujo `freq_id` coincide com o da execucao e, dentro dele, o
  registro numero `ii`; copia todas as variaveis de estado (tempo, site,
  neve e gelo, vegetacao baixa, dossel, nos do perfil, meteorologia
  inicial) para `state`, e devolve `sdensi` e `phie`.
* 2 -- grava o estado atual no arquivo de hot-start (unidade 17) e a
  data/hora do registro no arquivo de carimbos de tempo (unidade 18).
* 3 -- le os carimbos de tempo de uma execucao anterior (unidade 12) em
  `timeo(moverlap,4)` e anota em `moverlapr` quantos registros havia.

Assinatura original:
    subroutine read_old_data(iomode,ii,wstart,sdensi,phie,timeo)
    integer(ip),intent(in):: iomode,ii,wstart
    real(dp),intent(inout):: sdensi,phie
    real(dp),intent(out):: timeo(moverlap,4)

Nenhuma sub-rotina e chamada (o proprio .F90 diz "no subroutines
called").

Arquitetura
-------------
`read_old_data` recebe `state: FasstState` como primeiro argumento. As
constantes (`MAXCOL`, `MAXP`, `MOVERLAP`, `NCLAYERS`) vem de
`fasst.constants`. `sdensi` e `phie` sao `intent(inout)` escalares:
devolvidos em tupla, `sdensi, phie = read_old_data(...)`, igual aos
chamadores ja existentes (`initwater.py`, `initsurface.py`). `timeo` e
`intent(out)` e e um array que o chamador ja cria
(`[[0.0]*5 for _ in range(MOVERLAP+1)]`, 1-based): o modo 3 o preenche
IN-PLACE (unico parametro do projeto tratado assim, para manter a
chamada dos dois chamadores inalterada).

Os argumentos dos chamadores sao `d1i` (= `iomode`) e `state.mpos`
(= `ii`); os stubs que `initwater.py`/`initsurface.py` tinham foram
substituidos por import desta funcao. O driver chama com `d1i = 3` e
`state.dstart`.

Arquivos binarios (formato)
-----------------------------
O driver abre as quatro unidades com `form='unformatted',
access='stream'` (fasst_driver.F90:990-1010): os dados sao gravados em
sequencia, SEM marcadores de registro -- cada `write` so acrescenta os
bytes da lista. `integer(ip)` ocupa 4 bytes, `real(dp)` 8 e
`character(len=n)` n bytes, little-endian. Este modulo emula isso com
`_rd`/`_wr`, conferidos contra o gfortran (TDM-GCC) com um programa de
sondagem:

* leitura que passa do fim do arquivo -> as variaveis que couberam
  inteiras sao atribuidas; a que ficou pela metade recebe os bytes
  disponiveis SOBRE o valor anterior (um real de 8 bytes com so 4
  disponiveis troca os 4 primeiros bytes e mantem os outros 4 do valor
  antigo); as demais ficam como estavam; `iostat` vale -1;
* a execucao continua depois de uma leitura com erro, igual ao Fortran
  (que so olha `iostat` no teste `io /= -1` do laco externo).

As unidades sao os `_Unit` binarios do driver, ja abertos em
`state.unit11`, `state.unit12` (leitura, `'rb'`) e `state.unit17`,
`state.unit18` (escrita, `'wb'`) -- conferido contra `fasst_driver.py`.

Pontos de atencao (preservados do Fortran, nao corrigidos)
------------------------------------------------------------
* Modo 1, ramo "outro freq_id" INALCANCAVEL: depois do bloco
  `if(ijunk == freq_id)` vem `else if(ijunk == freq_id)` -- a mesma
  condicao, nunca verdadeira. Pela estrutura, era para ser um `else`
  que pulasse os blocos de outras execucoes ate `'endr'`. Na pratica,
  quando o primeiro registro lido nao e o `freq_id` procurado, o laco
  apenas le o registro seguinte e o interpreta como inteiro (lixo, com
  a primeira palavra de cada registro). O ramo morto foi mantido como
  `elif` equivalente, ligado a `_pular_bloco`.
* Modo 3: o ramo "outro freq_id" (pular blocos ate `'endr'`) e
  alcancavel, e `jline` NAO e reinicializado entre blocos: depois de
  pular um bloco, `jline` fica em `'endr'`, entao, se o bloco procurado
  for o segundo ou posterior do arquivo, o laco interno nao executa e
  `timeo` nao e lido. Quando o primeiro registro nao e o `freq_id`
  procurado (inclusive arquivo vazio) e o arquivo termina sem `'endr'`,
  o Fortran fica em laco infinito (leitura em fim de arquivo nao altera
  `jline`); aqui levanta `EOFError` (unico desvio de comportamento do
  modulo). Com `freq_id` casando e sem `'endr'` o laco acaba apos
  `MOVERLAP` registros, igual ao Fortran.
* Modo 2 grava a data/hora (unidade 18) em DOIS `write` -- quatro
  reais e depois `'loop'`/`'endr'` -- e o modo 3 le `timeo(i,1:4)` e
  `jline` num `read` so. Em `access='stream'` isso e consistente (36
  bytes em sequencia nos dois lados); so seria defeito em acesso
  sequencial, que o driver nao usa.
* Modo 2 grava `canopy_temp(j,ii)`; o modo 1 le para
  `canopy_temp(j,oldpos)` (posicoes diferentes).
* `avect(0:4,nclayers)`: acessado aqui como `state.avect[i, j]` com
  `j = 1..nclayers`, como em `module_canopy.py`. Exigia `state.avect` com
  forma `(5, NCLAYERS+1)`; o `state.py` declarava `(5, NCLAYERS)` e foi
  corrigido junto com esta integracao.
* `km` e `sphm` sao vetores de tamanho `EXTRAN` (3), lidos/gravados
  inteiros no meio da lista de variaveis de neve (17 reais no registro,
  e nao 13). ASSUNCAO: declarados `(1:extran)` em `fasst_global` (o
  `fasst_global.F90` nao foi fornecido); se forem `(0:extran)` o registro
  tem 19 reais e o layout do arquivo muda.
* DEFEITO DO ORIGINAL, o mais grave do arquivo: o ramo que PULA um
  registro (`kk /= ii`, linha 118, e o ramo morto, 155) le 13 reais de
  neve, mas o modo 2 grava 17 (`hi..vimelt` sao 8, mais `km(3)`,
  `sphm(3)`, `atopf`, `refreeze`, `phie`). O arquivo desalinha 32 bytes,
  `nnodes`/`ntot` saem de lixo e os lacos de nos leem no fim do arquivo
  indefinidamente. Como o laco `kk` continua DEPOIS do registro certo ate
  achar 'endr', isso pega todo hot start cujo bloco tenha mais de um
  registro -- ou seja, quase todos (`open_water` grava os ultimos
  `moverlap` passos). Conferido com o FASST compilado: o Fortran trava, e
  esta traducao trava no mesmo ponto (`_pular_bloco`). Preservado. Com um
  bloco de UM registro o hot start funciona e foi conferido byte a byte
  (`cmp_fortran_agua.py`, caso `um_hot`).
* Modo 1, `freq_id` que nao casa: sem marcadores, o laco externo le os
  4 bytes SEGUINTES como `ijunk` -- avanca pelo arquivo de 4 em 4 bytes
  ate achar uma palavra igual a `freq_id` (que pode estar no meio de um
  real) ou o fim. Preservado.
* Registro de tempo (36 bytes) lido incompleto no fim do arquivo:
  `nnodes`/`ntot` ficam com o valor antigo. Nao foi tratado.
* No ramo de "pular" os registros, `nnodes` e `ntot` sao lidos para o
  estado mesmo quando `kk /= ii` (o proprio Fortran os sobrescreve).
"""

import struct

from .constants import MAXCOL, MAXP, MOVERLAP, NCLAYERS
from .indices import MetCol
from .state import FasstState

_I = "i"      # integer(ip), 4 bytes
_D = "d"      # real(dp), 8 bytes
_A2 = "a2"    # character(len=2)
_A4 = "a4"    # character(len=4)
_SIZE = {_I: 4, _D: 8, _A2: 2, _A4: 4}


class _P:
    """Item lido so em parte: os bytes que sobraram no registro."""
    __slots__ = ("kind", "raw")

    def __init__(self, kind, raw):
        self.kind = kind
        self.raw = raw

    def merge(self, old):
        """Sobrepoe os bytes disponiveis ao valor anterior de `old`."""
        k = self.kind
        if k == _D:
            b = bytearray(struct.pack("<d", float(old)))
        elif k == _I:
            b = bytearray(struct.pack("<i", int(old)))
        else:
            n = _SIZE[k]
            b = bytearray(str(old).ljust(n)[:n].encode("latin-1"))
        b[:len(self.raw)] = self.raw
        if k == _D:
            return struct.unpack("<d", bytes(b))[0]
        if k == _I:
            return struct.unpack("<i", bytes(b))[0]
        return bytes(b).decode("latin-1")


def _m(v, old):
    """Valor a atribuir: `old` se nada foi lido, o item (ou a fusao dos
    bytes parciais com `old`) caso contrario."""
    if v is None:
        return old
    if isinstance(v, _P):
        return v.merge(old)
    return v


def _rd(f, spec):
    """`read(unit,iostat=io) lista` em `access='stream'` -- (valores, io).

    `valores` tem um item por codigo de `spec`: o valor lido, `_P` se o
    arquivo acabou no meio do item, ou None se nada foi lido.
    io: 0 ok, -1 fim de arquivo (mesmo com leitura parcial).
    """
    rec = f.read(sum(_SIZE[c] for c in spec))
    out = []
    pos = 0
    for k, c in enumerate(spec):
        sz = _SIZE[c]
        if pos + sz > len(rec):
            out.append(_P(c, rec[pos:]) if pos < len(rec) else None)
            out.extend([None] * (len(spec) - k - 1))
            return out, -1
        if c == _I:
            out.append(struct.unpack_from("<i", rec, pos)[0])
        elif c == _D:
            out.append(struct.unpack_from("<d", rec, pos)[0])
        else:
            out.append(rec[pos:pos + sz].decode("latin-1"))
        pos += sz
    return out, 0


def _wr(f, items):
    """`write(unit) lista` em `access='stream'` -- so os bytes, sem marcador."""
    parts = []
    for c, v in items:
        if c == _I:
            parts.append(struct.pack("<i", int(v)))
        elif c == _D:
            parts.append(struct.pack("<d", float(v)))
        else:
            n = _SIZE[c]
            parts.append(str(v).ljust(n)[:n].encode("latin-1"))
    f.write(b"".join(parts))


def _put(state, names, vals):
    """Atribui os valores lidos aos campos escalares de `state`."""
    for nm, v in zip(names, vals):
        if nm is not None:
            setattr(state, nm, _m(v, getattr(state, nm)))


def _put_arr(arr, idx, vals):
    for i, v in zip(idx, vals):
        arr[i] = _m(v, arr[i])


_SCAL_LOW = ["isigfl", "ihfol", "iepf", "ifola", "ilail", "ftemp", "storll",
             "storls", "hfol", "hfol_tot", "lail", "epf", "sigfl", "albf",
             "fola"]


def _pular_bloco(state, f):
    """Le (e descarta) um bloco de estado -- os `read` do ramo
    `kk /= ii` (linhas 113-145) e do ramo morto (linhas 150-184).
    Devolve (jline, io); `jline` e None se o ultimo registro nao foi lido.
    `nnodes` e `ntot` sao atribuidos ao estado, como no Fortran."""
    v, io = _rd(f, [_D, _I, _D, _I, _I, _I, _I])
    _put(state, [None, None, None, None, "nnodes", "ntot", None], v)
    _, io = _rd(f, [_D] * 8)
    _, io = _rd(f, [_D] * 13)
    _, io = _rd(f, [_D] * 15)
    _, io = _rd(f, [_I, _I])
    _, io = _rd(f, [_I, _I, _I] + [_D] * 15)
    _, io = _rd(f, [_I, _I, _D, _D])
    for _ in range(NCLAYERS):
        _, io = _rd(f, [_D] * 10)
        _, io = _rd(f, [_D] * 5)
    _, io = _rd(f, [_D] * (NCLAYERS + NCLAYERS))
    for _ in range(state.nnodes):
        _, io = _rd(f, [_A2])
        _, io = _rd(f, [_I] + [_D] * 11)
        _, io = _rd(f, [_D] * MAXP)
    for _ in range(state.ntot):
        _, io = _rd(f, [_D] * 12)
        _, io = _rd(f, [_D] * 5)
    _, io = _rd(f, [_D] * 13)
    _, io = _rd(f, [_D] * MAXCOL)
    v, io = _rd(f, [_A4])
    return v[0], io


def _ler_bloco(state, f, sdensi, phie):
    """Le o bloco de estado do registro escolhido (linhas 57-111).
    Devolve (sdensi, phie, jline, io)."""
    # time variables
    v, io = _rd(f, [_D, _I, _D, _I, _I, _I, _I])
    _put(state, ["deltati", "stepi", "deltat_fasst", "step", "nnodes",
                 "ntot", "refn"], v)

    # site variables, surface temperature
    v, io = _rd(f, [_D] * 8)
    _put(state, ["albedo_fasst", "emis", "sgralbedo", "sgremis", "hpond",
                 "gwl", "isurfoldg", "isurfoldf"], v)

    # snow and ice variables
    v, io = _rd(f, [_D] * 17)
    _put(state, ["hi", "hsaccum"], v[0:2])
    sdensi = _m(v[2], sdensi)
    _put(state, ["toptemp", "ptemp", "iheightn", "vsmelt", "vimelt"], v[3:8])
    _put_arr(state.km, (1, 2, 3), v[8:11])
    _put_arr(state.sphm, (1, 2, 3), v[11:14])
    _put(state, ["atopf", "refreeze"], v[14:16])
    phie = _m(v[16], phie)
    v, io = _rd(f, [_D] * 15)
    _put_arr(state.sn_stat, range(1, 16), v)
    v, io = _rd(f, [_I, _I])
    _put_arr(state.sn_istat, (1, 2), v)

    # low vegetation variables
    v, io = _rd(f, [_I, _I, _I] + [_D] * 15)
    _put(state, ["veg_flagl", "vegl_type", "iseason"] + _SCAL_LOW, v)

    # canopy variables
    v, io = _rd(f, [_I, _I, _D, _D])
    _put(state, ["veg_flagh", "vegh_type", "zh", "isigfh"], v)
    for j in range(1, NCLAYERS + 1):
        v, io = _rd(f, [_D] * 10)
        for arr, x in zip((state.ifoliage_type, state.iclump, state.irho,
                           state.itau, state.ialp, state.ieps, state.ilai,
                           state.laif, state.dzveg), v[0:9]):
            arr[j] = _m(x, arr[j])
        ct = state.canopy_temp
        ct[j, state.oldpos] = _m(v[9], ct[j, state.oldpos])
        v, io = _rd(f, [_D] * 5)
        for i in range(0, 5):
            state.avect[i, j] = _m(v[i], state.avect[i, j])
    v, io = _rd(f, [_D] * (NCLAYERS + NCLAYERS))
    _put_arr(state.stor, range(1, NCLAYERS + NCLAYERS + 1), v)

    # nodal properties
    for i in range(1, state.nnodes + 1):
        v, io = _rd(f, [_A2])
        state.node_type[i] = _m(v[0], state.node_type[i])
        v, io = _rd(f, [_I] + [_D] * 11)
        state.ntype[i] = _m(v[0], state.ntype[i])
        for arr, x in zip((state.soil_moist, state.nzi, state.nz,
                           state.delzsi, state.delzs, state.grthcond,
                           state.grspheat, state.frl, state.frh,
                           state.sinkr, state.pheadmin), v[1:12]):
            arr[i] = _m(x, arr[i])
        v, io = _rd(f, [_D] * MAXP)
        for j in range(1, MAXP + 1):
            state.nsoilp[i][j] = _m(v[j - 1], state.nsoilp[i][j])

    for i in range(1, state.ntot + 1):
        v, io = _rd(f, [_D] * 12)
        for arr, x in zip((state.stt, state.phead, state.ice, state.wvc,
                           state.vin, state.flowu, state.flowl, state.fv1,
                           state.source, state.sink, state.khl, state.khu), v):
            arr[i] = _m(x, arr[i])
        v, io = _rd(f, [_D] * 5)
        for arr, x in zip((state.too, state.phoo, state.ioo, state.woo,
                           state.smoo), v):
            arr[i] = _m(x, arr[i])

    # initial met information
    v, io = _rd(f, [_D] * 13)
    for j in range(1, 14):
        state.dmet1[state.oldpos][j] = _m(v[j - 1], state.dmet1[state.oldpos][j])
    v, io = _rd(f, [_D] * MAXCOL)
    for j in range(1, MAXCOL + 1):
        state.met[state.oldpos][j] = _m(v[j - 1], state.met[state.oldpos][j])

    v, io = _rd(f, [_A4])
    return sdensi, phie, v[0], io



def read_old_data(state: FasstState, iomode, ii, wstart, sdensi, phie, timeo):
    """Le/grava os arquivos de hot-start. Devolve `(sdensi, phie)`.

    `timeo` (1-based, `[MOVERLAP+1][5]`) e preenchido in-place no modo 3.
    """
    # ******************************************************************
    if iomode == 1:  # read in the old data
        f = state.unit11
        ijunk = -1
        io = 0
        jline = "    "

        while io != -1 and ijunk != state.freq_id:
            v, io = _rd(f, [_I])
            ijunk = _m(v[0], ijunk)

            if ijunk == state.freq_id:
                kk = 1

                while jline != "endr" and kk <= MOVERLAP:
                    if kk == ii:
                        sdensi, phie, nj, io = _ler_bloco(state, f, sdensi, phie)
                    else:  # not the correct location
                        nj, io = _pular_bloco(state, f)
                    jline = _m(nj, jline)
                    kk = kk + 1
            elif ijunk == state.freq_id:  # pragma: no cover -- INALCANCAVEL
                while jline[0:4] != "endr":
                    nj, io = _pular_bloco(state, f)
                    jline = _m(nj, jline)
                    if io == -1:
                        raise EOFError("read_old_data: fim de arquivo sem 'endr'")

    # ******************************************************************
    elif iomode == 2:  # write data to file used in "hot-start" runs
        f = state.unit17

        # write soil properties for hot run
        if ii == wstart:
            _wr(f, [(_I, state.freq_id)])

        # time variables
        _wr(f, [(_D, state.deltati), (_I, state.stepi),
                (_D, state.deltat_fasst), (_I, state.step),
                (_I, state.nnodes), (_I, state.ntot), (_I, state.refn)])

        # site variables, surface temperature
        _wr(f, [(_D, state.albedo_fasst), (_D, state.emis),
                (_D, state.sgralbedo), (_D, state.sgremis),
                (_D, state.hpond), (_D, state.gwl),
                (_D, state.isurfoldg), (_D, state.isurfoldf)])

        # snow and ice variables
        _wr(f, [(_D, state.hi), (_D, state.hsaccum), (_D, sdensi),
                (_D, state.toptemp), (_D, state.ptemp),
                (_D, state.iheightn), (_D, state.vsmelt),
                (_D, state.vimelt)]
               + [(_D, state.km[k]) for k in (1, 2, 3)]
               + [(_D, state.sphm[k]) for k in (1, 2, 3)]
               + [(_D, state.atopf), (_D, state.refreeze), (_D, phie)])
        _wr(f, [(_D, state.sn_stat[k]) for k in range(1, 16)])
        _wr(f, [(_I, state.sn_istat[1]), (_I, state.sn_istat[2])])

        # low vegetation variables
        _wr(f, [(_I, state.veg_flagl), (_I, state.vegl_type),
                (_I, state.iseason)]
               + [(_D, getattr(state, nm)) for nm in _SCAL_LOW])

        # canopy variables
        _wr(f, [(_I, state.veg_flagh), (_I, state.vegh_type),
                (_D, state.zh), (_D, state.isigfh)])
        for j in range(1, NCLAYERS + 1):
            _wr(f, [(_D, state.ifoliage_type[j]), (_D, state.iclump[j]),
                    (_D, state.irho[j]), (_D, state.itau[j]),
                    (_D, state.ialp[j]), (_D, state.ieps[j]),
                    (_D, state.ilai[j]), (_D, state.laif[j]),
                    (_D, state.dzveg[j]), (_D, state.canopy_temp[j, ii])])
            _wr(f, [(_D, state.avect[i, j]) for i in range(0, 5)])
        _wr(f, [(_D, state.stor[j]) for j in range(1, NCLAYERS + NCLAYERS + 1)])

        # node properties
        for i in range(1, state.nnodes + 1):
            _wr(f, [(_A2, state.node_type[i])])
            _wr(f, [(_I, state.ntype[i]), (_D, state.soil_moist[i]),
                    (_D, state.nzi[i]), (_D, state.nz[i]),
                    (_D, state.delzsi[i]), (_D, state.delzs[i]),
                    (_D, state.grthcond[i]), (_D, state.grspheat[i]),
                    (_D, state.frl[i]), (_D, state.frh[i]),
                    (_D, state.sinkr[i]), (_D, state.pheadmin[i])])
            _wr(f, [(_D, state.nsoilp[i][j]) for j in range(1, MAXP + 1)])

        for i in range(1, state.ntot + 1):
            _wr(f, [(_D, state.stt[i]), (_D, state.phead[i]),
                    (_D, state.ice[i]), (_D, state.wvc[i]),
                    (_D, state.vin[i]), (_D, state.flowu[i]),
                    (_D, state.flowl[i]), (_D, state.fv1[i]),
                    (_D, state.source[i]), (_D, state.sink[i]),
                    (_D, state.khl[i]), (_D, state.khu[i])])
            _wr(f, [(_D, state.too[i]), (_D, state.phoo[i]),
                    (_D, state.ioo[i]), (_D, state.woo[i]),
                    (_D, state.smoo[i])])

        # final met properties
        _wr(f, [(_D, state.dmet1[ii][j]) for j in range(1, 14)])
        _wr(f, [(_D, state.met[ii][j]) for j in range(1, MAXCOL + 1)])

        if ii != state.iend:
            _wr(f, [(_A4, "loop")])
        else:
            _wr(f, [(_A4, "endr")])

        # write to the time stamp file
        g = state.unit18
        if ii == wstart:
            _wr(g, [(_I, state.freq_id)])
        _wr(g, [(_D, state.met[ii][MetCol.YEAR]), (_D, state.met[ii][MetCol.DOY]),
                (_D, state.met[ii][MetCol.HR]), (_D, state.met[ii][MetCol.MIN])])

        if ii != state.iend:
            _wr(g, [(_A4, "loop")])
        else:
            _wr(g, [(_A4, "endr")])

    # ******************************************************************
    elif iomode == 3:  # read in time stamp of the old data
        f = state.unit12
        ijunk = -1
        io = 0
        jline = "    "

        while io != -1 and ijunk != state.freq_id:
            v, io = _rd(f, [_I])
            ijunk = _m(v[0], ijunk)
            if ijunk == state.freq_id:
                i = 1
                while jline[0:4] != "endr" and i <= MOVERLAP:
                    v, io = _rd(f, [_D, _D, _D, _D, _A4])
                    for k in range(4):
                        timeo[i][k + 1] = _m(v[k], timeo[i][k + 1])
                    jline = _m(v[4], jline)
                    if jline == "endr":
                        state.moverlapr = i
                    i = i + 1
            else:
                while jline[0:4] != "endr":
                    v, io = _rd(f, [_D, _D, _D, _D, _A4])
                    jline = _m(v[4], jline)
                    if io == -1:
                        raise EOFError("read_old_data: fim de arquivo sem 'endr' "
                                       "(o Fortran entraria em laco infinito)")

    return sdensi, phie
