"""
read_complex.py -- traducao de read_complex.F90

Fonte Fortran: read_complex.F90 (450 linhas), uma unica sub-rotina
(`read_complex`). Le UMA linha de entrada da unidade 5 (leitura em lista
livre) com a descricao de um poligono -- terreno, agua aberta, vegetacao
baixa e alta, rugosidade, temperatura de superficie, nivel d'agua e as
camadas de solo -- e converte isso nas variaveis de estado do modelo
(`slope_fasst`, `aspect`, `water_flag`, `vegl_type`, `soiltype`,
`lthick`, `rho_fac`, ...).

Assinatura original:
    subroutine read_complex(wtype,io,lcount,nm0,wvel,wdepth,mtest,sname)
    integer(ip),intent(inout):: wtype,io
    integer(ip),intent(out):: lcount,nm0
    real(dp),intent(inout):: wvel,wdepth
    character(len=1),intent(out):: mtest
    character(len=4),intent(out):: sname(maxl)

Arquitetura
-------------
`read_complex(state, wtype, io, wvel, wdepth, sname)` recebe `state:
FasstState` como primeiro argumento. Os argumentos `intent(inout)`
(`wtype`, `io`, `wvel`, `wdepth`) entram como parametros e saem na tupla
de retorno, junto com os `intent(out)` escalares:

    wtype, io, lcount, nm0, wvel, wdepth, mtest = read_complex(
        state, wtype, io, wvel, wdepth, sname)

`sname` (`intent(out)`, array de 4 caracteres, 1-based com a posicao 0
fantasma) e preenchido IN-PLACE: so os elementos das camadas de solo
definido pelo usuario (`soilt = -1`) sao atribuidos; os demais ficam
como o chamador os deixou (como no Fortran). `lcount` vale sempre 0 (o
Fortran so o zera). As variaveis `intent(out)` nao sao entrada.

Constantes (`EPS`, `PI`, `SPFLAG`, `TREF`, `MAXL`, `MAXN`) vem de
`fasst.constants`. `state.unit5` e um objeto de texto com `readline()`,
ja aberto (a convencao `state.unitNN` segue a dos demais modulos e nao foi
confirmada contra o driver).

Leitura em lista livre (`read(5,*,iostat=io) ...`)
----------------------------------------------------
A lista tem 20 itens inteiros seguidos de `nlayers` triplas (tipo de solo,
espessura em cm, nome de 4 caracteres). `nlayers` e lido NA PROPRIA lista
e controla o laco implicito da mesma instrucao. O Fortran nao confere
`io` depois da leitura: o que foi lido antes de um erro (ou de um fim de
arquivo) fica atribuido e o resto do codigo roda mesmo assim.

A classe `_ListaLivre` emula o `read(*,*)` do gfortran (comparada com o
executavel em milhares de entradas, ver o documento):

* itens separados por brancos, virgulas ou fins de linha (a lista continua
  na linha seguinte; uma virgula no comeco da linha seguinte e o
  separador do item anterior);
* valor nulo (`,,`: o item fica como esta); `/` encerra a lista (itens
  restantes ficam como estao); nomes com aspas simples ou duplas (aspas
  duplicadas escapam) ou sem aspas (termina em branco, virgula ou `/`),
  truncados/completados com brancos para 4 caracteres;
* ao terminar a instrucao o resto da linha e descartado;
* `iostat`: 0; -1 em fim de arquivo; 5010 para valor invalido (token nao
  inteiro como `1.5` ou `abc`, ou texto colado ao fecha-aspas), caso em
  que tambem se descarta o resto da linha; 5001 se a unidade ja estava em
  fim de arquivo e e lida de novo;
* um `/` como ULTIMO caractere do arquivo, sem fim de linha, vira fim de
  arquivo (-1), mesmo com todos os itens ja lidos.

Limites (nao emulados com precisao): repeticoes `r*c` (so quando todos os
itens repetidos tem o mesmo tipo), contagem de repeticao zero, inteiro
fora de 32 bits, aspas sem fechar, `+`/`-` soltos no fim do arquivo e uma
primeira linha em branco seguida de virgula. Esses casos nao aparecem em
uma linha de entrada bem formada.

sort2 (dependencia nao fornecida)
-----------------------------------
O cabecalho do .F90 diz "no subroutines called", mas o ultimo ramo
(`nlayers1 > 0`, linhas 397-447) chama `sort2(n, a, b, m, a_out, b_out)`,
que nao esta em nenhum arquivo recebido (`sort_locate.py` so tem `sort`
e `locate` provisorios). Aqui `sort2` e um stub explicito que levanta
`NotImplementedError`. A assinatura assumida -- `sort2(n, a, b, a_out,
b_out)` preenche `a_out`/`b_out` in-place e devolve `m` -- e so a forma de
chamada; nada sobre o que a funcao faz foi suposto.

Pontos de atencao (preservados do Fortran, nao corrigidos)
------------------------------------------------------------
* Linha 106: `lthick1(2) = wdepth - lthick(1)` usa `lthick(1)` (estado,
  possivelmente de uma leitura anterior) e nao `lthick1(1)` (0,02),
  diferente das linhas 143, 159, 195. O mesmo vale para a linha 163
  (`wdepth = lthick(1) + lthick1(2)`).
* Linha 154: no ramo `vegl_type == 14` atribui `nlayers = 3` (e nao
  `nlayers1`), mas deixa `soiltype`/`lthick` por preencher a partir de
  `soiltype1`/`lthick1`; o laco de camadas seguinte le `soilt`/`ithick`
  zerados.
* Linhas 276-281: os ramos `soiltype(i) == -2/-3/-4` testam o ESTADO
  `soiltype(i)` (valor anterior) e nao `soilt(i)`; so sao alcancaveis com
  `soiltype` de uma execucao anterior.
* Ramos inalcancaveis: `wflag2 == 1 .and. nlayers <= 1` (linha 389) nunca
  e verdadeiro porque `wflag2` so e atribuido quando `nlayers > 0`, e o
  primeiro ramo ja trata isso.
* Ultimo ramo (linhas 397-447): usa `lthick2(i)` sem inicializar (o
  Fortran le lixo de pilha); aqui e zero (assuncao). Os indices
  `soiltype(i-nlayers)` e `lthick(i-nlayers1)` da condicao e da
  atribuicao divergem; com `nlayers <= 0` essas condicoes nunca valem. O
  indice `k = i - j` pode ser 0 ou negativo (escrita fora dos limites no
  Fortran); aqui, 0 usa a posicao fantasma e negativo levanta
  `IndexError`.
* `nlayers > maxl` escreve fora de `soilt`/`ithick`/`isname`/`lthick` no
  Fortran; aqui levanta `IndexError`.
* A saida `sname` so e preenchida com `soilt = -1`.
"""

import math
import re
import weakref

from .constants import EPS, MAXL, MAXN, PI, SPFLAG, TREF
from .state import FasstState

_NT = MAXN + 3 + 1   # tamanho 1-based (posicao 0 fantasma) dos vetores soiltype2..4, lthick2..4


# ---------------------------------------------------------------------------
# Helpers numericos locais (mesma convencao de initwater.py)
# ---------------------------------------------------------------------------
def _aint(x):
    """Fortran AINT(x): trunca em direcao a zero."""
    if not math.isfinite(x):
        return x
    return math.copysign(math.floor(abs(x)), x)


def _anint(x):
    """Fortran ANINT(x): arredonda para longe de zero."""
    if not math.isfinite(x):
        return x
    return math.copysign(math.floor(abs(x) + 0.5), x)


# ---------------------------------------------------------------------------
# Leitura em lista livre, emulando read(5,*,iostat=io)
# ---------------------------------------------------------------------------
_NULO = object()          # item nulo: a variavel fica como esta
_ERRO = object()
_RE_REP = re.compile(r"(\d+)\*")
# unidades que ja chegaram ao fim de arquivo: novo read -> iostat 5001 (gfortran)
_APOS_EOF = weakref.WeakSet()
_INT32 = 2**31


class _ListaLivre:
    """Emula uma instrucao `read(unit,*,iostat=io) lista` do gfortran.

    Uso: `r.inteiro(velho)` / `r.texto4(velho)` devolvem o novo valor do
    item (ou `velho` se o item for nulo, se a lista foi encerrada por `/`
    ou se ja houve erro); `r.io` e 0, -1 (fim de arquivo), 5001 (novo read
    depois do fim de arquivo) ou 5010 (valor invalido); `r.fim()` descarta o resto da linha atual.
    """

    def __init__(self, unit):
        self.unit = unit
        self.linha = None      # registro corrente (sem o fim de linha)
        self.pos = 0
        self.io = 0
        self.fechada = False   # '/' encontrado
        self.rep = 0           # repeticoes pendentes
        self.rep_bruto = None
        self.rep_tipo = None
        self.sem_nl = False    # ultima linha sem fim de linha
        self.n = 0             # itens ja processados
        self.nulo_virgula = False   # item anterior foi nulo, encerrado por virgula
        try:
            if unit in _APOS_EOF:
                self.io = 5001   # "READ not allowed after EOF marker"
        except TypeError:
            pass

    # -- registros ---------------------------------------------------------
    def _proxima_linha(self):
        s = self.unit.readline()
        if s == "":
            return False
        self.sem_nl = not s.endswith("\n")
        self.linha = s.rstrip("\r\n")
        self.pos = 0
        return True

    def _pula_brancos(self):
        """Avanca ate o proximo caractere util; False em fim de arquivo."""
        while True:
            if self.linha is None or self.pos >= len(self.linha):
                if not self._proxima_linha():
                    return False
                continue
            if self.linha[self.pos] in " \t":
                self.pos += 1
                continue
            return True

    def _separador_antes(self, consome_virgula=True):
        """Separador que segue o item anterior: brancos e fins de linha
        (a lista continua nas linhas seguintes) e, no maximo, uma virgula."""
        while True:
            if self.linha is None or self.pos >= len(self.linha):
                if not self._proxima_linha():
                    return
                continue
            if self.linha[self.pos] in " \t":
                self.pos += 1
                continue
            break
        if consome_virgula and self.linha[self.pos] == ",":
            self.pos += 1

    def _token(self):
        n = len(self.linha)
        i = self.pos
        while i < n and self.linha[i] not in " \t,/":
            i += 1
        tok = self.linha[self.pos:i]
        self.pos = i
        return tok

    # -- itens ---------------------------------------------------------------
    def _bruto(self):
        """Le o proximo valor da lista como texto: (texto, entre_aspas),
        None (valor nulo) ou _ERRO. Trata o prefixo de repeticao `r*`."""
        m = _RE_REP.match(self.linha, self.pos)
        reps = 1
        if m:
            reps = int(m.group(1))
            if reps == 0:
                self.io = 5010
                return _ERRO
            self.pos = m.end()
            n = len(self.linha)
            if self.pos >= n or self.linha[self.pos] in " \t,/":
                bruto = None                         # `r*` = r valores nulos
            else:
                bruto = self._valor()
        else:
            bruto = self._valor()
        if bruto is _ERRO:
            return _ERRO
        self.rep = reps - 1
        self.rep_bruto = bruto
        return bruto

    def _valor(self):
        c = self.linha[self.pos]
        if c in "'\"":
            s = self._aspas(c)
            return _ERRO if s is None else (s, True)
        return (self._token(), False)

    def _converte(self, tipo, bruto):
        """Bruto -> valor do item (ou _NULO); registra erro em io."""
        if bruto is None:
            return _NULO
        txt, aspas = bruto
        if tipo == "a":
            return txt
        if aspas:
            self.io = 5010
            return _NULO
        t = txt[1:] if txt[:1] in "+-" else txt
        if not (t.isdigit() and t.isascii()):
            self.io = 5010
            return _NULO
        v = int(txt)
        if not (-_INT32 <= v < _INT32):
            self.io = 5010
            return _NULO
        return v

    def _item(self, tipo, velho):
        if self.io != 0 or self.fechada:
            return velho
        try:
            return self._item2(tipo, velho)
        finally:
            self.n += 1

    def _item2(self, tipo, velho):
        if self.rep > 0:
            self.rep -= 1
            if self.rep_bruto is not None and tipo != self.rep_tipo:
                self.io = 5010      # constante repetida de outro tipo
                return velho
            v = self._converte(tipo, self.rep_bruto)
        else:
            termina_virgula = self.nulo_virgula   # o item anterior ja consumiu a sua virgula
            self.nulo_virgula = False
            if self.n > 0:
                self._separador_antes(not termina_virgula)
            if not self._pula_brancos():
                self.io = -1
                try:
                    _APOS_EOF.add(self.unit)
                except TypeError:
                    pass
                return velho
            c = self.linha[self.pos]
            if c == ",":
                self.pos += 1
                self.nulo_virgula = True
                return velho
            if c == "/":
                if self._barra_no_fim():
                    return velho
                self.fechada = True
                return velho
            bruto = self._bruto()
            if bruto is _ERRO:
                return velho
            self.rep_tipo = tipo
            v = self._converte(tipo, bruto)
        if v is _NULO:
            return velho
        return v.ljust(4)[:4] if tipo == "a" else v

    def inteiro(self, velho):
        return self._item("i", velho)

    def texto4(self, velho):
        return self._item("a", velho)

    def _aspas(self, q):
        """String delimitada por aspas, com aspas duplicadas como escape."""
        self.pos += 1
        out = []
        while True:
            if self.pos >= len(self.linha):
                if not self._proxima_linha():
                    self.io = 5010          # aspas sem fechar ate o fim do arquivo
                    return None
                continue
            c = self.linha[self.pos]
            if c == q:
                if self.pos + 1 < len(self.linha) and self.linha[self.pos + 1] == q:
                    out.append(q)
                    self.pos += 2
                    continue
                self.pos += 1
                if self.pos < len(self.linha) and self.linha[self.pos] not in " \t,/":
                    self.io = 5010              # caractere colado ao fecha-aspas
                    return None
                return "".join(out)
            out.append(c)
            self.pos += 1

    def _barra_no_fim(self):
        """gfortran: '/' como ULTIMO caractere do arquivo, sem fim de linha,
        nao encerra a lista: vira fim de arquivo (iostat -1)."""
        if (self.sem_nl and self.linha is not None
                and self.linha[self.pos:].lstrip(" \t") == "/"):
            self.io = -1
            try:
                _APOS_EOF.add(self.unit)
            except TypeError:
                pass
            return True
        return False

    def fim(self):
        """Fim da instrucao: o resto do registro corrente e descartado."""
        if self.io == 0 and self.linha is not None:
            self._barra_no_fim()
        self.linha = None


# ---------------------------------------------------------------------------
# Dependencia externa nao fornecida
# ---------------------------------------------------------------------------
def sort2(n, a, b, a_out, b_out):
    """TODO: `sort2` nao esta em nenhum arquivo recebido (ver docstring)."""
    raise NotImplementedError(
        "sort2 nao foi fornecida (read_complex.F90, linhas 413 e 431)")


def _grava(arr, k, valor):
    """arr[k] = valor; indice negativo seria escrita fora dos limites no Fortran."""
    if k < 0:
        raise IndexError(
            f"read_complex: escrita fora dos limites (indice {k}); comportamento "
            "indefinido no Fortran")
    arr[k] = valor


def read_complex(state: FasstState, wtype, io, wvel, wdepth, sname):
    """Le uma linha da unidade 5 e monta o estado do poligono.

    Devolve `(wtype, io, lcount, nm0, wvel, wdepth, mtest)`; `sname` e
    preenchido in-place (ver docstring do modulo).
    """
    # ------------------------------------------------------------------
    # initialize variables (linhas 36-63)
    lcount = 0
    state.freq_id = 0
    state.vitd_index = 0
    islope = 0
    iaspect = 0
    ielev = 0
    wtype = 0
    iwdepth = 0
    iwvel = 0
    state.vegl_type = 0
    ivegldens = 0
    iveglht = 0
    state.vegh_type = 0
    iveghdens = 0
    iveghht = 0
    r_type = 0
    iswater = 0
    irough = 0
    stemp = 0
    igwl = 0
    state.nlayers = 0
    soilt = [0] * (MAXL + 1)
    ithick = [0] * (MAXL + 1)
    wflag2 = 0
    nlayers1 = 0
    isname = ["    "] * (MAXL + 1)

    soiltype1 = [0] * 4
    lthick1 = [0.0] * 4
    soiltype2 = [0] * _NT
    soiltype3 = [0] * _NT
    soiltype4 = [0] * _NT
    lthick2 = [0.0] * _NT   # no Fortran, sem inicializar (aqui: zeros)
    lthick3 = [0.0] * _NT
    lthick4 = [0.0] * _NT

    # ------------------------------------------------------------------
    # read(5,*,iostat=io) ... (linhas 65-69); atribui item a item
    r = _ListaLivre(state.unit5)
    state.freq_id = r.inteiro(state.freq_id)
    state.vitd_index = r.inteiro(state.vitd_index)
    islope = r.inteiro(islope)
    iaspect = r.inteiro(iaspect)
    ielev = r.inteiro(ielev)
    wtype = r.inteiro(wtype)
    iwdepth = r.inteiro(iwdepth)
    iwvel = r.inteiro(iwvel)
    state.vegl_type = r.inteiro(state.vegl_type)
    ivegldens = r.inteiro(ivegldens)
    iveglht = r.inteiro(iveglht)
    state.vegh_type = r.inteiro(state.vegh_type)
    iveghdens = r.inteiro(iveghdens)
    iveghht = r.inteiro(iveghht)
    r_type = r.inteiro(r_type)
    iswater = r.inteiro(iswater)
    irough = r.inteiro(irough)
    stemp = r.inteiro(stemp)
    igwl = r.inteiro(igwl)
    state.nlayers = r.inteiro(state.nlayers)
    for i in range(1, state.nlayers + 1):
        if r.io != 0 or r.fechada:
            break   # nada mais sera atribuido (e evita indices fora dos limites)
        if i > MAXL:
            raise IndexError("read_complex: nlayers > maxl (fora dos limites no Fortran)")
        soilt[i] = r.inteiro(soilt[i])
        ithick[i] = r.inteiro(ithick[i])
        isname[i] = r.texto4(isname[i])
    r.fim()
    io = r.io

    # ------------------------------------------------------------------
    # terrain orientation (linhas 71-81)
    state.slope_fasst = 0.0
    state.slope_fasst = max(0.0, min(9e1, float(islope)))
    state.sloper = max(0.0, min(1.57, state.slope_fasst * PI / 1.8e2))
    state.sloper = _anint(state.sloper * 1e20) * 1e-20

    state.aspect = 0.0
    state.aspect = max(0.0, min(3.6e2, float(iaspect)))

    state.elev = 0.0
    state.elev = float(ielev)

    # ------------------------------------------------------------------
    # open water (linhas 83-116)
    if wtype < 0 or wtype == 998:
        wtype = 0

    state.water_flag = 0
    if wtype > 0:
        state.water_flag = 1
        if wtype == 2:
            state.water_flag = 2

        wdepth = float(iwdepth)
        if wdepth < 0.0:
            wdepth = SPFLAG

        wvel = float(iwvel)
        if wvel < 0.0:
            wvel = SPFLAG
    else:
        state.water_flag = 0

    if state.water_flag == 1 or state.water_flag == 2:
        nlayers1 = 3
        soiltype1[1] = 26
        lthick1[1] = 2e-2
        soiltype1[2] = 26
        if abs(_aint((wdepth - SPFLAG) * 1e5) * 1e-5) > EPS:
            lthick1[2] = wdepth - state.lthick[1]      # sic: lthick(1), nao lthick1(1)
        else:
            if wtype == 1:
                lthick1[2] = 9.98
            else:
                lthick1[2] = 2.98
        soiltype1[3] = 7
        lthick1[3] = 1.0

    # ------------------------------------------------------------------
    # low vegetation (linhas 118-170)
    if state.vegl_type <= 0 or state.vegl_type == 998:
        state.vegl_type = 0

    state.isigfl = SPFLAG
    state.ihfol = SPFLAG
    if state.vegl_type > 0:
        state.veg_flagl = 1
        state.isigfl = float(ivegldens) * 1e-2
        if state.isigfl > 0.98:
            state.isigfl = 0.98
        if state.isigfl < 0.0:
            state.isigfl = SPFLAG

        state.ihfol = float(iveglht)
        if state.ihfol < 0.0:
            state.ihfol = SPFLAG

    if ((13 <= state.vegl_type <= 15)
            or (wtype > 0 and state.vegl_type != 0)):      # bog/marsh
        if state.vegl_type == 13:
            state.water_flag = 1
            if nlayers1 == 0:
                nlayers1 = 3
                soiltype1[1] = 26
                lthick1[1] = 2e-2
                soiltype1[2] = 26
                if abs(_aint((wdepth - SPFLAG) * 1e5) * 1e-5) > EPS:
                    lthick1[2] = wdepth - lthick1[1]
                else:
                    lthick1[2] = 1.5 - lthick1[1]
                wdepth = lthick1[1] + lthick1[2]
                soiltype1[3] = 7
                lthick1[3] = 1.0
        elif state.vegl_type == 14:                          # inland water
            state.water_flag = 1
            state.vegl_type = 0
            state.nlayers = 3
            soiltype1[1] = 26
            lthick1[1] = 2e-2
            soiltype1[2] = 26
            if abs(_aint((wdepth - SPFLAG) * 1e5) * 1e-5) > EPS:
                lthick1[2] = wdepth - lthick1[1]
            else:
                lthick1[2] = 10.0 - lthick1[1]
            wdepth = state.lthick[1] + lthick1[2]            # sic: lthick(1)
            soiltype1[3] = 7
            lthick1[3] = 1.0
        elif state.vegl_type == 15:
            state.water_flag = 3
            state.vegl_type = 0

    # ------------------------------------------------------------------
    # trees, i.e., high vegetation (linhas 172-203)
    if state.vegh_type <= 0 or state.vegh_type == 998:
        state.vegh_type = 0

    state.isigfh = SPFLAG
    state.izh = SPFLAG
    if state.vegh_type > 0:
        state.veg_flagh = 1
        state.isigfh = float(iveghdens) * 1e-2
        if state.isigfh > 0.98:
            state.isigfh = 0.98
        if state.isigfh < 0.0:
            state.isigfh = SPFLAG

        state.izh = float(iveghht)
        if state.izh < 0.0:
            state.izh = SPFLAG

    if wtype > 0 and state.vegh_type > 0:                    # swamp
        state.water_flag = 1
        if nlayers1 == 0:
            nlayers1 = 3
            soiltype1[1] = 26
            lthick1[1] = 2e-2
            soiltype1[2] = 26
            if abs(_aint((wdepth - SPFLAG) * 1e5) * 1e-5) > EPS:
                lthick1[2] = wdepth - lthick1[1]
            else:
                lthick1[2] = 1.5 - lthick1[1]
            wdepth = lthick1[1] + lthick1[2]
            soiltype1[3] = 7
            lthick1[3] = 1.0

    # ------------------------------------------------------------------
    # initial surface moisture factor (linhas 205-220)
    swater = 0.0
    if iswater == 0 or iswater == -998:
        swater = SPFLAG
    else:
        swater = float(iswater) * 1e-2

    mtest = "n"
    nm0 = 0
    if swater > 0.0 and abs(swater - SPFLAG) > EPS:
        mtest = "y"
        nm0 = 1
        state.sm[1] = swater
        state.zm[1] = 0.0

    # surface roughness (linhas 222-232)
    state.rough = 0.0
    if irough == 1:
        state.rough = 1e-3
    elif irough == 2:
        state.rough = 5e-3
    elif irough == 3:
        state.rough = 1e-2
    elif irough == 4:
        state.rough = 5e-2

    # surface temperature (linhas 234-240)
    state.toptemp = state.mflag
    if -101 < stemp < 101:
        state.toptemp = float(stemp) + TREF
    else:
        state.toptemp = state.mflag

    # groundwater level (linhas 242-244)
    state.gwl = -1.0
    if igwl > 0:
        state.gwl = float(igwl)

    # ------------------------------------------------------------------
    # determine the soil profile (linhas 246-448)
    if r_type == 998:
        r_type = 0

    nlayers = state.nlayers
    if nlayers > 0:
        for i in range(1, nlayers + 1):
            if i > MAXL:
                raise IndexError("read_complex: nlayers > maxl (fora dos limites no Fortran)")
            state.lthick[i] = float(ithick[i]) * 1e-2
            if state.lthick[i] <= 0.0:
                if nlayers == 1:
                    state.lthick[i] = 1.0
                else:
                    state.lthick[i] = 0.25

            if soilt[i] > 0 and soilt[i] != 998:
                state.soiltype[i] = soilt[i]

                if state.soiltype[i] == 26:
                    state.water_flag = 1
                    wflag2 = 1
                    wtype = 1
                    wvel = 0.0

                if r_type > 1:
                    state.rho_fac[i] = 1.2
            elif soilt[i] == -1:
                sname[i] = isname[i]
                if sname[i] in ("    ", "UNK_", "NA__"):
                    state.soiltype[i] = 0
            elif state.soiltype[i] == -2:       # sic: estado, nao soilt(i)
                if r_type > 1 and nlayers > 1:
                    state.rho_fac[i] = 1.2
            elif state.soiltype[i] == -3:
                if r_type > 1 and nlayers > 1:
                    state.rho_fac[i] = 1.2
            elif state.soiltype[i] == -4:
                state.rho_fac[i] = 0.8
            elif soilt[i] == 0 or soilt[i] == 998:
                state.soiltype[i] = 0
                if r_type > 1 and nlayers > 1:
                    state.rho_fac[i] = 1.2

        if nlayers == 1 and r_type >= 1:
            if r_type == 1:                                     # dirt
                nlayers = 2
                state.soiltype[2] = state.soiltype[1]
                state.lthick[2] = 0.5 * state.lthick[1]
                state.rho_fac[2] = 1.0

                state.lthick[1] = state.lthick[2]
            elif r_type == 2:                                   # asphalt
                nlayers = 3
                state.soiltype[3] = state.soiltype[1]
                state.lthick[3] = state.lthick[1]
                state.rho_fac[3] = 1.0

                state.soiltype[1] = 21                          # 'AS'
                state.lthick[1] = 0.08
                state.rho_fac[1] = 1.0

                state.soiltype[2] = 1                           # 'GW'
                state.lthick[2] = 0.6
                state.rho_fac[2] = 1.2
            elif r_type == 3:                                   # concrete
                nlayers = 3
                state.soiltype[3] = state.soiltype[1]
                state.lthick[3] = state.lthick[1]
                state.rho_fac[3] = 1.0

                state.soiltype[1] = 20                          # 'CO'
                state.lthick[1] = 0.08
                state.rho_fac[1] = 1.0

                state.soiltype[2] = 1                           # 'GW'
                state.lthick[2] = 0.6
                state.rho_fac[2] = 1.2
            elif r_type == 4:                                   # grass
                nlayers = 2
                state.soiltype[2] = state.soiltype[1]
                state.lthick[2] = 0.5 * state.lthick[1]
                state.rho_fac[2] = 1.0

                state.lthick[1] = state.lthick[2]

                if state.vegl_type == 0:
                    state.vegl_type = 2                         # short grass
    elif nlayers <= 0 and state.water_flag == 0:
        if r_type == 1:                                         # dirt
            nlayers = 2
            state.soiltype[1] = 7
            state.soiltype[2] = 7

            state.lthick[1] = 5e-1
            state.lthick[2] = 1.0

            state.rho_fac[1] = 1.2
            state.rho_fac[2] = 1.0
        elif r_type == 2:                                       # asphalt
            nlayers = 3
            state.soiltype[1] = 21
            state.lthick[1] = 0.08
            state.rho_fac[1] = 1.0

            state.soiltype[2] = 1
            state.lthick[2] = 0.6
            state.rho_fac[2] = 1.2

            state.soiltype[3] = 7
            state.lthick[3] = 0.4
            state.rho_fac[3] = 1.0
        elif r_type == 3:                                       # concrete
            nlayers = 3
            state.soiltype[1] = 20
            state.lthick[1] = 0.08
            state.rho_fac[1] = 1.0

            state.soiltype[2] = 7                               # sic: comentario diz 'GW'
            state.lthick[2] = 0.6
            state.rho_fac[2] = 1.2

            state.soiltype[3] = 7
            state.lthick[3] = 0.4
            state.rho_fac[3] = 1.0
        elif r_type == 4:                                       # grass
            if state.vegl_type == 0:
                state.vegl_type = 2
            nlayers = 2

            state.soiltype[1] = 7
            state.soiltype[2] = 7

            state.lthick[1] = 5e-1
            state.lthick[2] = 1.0

            state.rho_fac[1] = 1.2
            state.rho_fac[2] = 1.0
        else:                                                   # no road
            nlayers = 1
            state.lthick[1] = 1.0
            state.soiltype[1] = 7
            state.rho_fac[1] = 1.0
    elif wflag2 == 1 and nlayers <= 1:          # inalcancavel (ver docstring)
        nlayers = 3
        state.soiltype[2] = state.soiltype[1]
        state.lthick[2] = state.lthick[1] - 2e-2
        state.soiltype[1] = 26
        state.lthick[1] = 2e-2
        state.soiltype[3] = 7
        state.lthick[3] = 1.0
    elif nlayers1 > 0:
        kk = nlayers1 + nlayers
        j = 1
        d1 = 0.0
        d2 = 0.0
        for i in range(1, kk + 1):
            if i <= nlayers1 and soiltype1[i] == 26:
                lthick2[i] = d1 + lthick1[i]
                soiltype2[i] = soiltype1[i]
            elif i > nlayers1 and state.soiltype[i - nlayers] == 26:
                lthick2[i] = d2 + state.lthick[i - nlayers1]
                soiltype2[i] = state.soiltype[i - nlayers1]
            if lthick2[i] > EPS:
                j = i

        k = sort2(j, lthick2, soiltype2, lthick3, soiltype3)

        for i in range(1, k + 1):
            lthick3[i] = lthick3[i] - d1
            d1 = d1 + lthick3[i]
        d2 = d1

        for i in range(1, kk + 1):
            k = i - j
            if i <= nlayers1 and soiltype1[i] != 26:
                _grava(lthick2, k, d1 + lthick1[i])
                _grava(soiltype2, k, soiltype1[i])
            elif i > nlayers1 and state.soiltype[i - nlayers] != 26:
                _grava(lthick2, k, d2 + state.lthick[i - nlayers1])
                _grava(soiltype2, k, state.soiltype[i - nlayers1])
        jj = sort2(k, lthick2, soiltype2, lthick4, soiltype4)

        for i in range(1, jj + 1):
            lthick4[i] = lthick4[i] - d1
            d1 = d1 + lthick3[i]

        nlayers = k + jj
        for i in range(1, nlayers + 1):
            if i <= k:
                state.lthick[i] = lthick3[i]
                state.soiltype[i] = soiltype3[i]
            else:
                state.lthick[i] = lthick4[i - k]
                state.soiltype[i] = soiltype4[i - k]

    state.nlayers = nlayers
    return wtype, io, lcount, nm0, wvel, wdepth, mtest
