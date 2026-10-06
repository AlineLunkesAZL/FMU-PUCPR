"""
fasst_driver.py -- traducao de fasst_driver.F90

Fonte Fortran: fasst_driver.F90 (1535 linhas). E o `program` principal do
FASST -- nao uma sub-rotina. Integra as versoes de ponto unico e de
multiplos pontos (btra) numa unica casca, organizada pelos proprios
marcadores STEP do original:

  STEP 0  (109-176)    inicializa variaveis nao-alocaveis
  STEP 1  (178-699)    le o arquivo de controle (`fassti.inp`): tipo de
                       solo, numero de camadas, perfil inicial, vegetacao,
                       agua aberta, nomes dos arquivos de saida
  STEP 2  (702-822)    le o arquivo meteorologico (todos os sitios)
  STEP 3  (825-1153)   abre as 15 unidades de arquivo, le dados de execucao
                       anterior ("hot start"), conta lacunas de tempo no
                       met, calcula `maxlines` e ALOCA as ~40 arrays
                       dimensionadas em tempo de execucao
  STEP 4  (1156-1425)  os dois lacos principais de calculo
  STEP 6  (1433-1480)  fecha arquivos, desaloca, imprime o tempo de CPU
                       (nao existe STEP 5 no original)

Mais a sub-rotina `sort2` (1483-1534), apendada ao proprio arquivo.

Assinatura original:
    program fasst_driver        (sem argumentos; le `get_command_argument(1)`)

Arquitetura
-------------
Sendo um `program`, nao recebe `state`: e ele que CRIA o `FasstState` e o
passa para todo o resto. `main(argv)` ocupa esse lugar, e
`get_command_argument(1)` vira `sys.argv[1]`.

`main()` e uma funcao unica e longa, nao um conjunto de helpers por STEP.
Isso e deliberado: no Fortran e um `program`, ou seja UM escopo, e varias
variaveis atravessam fronteiras de STEP de maneira que importa. O caso
extremo e o ponto de atencao 5 (o `i` que sobra do laco do STEP 3 e e
reusado na contagem de lacunas): ele so e reproduzivel com fidelidade se
`i` for uma unica variavel local de `main`. Quebrar `main` em helpers
criaria escopos novos que o original nao tem e consertaria defeitos por
acidente.

O driver e o dono de tres coisas que nenhum outro modulo pode assumir:

1. As 15 unidades de arquivo (`open`/`close` nas linhas 225-1011 e
   1435-1450). Ficam em campos `state.unit*`, que `initsurface.py:56-60` ja
   assume existirem. Esses nomes entram em `_ALLOCATABLE_FIELDS`
   (state.py) -- nao porque sejam alocaveis, mas porque essa tupla e de
   fato "campos que `FasstState.reset()` nao mexe", e um handle de arquivo
   aberto nao pode ser zerado entre pontos do laco multiponto (o Fortran
   abre uma vez no STEP 3 e chama `zero_parameters` por registro VITD, na
   linha 1289).
2. A alocacao das ~40 arrays que `state.py:265-271` deixa explicitamente
   "a cargo de quem transcrever o driver" (linhas 1126-1148). ORDEM
   OBRIGATORIA: alocar primeiro, `zero_mxl_params` depois (Fortran:
   1126-1148, entao 1167) -- `zerovars._MXL_REAL_1D` itera sobre essas
   arrays e precisa que existam.
3. A leitura list-directed de `fassti.inp` (`read(4,*) a,b,c`), que
   consome varios valores por registro e continua no registro seguinte
   quando a lista nao cabe. `us_soil_tools._read_token` e o padrao de
   referencia da casa, mas le UM token por chamada e esta acoplado a
   unidade 90 -- nao serve aqui, e nao sera refatorado. A classe `_Unit`
   deste modulo cobre os dois `read` que o driver usa (`read(u,'(a)')` e
   `read(u,*)`) com a semantica de REGISTRO do Fortran, que e o que faz os
   comentarios `!` de fim de linha funcionarem -- ver `_Unit.read_star`.

Estado desta transcricao
--------------------------
EM ANDAMENTO (cronograma de 10 dias uteis).
  Dia 1  esqueleto, docstring, os 5 stubs externos
  Dia 2  `_Unit`/`_strip_comment`, secao de declaracoes, STEP 0 e o
         cabecalho de `fassti.inp` (linhas 109-266)  <-- aqui
`main()` levanta `NotImplementedError` no ponto exato onde a transcricao
para, de modo que o que ja existe e executavel e verificavel (ver
`demo()`).

Dependencias externas
-----------------------
Reais (reusadas, nao reescritas):
  `functions.met_date`, `functions.map_usda_soiltype_to_uscs`,
  `zerovars.zero_parameters` / `zero_mxl_params`,
  `initsurface.initsurface`, `initwater.initwater`.

Stubs -- as 5 sub-rotinas que o driver chama e que ainda nao foram
transcritas. Assinaturas e `intent` conferidos nos `.F90` reais (nao
inferidos das chamadas), para que a assinatura Python nao mude quando a
transcricao de verdade chegar:

  fasst_main.F90    (586 linhas)  o orquestrador por passo de tempo; e a
                                  que desbloqueia a execucao ponta a ponta
  write_outputs.F90 (375 linhas)
  open_water.F90    (361 linhas)
  read_complex.F90  (450 linhas)  so no modo multiponto
  read_old_data.F90 (286 linhas)  stub ja existente em `initsurface.py` --
                                  importado dali, para nao criar uma
                                  terceira copia (ha uma em `initwater.py`
                                  tambem). O driver passa `dstart` na
                                  posicao que o `.F90` chama `ii` e que o
                                  stub de `initsurface` chama `mpos`: e o
                                  mesmo parametro.

Divergencias deliberadas
--------------------------
* `status='unknown'` nos 15 `open` do original significa "crie o arquivo
  se nao existir". Ou seja: abrir um `fassti.inp` inexistente NAO falha no
  Fortran -- cria um arquivo vazio, e so entao o `read` seguinte morre com
  erro de runtime nao tratado (a linha 233 nao tem `iostat=`). Isso torna
  praticamente mortos os 15 ramos "Error opening ..." do original. Aqui os
  arquivos de ENTRADA sao abertos para leitura e a falha cai justamente no
  ramo de erro do Fortran, com a mesma mensagem: um arquivo de entrada
  que falta e um erro, nao um arquivo vazio recem-criado. Os arquivos de
  SAIDA continuam sendo criados. Precedente de divergencia documentada:
  `functions.map_usda_soiltype_to_uscs` (levanta `ValueError`).
* Varias variaveis locais do original nunca sao inicializadas antes do
  primeiro uso (`keep`, `lcount`, os arrays `soiltype1..4` / `lthick1..4`,
  `nstr_flag`). Em Fortran ficam indefinidas -- o valor real depende do
  compilador, e na pratica muitos zeram. Aqui sao zeros explicitos. Isso
  importa no ponto de atencao 1, onde o bloco `nlayers1` LE posicoes que
  nunca foram atribuidas.
* `modis_type(new_vtl)` / `umd_type(new_vth)` (linhas 404, 411, 483, 485)
  nao tem checagem de limite. Os arrays tem 20 e 14 posicoes, mas
  `new_vtl`/`new_vth` podem chegar a 999 (ex. `vegl_type = 1025` ou
  `2999`). Em Fortran isso le memoria fora do array em silencio; aqui
  levanta `IndexError`. Falhar alto e melhor que propagar lixo, e segue o
  precedente de `map_usda_soiltype_to_uscs`.

Pontos de atencao
-------------------
O repo ja fixou a regra (`indices.py:112`): "Preserve literally when
transcribing; do not 'fix'". Os 13 pontos abaixo sao exatamente onde a
transcricao vai QUERER consertar. Cada um leva um comentario no ponto
correspondente do codigo, conforme for preenchido.

1.  Bloco `nlayers1 > 0` (573-626) -- a regiao mais perigosa do arquivo.
    A linha 581 testa `soiltype(i-nlayers)` mas a 582 atribui de
    `soiltype(i-nlayers1)`: variaveis diferentes. `d1`/`d2` ainda valem 0
    (linha 170) na primeira passagem. A linha 613 avanca `d1` com
    `lthick3(i)` enquanto ajusta `lthick4`. `sort2` declara saidas
    `b(maxn)` com entradas `arr(maxn+3)`. Alcancavel so via
    `vegl_type` 13/14. Transcrever indice por indice, nunca "pelo que faz
    sentido".
2.  Linha 217 -- tautologia. A condicao do `do while` compara `fp` com o
    mesmo literal duas vezes (`.or.` entre operandos identicos), logo e
    sempre verdadeira: o laco sempre roda ate `dirpos == 1`, o
    `if(dirpos == 0)` da linha 221 e inalcancavel, e `dirpos` nunca e lido
    depois. Codigo morto que calcula um valor nao usado.
3.  Linhas 492-512 -- ramo morto. `vegh_type` e zerado na linha 494 e so
    entao testado (`if(vegh_type == 14)`, linha 497): o ramo "inland
    water" nunca executa, sempre cai em "ocean" (`water_flag = 3`), e as
    linhas 510-511 sobrescrevem `wdepth = 3`/`wvel = 0` nos dois casos --
    o que anularia o `wdepth = 10` mesmo se o ramo fosse alcancavel. O
    bloco equivalente da vegetacao baixa (422-452) NAO tem esse defeito:
    la `vegl_type` e testado antes de ser zerado, e os tres ramos
    (swamp/bog, inland water, ocean) funcionam. Entradas que alcancam
    cada um: MODIS 11 -> 13, MODIS 17 -> 14, UMD 14 -> 15.
4.  Assimetrias entre o bloco de vegetacao baixa (391-467) e o de alta
    (471-534), que no mais sao quase identicos:
    a) linha 399 `new_vtl = vegl_type - biome_source` NAO tem o `max(0,.)`
       que a linha 479 tem para `new_vth`;
    b) a vegetacao baixa exclui explicitamente os tipos MODIS 13 e 16
       (linha 403) e UMD 13 (linha 410), zerando `veg_flagl` junto; a
       vegetacao alta nao exclui nada (linha 483), entao MODIS 13 ->
       `vegh_type = 0` COM `veg_flagh` continuando 1 -- um estado que o
       lado da vegetacao baixa nunca produz.
14. `biome_source`, `new_vtl` e `new_vth` sao campos de estado que NINGUEM
    rezera: nao estao em `zero_parameters` (conferido em
    `fasst/zerovars.py`), e as linhas 394-398 / 474-478 so os ESCREVEM
    quando o tipo lido e >= 1000. O valor vaza em dois niveis:
    a) da vegetacao baixa para a alta, dentro do mesmo ponto. Um arquivo
       com baixa MODIS/UMD (ex. 1003) e alta no padrao FASST (ex. 5)
       deixa `biome_source = 1000` valendo para a alta, fazendo
       `new_vth = max(0, 5 - 1000) = 0`, e a vegetacao alta e descartada
       em silencio (`vegh_type = 0`, `veg_flagh = 0`) -- embora os 4
       registros do dossel continuem sendo consumidos do arquivo. Na
       ordem inversa nao acontece, porque a baixa e lida primeiro.
    b) de um ponto de grade para o SEGUINTE, no modo multiponto. O STEP 0
       (linhas 147-149) zera os tres UMA vez, no inicio do programa; o
       laco da linha 1289 chama `zero_parameters` por registro VITD, e ela
       nao toca neles. Um ponto com vegetacao MODIS contamina todos os
       pontos seguintes que usem o padrao FASST, perdendo tanto a
       vegetacao baixa quanto a alta.
    E so por (b) que o ponto 4a tem efeito observavel: como desvio,
    `max(0,x) > 0` e `x > 0` sao equivalentes. O que o `max` muda e o
    VALOR guardado em `new_vtl`/`new_vth`, que `initsurface.py:269,272`
    repassa a `veg_propl`/`veg_proph`: sem o `max`, `new_vtl` chega la
    negativo. Como o STEP 0 zera `biome_source` no inicio, (b) e o 4a so
    sao demonstraveis a partir da SEGUNDA iteracao do laco multiponto --
    a checagem dos dois fica para o Dia 9, quando esse laco existir.
5.  Vazamento do `i` (1063 + 1095-1099). Depois de `i = i - 1`, o laco de
    contagem de lacunas indexa `metm(i,k,...)`. No ramo ponto-unico `i` e o
    sitio casado; no multiponto e um contador de linhas da unidade 5, sem
    relacao nenhuma. Python vai acidentalmente consertar OU acidentalmente
    quebrar isso (ver "Arquitetura": e por isso que `main` e um escopo so).
6.  STEP 4: os dois ramos sao ~90% iguais e NAO sao identicos. Linha 1185:
    `met(k,ip_upsol)-met(k,ip_tsol) > eps`; linha 1324:
    `met(k,ip_upsol) > met(k,ip_tsol)`. Linha 1275 usa `i <= met_count`;
    linha 1418 usa `i >= met_count`. NAO fatorar os dois ramos num helper
    comum: e exatamente esse diff que uma fatoracao destroi.
7.  Reporte de erro morto no multiponto. As linhas 1220 e 1358 tem a mesma
    guarda `if(error_code == 1.and.single_multi_flag == 0)`, mas a 1358
    esta DENTRO do ramo `single_multi_flag == 1` -- nunca dispara.
8.  `met_col = 32` vs `MAXCOL = 35` (linha 97, lido na 796): a coluna 33
    (`MetCol.HI`) nunca e lida do arquivo e fica com `mflag`.
9.  `tstps == 0` quando `timstep <= 0` (linha 1092): `rtstps` fica 0, logo
    `msteps1 = 0` sempre e nenhuma lacuna e detectada.
10. Linhas 565-572: `soiltype(2) = soiltype(1)` acontece ANTES de
    `soiltype(1) = 26`. A ordem e o comportamento.
11. Linha 1006: `ro2o = ro1o(1:len_trim(ro1o)-9)//'_time.dat'` assume que o
    nome termina em `_soil.dat` (9 caracteres). Sem validacao no original.
12. Linhas 1439-1440: `close(11/12, status='delete')` APAGA os arquivos de
    hot-start antigos ao final de uma execucao `infer_test == 1`.
13. Strip de comentario: o idioma `s(1:len_trim(s(1:index(s,'!')-1)))`
    (linhas 234, 257, 275, ...) da string VAZIA quando nao ha `!` no
    registro, porque `index` devolve 0 e `0-1 = -1`. A linha 878
    (`fasstusersoil = 'fasstusersoil.inp'` quando `lcount == 0`) existe
    justamente para consertar esse efeito. Ver `_strip_comment`.

Pendencias
------------
* As 5 transcricoes acima. `fasst_main` e a maior e a unica que desbloqueia
  execucao ponta a ponta.
* Convencao de unidade de arquivo dividida: `initsurface` assume
  `state.unit*`, `us_soil_tools.get_user_soil_params` recebe `unit90` por
  parametro. O driver passa a preencher `state.unit*`; unificar as duas
  convencoes fica para depois.
* `_read_table_lines` (`initsurface.py:181`) e `_read_veg_table_row`
  (`module_canopy.py:2039`) leem os mesmos arquivos (units 31/32/33) com
  contratos incompativeis. O driver abre as unidades; reconciliar os dois
  leitores e trabalho de `fasst_main`.
* `makemtagroundout` / `makemtafasstout` (1428-1430) estao comentadas no
  proprio Fortran. Nao transcrever.
"""

import math
import sys
import time

from .constants import EPS, MAXCOL, MAXL, MAXN, MOVERLAP, NCLAYERS, PI, SPFLAG, TREF
from .indices import MetCol
from .state import FasstState
from .functions import map_usda_soiltype_to_uscs, met_date
from .zerovars import zero_mxl_params, zero_parameters
from .initsurface import initsurface
from .initwater import initwater

# read_old_data.F90 ainda nao foi transcrita; o stub ja existe em
# initsurface.py (e tambem em initwater.py). Importado dali para nao criar
# uma terceira copia -- ver "Dependencias externas" na docstring.
from .initsurface import read_old_data


# ---------------------------------------------------------------------------
# `data` do original (linhas 97-106)
# ---------------------------------------------------------------------------
# data met_col/32/ -- numero de colunas lidas do arquivo met. NAO e o mesmo
# que MAXCOL (35): ver ponto de atencao 8.
MET_COL = 32

# Conversao dos tipos de vegetacao MODIS e UMD para os defaults do FASST.
# biome_source: FASST default = 0, MODIS/IGBP = 1000, UMD/LDAS = 2000.
# Tabelas 1-based: slot fantasma no indice 0, como no resto do pacote.
#                1  2  3  4   5   6   7   8  9 10  11  12 13 14  15 16
_MODIS_TYPE = [0, 3, 6, 4, 5, 18, 16, 16, 18, 7, 2, 13, 10, 0, 1, 12, 8,
#               17 18 19 20
                14, 9, 9, 9]
#              1  2  3  4   5  6  7   8   9 10 11 12 13  14
_UMD_TYPE = [0, 3, 6, 5, 4, 18, 3, 2, 11, 11, 2, 1, 8, 0, 15]

assert len(_MODIS_TYPE) == 21, "modis_type(20) + slot fantasma"
assert len(_UMD_TYPE) == 15, "umd_type(14) + slot fantasma"


# ---------------------------------------------------------------------------
# helpers numericos (identicos aos de initsurface.py:157-170)
# ---------------------------------------------------------------------------
def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _aint(x, p=0):
    if not math.isfinite(x):
        return x
    return math.trunc(x * (10.0**p)) * (10.0**-p)


def _anint(x, p=0):
    if not math.isfinite(x):
        return x
    xs = x * (10.0**p)
    ai = math.copysign(math.floor(abs(xs) + 0.5), xs)
    return ai * (10.0**-p)


# ---------------------------------------------------------------------------
# helpers de string/arquivo no estilo Fortran
# ---------------------------------------------------------------------------
def _strip_comment(s):
    """O idioma `s(1:len_trim(s(1:index(s,'!')-1)))` do original.

    Ponto de atencao 13: quando NAO ha `!` no registro, `index` devolve 0,
    a substring fica `s(1:-1)` (comprimento zero) e o resultado e uma
    string VAZIA -- nao a linha inteira. Aparece nas linhas 234, 257, 275,
    284, 293, 302, 310, 629, 674 e 684 do original, e e a razao de existir
    a linha 878.
    """
    p = s.find("!")
    if p < 0:
        return ""
    return s[:p].rstrip()


def _ch(s, i):
    """`s(i:i)` num `character(len=N)`.

    A variavel Fortran e preenchida com brancos ate o comprimento
    declarado, entao um indice alem do texto devolve ' ' em vez de estourar.
    """
    return s[i - 1] if 1 <= i <= len(s) else " "


def _stop(msg):
    """`write(*,'(/,''...'',/)')` seguido de `stop`.

    O FORMAT `(/,...,/)` do original imprime uma linha em branco antes e
    depois da mensagem. Usado nos 8 pontos de validacao do STEP 1.
    """
    print()
    print(msg)
    print()
    sys.exit(1)


def _to_int(tok, current=0):
    """Item list-directed inteiro. `None` (valor nulo) deixa a variavel
    INALTERADA, como no Fortran."""
    if tok is None:
        return current
    return int(tok)


def _to_real(tok, current=0.0):
    """Item list-directed real. `None` (valor nulo) deixa a variavel
    INALTERADA. Aceita o expoente `d`/`D` do Fortran."""
    if tok is None:
        return current
    return float(tok.replace("d", "e").replace("D", "E"))


def _list_items(rec):
    """Itens list-directed de UM registro; `None` = valor nulo explicito.

    Separador Fortran: uma virgula (com brancos opcionais em volta) ou um
    ou mais brancos. Duas virgulas seguidas, ou uma virgula inicial, sao um
    valor nulo -- a variavel correspondente fica inalterada. Uma virgula no
    FIM do registro nao e valor nulo: a leitura simplesmente continua no
    registro seguinte. Um registro todo em branco nao tem item nenhum (nao
    e um nulo), e a leitura segue adiante -- igual ao que
    `us_soil_tools._read_token` ja documenta.
    """
    fields = rec.split(",")
    items = []
    for k, field in enumerate(fields):
        toks = field.split()
        if toks:
            items.extend(toks)
        elif k < len(fields) - 1:
            items.append(None)
    return items


class _Unit:
    """Um arquivo texto aberto, com a semantica de REGISTRO do Fortran.

    Cobre os dois `read` que o driver usa. O ponto importante e que cada
    `read` consome registros INTEIROS: o que sobra no ultimo registro usado
    e descartado, e nunca aproveitado pelo `read` seguinte.

    E exatamente isso que faz os comentarios `!` de fim de linha de
    `fassti.inp` funcionarem com `read(4,*)`: a lista de variaveis e
    satisfeita antes de a leitura chegar ao `!`, e o resto do registro vai
    embora. Note que isso e um acidente feliz, nao tratamento de
    comentario -- um registro com MENOS valores do que a lista pede faz o
    Fortran continuar no registro seguinte, e ai um `!` seria lido como
    valor e quebraria a leitura.
    """

    def __init__(self, path, mode="r"):
        self.path = path
        self.f = open(path, mode)

    def read_a(self):
        """`read(u,'(a)') s` -- um registro inteiro, sem o fim de linha."""
        line = self.f.readline()
        if line == "":
            raise EOFError("fim do arquivo em read('(a)'): " + str(self.path))
        return line.rstrip("\n").rstrip("\r")

    def read_star(self, n):
        """`read(u,*) v1,...,vn` -- leitura list-directed de `n` itens.

        Devolve `n` tokens. `None` num slot significa valor nulo: no
        Fortran a variavel correspondente fica inalterada (e por isso que
        `_to_int`/`_to_real` recebem o valor atual). Um `/` encerra a
        leitura e os itens restantes voltam como `None`.
        """
        out = []
        terminated = False
        while len(out) < n and not terminated:
            line = self.f.readline()
            if line == "":
                raise EOFError("fim do arquivo em read(*): " + str(self.path))
            for t in _list_items(line.rstrip("\n").rstrip("\r")):
                if t is not None and t.startswith("/"):
                    terminated = True
                    break
                out.append(t)
                if len(out) == n:
                    break
        return (out + [None] * n)[:n]

    def close(self):
        self.f.close()


# ---------------------------------------------------------------------------
# Dependencias externas ainda sem implementacao (ver docstring).
# As assinaturas vem dos `.F90` reais, com os `intent` conferidos: os
# parametros `intent(out)`/`intent(inout)` entram na tupla de retorno, no
# padrao do resto do repo (cf. `initsurface` -> 10-tupla).
# ---------------------------------------------------------------------------
def fasst_main(state: FasstState, nprint, sprint, frozen, nstr_flag, sdensi,
               surfdepth, code2, code3, code4, code5, first_loop, phie,
               pdens, oldsd, oldhi):
    """TODO: fasst_main.F90 (586 linhas) -- orquestrador por passo de tempo.

    subroutine fasst_main(nprint,sprint,frozen,nstr_flag,sdensi,surfdepth,
                          code2,code3,code4,code5,first_loop,phie,pdens,
                          oldsd,oldhi)
      intent(in)    :: nprint,sprint,frozen,nstr_flag(maxn),sdensi,surfdepth
      intent(inout) :: code2,code3,code4,code5,first_loop,phie,pdens,oldsd,oldhi

    Retorno esperado quando transcrita:
        (code2, code3, code4, code5, first_loop, phie, pdens, oldsd, oldhi)
    """
    raise NotImplementedError(
        "fasst_main ainda nao foi traduzida (fasst_main.F90)"
    )


def write_outputs(state: FasstState, flprint, keep, vegint):
    """TODO: write_outputs.F90 (375 linhas) -- escreve as saidas da execucao.

    subroutine write_outputs(flprint,keep,vegint)
      intent(in) :: flprint,keep,vegint

    Sem `intent(out)`: nao retorna nada.
    """
    raise NotImplementedError(
        "write_outputs ainda nao foi traduzida (write_outputs.F90)"
    )


def open_water(state: FasstState, wtype, sprint, wvel, wdepth, sdensi,
               first_loop, code2, code4, oldsd):
    """TODO: open_water.F90 (361 linhas) -- caminho de agua aberta.

    subroutine open_water(wtype,sprint,wvel,wdepth,sdensi,first_loop,
                          code2,code4,oldsd)
      intent(in)    :: wtype,sprint,wvel,wdepth,sdensi
      intent(inout) :: first_loop,code2,code4,oldsd

    Retorno esperado quando transcrita:
        (first_loop, code2, code4, oldsd)
    """
    raise NotImplementedError(
        "open_water ainda nao foi traduzida (open_water.F90)"
    )


def read_complex(state: FasstState, wtype, io, lcount, nm0, wvel, wdepth,
                 mtest, sname):
    """TODO: read_complex.F90 (450 linhas) -- le um registro VITD (multiponto).

    subroutine read_complex(wtype,io,lcount,nm0,wvel,wdepth,mtest,sname)
      intent(inout) :: wtype,io,wvel,wdepth
      intent(out)   :: lcount,nm0,mtest,sname(maxl)

    Retorno esperado quando transcrita:
        (wtype, io, lcount, nm0, wvel, wdepth, mtest, sname)
    """
    raise NotImplementedError(
        "read_complex ainda nao foi traduzida (read_complex.F90)"
    )


# ---------------------------------------------------------------------------
# program fasst_driver
# ---------------------------------------------------------------------------
def main(argv=None, _state=None):
    """`program fasst_driver` (fasst_driver.F90, linhas 1-1480).

    Um escopo unico e longo, de proposito -- ver "Arquitetura" na docstring
    do modulo.
    """
    if argv is None:
        argv = sys.argv

    # O `program` nao recebe estado: e ele que o cria.
    # ponytail: `_state` e gancho de teste -- deixa o `demo()` inspecionar o
    # FasstState pelo caminho de codigo REAL enquanto `main` ainda para no
    # meio da transcricao. Sem ele, os ~590 linhas de parsing dos dias 3-5
    # (a classe de bug dominante) so teriam verificacao no Dia 9. Remover
    # quando `main` rodar ate o fim e o estado puder ser conferido pelas
    # saidas.
    state = _state if _state is not None else FasstState()

    # ---------------- declaracoes (linhas 61-97) ----------------
    # Apenas os arrays: escalares nao precisam ser declarados em Python.
    # O Fortran NAO inicializa nenhum destes (nem no STEP 0); aqui sao
    # zeros explicitos -- ver "Divergencias deliberadas" e o ponto de
    # atencao 1, onde o bloco `nlayers1` le posicoes nunca atribuidas.
    nstr_flag = _new_array(MAXN, 0)              # integer nstr_flag(maxn)
    soiltype1 = _new_array(3, 0)                 # integer soiltype1(3)
    soiltype2 = _new_array(MAXN + 3, 0)          # integer soiltype2(maxn+3)
    soiltype3 = _new_array(MAXN + 3, 0)
    soiltype4 = _new_array(MAXN + 3, 0)
    lthick1 = _new_array(3)                      # real lthick1(3)
    lthick2 = _new_array(MAXN + 3)               # real lthick2(maxn+3)
    lthick3 = _new_array(MAXN + 3)
    lthick4 = _new_array(MAXN + 3)
    sname = _new_array(MAXL, "    ")             # character(len=4) sname(maxl)
    t1o = _new_array(MOVERLAP)                   # real t1o(moverlap)
    timeo = [[0.0] * 5 for _ in range(MOVERLAP + 1)]   # real timeo(moverlap,4)

    # Escalares que o original declara mas nao inicializa antes do primeiro
    # uso (o Fortran os deixa indefinidos) -- ver "Divergencias". Os sete
    # ultimos sao `intent(out)` de `zero_parameters`, que so e chamada
    # dentro dos ramos de ponto unico/multiponto: no ramo de
    # `single_multi_flag` invalido eles ficam sem valor nenhum no original.
    keep = 0
    tflag = 0
    nt0 = 0
    nm0 = 0
    stflag = 0
    lcount = 0
    wtype = 0
    wvel = 0.0
    wdepth = 0.0
    ttest = " "
    mtest = " "

    it = 0.0
    et = 0.0
    it = time.process_time()                                   # call cpu_time(it)

    # ---------------- STEP 0 (linhas 114-174) ----------------
    # inicializa algumas variaveis
    state.single_multi_flag = 0
    state.infer_test = 0
    vegint = 0.0
    nprint = 0
    flprint = 0
    sprint = 0
    frozen = 100
    state.istart = 0
    state.iend = 0
    state.oldpos = 0
    state.mpos = 0
    iendi = 0
    wstart = 0
    wflag2 = 0
    nlayers1 = 0
    fasst_input = " "
    met_file_name = " "
    fasstusersoil = "fasstusersoil.inp"
    fasstoutmta = "fasstout.mta"                 # so para makemtafasstout, que esta
    groundoutmta = "groundout.mta"               # comentada no original (1428-1430)
    vitd_input = " "
    output1 = " "
    output2 = " "
    output3 = " "
    output4 = " "
    output5 = " "
    ro1 = " "
    ro2 = " "
    ro1o = " "
    ro2o = " "

    state.biome_source = 0
    state.new_vtl = 0
    state.new_vth = 0
    icnti = 0
    state.ncols = 0
    pindex = 0
    state.met_count = 0
    state.mstflag = 0
    npos = 0
    ijunk = 0
    lis = 0
    fnlen = 0
    mlen = 0
    state.elev = 0.0
    state.lat = 0.0
    state.mlong = 0.0
    state.mflag = 0.0
    state.timeoffset = 0.0
    state.timstep = 0.0
    state.iheight = 0.0
    junk = 0.0
    rtstps = 0.0
    surfdepth = 0.0
    d1 = 0.0
    d2 = 0.0
    mody = 0.0
    modyo = 0.0
    daylim = 0.0

    # ---------------- STEP 1 (linhas 178-699) ----------------
    # Este arquivo contem a informacao necessaria para rodar o modelo: por
    # exemplo o tipo de solo, o numero de camadas, dados de perfil inicial
    # e informacao de profundidade de neve.

    # Se existir argumento de linha de comando, usa os nomes de arquivo
    # dele (linha 184).
    if len(argv) > 1:
        fasst_input = argv[1]
    fasst_input = fasst_input.rstrip()           # fasst_input(1:len_trim(fasst_input))

    # Linhas 214-221 -- ponto de atencao 2: a condicao do `do while` compara
    # `fp` com o MESMO literal nos dois lados do `.or.`, logo e sempre
    # verdadeira e o laco so para por `dirpos > 1`. Consequencias: o
    # `if(dirpos == 0)` da linha 221 e inalcancavel, e `dirpos`/`fp` nunca
    # sao lidos depois daqui. Codigo morto preservado como no original.
    dirpos = len(fasst_input.rstrip()) + 1
    fp = _ch(fasst_input, dirpos)
    while dirpos > 1:                            # a tautologia do `.or.` sempre vale
        dirpos = dirpos - 1
        fp = _ch(fasst_input, dirpos)
    if dirpos == 0:                              # inalcancavel (ver acima)
        dirpos = 1

    # abre o arquivo de entrada (linhas 224-230)
    io = 0
    try:
        unit4 = _Unit(fasst_input)
    except OSError as err4:
        # Ver "Divergencias deliberadas": no original `status='unknown'`
        # criaria o arquivo e a falha so apareceria no `read` seguinte.
        io = err4.errno or 1
        print(" Error opening FASST input file = %s error=%3d" % (fasst_input, io))
        sys.exit(1)

    # le o arquivo de entrada (linhas 233-266)
    met_file_name = unit4.read_a()
    met_file_name = _strip_comment(met_file_name)        # ponto de atencao 13
    mlen = len(met_file_name.rstrip())

    tok = unit4.read_star(4)
    state.infer_test = _to_int(tok[0], state.infer_test)
    state.single_multi_flag = _to_int(tok[1], state.single_multi_flag)
    state.mgap = _to_real(tok[2], state.mgap)
    keep = _to_int(tok[3], keep)
    if state.infer_test != 1 and state.infer_test != 0:
        _stop("!!! The infer test flag must be 0 or 1, stopping. !!!")

    state.oldpos = 1
    state.istart = 1
    if state.infer_test == 1:
        state.istart = 2
        ro1o = unit4.read_a()                    # dados antigos que inicializam a nova execucao
        ro1o = _strip_comment(ro1o)

    # mstflag: 0 = produto AFWA, nao usar temperatura de solo medida
    #          1 = outro
    tok = unit4.read_star(7)
    state.gwl = _to_real(tok[0], state.gwl)
    vegint = _to_real(tok[1], vegint)
    nprint = _to_int(tok[2], nprint)
    flprint = _to_int(tok[3], flprint)
    sprint = _to_int(tok[4], sprint)
    frozen = _to_int(tok[5], frozen)
    state.mstflag = _to_int(tok[6], state.mstflag)

    # ---- ramo de ponto unico (linhas 268-670) ----
    if state.single_multi_flag == 0:
        zp = zero_parameters(state)
        # `zero_parameters` preserva tudo que foi lido acima (infer_test,
        # single_multi_flag, mgap, oldpos, istart, gwl, mstflag) e zera
        # exatamente os campos que este bloco vai ler em seguida.
        tflag = zp.tflag
        nt0 = zp.nt0
        nm0 = zp.nm0
        stflag = zp.stflag
        lcount = zp.lcount
        wtype = zp.wtype                             # = int(SPFLAG) = 999
        wvel = zp.wvel                               # = SPFLAG
        wdepth = zp.wdepth                           # = SPFLAG (ver sentinela abaixo)
        ttest = zp.ttest
        mtest = zp.mtest
        sname = zp.sname

        tok = unit4.read_star(1)
        state.vitd_index = _to_int(tok[0], state.vitd_index)

        # Os cinco nomes de arquivo de saida passam pelo mesmo strip do
        # ponto de atencao 13: uma linha SEM `!` deixa o nome VAZIO. Para
        # `fasstusersoil` a linha 878 conserta isso depois; para estes cinco
        # nao ha conserto nenhum -- o `open` do STEP 3 recebe nome vazio.
        output1 = _strip_comment(unit4.read_a())     # met and surface data
        output2 = _strip_comment(unit4.read_a())     # vertical profile data
        output3 = _strip_comment(unit4.read_a())     # surface flux data
        output4 = _strip_comment(unit4.read_a())     # vegetation temperature output data
        output5 = _strip_comment(unit4.read_a())     # snow processes output data

        # graus a partir da horizontal, graus a partir do N
        tok = unit4.read_star(2)
        state.slope_fasst = _to_real(tok[0], state.slope_fasst)
        state.aspect = _to_real(tok[1], state.aspect)
        # PI aqui e a constante do fasst_global (3.141592654, truncada), nao
        # `math.pi` -- ver fasst/constants.py.
        state.sloper = max(0.0, min(1.57, state.slope_fasst * PI / 1.8e2))
        state.sloper = _anint(state.sloper, 20)
        # `sloper` e calculado ANTES da validacao, como no original: um
        # slope invalido ainda e convertido e so entao o programa para.
        if state.slope_fasst > 90.0:
            _stop("!!! Slope can not be greater than 90, stopping. !!!")

        if abs(state.aspect) > 360.0:
            _stop("!!! Aspect can not be greater than +/-360, stopping. !!!")

        tok = unit4.read_star(1)
        state.rough = _to_real(tok[0], state.rough)  # rugosidade da superficie (m)

        # determina se o poligono e agua aberta (1) ou terra (0)
        tok = unit4.read_star(1)
        state.water_flag = _to_int(tok[0], state.water_flag)
        if state.water_flag != 1 and state.water_flag != 0:
            _stop("!!! Open water flag must be 0 or 1, stopping. !!!")

        if state.water_flag == 1:
            tok = unit4.read_star(1)
            wtype = _to_int(tok[0], wtype)           # 1 = lagos, lagoas; 2 = rios
            if wtype == 2:
                state.water_flag = 2
            tok = unit4.read_star(1)
            wdepth = _to_real(tok[0], wdepth)
            tok = unit4.read_star(1)
            wvel = _to_real(tok[0], wvel)

            state.soiltype[1] = 26                   # WAter
            state.lthick[1] = 2e-2
            state.soiltype[2] = 26
            state.soiltype[3] = 7                    # SM
            # `aint((wdepth-spflag)*1d5)*1d-5`: testa se `wdepth` difere do
            # sentinela. `zero_parameters` o deixa em SPFLAG, entao so um
            # valor realmente lido do arquivo passa aqui -- um item nulo
            # (virgula dupla) mantem o sentinela e cai no `else`.
            if abs(_aint(wdepth - SPFLAG, 5)) > EPS:
                state.lthick[2] = wdepth - state.lthick[1]
            else:
                if wtype == 1:
                    state.lthick[2] = 1e1 - state.lthick[1]
                else:
                    state.lthick[2] = 3e0 - state.lthick[1]
            state.lthick[3] = 1e0
            state.nlayers = 3

        if state.infer_test == 0:
            # profundidade inicial de neve (m), equivalente em agua (m)
            tok = unit4.read_star(2)
            state.hsaccum = _to_real(tok[0], state.hsaccum)
            state.iswe = _to_real(tok[1], state.iswe)
            if state.hsaccum < 0.0:
                _stop("!!! Initial snow depth can not be less than 0, stopping. !!!")

            if state.iswe < 0.0:
                _stop("!!! Initial snow water equivalent can not be less than 0, "
                      "stopping. !!!")

            tok = unit4.read_star(1)
            state.hi = _to_real(tok[0], state.hi)    # espessura inicial de gelo (m)
            if state.hi < 0.0:
                _stop("!!! Initial ice thickness can not be less than 0, stopping. !!!")

            # ---- vegetacao e camadas de solo (linhas 388-669) ----
            # Note a aninhagem: todo este bloco esta DENTRO de
            # `infer_test == 0`. Numa execucao de hot start (infer_test==1)
            # nada disso e lido -- o perfil vem do arquivo da execucao
            # anterior.
            if state.water_flag == 0:
                tok = unit4.read_star(1)
                # tipo de vegetacao baixa na superficie
                state.vegl_type = _to_int(tok[0], state.vegl_type)

                if state.vegl_type != 0:
                    state.veg_flagl = 1

                    if 1000 <= state.vegl_type < 2000:
                        state.biome_source = 1000                # MODIS/IGBP
                    elif 2000 <= state.vegl_type < 3000:
                        state.biome_source = 2000                # UMD/LDAS
                    # Ponto de atencao 4a: aqui NAO ha `max(0,...)` -- a
                    # linha 479, do lado da vegetacao alta, tem.
                    state.new_vtl = state.vegl_type - state.biome_source

                    if state.new_vtl > 0:
                        if state.biome_source == 1000:           # MODIS
                            # Ponto de atencao 4b: os tipos 13 e 16 sao
                            # excluidos aqui; do lado da vegetacao alta nao
                            # ha exclusao nenhuma.
                            if state.new_vtl != 13 and state.new_vtl != 16:
                                state.vegl_type = _MODIS_TYPE[state.new_vtl]
                            else:
                                state.vegl_type = 0
                                state.veg_flagl = 0
                        elif state.biome_source == 2000:         # UMD
                            if state.new_vtl != 13:
                                state.vegl_type = _UMD_TYPE[state.new_vtl]
                            else:
                                state.vegl_type = 0
                                state.veg_flagl = 0
                        # Com `biome_source == 0` (padrao FASST) nenhum dos
                        # dois ramos roda e `vegl_type` segue com o valor
                        # lido do arquivo.
                    else:
                        state.vegl_type = 0
                        state.veg_flagl = 0

                    # Ao contrario do bloco equivalente da vegetacao alta
                    # (ponto de atencao 3), aqui `vegl_type` e testado ANTES
                    # de ser zerado, e os tres ramos funcionam.
                    if 13 <= state.vegl_type <= 15:
                        wtype = 1
                        wvel = 0.0
                        if state.vegl_type == 13:                # swamp/bog
                            state.water_flag = 1
                            wdepth = 1.5
                            nlayers1 = 3
                            soiltype1[1] = 26
                            lthick1[1] = 2e-2
                            soiltype1[2] = 26
                            lthick1[2] = wdepth - lthick1[1]
                            soiltype1[3] = 7
                            lthick1[3] = 1e0
                        elif state.vegl_type == 14:              # inland water
                            state.water_flag = 1
                            state.veg_flagl = 0
                            state.vegl_type = 0
                            wdepth = 10e0
                            nlayers1 = 3
                            soiltype1[1] = 26
                            lthick1[1] = 2e-2
                            soiltype1[2] = 26
                            lthick1[2] = wdepth - lthick1[1]
                            soiltype1[3] = 7
                            lthick1[3] = 1e0
                        else:                                    # ocean
                            state.water_flag = 3
                            state.vegl_type = 0
                            state.veg_flagl = 0
                            # `wdepth` NAO e tocado neste ramo: segue com o
                            # sentinela SPFLAG, porque o bloco de agua da
                            # linha 344 nao rodou (water_flag era 0). Isso
                            # chega ao bloco de fusao de camadas (ponto 1).

                    tok = unit4.read_star(1)
                    # densidade inicial de folhagem; limitada a 0.98, exceto
                    # se for o sentinela SPFLAG (que e > 0.98) -- mesmo
                    # idioma do `wdepth` na linha 354
                    state.isigfl = _to_real(tok[0], state.isigfl)
                    if state.isigfl > 0.98 and abs(_aint(state.isigfl - SPFLAG, 5)) > EPS:
                        state.isigfl = 0.98
                    tok = unit4.read_star(1)
                    state.iepf = _to_real(tok[0], state.iepf)    # emissividade inicial
                    tok = unit4.read_star(1)
                    state.ifola = _to_real(tok[0], state.ifola)  # absortividade (1 - albedo)
                    tok = unit4.read_star(1)
                    state.ihfol = _to_real(tok[0], state.ihfol)  # altura inicial (cm)

                    if abs(state.isigfl) <= EPS:
                        state.vegl_type = 0
                        state.veg_flagl = 0

                    if state.vegl_type == 12:                    # ice cap / neve permanente
                        state.veg_flagl = 0

                tok = unit4.read_star(1)
                # tipo de vegetacao alta (arvores)
                state.vegh_type = _to_int(tok[0], state.vegh_type)

                if state.vegh_type != 0:
                    state.veg_flagh = 1

                    if 1000 <= state.vegh_type < 2000:
                        state.biome_source = 1000
                    elif 2000 <= state.vegh_type < 3000:
                        state.biome_source = 2000
                    # Ponto de atencao 4a: COM `max(0,...)`, ao contrario da
                    # linha 399. Ponto de atencao 14: `biome_source` nao foi
                    # rezerado depois da vegetacao baixa, entao um arquivo
                    # com baixa MODIS/UMD e alta no padrao FASST cai aqui em
                    # `max(0, negativo) = 0` e perde a vegetacao alta.
                    state.new_vth = max(0, state.vegh_type - state.biome_source)

                    if state.new_vth > 0:
                        # Ponto de atencao 4b: sem a exclusao de 13/16 que a
                        # vegetacao baixa tem. MODIS 13 da `vegh_type = 0`
                        # com `veg_flagh` continuando 1.
                        if state.biome_source == 1000:
                            state.vegh_type = _MODIS_TYPE[state.new_vth]
                        elif state.biome_source == 2000:
                            state.vegh_type = _UMD_TYPE[state.new_vth]
                    else:
                        state.vegh_type = 0
                        state.veg_flagh = 0

                    # Ponto de atencao 3 -- ramo morto. `vegh_type` e zerado
                    # duas linhas abaixo e SO ENTAO testado contra 14, entao
                    # o ramo "inland water" nunca executa: cai sempre no
                    # `else` (ocean). E as duas ultimas atribuicoes
                    # sobrescrevem `wdepth`/`wvel` nos dois casos, anulando
                    # o `wdepth = 10` mesmo se o ramo fosse alcancavel.
                    if state.vegh_type == 14 or state.vegh_type == 15:  # inland water, ocean
                        wtype = 1
                        state.vegh_type = 0
                        state.veg_flagh = 0

                        if state.vegh_type == 14:                # inalcancavel
                            state.water_flag = 1
                            wdepth = 10e0
                            nlayers1 = 3
                            soiltype1[1] = 26
                            lthick1[1] = 2e-2
                            soiltype1[2] = 26
                            lthick1[2] = wdepth - lthick1[1]
                            soiltype1[3] = 7
                            lthick1[3] = 1e0
                        else:                                    # ocean (sempre)
                            state.water_flag = 3
                        wdepth = 3e0
                        wvel = 0e0

                    tok = unit4.read_star(2)
                    state.izh = _to_real(tok[0], state.izh)      # altura do dossel (m)
                    state.isigfh = _to_real(tok[1], state.isigfh)  # densidade do dossel (0 - 1)

                    if state.isigfh > 0.98 and abs(_aint(state.isigfh - SPFLAG, 5)) > EPS:
                        state.isigfh = 0.98

                    # para cada camada: indice de area foliar, fator de
                    # agrupamento de Markov ("Markox" no original,
                    # tipicamente 1.0), espessura, reflectancia,
                    # transmissao, absorcao de onda curta (1 - albedo),
                    # emissividade de onda longa
                    for j in range(1, NCLAYERS + 1):
                        tok = unit4.read_star(8)
                        state.ifoliage_type[j] = _to_real(tok[0], state.ifoliage_type[j])
                        state.ilai[j] = _to_real(tok[1], state.ilai[j])
                        state.iclump[j] = _to_real(tok[2], state.iclump[j])
                        state.idzveg[j] = _to_real(tok[3], state.idzveg[j])
                        state.irho[j] = _to_real(tok[4], state.irho[j])
                        state.itau[j] = _to_real(tok[5], state.itau[j])
                        state.ialp[j] = _to_real(tok[6], state.ialp[j])
                        state.ieps[j] = _to_real(tok[7], state.ieps[j])
                        if abs(state.ilai[j]) <= EPS:
                            state.ilai[j] = 0.5
                        if abs(state.iclump[j]) <= EPS:
                            state.iclump[j] = 1e-2

                    if abs(state.isigfh) <= EPS:
                        state.vegh_type = 0
                        state.veg_flagh = 0

                    if state.vegh_type == 12:                    # ice cap / neve permanente
                        state.veg_flagh = 0

    # ------------------------------------------------------------------
    # A transcricao para aqui: linha 534 do original (fim do bloco de
    # vegetacao alta). Falta:
    #   STEP 1   (536-669)   camadas de solo, USDA->USCS, bloco `nlayers1`
    #                        e perfil inicial de temperatura/umidade
    #   STEP 1   (672-696)   ramo multiponto (single_multi_flag == 1)
    #   STEP 2   (702-822)   le o arquivo meteorologico
    #   STEP 3   (825-1153)  abre unidades, hot start, lacunas, aloca
    #   STEP 4   (1156-1425) lacos principais
    #   STEP 6   (1433-1480) fecha arquivos, desaloca, tempo de CPU
    # ------------------------------------------------------------------
    unit4.close()
    raise NotImplementedError(
        "transcricao de fasst_driver.F90 em andamento: parou na linha 534 "
        "do original (fim do bloco de vegetacao alta)"
    )

    et = time.process_time()                                   # call cpu_time(et)
    print()
    print(" total run time: ", et - it)


# ---------------------------------------------------------------------------
# Verificacao executavel (cresce a cada dia do cronograma)
# ---------------------------------------------------------------------------
def demo():
    """Checagens minimas do que ja esta transcrito. `python -m fasst.fasst_driver --demo`"""
    import tempfile
    import os

    # ---- _strip_comment: ponto de atencao 13 ----
    assert _strip_comment("met.dat   ! o arquivo met") == "met.dat"
    assert _strip_comment("  met.dat !x") == "  met.dat"      # brancos a esquerda ficam
    assert _strip_comment("met.dat") == ""                    # SEM `!` -> vazio, nao a linha
    assert _strip_comment("!tudo comentario") == ""
    assert _strip_comment("") == ""

    # ---- _ch: `s(i:i)` com preenchimento de brancos ----
    assert _ch("abc", 1) == "a"
    assert _ch("abc", 3) == "c"
    assert _ch("abc", 4) == " "                               # alem do texto
    assert _ch("abc", 0) == " "

    # ---- _list_items: separadores e valores nulos ----
    assert _list_items("1 2 3") == ["1", "2", "3"]
    assert _list_items("1,2,3") == ["1", "2", "3"]
    assert _list_items("1, ,3") == ["1", None, "3"]            # virgula dupla -> nulo
    assert _list_items("1,,3") == ["1", None, "3"]
    assert _list_items(",5") == [None, "5"]                    # virgula inicial -> nulo
    assert _list_items("1 2,") == ["1", "2"]                   # virgula final NAO e nulo
    assert _list_items("   ") == []                            # registro em branco: sem itens
    assert _list_items("") == []

    # ---- _Unit.read_star: semantica de registro ----
    d = tempfile.mkdtemp()
    p = os.path.join(d, "rec.inp")
    with open(p, "w") as f:
        # registro 1 tem so 2 itens para uma lista de 3: a leitura continua
        # no registro 2, e o que sobra DELE ("4 5") e descartado.
        f.write("1 2\n3 4 5\n6   ! comentario\n7 / 9\n")
    u = _Unit(p)
    assert u.read_star(3) == ["1", "2", "3"]
    assert u.read_star(1) == ["6"]                             # "4 5" foi descartado
    assert u.read_star(3) == ["7", None, None]                 # `/` encerra a leitura
    u.close()

    # ---- conversores list-directed ----
    assert _to_int("42") == 42
    assert _to_int(None, 7) == 7                               # valor nulo: inalterado
    assert _to_real("1.5d2") == 150.0                          # expoente `d` do Fortran
    assert _to_real(None, 2.5) == 2.5

    # ---- cabecalho de fassti.inp: STEP 0 + linhas 178-266 ----
    p2 = os.path.join(d, "fassti.inp")
    with open(p2, "w") as f:
        f.write("c:/met/yuma.met              ! arquivo met\n")
        f.write("0 0 0.5 1                    ! infer_test, single_multi, mgap, keep\n")
        f.write("1.5 0.0 1 1 0 100 0          ! gwl, vegint, nprint, flprint, sprint, frozen, mstflag\n")
    st = FasstState()
    u = _Unit(p2)
    assert _strip_comment(u.read_a()) == "c:/met/yuma.met"
    t = u.read_star(4)
    st.infer_test = _to_int(t[0], st.infer_test)
    st.single_multi_flag = _to_int(t[1], st.single_multi_flag)
    st.mgap = _to_real(t[2], st.mgap)
    assert (st.infer_test, st.single_multi_flag, st.mgap, _to_int(t[3])) == (0, 0, 0.5, 1)
    t = u.read_star(7)
    st.gwl = _to_real(t[0], st.gwl)
    st.mstflag = _to_int(t[6], st.mstflag)
    assert st.gwl == 1.5
    assert [_to_int(t[i]) for i in (2, 3, 4, 5)] == [1, 1, 0, 100]
    assert st.mstflag == 0
    u.close()

    # ---- `data` do original ----
    assert _MODIS_TYPE[14] == 1 and _MODIS_TYPE[20] == 9       # linhas 102-104
    assert _UMD_TYPE[1] == 3 and _UMD_TYPE[14] == 15           # linha 106
    assert MET_COL == 32 and MAXCOL == 35                      # ponto de atencao 8

    # ---- STEP 1 ponto unico (linhas 268-386), pelo caminho de main() ----
    def _run(linhas):
        """Roda `main` ate onde a transcricao chegou e devolve o estado."""
        import contextlib
        import io as _io
        caminho = os.path.join(d, "caso.inp")
        with open(caminho, "w") as f:
            f.write(linhas)
        st = FasstState()
        saida = _io.StringIO()
        with contextlib.redirect_stdout(saida):
            try:
                main(["fasst_driver", caminho], _state=st)
            except NotImplementedError:
                pass                                  # parou na linha 534
        return st, saida.getvalue()

    _CABECALHO = (
        "yuma.met                     ! met\n"
        "0 0 0.5 1                    ! infer_test, single_multi, mgap, keep\n"
        "1.5 0.0 1 1 0 100 0          ! gwl, vegint, nprint..., mstflag\n"
        "1234                         ! vitd_index\n"
        "saida.fst                    ! output1\n"
        "saida.grd                    ! output2\n"
        "saida.flx                    ! output3\n"
        "saida.vgt                    ! output4\n"
        "saida.snw                    ! output5\n"
    )

    # terra (water_flag = 0). Note as duas ultimas linhas: com water_flag
    # = 0 a leitura CONTINUA no bloco de vegetacao (linhas 388-534).
    st, _ = _run(_CABECALHO + (
        "30.0 180.0                   ! slope, aspect\n"
        "0.01                         ! rough\n"
        "0                            ! water_flag\n"
        "0.3 0.05                     ! hsaccum, iswe\n"
        "0.02                         ! hi\n"
        "0                            ! vegl_type\n"
        "0                            ! vegh_type\n"
    ))
    assert st.vitd_index == 1234
    assert st.slope_fasst == 30.0 and st.aspect == 180.0
    assert abs(st.sloper - 0.5235987756666667) < 1e-15        # 30 graus com PI truncado
    assert st.rough == 0.01
    assert st.water_flag == 0
    assert st.hsaccum == 0.3 and st.iswe == 0.05 and st.hi == 0.02
    assert st.nlayers == 0                                     # bloco de agua nao rodou

    # agua aberta com profundidade informada (linhas 344-365)
    st, _ = _run(_CABECALHO + (
        "0.0 0.0\n"
        "0.01\n"
        "1                            ! water_flag\n"
        "1                            ! wtype: lagos/lagoas\n"
        "5.0                          ! wdepth\n"
        "0.0                          ! wvel\n"
        "0.0 0.0\n"
        "0.0\n"
    ))
    assert st.water_flag == 1 and st.nlayers == 3
    assert list(st.soiltype[1:4]) == [26, 26, 7]
    assert st.lthick[1] == 2e-2 and st.lthick[2] == 5.0 - 2e-2 and st.lthick[3] == 1.0

    # wtype = 2 (rios) muda water_flag para 2
    st, _ = _run(_CABECALHO + "0.0 0.0\n0.01\n1\n2\n5.0\n0.5\n0.0 0.0\n0.0\n")
    assert st.water_flag == 2

    # `wdepth` nulo (virgula sozinha) mantem o sentinela SPFLAG e cai no
    # `else`: lagos -> 10 m, rios -> 3 m
    st, _ = _run(_CABECALHO + "0.0 0.0\n0.01\n1\n1\n,\n0.0\n0.0 0.0\n0.0\n")
    assert st.lthick[2] == 1e1 - 2e-2
    st, _ = _run(_CABECALHO + "0.0 0.0\n0.01\n1\n2\n,\n0.0\n0.0 0.0\n0.0\n")
    assert st.lthick[2] == 3e0 - 2e-2

    # infer_test = 1 pula o bloco de neve/gelo e consome a linha de ro1o
    st, _ = _run(
        "yuma.met    ! met\n"
        "1 0 0.5 1   ! infer_test = 1\n"
        "velho_ro_soil.dat ! ro1o\n"
        "1.5 0.0 1 1 0 100 0\n"
        "1234\n" + "saida.fst !\n" * 5 +
        "0.0 0.0\n0.01\n0\n"
    )
    assert st.istart == 2 and st.hsaccum == 0.0 and st.hi == 0.0

    # ---- as validacoes que param o programa (6 pontos do STEP 1) ----
    def _para(linhas, trecho):
        import contextlib
        import io as _io
        caminho = os.path.join(d, "ruim.inp")
        with open(caminho, "w") as f:
            f.write(linhas)
        saida = _io.StringIO()
        try:
            with contextlib.redirect_stdout(saida):
                main(["fasst_driver", caminho], _state=FasstState())
        except SystemExit as e:
            assert e.code == 1
            assert trecho in saida.getvalue(), saida.getvalue()
            return
        raise AssertionError("deveria ter parado: " + trecho)

    _para("m.met !\n2 0 0.5 1 !\n", "infer test flag must be 0 or 1")
    _para(_CABECALHO + "91.0 0.0\n", "Slope can not be greater than 90")
    _para(_CABECALHO + "0.0 361.0\n", "Aspect can not be greater than +/-360")
    _para(_CABECALHO + "0.0 0.0\n0.01\n2\n", "Open water flag must be 0 or 1")
    _para(_CABECALHO + "0.0 0.0\n0.01\n0\n-1.0 0.0\n", "snow depth can not be less than 0")
    _para(_CABECALHO + "0.0 0.0\n0.01\n0\n0.0 -1.0\n", "snow water equivalent can not")
    _para(_CABECALHO + "0.0 0.0\n0.01\n0\n0.0 0.0\n-1.0\n", "ice thickness can not be less")

    # ---- STEP 1 vegetacao (linhas 388-534) ----
    _TERRA = (
        "0.0 0.0                      ! slope, aspect\n"
        "0.01                         ! rough\n"
        "0                            ! water_flag\n"
        "0.0 0.0                      ! hsaccum, iswe\n"
        "0.0                          ! hi\n"
    )
    _FOLHAGEM = (
        "0.5    ! isigfl\n"
        "0.95   ! iepf\n"
        "0.8    ! ifola\n"
        "10.0   ! ihfol\n"
    )
    _DOSSEL = (
        "2.0 0.9   ! izh, isigfh\n"
        + "1 2.0 1.0 0.5 0.1 0.1 0.8 0.95   ! camada do dossel\n" * 3
    )

    def _veg(baixa, alta, folhagem=_FOLHAGEM, dossel=_DOSSEL):
        """Monta um caso com `vegl_type` = baixa e `vegh_type` = alta."""
        texto = _CABECALHO + _TERRA + str(baixa) + " ! vegl_type\n"
        if baixa != 0:
            texto += folhagem
        texto += str(alta) + " ! vegh_type\n"
        if alta != 0:
            texto += dossel
        return _run(texto)[0]

    # padrao FASST: passa sem mapeamento, biome_source fica 0
    st = _veg(5, 0)
    assert st.vegl_type == 5 and st.veg_flagl == 1 and st.biome_source == 0
    assert st.isigfl == 0.5 and st.iepf == 0.95 and st.ifola == 0.8 and st.ihfol == 10.0

    # MODIS/IGBP e UMD/LDAS
    st = _veg(1003, 0)
    assert st.vegl_type == _MODIS_TYPE[3] == 4
    assert st.biome_source == 1000
    # `new_vtl` nao e so uma variavel de desvio: initsurface.py:269 o
    # repassa a `veg_propl`, entao o valor guardado importa (ponto 4a)
    assert st.new_vtl == 3
    st = _veg(2001, 0)
    assert st.vegl_type == _UMD_TYPE[1] == 3
    assert st.biome_source == 2000 and st.new_vtl == 1
    # ponto 4b: a vegetacao baixa exclui MODIS 13/16 e UMD 13
    for entrada in (1013, 1016, 2013):
        st = _veg(entrada, 0)
        assert st.vegl_type == 0 and st.veg_flagl == 0, entrada

    # os tres ramos de agua da vegetacao baixa funcionam (ao contrario do
    # equivalente da vegetacao alta -- ponto de atencao 3)
    st = _veg(1011, 0)                                   # MODIS 11 -> 13 swamp/bog
    assert st.water_flag == 1 and st.vegl_type == 13 and st.veg_flagl == 1
    st = _veg(1017, 0)                                   # MODIS 17 -> 14 inland water
    assert st.water_flag == 1 and st.vegl_type == 0 and st.veg_flagl == 0
    st = _veg(2014, 0)                                   # UMD 14 -> 15 ocean
    assert st.water_flag == 3 and st.vegl_type == 0 and st.veg_flagl == 0

    assert _veg(1015, 0).veg_flagl == 0                  # MODIS 15 -> 12 ice cap
    assert _veg(1015, 0).vegl_type == 12                 # ...mas o tipo fica

    # limite de 0.98 na densidade de folhagem, e o sentinela que escapa dele
    assert _veg(5, 0, folhagem="0.99 !\n0.9 !\n0.8 !\n1.0 !\n").isigfl == 0.98
    assert _veg(5, 0, folhagem=", !\n0.9 !\n0.8 !\n1.0 !\n").isigfl == SPFLAG
    # isigfl = 0 desliga a vegetacao baixa depois de ja ter sido mapeada
    st = _veg(5, 0, folhagem="0.0 !\n0.9 !\n0.8 !\n1.0 !\n")
    assert st.vegl_type == 0 and st.veg_flagl == 0

    # ---- vegetacao alta ----
    st = _veg(0, 7)
    assert st.vegh_type == 7 and st.veg_flagh == 1
    assert st.izh == 2.0 and st.isigfh == 0.9
    assert list(st.ilai[1:4]) == [2.0, 2.0, 2.0]
    assert list(st.ieps[1:4]) == [0.95, 0.95, 0.95]

    # os limites de ilai e iclump quando vem zero
    st = _veg(0, 7, dossel="2.0 0.9 !\n" + "1 0.0 0.0 0.5 0.1 0.1 0.8 0.95 !\n" * 3)
    assert list(st.ilai[1:4]) == [0.5, 0.5, 0.5]
    assert list(st.iclump[1:4]) == [1e-2, 1e-2, 1e-2]

    # PONTO DE ATENCAO 3 -- ramo morto: MODIS 17 -> 14 ("inland water")
    # NAO produz water_flag = 1; cai no `else` do ocean e termina com
    # water_flag = 3. E o `wdepth = 10` do ramo morto seria sobrescrito
    # por 3 de qualquer forma.
    st = _veg(0, 1017)
    assert st.water_flag == 3, "ponto 3: o ramo inland water deveria estar morto"
    assert st.vegh_type == 0 and st.veg_flagh == 0
    st = _veg(0, 2014)                                   # UMD 14 -> 15 ocean
    assert st.water_flag == 3 and st.vegh_type == 0 and st.veg_flagh == 0

    # PONTO DE ATENCAO 4b: sem exclusao de 13/16 do lado da vegetacao alta.
    # MODIS 13 -> 0, mas veg_flagh CONTINUA 1 -- estado que o lado da
    # vegetacao baixa nunca produz.
    st = _veg(0, 1013)
    assert st.vegh_type == _MODIS_TYPE[13] == 0
    assert st.veg_flagh == 1, "ponto 4b: veg_flagh nao e zerado aqui"
    # O 16 e o caso que distingue de verdade os dois lados: a vegetacao
    # baixa zera (testado acima, com 1016), a alta MAPEIA para 8.
    st = _veg(0, 1016)
    assert st.vegh_type == _MODIS_TYPE[16] == 8, "ponto 4b: a alta nao exclui o 16"
    assert st.veg_flagh == 1

    # PONTO DE ATENCAO 14: `biome_source` vaza da vegetacao baixa para a
    # alta. Baixa MODIS (1003) + alta no padrao FASST (5) ->
    # max(0, 5 - 1000) = 0 -> a vegetacao alta e descartada em silencio.
    st = _veg(1003, 5)
    assert st.vegl_type == 4 and st.veg_flagl == 1
    assert st.biome_source == 1000
    assert st.new_vth == 0
    assert st.vegh_type == 0 and st.veg_flagh == 0, "ponto 14: alta descartada"
    # ...mas os 4 registros do dossel continuam sendo CONSUMIDOS do
    # arquivo: o `if(vegh_type /= 0)` da linha 471 foi avaliado com o valor
    # original, e so o conteudo interno zerou o tipo. Quem monta o .inp
    # ainda precisa fornece-los, e os valores vao parar em ilai/ieps.
    assert list(st.ilai[1:4]) == [2.0, 2.0, 2.0]
    assert st.izh == 2.0
    # na ordem inversa nao acontece: a baixa e processada primeiro
    st = _veg(5, 1003)
    assert st.vegl_type == 5 and st.vegh_type == _MODIS_TYPE[3] == 4
    assert st.veg_flagh == 1

    assert _veg(2000, 0).vegl_type == 0                  # new_vtl = 0 -> ramo `else`

    # SEM CHECAGEM AINDA: os pontos de atencao 4a e 14b. O STEP 0 zera
    # `biome_source` no inicio do programa, entao um `biome_source` herdado
    # so aparece a partir da SEGUNDA iteracao do laco multiponto -- que so
    # existe no Dia 9. Nao da para demonstrar nenhum dos dois aqui sem
    # simular um caminho que o codigo ainda nao tem, e um teste assim
    # passaria sem provar nada.

    # divergencia documentada: indice fora da tabela levanta IndexError
    # (no Fortran seria leitura de memoria fora do array)
    try:
        _veg(1025, 0)                                    # new_vtl = 25 > 20
    except IndexError:
        pass
    else:
        raise AssertionError("esperado IndexError para new_vtl fora da tabela")

    print("demo OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
    else:
        main(sys.argv)
