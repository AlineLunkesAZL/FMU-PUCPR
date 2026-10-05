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
   unidade 90 -- nao serve aqui, e nao sera refatorado.

Estado desta transcricao
--------------------------
EM ANDAMENTO (cronograma de 10 dias uteis). Hoje: esqueleto, docstring e
os 5 stubs de dependencia externa. `main()` levanta `NotImplementedError`
-- os STEPs sao preenchidos nos dias 2-9. O modulo importa limpo de
proposito, para que a verificacao do pacote possa rodar desde o Dia 1.

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
    linhas 510-511 sobrescrevem `wdepth = 3`/`wvel = 0` nos dois casos.
4.  Assimetria `new_vtl` vs `new_vth`. Linha 399:
    `new_vtl = vegl_type - biome_source`. Linha 479:
    `new_vth = max(0, vegh_type - biome_source)`. So um tem o `max`.
5.  Vazamento do `i` (1063 + 1095-1099). Depois de `i = i - 1`, o laco de
    contagem de lacunas indexa `metm(i,k,...)`. No ramo ponto-unico `i` e o
    sitio casado; no multiponto e um contador de linhas da unidade 5, sem
    relacao nenhuma. Python vai acidentalmente consertar OU acidentalmente
    quebrar isso.
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
    justamente para consertar esse efeito.

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
def main(argv=None):
    """`program fasst_driver` (fasst_driver.F90, linhas 1-1480).

    EM ANDAMENTO -- so o esqueleto existe. Os STEPs abaixo marcam as
    fronteiras do original e a ordem em que serao preenchidos.
    """
    it = time.process_time()                                   # call cpu_time(it)

    # STEP 0   (109-176)   inicializa variaveis
    # STEP 1   (178-699)   le o arquivo de controle
    # STEP 2   (702-822)   le o arquivo meteorologico
    # STEP 3   (825-1153)  abre unidades, hot start, lacunas, aloca
    # STEP 4   (1156-1425) lacos principais (ponto unico / multiponto)
    # STEP 6   (1433-1480) fecha arquivos, desaloca, tempo de CPU
    raise NotImplementedError(
        "transcricao de fasst_driver.F90 em andamento: STEP 0-6 ainda nao "
        "preenchidos (ver 'Estado desta transcricao' na docstring do modulo)"
    )

    et = time.process_time()                                   # call cpu_time(et)
    print()
    print(" total run time: ", et - it)


if __name__ == "__main__":
    main(sys.argv)
