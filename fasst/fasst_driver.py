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
         cabecalho de `fassti.inp` (linhas 109-266)
  Dia 3  saidas, declive, rugosidade, agua aberta, neve/gelo (268-386)
  Dia 4  vegetacao baixa e alta, mapeamento MODIS/UMD (388-534)
  Dia 5  camadas de solo, bloco de fusao, perfis iniciais, ramo
         multiponto, `sort2` (536-699 e 1483-1534)
  Dia 6  STEP 2 -- leitura do arquivo meteorologico (702-822)
  Dia 7  STEP 3 -- as 16 unidades de arquivo, os cabecalhos das
         saidas e o bloco de hot start (825-1084)
  Dia 8  STEP 3 -- lacunas de tempo, `maxlines`, `allocate_run_arrays`
         (1086-1153)  <-- aqui
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

CORRECAO ao paragrafo acima: a assinatura de `read_old_data` E a unica das
cinco que VAI mudar quando a transcricao chegar. O `.F90` dela faz
`use fasst_global` de mais de 100 variaveis e ESCREVE em boa parte --
inclusive `moverlapr`, que o proprio driver le tres linhas depois da
chamada (`if(jline == 'endr') moverlapr = i`). Como nada disso passa por
argumento, a versao Python vai precisar de `state` como primeiro
parametro, no padrao do resto do pacote. O stub herdado de
`initsurface.py` nao tem esse parametro, e trocar isso mexeria nos
chamadores em `initsurface.py` e `initwater.py` -- fica para quando
`read_old_data` for transcrita de verdade. As outras quatro ja recebem
`state`.

Nota relacionada: `moverlapr` nao e zerado por ninguem -- nem pelo STEP 0,
nem por `zero_parameters`. Sai de `read_old_data` e fica.

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
    Alcancavel so quando a vegetacao baixa deu swamp/bog (13) ou inland
    water (14), o que preencheu `nlayers1`/`soiltype1`/`lthick1`, E o
    primeiro ramo (`wflag2 == 1 .and. nlayers <= 1`) nao pegou. Cinco
    defeitos distintos, todos tracados a mao indice por indice e
    confirmados contra o codigo (tres casos no `demo()`):
    a) linhas 581-583 e 603-605 TESTAM `soiltype(i-nlayers)` e ATRIBUEM de
       `lthick(i-nlayers1)`/`soiltype(i-nlayers1)` -- variaveis
       diferentes. Com `nlayers == nlayers1` o defeito e invisivel; fora
       disso o bloco detecta agua numa camada e copia outra. Com
       `nlayers > nlayers1 + 1` o indice do teste fica NEGATIVO, e com
       `nlayers == nlayers1 + 1` fica 0 (slot fantasma aqui, memoria antes
       do array no Fortran).
    b) `k` sai de `sort2` como a contagem de camadas deduplicadas e e
       SOBRESCRITO pelo `k = i - j` da linha 598. O valor original se
       perde, e e a versao nova (`kk - j`) que as linhas 609, 616 e 618
       usam.
    c) a linha 613 avanca `d1` com `lthick3(i)` enquanto a 612 ajusta
       `lthick4(i)` -- o laco analogo de 591-594 usa o mesmo array nas
       duas linhas.
    d) o laco final (617-625) le `lthick3`/`soiltype3` ate o `k` NOVO, mas
       esses arrays so tem entradas validas ate o `k` antigo. As posicoes
       extras valem 0 -- definidas, nao lixo, porque `sort2` zera `b`/`b1`
       de 1 a maxn antes de tudo. O `nlayers` resultante tambem pode
       passar de `maxl`.
    e) consequencia das anteriores: os VALORES que o `else if` do primeiro
       laco grava nunca chegam a saida. Quando ele dispara, `j` passa a ser
       aquele `i`, o `k` final encolhe e o laco 617-625 le so as primeiras
       `kk-j` posicoes de `lthick3` (as menores, tipicamente as camadas do
       swamp); e o offset `d2` se cancela em 611-614, porque ali
       `d1 == d2`. So o TESTE do `else if` tem efeito observavel, pelo `j`
       que ele move. Verificado por mutacao: trocar o indice das duas
       atribuicoes nao muda nenhuma saida.
    Ainda: `sort2` declara as saidas `b(maxn)`/`b1(maxn)` com as entradas
    `arr(maxn+3)`/`arr1(maxn+3)`. Transcrever indice por indice, nunca
    "pelo que faz sentido".
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
    (`MetCol.HI`) nunca e lida do arquivo e fica com `mflag`. Com
    `ncols <= 32` as colunas 34 e 35 tambem ficam.
15. Os conjuntos "lido" e "quantizado" nao batem no ramo `ncols > met_col`
    (796-800): a leitura cobre as colunas 1..32, 34 e 35, mas o laco de
    quantizacao vai de 1 a `met_col+2` = 34. Resultado, com as tres pontas
    observaveis:
    * a coluna 34 (VIS) e lida E quantizada;
    * a coluna 35 (AER) e lida e NAO quantizada -- fica com o valor cru;
    * a coluna 33 (HI) nao e lida mas E quantizada, e por causa do ponto 16
      isso tira o sentinela da igualdade: `mflag = -999.0` sai de la como
      -999.0000000000001. E provavelmente por isso que todo teste de dado
      faltante no codigo usa TOLERANCIA (linha 1176:
      `aint(dabs(met(k,ip_zen)-mflag)*1d5)*1d-5 > eps`) e nunca igualdade.
16. `anint(x*1d15)*1d-15` (linhas 799 e 804) NAO e identidade. Nem `1d15`
    nem `1d-15` sao exatos em binario, e para x acima de ~9 o produto
    `x*1d15` ja passa de 2^53, onde o espacamento entre doubles e maior
    que 1 -- entao a ida e volta perde bits. 7.0 vira 7.000000000000001,
    12.0 vira 12.000000000000002, 1013.0 vira 1013.0000000000001; outros
    (0.2, 180.0) sobrevivem por sorte de arredondamento. Ou seja: TODO
    valor meteorologico entra no modelo levemente diferente do que esta no
    arquivo. O mesmo vale para o `anint(.,10)` do cabecalho e para o
    `anint(.,20)` do `sloper` (onde 1.57 volta 1.5699999999999998). Nao e
    defeito de transcricao -- e a aritmetica IEEE do original, preservada.
9.  `tstps == 0` quando `timstep <= 0` (linha 1092): `rtstps` fica 0, logo
    `msteps1 = 0` sempre e nenhuma lacuna e detectada. Note que `tstps` sai
    de `int(rtstps)`, truncando: com `timstep = 9` o quociente 2.666... vira
    2, nao 3.
17. `modyo` (linha 1102) e calculado e NUNCA lido. A unica linha que o
    usaria -- `if(mody /= modyo) daylim = 1.07d0`, 1105 -- esta comentada no
    proprio Fortran. Custo de calculo sem efeito.
19. O ramo bissexto do laco de lacunas (fator 367, linhas 1109-1111) e
    INALCANCAVEL para qualquer ano real. `mody` e calculado sobre o ano JA
    QUANTIZADO pelo ponto 16, e `anint(2012*1d15)*1d-15` devolve
    2012.0000000000002274, logo `mody = 2.27e-13` -- muito acima de
    `eps = 2.22e-16`. Os unicos anos que passam o teste sao os multiplos de
    4 pequenos o bastante para `ano*1d15` caber exato num double (< 2^53):
    0, 4, 8 e 16 (o 12 NAO passa, da 1.78e-15). Em dado de verdade, o laco
    sempre usa 366 -- de acordo com `met_date`, que faz o mesmo teste sobre
    o mesmo valor quantizado, entao os dois lados sao consistentes.
20. `mtsteps` nao e uma contagem de registros faltantes: e uma aproximacao
    cujo arredondamento depende de QUAIS horas a lacuna cai entre. Com dado
    horario, uma lacuna de um registro da produto 1.999999974 (-> 1) em 21
    dos 22 pares de horas possiveis, e 2.000000852 (-> 2) entre 20h e 22h,
    porque o erro do ponto 16 muda de sinal. Na virada de ano o vies e
    maior: uma hora faltando entre 365/22h e 001/00h conta 25 passos,
    porque `met_date` divide o dia-do-ano por 366 e nao por 365. E
    justamente para mascarar isso que `daylim` sobe para 25.07 (linha 1106).
21. Os dois `msteps1 = 0` (linhas 1088 e 1108) sao armazenamentos MORTOS:
    os dois ramos seguintes atribuem `msteps1` incondicionalmente, e depois
    do laco ninguem le `msteps1` -- so `mtsteps` sobrevive.
10. Linhas 565-572: `soiltype(2) = soiltype(1)` acontece ANTES de
    `soiltype(1) = 26`, o que parece dependencia de ordem -- mas NAO e.
    Para chegar a este ramo e preciso `wflag2 == 1` com `nlayers <= 1`, e
    `wflag2` so vira 1 quando alguma camada lida tem tipo 26: com uma
    camada so, essa camada E a de tipo 26, entao `soiltype(2)` recebe 26
    de qualquer jeito. Mantido na lista porque a leitura inicial sugere o
    contrario; conferido, nao ha defeito aqui.
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
* VALOR FINAL DE VARIAVEL DE LACO -- conferir em cada dia que falta. Um
  `do i=1,n` do Fortran termina com `i == n+1` (ou com o valor inicial, se
  o laco nao rodar nenhuma vez); um `for i in range(1, n+1)` do Python
  termina com `i == n`. Qualquer leitura de `i`/`j`/`k`/`kk` DEPOIS de um
  laco difere em 1 entre as duas linguagens. Auditado nas linhas 109-822:
  nenhum `for` transcrito ate agora tem a variavel lida depois do laco, e
  os pontos onde o original depende disso (o laco de sitios do STEP 2, o
  `i = i - 1` da linha 1063) sao `do while` com incremento explicito nos
  dois lados, entao sao fieis. CONFERIDO ate a linha 1153: o ponto de
  atencao 5 era justamente essa fronteira -- o laco de lacunas (1095-1117)
  le `i` de fora e e um `for` que nao vaza `k`, porque `k` nao e lido depois
  dele. Falta conferir o STEP 4 (Dia 9), que tem `do k=istart,iend` e
  `do iw=istart,iend` com `iw` em `fasst_global`.
"""

import math
import sys
import time

import numpy as np

from .constants import EPS, MAXCOL, MAXL, MAXN, MOVERLAP, NCLAYERS, PI, SPFLAG, TREF
from .indices import MetCol
from .state import FasstState, _chars, _ints, _reals, _reals2
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


def _open_or_stop(path, mode, msg):
    """`open(unit=N,file=...,iostat=io)` mais o bloco de erro que sempre vem
    depois dele: mensagem com o `iostat`, uma linha em branco, o nome do
    arquivo, e `stop`.

    Os 12 blocos de abertura dos STEPs 2 e 3 tem exatamente essa forma,
    variando so o texto -- inclusive os espacos em volta do "error", que sao
    inconsistentes no original (`error = `, `error=`, com e sem espaco
    inicial) e ficam preservados no `msg` de cada chamada.

    Ver divergencia D1: no original, `status='unknown'` CRIA o arquivo que
    falta, o que torna estes ramos de erro praticamente mortos para arquivos
    de ENTRADA. Aqui os de entrada ('r'/'rb') sao abertos para leitura e a
    falha cai justamente neste ramo; os de SAIDA ('w'/'wb') continuam sendo
    criados.
    """
    try:
        return _Unit(path, mode)
    except OSError as e:
        print(msg % (e.errno or 1))
        print(" ")
        print(path)
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


def _to_char(tok, current=" ", width=1):
    """Item list-directed de caractere.

    `None` (valor nulo) deixa a variavel INALTERADA. Aspas simples ou
    duplas em volta sao removidas, e o resultado e truncado no comprimento
    declarado da variavel Fortran (`character(len=width)`).
    """
    if tok is None:
        return current
    if len(tok) >= 2 and tok[0] in "'\"" and tok[-1] == tok[0]:
        tok = tok[1:-1]
    return tok[:width]


def _ix(i):
    """Indice 1-based calculado exatamente como o original calcula.

    Divergencia deliberada D4, usada nos dois pontos em que o original
    indexa com uma conta que pode dar negativo: o bloco de fusao de
    camadas (ponto de atencao 1) e o laco de lacunas de tempo do STEP 3,
    que indexa `metm` com o `i` que sobrou do laco anterior (ponto de
    atencao 5). O indice 0 existe aqui (slot fantasma, valor 0) e e
    inofensivo, mas um indice NEGATIVO faria o Python pegar o fim do array
    por wraparound, enquanto o Fortran leria memoria antes do array. Os
    dois sao lixo; falhar alto e mais util do que produzir numeros
    plausiveis e diferentes dos do Fortran.
    """
    if i < 0:
        raise IndexError(
            "indice %d < 0: no bloco de fusao de camadas (ponto de atencao "
            "1) `i - nlayers` fica negativo com nlayers > nlayers1 + 1; no "
            "laco de lacunas do STEP 3 (ponto de atencao 5) o `i` que sobra "
            "do laco anterior vale -1 quando num_frq == 0" % i
        )
    return i


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
        """`read(u,'(a)') s` -- um registro inteiro, sem o fim de linha.

        Sem `iostat=`: o fim de arquivo e erro de runtime, como no Fortran.
        """
        line = self.f.readline()
        if line == "":
            raise EOFError("fim do arquivo em read('(a)'): " + str(self.path))
        return line.rstrip("\n").rstrip("\r")

    def read_a_io(self, current=""):
        """`read(u,'(a)',iostat=io) s` -- devolve `(s, io)`.

        Com `iostat=`, o fim de arquivo nao e erro: `io` volta -1 e a
        variavel fica INALTERADA, por isso o `current` de volta. Isso nao e
        detalhe: o laco de contagem de linhas do STEP 2 (738-754) continua
        escaneando `a` depois do EOF, e como `a` ainda tem a ULTIMA linha de
        dados, ela e contada DUAS vezes. O `ntlines - 1` da linha 759 e o
        que compensa esse duplo.
        """
        line = self.f.readline()
        if line == "":
            return current, -1
        return line.rstrip("\n").rstrip("\r"), 0

    def read_star_io(self, n):
        """`read(u,*,iostat=io) v1,...,vn` -- devolve `(tokens, io)`.

        `io` = -1 no fim de arquivo; nesse caso os itens que faltaram voltam
        como `None`, o que faz `_to_int`/`_to_real` manterem o valor atual
        (= variavel inalterada, como no Fortran).
        """
        try:
            return self.read_star(n), 0
        except EOFError:
            return [None] * n, -1

    def rewind(self):
        """`rewind(u)` -- volta ao primeiro registro."""
        self.f.seek(0)

    # Passagens diretas, para os modulos que recebem a unidade e leem/escrevem
    # por conta propria: `us_soil_tools.get_user_soil_params` chama
    # `.readline()` na unidade 90, e as unidades de saida recebem `.write()`.
    def readline(self):
        return self.f.readline()

    def write(self, s):
        return self.f.write(s)

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
# allocate(...) do STEP 3 (fasst_driver.F90, linhas 1126-1148)
# ---------------------------------------------------------------------------
# As 34 alocaveis `real(dp)` de uma dimensao, na ordem em que o original as
# aloca. As outras seis sao casos a parte (tipo ou forma diferentes) e estao
# escritas uma a uma em `allocate_run_arrays`.
_RUN_REAL_1D = (
    "sdens",
    "ft", "tt",
    "surfice", "surfmoist",
    "surficep", "surfmoistp",
    "surfci", "surfrci",
    "surfcbr", "surfd",
    "frthick", "twthick",
    "surfemis", "surfemisf",
    "surfemisc",
    "tot_moist",
    "tot_thick",
    # fluxos de energia de superficie (1140-1148)
    "cheat", "cheat1",
    "sdown", "sup",
    "irdown", "irup",
    "pheat1", "lheat",
    "evap_heat", "sheat",
    "melt", "tmelt",
    "lheatf", "lhes",
    "cevap", "ponding",
)


def allocate_run_arrays(state):
    """Aloca as ~40 arrays alocaveis de `fasst_global` (linhas 1126-1148).

    Todas sao dimensionadas por `state.maxlines`, que o bloco imediatamente
    anterior no STEP 3 acabou de calcular. `fasst/state.py` as deixa em
    `None` de proposito e delega a alocacao a quem transcreve o driver --
    ou seja, a esta funcao (ver o comentario da secao "allocatable arrays"
    em `state.py`).

    A ORDEM IMPORTA: isto roda ANTES de `zero_mxl_params` (Fortran: 1126-1148
    e so depois 1167/1292), que itera sobre estas arrays e precisa que elas
    existam. As formas tambem importam -- sao elas que o slicing de
    `zerovars.zero_mxl_params` assume.

    Dois tipos fora do padrao, lidos de `fasst_global.F90`: `slushy` e
    `integer(ip)` (linha 179) e `sstate` e `character(len=1)` (linha 191).
    As outras 38 sao `real(dp)`.

    O original passa `stat=err` nos 23 `allocate` e NUNCA le `err`: um
    erro de alocacao passaria em silencio. Aqui nao ha o que preservar --
    falta de memoria levanta `MemoryError` sozinho.
    """
    n = state.maxlines
    for nome in _RUN_REAL_1D:
        setattr(state, nome, _reals(n))

    state.slushy = _ints(n)                      # integer(ip), nao real
    state.sstate = _chars(n, 1)                  # character(len=1)
    state.met = _reals2(n, MAXCOL)               # met(maxlines,maxcol)
    state.dmet1 = _reals2(n, 13)                 # dmet1(maxlines,13)
    state.airt = _reals2(2, n)                   # airt(2,maxlines)
    state.canopy_temp = _reals2(NCLAYERS, n)     # canopy_temp(nclayers,maxlines)


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
def main(argv=None, _state=None, _locals=None):
    """`program fasst_driver` (fasst_driver.F90, linhas 1-1480).

    Um escopo unico e longo, de proposito -- ver "Arquitetura" na docstring
    do modulo.
    """
    if argv is None:
        argv = sys.argv

    # O `program` nao recebe estado: e ele que o cria.
    # ponytail: `_state` e `_locals` sao ganchos de teste -- deixam o
    # `demo()` inspecionar o estado E as variaveis locais pelo caminho de
    # codigo REAL enquanto `main` ainda para no meio da transcricao. Sem
    # eles, as ~700 linhas de parsing dos dias 2-6 (a classe de bug
    # dominante) ficariam quase sem verificacao: metade do que o STEP 1 e o
    # STEP 2 produzem -- `metm`, `latm`, `icnti`, `output1..5`,
    # `fasstusersoil`, `sname`, `surfdepth` -- sao LOCAIS do `program`, nao
    # campos de `FasstState`. Remover no Dia 9, quando `main` rodar ate o
    # fim e essas coisas ficarem observaveis pelos arquivos de saida.
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
    # idem, usados a partir do STEP 3. `sdensi`/`phie` sao passados a
    # `read_old_data` na linha 1069 sem nunca terem sido escritos (o
    # original so os inicializa na linha 1196, dentro do STEP 4); `kk` so e
    # atribuido quando ha casamento de indice VITD, mas e lido na 1076.
    fdi = 0
    kk = 0
    d1i = 0
    num_frq = 0
    vitd_att = " "
    t1n = 0.0
    t1oo = 0.0
    sdensi = 0.0
    phie = 0.0
    # `rtstps` e o caso que IMPORTA desta lista: o original so o atribui
    # dentro de `if(timstep > 0d0)` (linha 1092) e o le na linha seguinte
    # sem guarda nenhuma. E o ponto de atencao 9.
    rtstps = 0.0

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

                # ---- tipo de solo, numero de camadas, espessura (536-563) ----
                tok = unit4.read_star(1)
                state.nlayers = _to_int(tok[0], state.nlayers)   # numero de camadas (max = 10)
                for i in range(1, state.nlayers + 1):
                    tok = unit4.read_star(2)
                    # tipo da camada, espessura da camada (m)
                    state.soiltype[i] = _to_int(tok[0], state.soiltype[i])
                    state.lthick[i] = _to_real(tok[1], state.lthick[i])
                    if abs(state.soiltype[i]) < 100:             # tipo de solo USCS
                        state.sclass[i] = "USCS"
                    elif abs(state.soiltype[i]) >= 100:          # tipo de solo USDA
                        state.sclass[i] = "USDA"
                        if state.soiltype[i] > 0:
                            lis = state.soiltype[i] - 100
                            state.soiltype[i] = map_usda_soiltype_to_uscs(lis)
                        else:
                            state.soiltype[i] = -(abs(state.soiltype[i]) - 100)

                    if state.soiltype[i] == -1:
                        lcount = i
                        # `sname` e character(len=4): o Fortran trunca em 4.
                        # Aqui NAO ha strip de `!` -- ao contrario de todos
                        # os outros `read(4,'(a)')` deste arquivo.
                        sname[i] = unit4.read_a()[:4]

                    if state.soiltype[i] == 26:
                        state.water_flag = 1
                        wflag2 = 1
                        wtype = 1
                        wvel = 0.0

                # ---- fusao de camadas (565-626) ----
                if wflag2 == 1 and state.nlayers <= 1:
                    state.nlayers = 3
                    # `soiltype(2) = soiltype(1)` acontece antes de
                    # `soiltype(1) = 26`, mas a ordem aqui e indiferente:
                    # para chegar neste ramo e preciso wflag2 == 1 com
                    # nlayers <= 1, ou seja a UNICA camada ja e a de tipo 26
                    # -- ver ponto de atencao 10.
                    state.soiltype[2] = state.soiltype[1]
                    state.lthick[2] = state.lthick[1] - 2e-2
                    state.soiltype[1] = 26
                    state.lthick[1] = 2e-2
                    state.soiltype[3] = 7
                    state.lthick[3] = 1e0
                elif nlayers1 > 0:
                    # ==================================================
                    # PONTO DE ATENCAO 1 -- a regiao mais perigosa do
                    # arquivo. Alcancavel so quando a vegetacao baixa deu
                    # swamp/bog (13) ou inland water (14), o que preencheu
                    # nlayers1/soiltype1/lthick1. Transcrito indice por
                    # indice. Os quatro defeitos distintos estao marcados
                    # (1a) a (1d) abaixo.
                    # ==================================================
                    kk = 0
                    kk = nlayers1 + state.nlayers
                    j = 1
                    for i in range(1, kk + 1):
                        if i <= nlayers1 and soiltype1[i] == 26:
                            lthick2[i] = d1 + lthick1[i]         # d1 ainda e 0
                            soiltype2[i] = soiltype1[i]
                        elif i > nlayers1 and state.soiltype[_ix(i - state.nlayers)] == 26:
                            # (1a) O TESTE acima usa `i - nlayers`; as duas
                            # atribuicoes abaixo usam `i - nlayers1`. Sao
                            # variaveis diferentes. Quando nlayers ==
                            # nlayers1 o defeito e invisivel; fora disso o
                            # bloco testa uma camada e copia outra.
                            #
                            # (1e) Estes dois VALORES, porem, nunca chegam a
                            # saida -- conferido por mutacao e por traco
                            # manual. Quando este ramo dispara, `j` passa a
                            # ser este `i`, logo o `k` final (kk - j) encolhe
                            # e o laco de 617-625 le so as primeiras kk-j
                            # posicoes de lthick3, que sao as MENORES
                            # (tipicamente as camadas do swamp). E o offset
                            # `d2` se cancela no laco de 611-614, porque ali
                            # `d1 == d2`. Ou seja: o trabalho deste ramo e
                            # descartado. So o TESTE tem efeito observavel,
                            # via o `j` que ele move.
                            lthick2[i] = d2 + state.lthick[i - nlayers1]
                            soiltype2[i] = state.soiltype[i - nlayers1]
                        # Posicoes que nenhum dos dois ramos atingiu ficam
                        # com o valor anterior (zero aqui, indefinido no
                        # Fortran) e ainda assim sao lidas nesta linha:
                        if lthick2[i] > EPS:
                            j = i

                    k, lthick3, soiltype3 = sort2(j, lthick2, soiltype2)

                    # profundidades acumuladas -> espessuras por camada
                    d1 = 0
                    for i in range(1, k + 1):
                        lthick3[i] = lthick3[i] - d1
                        d1 = d1 + lthick3[i]
                    d2 = d1

                    for i in range(1, kk + 1):
                        # (1b) `k` acabou de sair de sort2 como a contagem de
                        # camadas deduplicadas e e SOBRESCRITO aqui, dentro
                        # do laco. O valor original se perde; e esta nova
                        # versao (no fim, kk - j) que as linhas 609, 616 e
                        # 618 usam.
                        k = i - j
                        if k > 0:
                            if i <= nlayers1 and soiltype1[i] != 26:
                                lthick2[k] = d1 + lthick1[i]
                                soiltype2[k] = soiltype1[i]
                            elif i > nlayers1 and state.soiltype[_ix(i - state.nlayers)] != 26:
                                # (1a) de novo: testa `i - nlayers`, copia
                                # `i - nlayers1`.
                                lthick2[k] = d2 + state.lthick[i - nlayers1]
                                soiltype2[k] = state.soiltype[i - nlayers1]

                    jj, lthick4, soiltype4 = sort2(k, lthick2, soiltype2)

                    for i in range(1, jj + 1):
                        lthick4[i] = lthick4[i] - d1
                        # (1c) avanca `d1` com `lthick3`, nao com `lthick4`
                        # -- ao contrario do laco analogo de 591-594, que usa
                        # o mesmo array nas duas linhas.
                        d1 = d1 + lthick3[i]

                    state.nlayers = k + jj
                    for i in range(1, state.nlayers + 1):
                        # (1d) `lthick3`/`soiltype3` so tem entradas validas
                        # ate o k ANTIGO (o de sort2), mas sao lidos ate o k
                        # novo. As posicoes extras valem 0 -- definidas, nao
                        # lixo, porque sort2 zera `b`/`b1` de 1 a maxn antes
                        # de qualquer coisa. `nlayers` tambem pode passar de
                        # maxl aqui: no Fortran isso escreve fora do array,
                        # em Python levanta IndexError.
                        if i <= k:
                            state.lthick[i] = lthick3[i]
                            state.soiltype[i] = soiltype3[i]
                        else:
                            state.lthick[i] = lthick4[i - k]
                            state.soiltype[i] = soiltype4[i - k]

                if lcount != 0:
                    fasstusersoil = unit4.read_a()
                # Ponto de atencao 13: o strip roda SEMPRE, inclusive quando
                # nada foi lido acima -- e ai o default 'fasstusersoil.inp'
                # (que nao tem `!`) vira string VAZIA. A linha 878 existe
                # exatamente para consertar isso.
                fasstusersoil = _strip_comment(fasstusersoil)

                # ---- temperaturas iniciais de solo (641-654) ----
                tok = unit4.read_star(1)
                # "y" = ha temperaturas iniciais de solo, em C
                ttest = _to_char(tok[0], ttest)
                if ttest == "y" or ttest == "Y":
                    tok = unit4.read_star(1)
                    nt0 = _to_int(tok[0], nt0)               # quantas temperaturas
                    for j in range(1, nt0 + 1):
                        tok = unit4.read_star(2)
                        state.zti[j] = _to_real(tok[0], state.zti[j])  # profundidade (m)
                        state.tm[j] = _to_real(tok[1], state.tm[j])    # temperatura (C)
                        if abs(state.tm[j]) > 100.0:
                            _stop("!!! Measured soil temperatures out of range, "
                                  "stopping. !!!")
                        state.tm[j] = state.tm[j] + TREF     # converte para Kelvin

                # ---- umidades iniciais de solo (656-668) ----
                tok = unit4.read_star(1)
                # "y" = ha umidades iniciais de solo (fracao)
                mtest = _to_char(tok[0], mtest)
                if mtest == "y" or mtest == "Y":
                    tok = unit4.read_star(1)
                    # o comentario do original diz "number of initial temps"
                    # tambem aqui (copia e cola); sao umidades
                    nm0 = _to_int(tok[0], nm0)
                    for j in range(1, nm0 + 1):
                        tok = unit4.read_star(2)
                        state.zm[j] = _to_real(tok[0], state.zm[j])    # profundidade (m)
                        state.sm[j] = _to_real(tok[1], state.sm[j])    # umidade (0.0 - 1.0)
                        if state.sm[j] > 1.0:
                            _stop("!!! Measured soil moistures out of range, "
                                  "stopping. !!!")

    elif state.single_multi_flag == 1:                       # calculo multiponto
        vitd_input = _strip_comment(unit4.read_a())
        fasstusersoil = _strip_comment(unit4.read_a())
    else:
        # Nao ha `stop` aqui: o programa segue para o STEP 2 sem ter lido
        # nada do arquivo de controle alem do cabecalho.
        print(" Model does not recognize single_multi_flag input of = %3d"
              % state.single_multi_flag)

    surfdepth = 5e-2

    unit4.close()                                            # fecha o arquivo de entrada

    # ---------------- STEP 2 (linhas 702-822) ----------------
    # Le o meteorologico do arquivo apropriado. TODOS os valores sao lidos
    # de uma vez; o modelo depois itera sobre eles. As condicoes
    # meteorologicas sao as condicoes de contorno do sistema de equacoes.

    # abre o arquivo met (linhas 708-715)
    io = 0
    unit1 = _open_or_stop(met_file_name, "r",
                          " Error opening met file, error=%3d")

    # Cabecalho global: 11 valores. ATENCAO -- esta mesma linha e lida DE
    # NOVO adiante (linha 786) como o registro do PRIMEIRO sitio, com um
    # alinhamento de campos diferente: posicao 3 -> latm(1), 4 -> longm(1),
    # 6 -> vitd_indexm(1). Ou seja o `pindex` daqui E o indice VITD do
    # primeiro sitio. Nao e defeito, e o formato do arquivo.
    tok = unit1.read_star(11)
    icnti = _to_int(tok[0], icnti)                       # numero de linhas
    state.ncols = _to_int(tok[1], state.ncols)           # numero de colunas
    state.lat = _to_real(tok[2], state.lat)              # latitude
    state.mlong = _to_real(tok[3], state.mlong)          # longitude
    state.elev = _to_real(tok[4], state.elev)            # elevacao (m)
    pindex = _to_int(tok[5], pindex)                     # ponteiro
    state.met_count = _to_int(tok[6], state.met_count)   # numero de sitios met
    state.timeoffset = _to_real(tok[7], state.timeoffset)  # defasagem do zulu
    state.timstep = _to_real(tok[8], state.timstep)      # passo de tempo do dado
    state.mflag = _to_real(tok[9], state.mflag)          # sentinela de dado faltante
    state.iheight = _to_real(tok[10], state.iheight)     # altura da instrumentacao

    state.lat = _anint(state.lat, 10)
    state.mlong = _anint(state.mlong, 10)
    state.elev = _anint(state.elev, 10)
    state.timstep = _anint(state.timstep, 10)
    state.mflag = _anint(state.mflag, 10)

    if state.iheight < EPS:                              # na pratica: <= 0
        state.iheight = 2e0
    state.iheight = _anint(state.iheight, 10)

    # confere se o numero de linhas de met esta correto (linhas 730-754)
    ntlines = 0
    header = unit1.read_a()
    header = unit1.read_a()

    ltest = 0
    io = 0
    a = " "          # character(len=512); o original nao inicializa (ver D2)

    while io != -1 and ltest == 0:
        a, io = unit1.read_a_io(a)
        i = 1
        itest = 0
        while i <= len(a.rstrip()) and itest == 0:
            if ord(_ch(a, i)) != 32:                     # nao e branco
                # 45..57 = '-' '.' '/' e os digitos '0'..'9'
                if 45 <= ord(_ch(a, i)) <= 57:
                    ntlines = ntlines + 1
                    itest = 1
                else:
                    break                                # `exit` do original
            i = i + 1
        if itest == 0:
            ltest = 1

    # Conta linhas de dados consecutivas ate achar uma linha cujo primeiro
    # caractere nao-branco nao seja `-./0-9` (ou uma linha em branco, ou o
    # fim do arquivo). Como o `read` com `iostat` deixa `a` INALTERADA no
    # EOF, a ultima linha de dados e contada duas vezes -- e e isso que o
    # `- 1` abaixo compensa. Se o arquivo terminar com uma linha em branco
    # ou de texto, nao ha duplo, `ntlines` sai igual a `icnti` e a correcao
    # nao dispara. As linhas 756-758 do original mostram o autor iterando
    # nessa conta.
    if icnti != ntlines and ntlines > 1:
        icnti = ntlines - 1
    unit1.rewind()

    state.dstart = 0
    if state.timstep == 3e0 and state.met_count > 1:
        state.dstart = min(14, icnti)
    else:
        state.dstart = 2

    # dimensiona os arrays (linhas 769-771)
    vitd_indexm = _new_array(state.met_count, 0)
    latm = _new_array(state.met_count)
    longm = _new_array(state.met_count)
    metm = np.zeros((state.met_count + 1, icnti + 1, MAXCOL + 1), dtype=np.float64)

    # zera o array de met -- com o sentinela, nao com 0. O laco triplo do
    # original cobre i=1..met_count, j=1..icnti, k=1..maxcol; os slots
    # fantasma ficam em 0.
    metm[1:, 1:, 1:] = state.mflag

    # le o array de met e a informacao de lat/long/indice (linhas 782-813)
    i = 1
    io = 0
    while io != -1 and i <= state.met_count:
        # 10 campos; na primeira volta este e o cabecalho global relido
        tok, io = unit1.read_star_io(10)
        ijunk = _to_int(tok[0], ijunk)
        ijunk = _to_int(tok[1], ijunk)
        latm[i] = _to_real(tok[2], latm[i])
        longm[i] = _to_real(tok[3], longm[i])
        junk = _to_real(tok[4], junk)
        vitd_indexm[i] = _to_int(tok[5], vitd_indexm[i])
        ijunk = _to_int(tok[6], ijunk)
        junk = _to_real(tok[7], junk)
        junk = _to_real(tok[8], junk)
        junk = _to_real(tok[9], junk)
        latm[i] = _anint(latm[i], 10)
        longm[i] = _anint(longm[i], 10)

        header, io = unit1.read_a_io(header)
        header, io = unit1.read_a_io(header)

        for j in range(1, icnti + 1):
            # ponytail: laco elemento a elemento em Python, como no
            # original. Se a leitura de arquivos multiponto grandes ficar
            # lenta, vetorizar com numpy -- mas o `_anint` tem que continuar
            # multiplicando por `10**-p`, nao dividindo por `10**p`: os dois
            # nao dao o mesmo double.
            if state.ncols > MET_COL:
                tok, io = unit1.read_star_io(MET_COL + 2)
                for k in range(1, MET_COL + 1):
                    metm[i, j, k] = _to_real(tok[k - 1], metm[i, j, k])
                metm[i, j, 34] = _to_real(tok[MET_COL], metm[i, j, 34])
                metm[i, j, 35] = _to_real(tok[MET_COL + 1], metm[i, j, 35])
                # Ponto de atencao 15: le as colunas 1..32, 34 e 35, mas
                # quantiza 1..34 -- quantiza a 33 (que ninguem leu, ainda
                # com o sentinela) e NAO quantiza a 35.
                for k in range(1, MET_COL + 3):
                    metm[i, j, k] = _anint(metm[i, j, k], 15)
            else:
                tok, io = unit1.read_star_io(MET_COL)
                for k in range(1, MET_COL + 1):
                    metm[i, j, k] = _to_real(tok[k - 1], metm[i, j, k])
                for k in range(1, MET_COL + 1):
                    metm[i, j, k] = _anint(metm[i, j, k], 15)

            # sem precipitacao mensuravel -> tipo de precipitacao = 1
            if metm[i, j, MetCol.PREC] <= 1e-4:
                metm[i, j, MetCol.PT] = 1
            if metm[i, j, MetCol.PREC2] <= 1e-4:
                metm[i, j, MetCol.PT2] = 1

        i = i + 1

    # Descarta a ultima linha se o ano pular mais de 1 em relacao a
    # penultima (linha parcial ou lixo no fim do arquivo). Usa `metm(1,...)`
    # fixo no sitio 1, independente de `met_count`.
    if icnti > 1:
        if metm[1, icnti, MetCol.YEAR] - metm[1, icnti - 1, MetCol.YEAR] > 1:
            icnti = icnti - 1

    state.iend = state.istart + icnti - 1
    iendi = state.iend
    unit1.close()                                            # fecha o arquivo met

    # ---------------- STEP 3 (linhas 825-1153) ----------------
    # Le a informacao antiga, se existir; abre arquivos e escreve cabecalhos.

    # arquivo de propriedades padrao de camada (829-839)
    io = 0
    state.unit30 = _open_or_stop(
        "FASST_defaultlayerprop.inp", "r",
        " Error opening default layerprop file, error = %3d")

    # arquivo de propriedades padrao de vegetacao (841-851)
    io = 0
    state.unit31 = _open_or_stop(
        "FASST_defaultvegprop.inp", "r",
        " Error opening default vegprop file, error=%3d")

    # arquivo de propriedades de vegetacao MODIS (853-863)
    io = 0
    state.unit32 = _open_or_stop(
        "MODIS_IGBP_vegprop.inp", "r",
        "Error opening MODIS vegprop file, error = %3d")

    # arquivo de propriedades de vegetacao UMD (865-874)
    io = 0
    state.unit33 = _open_or_stop(
        "UMD_LDAS_vegprop.inp", "r",
        "Error opening UMD vegprop file, error = %3d")

    # Arquivo de propriedades de solo do usuario; o padrao e
    # 'fasstusersoil.inp' (876-887).
    # ESTA e a linha 878 que conserta o ponto de atencao 13 -- mas conserta
    # so UM dos tres caminhos. Com `lcount == 0` e ponto unico, o nome que o
    # strip da linha 629 deixou vazio e restaurado aqui. Nos outros dois
    # casos nao ha conserto: com `lcount != 0` o nome veio do `read` da linha
    # 628, e no modo multiponto veio do `read` da linha 683 -- se nenhuma
    # dessas linhas tiver `!`, o nome chega vazio ate aqui e o `open` falha.
    io = 0
    if lcount == 0 and state.single_multi_flag == 0:
        fasstusersoil = "fasstusersoil.inp"
    state.unit90 = _open_or_stop(
        fasstusersoil, "r",
        " Error opening user soil property file, error = %3d")

    # Base para varios nomes de saida (889-896): posicao do ultimo '.' do
    # nome do arquivo met. `_ch` devolve ' ' para indice 0, que e o que
    # impede o `met_file_name(0:0)` fora de faixa do original -- ali a
    # condicao `npos > 0` so e avaliada DEPOIS da substring.
    npos = mlen
    while _ch(met_file_name, npos) != "." and npos > 0:
        npos = npos - 1
    if npos == 0:
        npos = mlen
    fnlen = 0
    # `fnlen` so era usado pelos `allocate(character(len=fnlen)...)` das
    # linhas 909-910, que estao comentados: atribuido e nunca lido.
    fnlen = npos + 3

    # nomeia os arquivos de saida no caso multiponto (898-916)
    if state.single_multi_flag == 1:
        io = 0
        state.unit5 = _open_or_stop(
            vitd_input, "r", " Error opening complex file, error = %3d")

        # `met_file_name(1:npos)` INCLUI o ponto, entao "yuma.met" com
        # npos = 5 da "yuma." + "fst" = "yuma.fst".
        output1 = met_file_name[:npos] + "fst"       # met e dados de superficie
        output2 = met_file_name[:npos] + "grd"       # informacao nodal
        output3 = met_file_name[:npos] + "flx"       # fluxos de superficie
        output4 = met_file_name[:npos] + "vgt"       # dados de vegetacao
        output5 = met_file_name[:npos] + "snw"       # processos de neve

    io = 0
    state.unit2 = _open_or_stop(                     # met e info de superficie
        output1, "w", " Error opening main output, error = %3d")

    if nprint == 1:
        io = 0
        state.unit3 = _open_or_stop(                 # saida nodal
            output2, "w", " Error opening node output, error = %3d")

    if flprint == 1:
        io = 0
        state.unit24 = _open_or_stop(                # saida de fluxos
            output3, "w", " Error opening flux output, error=%3d")

        state.unit24.write("Surface Flux File\n")
        state.unit24.write("Conventions used: \n")
        state.unit24.write("     + => flux towards surface\n")
        state.unit24.write("     sheat > 0 -> Tair > Tgr => surface warms, "
                           "air cools\n")
        state.unit24.write("     sheat < 0 -> Tair < Tgr => surface cools, "
                           "air warms\n")
        state.unit24.write("     lheat > 0 -> qair > qgr => condensation; "
                           "surface warms, air cools\n")
        state.unit24.write("     lheat < 0 -> qair < qgr => evaporation; "
                           "surface cools, air warms\n")
        state.unit24.write(" \n")

    if vegint > EPS:
        io = 0
        state.unit28 = _open_or_stop(                # saida de vegetacao
            output4, "w", " Error opening veg output, error = %3d")

    if sprint == 1:
        io = 0
        state.unit55 = _open_or_stop(                # saida de neve
            output5, "w", " Error opening snow output, error = %3d")

    # Abre o arquivo de flags de inferencia, apagando o anterior se existir
    # (984-985). Sem `iostat=` e sem ramo de erro no original.
    state.unit10 = _Unit("inferred.out", "w")

    # arquivos usados nas execucoes de "hot start" (987-1011)
    ro1 = met_file_name[:npos - 1] + "_ro_soil.dat"
    state.unit17 = _Unit(ro1, "wb")                  # unformatted, access=stream
    ro2 = met_file_name[:npos - 1] + "_ro_time.dat"
    state.unit18 = _Unit(ro2, "wb")

    if state.infer_test == 1:
        state.unit11 = _Unit(ro1o, "rb")
        # Ponto de atencao 11: assume que `ro1o` termina em `_soil.dat`
        # (9 caracteres), sem validar. Se o `read` da linha 256 nao tiver
        # `!`, o strip deixou `ro1o` vazio e isto da so "_time.dat" -- outra
        # consequencia do ponto 13.
        ro2o = ro1o[:len(ro1o.rstrip()) - 9] + "_time.dat"
        state.unit12 = _Unit(ro2o, "rb")

    # Determina o numero de linhas faltantes no novo arquivo met e
    # dimensiona os arrays de acordo (1013-1063).
    if state.single_multi_flag == 0:
        i = 1
        tflag = 0
        while i <= state.met_count and tflag == 0:
            if state.vitd_index == vitd_indexm[i]:
                fdi = state.freq_id
                tflag = 1
                kk = i
            i = i + 1
        if tflag == 0:
            print(" Something wrong with input meteorological file, stopping.")
            sys.exit(1)
    elif state.single_multi_flag == 1:
        io = 0
        i = 0
        vitd_att = " "
        while io != -1:
            vitd_att, io = state.unit5.read_a_io(vitd_att)
            i = i + 1
        i = i - 2
        num_frq = i
        state.unit5.rewind()

        io = 0
        vitd_att = " "
        vitd_att, io = state.unit5.read_a_io(vitd_att)    # linha de cabecalho

        j = 1
        while j <= num_frq:
            (wtype, io, lcount, nm0, wvel, wdepth, mtest,
             sname) = read_complex(state, wtype, io, lcount, nm0, wvel,
                                   wdepth, mtest, sname)
            i = 1
            tflag = 0
            while i <= state.met_count and tflag == 0:
                if state.vitd_index == vitd_indexm[i]:
                    if j == 1:
                        fdi = state.freq_id
                    tflag = 1
                    kk = i
                i = i + 1
            j = j + 1
        state.unit5.rewind()
    # Ponto de atencao 5, primeira metade: este `i = i - 1` esta FORA dos
    # dois ramos. No ponto unico `i` saiu do `do while` valendo kk+1 (a
    # execucao so chega aqui com tflag == 1, senao teria parado), entao
    # `i` passa a ser exatamente `kk`. No multiponto `i` e o contador do
    # laco interno da ULTIMA chamada de read_complex, sem relacao nenhuma
    # com sitio de met -- e e esse `i` que o STEP 3 usa depois para indexar
    # `metm` (linhas 1095-1099, Dia 8).
    i = i - 1

    if state.infer_test == 1:
        d1i = 3
        # No ponto unico isto e uma ida e volta sem efeito: `fdi` foi
        # atribuido DE `freq_id` na linha 1019. No multiponto tem sentido --
        # restaura o `freq_id` do primeiro registro, depois de o laco de
        # read_complex te-lo avancado.
        state.freq_id = fdi
        # `max(moverlap, istart)`: MOVERLAP e a constante 15 e `istart` vale
        # 2 quando infer_test == 1, entao o resultado e sempre 15.
        wstart = max(MOVERLAP, state.istart)
        (sdensi, phie, timeo) = read_old_data(d1i, state.dstart, wstart,
                                              sdensi, phie, timeo)

        for k in range(1, state.moverlapr + 1):
            t1o[k] = 0e0
            t1o[k] = met_date(timeo[k][1], timeo[k][2], timeo[k][3],
                              timeo[k][4])
        t1n = 0e0
        t1n = met_date(metm[kk, 1, 1], metm[kk, 1, 2], metm[kk, 1, 3],
                       metm[kk, 1, 4])

        state.mpos = 0
        for k in range(1, state.moverlapr):
            if t1o[k] < t1n and t1o[k + 1] >= t1n:
                state.mpos = k
        if state.mpos == 0:
            state.mpos = state.moverlapr

    # ponytail: primeira das capturas do gancho `_locals` -- aqui,
    # no fim do bloco de hot start (linha 1084). Ela existe porque o laco
    # logo abaixo LEVANTA de proposito no modo multiponto com
    # `num_frq == 0` (ponto de atencao 5 + divergencia D4): sem uma captura
    # ANTES dele, os fixtures multiponto dos Dias 5-7 voltariam com `loc`
    # vazio. As capturas seguintes sobrescrevem esta quando sao alcancadas.
    if _locals is not None:
        _locals.update(locals())
        _locals["_hot"] = dict(locals())

    # ---- lacunas de tempo no arquivo met (1086-1117) ----
    # Conta quantos passos de tempo FALTAM no met, somando as lacunas entre
    # registros consecutivos. O total (`mtsteps`) entra no dimensionamento
    # dos arrays de execucao: `missing_met` depois preenche essas posicoes.
    #
    # Ponto de atencao 5, segunda metade: o laco indexa `metm` com o `i`
    # que SOBROU do bloco anterior. No ponto unico isso e o sitio casado
    # (`i == kk`, pelo `i = i - 1` da 1063). No multiponto `i` e o contador
    # do laco interno da ULTIMA chamada de `read_complex`, sem relacao
    # nenhuma com sitio de met -- e com `num_frq == 0` o laco nem roda e
    # `i` chega aqui valendo -1. Por isso os indices passam por `_ix`
    # (divergencia D4): em NumPy, -1 pegaria o ULTIMO sitio em silencio.
    mtsteps = 0
    msteps1 = 0
    t1oo = 0e0
    t1n = 0e0

    # Ponto de atencao 9: com `timstep <= 0` o `if` nao dispara, `rtstps`
    # fica no zero com que foi declarado (D2), `tstps` sai 0 e TODO
    # `msteps1` sai 0 -- nenhuma lacuna e detectada e o arquivo passa como
    # se estivesse completo.
    if state.timstep > 0e0:
        rtstps = 24e0 / state.timstep
    tstps = int(rtstps)

    for k in range(2, icnti + 1):
        t1oo = met_date(metm[_ix(i), k - 1, 1], metm[_ix(i), k - 1, 2],
                        metm[_ix(i), k - 1, 3], metm[_ix(i), k - 1, 4])
        t1n = met_date(metm[_ix(i), k, 1], metm[_ix(i), k, 2],
                       metm[_ix(i), k, 3], metm[_ix(i), k, 4])

        # ano % 4 == 0, escrito com `aint`: o fator 367/366 abaixo sai
        # daqui. Note que o teste e sobre o ano do registro NOVO.
        mody = metm[_ix(i), k, 1] - _aint(metm[_ix(i), k, 1] * 2.5e-1) * 4e0
        # Ponto de atencao 17: `modyo` e calculado e NUNCA lido. A unica
        # linha que o usaria -- `if(mody /= modyo) daylim = 1.07d0`, 1105 --
        # esta comentada no proprio Fortran. Transcrito para preservar o
        # custo e o fato.
        modyo = (metm[_ix(i), k - 1, 1]
                 - _aint(metm[_ix(i), k - 1, 1] * 2.5e-1) * 4e0)

        daylim = 1.042e0                             # `!1.003d0` no original
        if metm[_ix(i), k - 1, 1] != metm[_ix(i), k, 1]:
            daylim = 25.07e0                         # virada de ano

        # PONTO DE ATENCAO 21: este `msteps1 = 0` (1108) e MORTO, e o
        # da 1088 tambem. Os dois ramos logo abaixo atribuem `msteps1`
        # incondicionalmente, e depois do laco ninguem mais o le -- so
        # `mtsteps` sobrevive. Preservado; e o unico mutante do Dia 8 que
        # escapa por ser provadamente equivalente.
        #
        # ATENCAO: a guarda compara o PRODUTO REAL com `daylim`, nao o
        # `msteps1` truncado. Sao duas expressoes distintas e e por isso
        # que a linha aparece duplicada: `int(x)` e `x < daylim`. O fator
        # tem de ser o MESMO nas duas (ha fixture para isso).
        msteps1 = 0
        if abs(mody) <= EPS:
            msteps1 = int((t1n - t1oo) * 367e0 * tstps)
            if (t1n - t1oo) * 367e0 * tstps < daylim:
                msteps1 = 0
        else:
            msteps1 = int((t1n - t1oo) * 366e0 * tstps)
            if (t1n - t1oo) * 366e0 * tstps < daylim:
                msteps1 = 0
        # PONTO DE ATENCAO 20: `mtsteps` nao e uma contagem de registros
        # faltantes, e uma aproximacao cujo arredondamento depende de QUAIS
        # horas a lacuna cai entre. Com dado horario, uma lacuna de um
        # registro da produto 1.999999974 (-> 1) em 21 dos 22 pares de
        # horas possiveis, mas 2.000000852 (-> 2) entre 20h e 22h: o erro
        # do ponto 16 muda de sinal. Como `maxlines` e so um limite de
        # alocacao, o efeito e alocar uma linha a mais ou a menos.
        mtsteps = mtsteps + msteps1

    # ---- maxlines e alocacao (1119-1148) ----
    # `maxlines` e campo de `fasst_global`, nao local: e de `state.maxlines`
    # que `zero_mxl_params` tira o limite do zeramento.
    if state.infer_test == 0:
        state.maxlines = state.iend + mtsteps
    else:
        # A variante comentada da 1123 -- `iend + (dstart - oldpos) +
        # mtsteps` -- nao e transcrita.
        state.maxlines = state.iend + mtsteps + state.moverlapr + 1

    allocate_run_arrays(state)

    # ponytail: segunda captura, no fim do STEP 3. A do fim de `main` nao
    # substitui esta: o STEP 4 REUSA `i`, `k`, `kk`, `j` e `tflag`, entao
    # quem quiser verificar o que os STEPs 1-3 produziram precisa olhar
    # aqui, antes de os lacos principais sobrescreverem.
    if _locals is not None:
        _locals.update(locals())
        _locals["_step3"] = dict(locals())

    # ------------------------------------------------------------------
    # A transcricao para aqui: linha 1153 do original. Falta:
    #   STEP 4   (1156-1425) os dois lacos principais
    #   STEP 6   (1433-1480) fecha arquivos, desaloca, tempo de CPU
    # ------------------------------------------------------------------
    unit4.close()
    raise NotImplementedError(
        "transcricao de fasst_driver.F90 em andamento: parou na linha 1153 "
        "do original (fim do STEP 3)"
    )

    et = time.process_time()                                   # call cpu_time(et)
    print()
    print(" total run time: ", et - it)


# ---------------------------------------------------------------------------
# subroutine sort2 (fasst_driver.F90, linhas 1483-1534)
# ---------------------------------------------------------------------------
def sort2(n, arr, arr1):
    """Ordena os nos em ordem crescente, carregando um array inteiro junto.

    Assinatura original:
        subroutine sort2(n,arr,arr1,ncount,b,b1)
          intent(in)    :: n
          intent(inout) :: arr(maxn+3), arr1(maxn+3)
          intent(out)   :: ncount, b(maxn), b1(maxn)

    E a versao de `initprofile.sort` (initprofile.py:1116) com um array
    inteiro paralelo: insertion sort IN-PLACE sobre `arr[1..n]` (e por isso
    que `arr`/`arr1` sao `intent(inout)` -- o bloco de fusao de camadas usa
    o array JA ORDENADO depois da chamada), seguido da deduplicacao de
    valores adjacentes cuja diferenca, escalada por 1000, fique abaixo de
    1e-4. A comparacao e sempre contra o elemento anterior do array ja
    ordenado (`arr(j-1)`), nunca contra o ultimo valor mantido em `b`.

    Devolve `(ncount, b, b1)` e modifica `arr`/`arr1` no lugar.

    Dois detalhes do original preservados:
    * `if(i == 1 .and. arr(i) > a) i = 0` nunca dispara -- se i==1 e
      arr(i)>a, o `while` anterior ainda estaria rodando. Igual em
      `initprofile.sort`.
    * `b`/`b1` sao zerados de 1 a maxn ANTES de tudo. Isso torna DEFINIDA
      (= 0) a leitura de posicoes nunca preenchidas que o bloco de fusao
      faz depois -- ver ponto de atencao 1d.

    Oddity do original: as saidas sao declaradas `b(maxn)`/`b1(maxn)`
    enquanto as entradas sao `arr(maxn+3)`/`arr1(maxn+3)`. Com n > maxn as
    escritas em `b` estourariam; nao acontece a partir do driver, onde
    n <= nlayers1 + maxl = 13.
    """
    b = _new_array(MAXN)
    b1 = _new_array(MAXN, 0)

    for j in range(2, n + 1):
        a = arr[j]
        nm = arr1[j]
        i = j - 1
        while i >= 1 and arr[i] > a:
            arr[i + 1] = arr[i]
            arr1[i + 1] = arr1[i]
            i = i - 1
        if i == 1 and arr[i] > a:        # nunca dispara; preservado por fidelidade
            i = 0
        arr[i + 1] = a
        arr1[i + 1] = nm

    ncount = 1
    b[ncount] = arr[1]
    b1[ncount] = arr1[1]
    for j in range(2, n + 1):
        a = arr[j] * 1000.0
        if abs(a - arr[j - 1] * 1e3) > 1e-4:
            ncount = ncount + 1
            b[ncount] = arr[j]
            b1[ncount] = arr1[j]

    return ncount, b, b1


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
    # Uma linha de dados met com as 32 colunas que o original le
    # (MetCol.YEAR .. MetCol.TSOIL). `hr` e parametrizada so para dar uma
    # marca distinguivel por linha.
    def _linha_met(hr, ano=2010, prec=0.0, extra="", doy=100, minuto=0):
        return (
            "%d %d %d %s 1013.0 25.0 50.0 3.0 180.0 %s 0 0.0 0 "
            "0.2 1.0 1 0.0 0.0 0 0.0 0.0 0 "
            "800.0 600.0 200.0 150.0 350.0 400.0 30.0 180.0 0.0 20.0%s\n"
            % (ano, doy, hr, minuto, prec, extra)
        )

    def _met(nome="yuma.met", cabecalho=None, corpo=None, cauda=""):
        """Escreve um arquivo met e devolve o caminho.

        Formato: linha 1 = cabecalho global de 11 campos (que o STEP 2 le
        DUAS vezes -- ver o comentario na transcricao), linhas 2-3 = texto,
        depois as linhas de dados.
        """
        if cabecalho is None:
            # icnti ncols lat mlong elev pindex met_count offset timstep mflag iheight
            cabecalho = "3 32 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 2.0"
        if corpo is None:
            corpo = "".join(_linha_met(h) for h in (12, 13, 14))
        caminho = os.path.join(d, nome)
        with open(caminho, "w") as f:
            f.write(cabecalho + "\n")
            f.write("Year DOY Hr Min AP Tmp RH WS WDir Prec ...\n")
            f.write("                 mbar  C   %  m/s deg mm/hr\n")
            f.write(corpo)
            f.write(cauda)
        return caminho

    _MET = _met()

    # O STEP 3 abre quatro arquivos de propriedades por nome FIXO, relativos
    # ao diretorio corrente: FASST_defaultlayerprop.inp,
    # FASST_defaultvegprop.inp, MODIS_IGBP_vegprop.inp e
    # UMD_LDAS_vegprop.inp. Aqui bastam stubs vazios -- quem le o conteudo
    # (`veg_propl`/`veg_proph`) nao e chamado pelo driver. Tambem e por isso
    # que `_run` troca de diretorio: assim `inferred.out` e os `output1..5`
    # de nome relativo caem no temporario, nao no repo.
    _PROPS = ("FASST_defaultlayerprop.inp", "FASST_defaultvegprop.inp",
              "MODIS_IGBP_vegprop.inp", "UMD_LDAS_vegprop.inp")
    for _n in _PROPS:
        with open(os.path.join(d, _n), "w") as f:
            f.write("")
    with open(os.path.join(d, "fasstusersoil.inp"), "w") as f:
        f.write("")
    # Os dois arquivos binarios da execucao ANTERIOR, para o caminho de hot
    # start (infer_test == 1). No original o `status='unknown'` os criaria
    # vazios e `read_old_data` leria de um arquivo vazio; aqui as unidades
    # 11 e 12 abrem para leitura e a falta e erro -- ver divergencia D1.
    for _n in ("velho_ro_soil.dat", "velho_ro_time.dat"):
        with open(os.path.join(d, _n), "wb") as f:
            f.write(b"")
    # Nomes que alguns fixtures informam no proprio `.inp` (unidades 90 e 5).
    # `complexo.vitd` tem UMA linha so de proposito: o laco de contagem da
    # linha 1034 faz `num_frq = (linhas + 1) - 2`, entao com 1 linha
    # `num_frq` sai 0 e o laco de `read_complex` nao roda -- e o unico jeito
    # de atravessar o STEP 3 multiponto sem bater no stub.
    for _n in ("meuarquivo.inp", "meusolo.inp"):
        with open(os.path.join(d, _n), "w") as f:
            f.write("")
    with open(os.path.join(d, "complexo.vitd"), "w") as f:
        f.write("cabecalho do arquivo complexo\n")

    def _exec(linhas, met, st, loc, tolera=()):
        """Escreve o `.inp`, troca para o diretorio temporario e roda `main`.

        `tolera` lista tipos de excecao a engolir ALEM de
        `NotImplementedError`. Serve para um unico caso: o ponto de atencao
        5 faz o laco de lacunas do STEP 3 indexar `metm` com `i == -1` no
        modo multiponto com `num_frq == 0`, e `_ix` (divergencia D4)
        levanta `IndexError` ali de proposito. Nao esta no `except` padrao
        porque esse `IndexError` E o sinal que D4 existe para dar -- quem
        quiser atravessa-lo tem de pedir explicitamente.
        """
        import contextlib
        import io as _io
        caminho = os.path.join(d, "caso.inp")
        with open(caminho, "w") as f:
            f.write(linhas.replace("@MET@", met if met else _MET))
        saida = _io.StringIO()
        antes = os.getcwd()
        try:
            os.chdir(d)
            with contextlib.redirect_stdout(saida):
                try:
                    main(["fasst_driver", caminho], _state=st, _locals=loc)
                except NotImplementedError:
                    pass                   # parou no marcador da transcricao
                except tolera:
                    pass
        finally:
            os.chdir(antes)
        return saida.getvalue()

    def _run(linhas, met=None, tolera=()):
        """Roda `main` ate onde a transcricao chegou e devolve (estado, saida)."""
        st = FasstState()
        return st, _exec(linhas, met, st, None, tolera)

    def _run_loc(linhas, met=None, tolera=()):
        """Como `_run`, mas tambem devolve as variaveis LOCAIS de `main`.

        Metade do que os STEPs 1 e 2 produzem (`icnti`, `metm`, `latm`,
        `output1..5`, `fasstusersoil`, `sname`, `surfdepth`, `npos`) sao
        locais do `program`, nao campos de `FasstState` -- sem isto ficariam
        sem verificacao nenhuma.

        Ha DUAS capturas em `main`: uma no fim do STEP 3 e outra no
        marcador de fim da transcricao. A segunda sobrescreve a primeira,
        entao o que volta aqui e o estado mais avancado que a execucao
        alcancou. Se um stub levantar antes das duas (e `read_complex`
        levanta, no modo multiponto), `loc` volta VAZIO -- nesse caso usar
        `_run` e checar so o estado.
        """
        st, loc = FasstState(), {}
        _exec(linhas, met, st, loc, tolera)
        return st, loc

    _CABECALHO = (
        "@MET@                        ! met\n"
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
        "1                            ! nlayers\n"
        "7 1.0                        ! soiltype, lthick\n"
        "n                            ! ttest\n"
        "n                            ! mtest\n"
    ))
    assert st.vitd_index == 1234
    assert st.slope_fasst == 30.0 and st.aspect == 180.0
    assert abs(st.sloper - 0.5235987756666667) < 1e-15        # 30 graus com PI truncado
    assert st.rough == 0.01
    assert st.water_flag == 0
    assert st.hsaccum == 0.3 and st.iswe == 0.05 and st.hi == 0.02
    # o bloco de agua (linhas 344-365) nao rodou: as camadas vieram do
    # arquivo, nao do [26, 26, 7] que ele montaria
    assert st.nlayers == 1 and st.soiltype[1] == 7

    # o limite de 1.57 rad do `dmin1`: 90 graus dao 1.570796... com o PI
    # truncado do fasst_global, acima do limite, entao o `dmin1` trava em
    # 1.57. 90 e o maior declive aceito (91 ja para o programa). O valor
    # final NAO e 1.57 exato: o `anint(x*1d20)*1d-20` da linha 320 perde o
    # ultimo bit no caminho de volta. Preservado -- e a mesma operacao em
    # IEEE double que o Fortran faz.
    st, _ = _run(_CABECALHO + "90.0 0.0\n0.01\n0\n0.0 0.0\n0.0\n0\n0\n"
                              "1\n7 1.0\nn\nn\n")
    assert st.slope_fasst == 90.0
    assert st.sloper == 1.5699999999999998 != 1.57

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
        "@MET@    ! met\n"
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
            f.write(linhas.replace("@MET@", _MET))
        saida = _io.StringIO()
        antes = os.getcwd()
        try:
            os.chdir(d)                    # como `_exec`: nao sujar o repo
            with contextlib.redirect_stdout(saida):
                main(["fasst_driver", caminho], _state=FasstState())
        except SystemExit as e:
            assert e.code == 1
            assert trecho in saida.getvalue(), saida.getvalue()
            return
        finally:
            os.chdir(antes)
        raise AssertionError("deveria ter parado: " + trecho)

    _para("@MET@ !\n2 0 0.5 1 !\n", "infer test flag must be 0 or 1")
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

    # bloco de solo minimo: uma camada e nenhum perfil inicial. Necessario
    # porque a leitura continua para as linhas 536-668 depois da vegetacao.
    _SOLO = (
        "1      ! nlayers\n"
        "7 1.0  ! soiltype, lthick\n"
        "n      ! ttest\n"
        "n      ! mtest\n"
    )

    def _veg(baixa, alta, folhagem=_FOLHAGEM, dossel=_DOSSEL, solo=_SOLO):
        """Monta um caso com `vegl_type` = baixa e `vegh_type` = alta."""
        texto = _CABECALHO + _TERRA + str(baixa) + " ! vegl_type\n"
        if baixa != 0:
            texto += folhagem
        texto += str(alta) + " ! vegh_type\n"
        if alta != 0:
            texto += dossel
        return _run(texto + solo)[0]

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

    # ---- STEP 1 camadas de solo e perfil inicial (linhas 536-699) ----
    def _solo(corpo, baixa=0, alta=0):
        """Caso de terra com o bloco de solo `corpo`."""
        texto = _CABECALHO + _TERRA + str(baixa) + " ! vegl\n"
        if baixa:
            texto += _FOLHAGEM
        texto += str(alta) + " ! vegh\n"
        if alta:
            texto += _DOSSEL
        return _run(texto + corpo)[0]

    # tipos USCS e USDA (positivo pela tabela, negativo so deslocado)
    st = _solo("3 !\n7 0.5 !\n103 1.0 !\n-105 2.0 !\nn !\nn !\n")
    assert st.nlayers == 3
    assert list(st.soiltype[1:4]) == [7, map_usda_soiltype_to_uscs(3), -5]
    assert [str(s) for s in st.sclass[1:4]] == ["USCS", "USDA", "USDA"]
    assert list(st.lthick[1:4]) == [0.5, 1.0, 2.0]

    # soiltype == -1: le uma linha extra com o nome do solo (truncado em 4,
    # e SEM strip de `!`), e depois lcount != 0 faz ler o fasstusersoil.
    # Nenhum dos dois e campo de estado -- a prova de que as duas linhas
    # foram consumidas e que o perfil de temperatura depois delas alinha.
    st = _solo("1 !\n-1 1.0 !\nMEUSOLO\nmeuarquivo.inp !\ny !\n1 !\n0.5 20.0 !\nn !\n")
    assert st.zti[1] == 0.5 and abs(st.tm[1] - (20.0 + TREF)) < 1e-12

    # uma camada de tipo 26 liga water_flag/wflag2
    st = _solo("2 !\n26 1.0 !\n7 2.0 !\nn !\nn !\n")
    assert st.water_flag == 1 and st.nlayers == 2   # fusao nao roda: nlayers > 1

    # ramo wflag2 (565-572): uma unica camada de agua vira tres
    st = _solo("1 !\n26 5.0 !\nn !\nn !\n")
    assert st.nlayers == 3
    assert list(st.soiltype[1:4]) == [26, 26, 7]
    assert st.lthick[1] == 2e-2 and abs(st.lthick[2] - 4.98) < 1e-12
    assert st.lthick[3] == 1.0

    # ===== PONTO DE ATENCAO 1: o bloco de fusao de camadas =====
    # vegl 1011 -> MODIS 11 -> 13 (swamp/bog) monta nlayers1 = 3,
    # soiltype1 = [26, 26, 7], lthick1 = [0.02, 1.48, 1.0] (wdepth = 1.5).
    # O solo do arquivo tem 2 camadas: [26, 1.0] e [7, 2.0].
    # O resultado abaixo foi tracado a mao, indice por indice, e e LIXO de
    # proposito -- e o que o Fortran produz:
    #   * com nlayers(2) != nlayers1(3), o teste `soiltype(i-nlayers)` olha
    #     a camada errada e a camada de agua do arquivo e tratada como se
    #     nao fosse agua (1a);
    #   * sobra uma camada de espessura 0.0 e tipo 0 no meio, vinda de
    #     posicoes de lthick3/soiltype3 que sort2 zerou e o laco final le
    #     mesmo assim (1b + 1d).
    st = _solo("2 !\n26 1.0 !\n7 2.0 !\nn !\nn !\n", baixa=1011)
    assert st.nlayers == 5
    assert list(st.soiltype[1:6]) == [26, 26, 0, 7, 7]
    for i, v in enumerate([0.02, 1.46, 0.0, 1.0, 1.98], start=1):
        assert abs(st.lthick[i] - v) < 1e-12, (i, st.lthick[i], v)

    # Segundo caso de fusao, tambem tracado a mao: aqui a camada de agua do
    # arquivo e a SEGUNDA, o que faz o ramo `elif` do primeiro laco ser
    # realmente tomado -- no caso acima ele nunca e. E a diferenca entre
    # testar `i - nlayers` e copiar `i - nlayers1` fica visivel nos dois
    # lacos ao mesmo tempo: em i=4 o bloco detecta agua na camada 2 e copia
    # a camada 1 (tipo 7, que nao e agua) para a lista de agua; em i=5
    # detecta "nao-agua" na camada 3 (que nem existe) e copia a camada 2
    # (tipo 26, agua) para a lista de nao-agua. Resultado: as tres camadas
    # do swamp desaparecem e sobram uma camada de espessura zero e tipo 0 e
    # uma de agua de 2 m.
    st = _solo("2 !\n7 1.0 !\n26 2.0 !\nn !\nn !\n", baixa=1011)
    assert st.nlayers == 2
    assert list(st.soiltype[1:3]) == [0, 26]
    for i, v in enumerate([0.0, 2.0], start=1):
        assert abs(st.lthick[i] - v) < 1e-12, (i, st.lthick[i], v)

    # Terceiro caso, tambem tracado a mao: nlayers = 4 com nlayers1 = 3 faz
    # `i - nlayers` valer 0 em i=4, lendo o SLOT FANTASMA `soiltype[0]`. No
    # Fortran isso le memoria imediatamente antes do array; aqui le um zero
    # definido, que nao e 26, entao o ramo nao dispara -- o mesmo efeito
    # pratico, por acidente da convencao de slot fantasma do pacote.
    st = _solo("4 !\n26 1.0 !\n7 2.0 !\n7 3.0 !\n7 4.0 !\nn !\nn !\n", baixa=1011)
    assert st.nlayers == 4
    assert list(st.soiltype[1:5]) == [0, 26, 7, 7]
    for i, v in enumerate([0.0, 0.02, 3.0, 4.0], start=1):
        assert abs(st.lthick[i] - v) < 1e-12, (i, st.lthick[i], v)

    # NOTA de verificacao (ponto 1e): as ATRIBUICOES do `elif` do primeiro
    # laco nao tem teste aqui, e nao e por falta de cenario -- os valores
    # que elas gravam sao descartados antes de chegar a `state`. Conferido
    # por mutacao (trocar `i - nlayers1` por `i - nlayers` nas duas linhas
    # nao muda nenhuma saida) e pelos tres tracos manuais acima. So o
    # TESTE do `elif` e observavel, pelo `j` que ele move.

    # indice negativo no mesmo bloco: nlayers > nlayers1 + 1 leva
    # `i - nlayers` abaixo de zero. No Fortran e leitura antes do array;
    # aqui `_ix` levanta IndexError em vez de dar wraparound silencioso.
    try:
        _solo("5 !\n26 1.0 !\n7 1.0 !\n7 1.0 !\n7 1.0 !\n7 1.0 !\nn !\nn !\n",
              baixa=1011)
    except IndexError as e:
        assert "ponto de atencao 1" in str(e)
    else:
        raise AssertionError("esperado IndexError com indice negativo na fusao")

    # ---- perfis iniciais de temperatura e umidade (641-668) ----
    st = _solo("1 !\n7 1.0 !\n"
               "y !\n2 !\n0.1 10.0 !\n0.5 -5.0 !\n"
               "y !\n2 !\n0.1 0.3 !\n0.5 0.4 !\n")
    assert list(st.zti[1:3]) == [0.1, 0.5]
    assert abs(st.tm[1] - (10.0 + TREF)) < 1e-12        # convertido para Kelvin
    assert abs(st.tm[2] - (-5.0 + TREF)) < 1e-12
    assert list(st.zm[1:3]) == [0.1, 0.5]
    assert list(st.sm[1:3]) == [0.3, 0.4]               # umidade NAO e convertida
    # "Y" maiusculo tambem vale; qualquer outra coisa pula o bloco
    st = _solo("1 !\n7 1.0 !\nY !\n1 !\n0.1 10.0 !\nN !\n")
    assert abs(st.tm[1] - (10.0 + TREF)) < 1e-12 and st.sm[1] == 0.0

    _CAB_TERRA = _CABECALHO + _TERRA + "0 ! vegl\n0 ! vegh\n"
    _para(_CAB_TERRA + "1 !\n7 1.0 !\ny !\n1 !\n0.1 150.0 !\n",
          "Measured soil temperatures out of range")
    _para(_CAB_TERRA + "1 !\n7 1.0 !\nn !\ny !\n1 !\n0.1 1.5 !\n",
          "Measured soil moistures out of range")

    # ---- ramo multiponto e flag invalido (672-696) ----
    # `vitd_input`/`fasstusersoil` sao locais de main(), nao campos de
    # estado: aqui so da para provar que os dois registros foram
    # consumidos e que o STEP 1 terminou. Vira verificavel no Dia 7,
    # quando essas unidades forem abertas.
    st, _ = _run("@MET@ !\n0 1 0.5 1 !\n1.5 0.0 1 1 0 100 0 !\n"
                 "complexo.vitd !\nmeusolo.inp !\n", tolera=(IndexError,))
    assert st.single_multi_flag == 1 and st.vitd_index == 0

    # `single_multi_flag` invalido: o `else` da linha 693 so imprime, sem
    # `stop` -- mas o programa morre alguns passos depois de qualquer forma,
    # no STEP 3, porque nenhum dos dois ramos definiu `output1` e o `open`
    # recebe o nome em branco que o STEP 0 deixou. Ou seja, o "nao para" e
    # so aparente.
    import contextlib
    import io as _io
    caminho_ruim = os.path.join(d, "flag.inp")
    with open(caminho_ruim, "w") as f:
        f.write("@MET@ !\n0 2 0.5 1 !\n1.5 0.0 1 1 0 100 0 !\n"
                .replace("@MET@", _MET))
    saida = _io.StringIO()
    antes = os.getcwd()
    try:
        os.chdir(d)
        with contextlib.redirect_stdout(saida):
            main(["fasst_driver", caminho_ruim], _state=FasstState())
    except SystemExit as e:
        assert e.code == 1
        texto = saida.getvalue()
        assert "Model does not recognize single_multi_flag input of =   2" in texto
        assert "Error opening main output" in texto
    else:
        raise AssertionError("esperado SystemExit com single_multi_flag = 2")
    finally:
        os.chdir(antes)

    # ---- sort2 (linhas 1483-1534) ----
    def _arr(vals, fill=0.0):
        a = _new_array(MAXN + 3, fill)
        for i, v in enumerate(vals, start=1):
            a[i] = v
        return a

    arr, arr1 = _arr([3.0, 1.0, 3.0, 2.0]), _arr([30, 10, 33, 20], 0)
    n, b, b1 = sort2(4, arr, arr1)
    assert n == 3                                        # a duplicata 3.0 sai
    assert b[1:4] == [1.0, 2.0, 3.0]
    assert b1[1:4] == [10, 20, 30]
    assert arr[1:5] == [1.0, 2.0, 3.0, 3.0]              # ordenado IN-PLACE
    assert arr1[1:5] == [10, 20, 30, 33]                 # inteiros acompanham
    # limiar de deduplicacao: |a - anterior| * 1000 > 1e-4, ou seja 1e-7 na
    # escala original. Os tres casos cercam o limiar de perto o suficiente
    # para que mexer no 1e-4 mude o resultado.
    assert sort2(2, _arr([1.0, 1.00000005]), _arr([1, 2], 0))[0] == 1   # 5e-5
    assert sort2(2, _arr([1.0, 1.000001]), _arr([1, 2], 0))[0] == 2     # 1e-3
    assert sort2(2, _arr([1.0, 1.001]), _arr([1, 2], 0))[0] == 2        # 1e0
    # n <= 1: nenhum laco roda, mas ncount = 1 com arr(1) mesmo assim
    n, b, b1 = sort2(1, _arr([7.0]), _arr([5], 0))
    assert n == 1 and b[1] == 7.0 and b1[1] == 5
    # `b`/`b1` zerados de 1 a maxn antes de tudo (base do ponto 1d)
    assert b[2:] == [0.0] * (MAXN - 1) and b1[2:] == [0] * (MAXN - 1)

    # ---- STEP 2: leitura do arquivo meteorologico (linhas 702-822) ----
    _LAND = _CABECALHO + _TERRA + "0 ! vegl\n0 ! vegh\n" + _SOLO

    # cabecalho global (717-728)
    st, loc = _run_loc(_LAND)
    assert st.ncols == 32
    assert st.lat == 32.5 and st.mlong == -114.0 and st.elev == 50.0
    assert st.met_count == 1
    assert st.timeoffset == -7.0 and st.timstep == 1.0
    assert st.mflag == -999.0 and st.iheight == 2.0
    assert loc["pindex"] == 1234
    assert st.dstart == 2                            # timstep != 3
    assert loc["icnti"] == 3
    assert st.iend == st.istart + 3 - 1 == 3
    assert loc["iendi"] == st.iend

    # A ultima linha de dados e contada DUAS vezes: com `iostat=`, o `read`
    # no fim do arquivo deixa `a` inalterada e o laco escaneia de novo o
    # conteudo anterior. Por isso ntlines = icnti + 1, e por isso a linha
    # 759 faz `ntlines - 1`.
    assert loc["ntlines"] == 4

    # Com uma linha em branco no fim nao ha duplo: o laco para nela
    # (len_trim == 0 -> itest == 0), ntlines sai igual a icnti e a correcao
    # da linha 759 nao dispara.
    st, loc = _run_loc(_LAND, met=_met("branco.met", cauda="\n"))
    assert loc["ntlines"] == 3 and loc["icnti"] == 3
    # uma linha de texto no fim para a contagem do mesmo jeito (primeiro
    # caractere nao-branco fora da faixa 45..57)
    st, loc = _run_loc(_LAND, met=_met("texto.met", cauda="FIM DOS DADOS\n"))
    assert loc["ntlines"] == 3 and loc["icnti"] == 3

    # Cabecalho mentindo sobre o numero de linhas: vale o que foi contado,
    # menos 1. Aqui ha 5 linhas de dados -> ntlines = 6 -> icnti = 5.
    st, loc = _run_loc(_LAND, met=_met(
        "cinco.met",
        cabecalho="99 32 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 2.0",
        corpo="".join(_linha_met(h) for h in (12, 13, 14, 15, 16))))
    assert loc["ntlines"] == 6 and loc["icnti"] == 5 and st.iend == 5

    # `ntlines > 1` na guarda da linha 759: com UMA linha de dados e uma
    # linha de texto no fim, ntlines sai 1 (nao ha duplo de EOF) e a
    # correcao NAO dispara -- o `icnti` do cabecalho prevalece, mesmo
    # estando errado.
    st, loc = _run_loc(_LAND, met=_met(
        "uma.met",
        cabecalho="5 32 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 2.0",
        corpo=_linha_met(12), cauda="FIM\n"))
    assert loc["ntlines"] == 1 and loc["icnti"] == 5 and st.iend == 5

    # A faixa 45..57 do scan inclui '-', '.' e '/', nao so os digitos: uma
    # linha de dados que comeca com `-999` (ano faltante) ainda e contada.
    st, loc = _run_loc(_LAND, met=_met(
        "neg.met",
        cabecalho="99 32 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 2.0",
        corpo="".join(_linha_met(h, ano=-999) for h in (12, 13, 14))))
    assert loc["ntlines"] == 4 and loc["icnti"] == 3 and st.iend == 3

    # `iheight` <= 0 cai no padrao de 2 m (linha 727)
    st, _ = _run_loc(_LAND, met=_met(
        "alt0.met",
        cabecalho="3 32 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 0.0"))
    assert st.iheight == 2.0

    # `dstart` exige as DUAS condicoes: com timstep == 3 mas um unico
    # sitio, continua 2.
    st, _ = _run_loc(_LAND, met=_met(
        "ts3um.met",
        cabecalho="3 32 32.5 -114.0 50.0 1234 1 -7.0 3.0 -999.0 2.0"))
    assert st.timstep == 3.0 and st.met_count == 1 and st.dstart == 2

    # `dstart` (762-767): 14 so com timstep == 3 E mais de um sitio.
    # Aproveita para exercitar dois sitios: o segundo tem registro proprio
    # de 10 campos, seguido de duas linhas de texto.
    st, loc = _run_loc(_LAND, met=_met(
        "ts3.met",
        cabecalho="3 32 32.5 -114.0 50.0 1234 2 -7.0 3.0 -999.0 2.0",
        corpo="".join(_linha_met(h) for h in (12, 13, 14))
              + "0 0 40.0 -80.0 100.0 5678 0 0.0 0.0 0.0\n"
              + "Year DOY Hr Min ...\n  unidades\n"
              + "".join(_linha_met(h, ano=2011) for h in (12, 13, 14))))
    assert st.timstep == 3.0 and st.met_count == 2
    assert st.dstart == min(14, loc["icnti"])
    assert loc["latm"][1] == 32.5 and loc["latm"][2] == 40.0
    assert loc["longm"][1] == -114.0 and loc["longm"][2] == -80.0
    assert loc["vitd_indexm"][1] == 1234 and loc["vitd_indexm"][2] == 5678
    assert loc["metm"][2, 1, MetCol.YEAR] == _anint(2011.0, 15)

    # O cabecalho global E o registro do primeiro sitio: `latm(1)`/`longm(1)`
    # saem das posicoes 3 e 4, e `vitd_indexm(1)` da posicao 6 -- que no
    # cabecalho e o `pindex`. Nao e defeito, e o formato do arquivo.
    st, loc = _run_loc(_LAND)
    assert loc["latm"][1] == st.lat == 32.5
    assert loc["longm"][1] == st.mlong == -114.0
    assert loc["vitd_indexm"][1] == loc["pindex"] == 1234

    # conteudo do array de met.
    # PONTO DE ATENCAO 16: todo valor lido passa por
    # `anint(x*1d15)*1d-15`, que NAO e identidade. Nem o `1d15` nem o
    # `1d-15` sao exatos em binario, e para x acima de ~9 o produto `x*1d15`
    # ja passa de 2^53, onde o espacamento entre doubles e maior que 1 --
    # entao o caminho de ida e volta perde bits. 7.0 vira 7.000000000000001,
    # 12.0 vira 12.000000000000002, 1013.0 vira 1013.0000000000001; outros
    # (0.2, 180.0) sobrevivem por sorte de arredondamento. Nao e erro de
    # transcricao: e a mesma aritmetica IEEE que o Fortran faz. Por isso os
    # valores esperados aqui passam pelo mesmo `_anint`.
    def q(x):
        return _anint(x, 15)

    assert q(7.0) != 7.0 and q(12.0) != 12.0 and q(1013.0) != 1013.0
    assert q(0.2) == 0.2 and q(180.0) == 180.0         # nao atinge todos

    metm = loc["metm"]
    assert metm.shape == (st.met_count + 1, loc["icnti"] + 1, MAXCOL + 1)
    assert [metm[1, j, MetCol.HR] for j in (1, 2, 3)] == [q(12.0), q(13.0), q(14.0)]
    assert metm[1, 1, MetCol.TSOL] == q(800.0)
    assert metm[1, 1, MetCol.TSOIL] == q(20.0)
    assert metm[1, 1, MetCol.AP] == q(1013.0)

    # PONTO DE ATENCAO 8: a coluna 33 (MetCol.HI) nunca e lida do arquivo
    # -- `met_col` e 32, MAXCOL e 35 -- e fica com o sentinela. Com
    # ncols <= 32 as colunas 34 e 35 tambem ficam.
    assert metm[1, 1, MetCol.HI] == st.mflag == -999.0
    assert metm[1, 1, MetCol.VIS] == st.mflag
    assert metm[1, 1, MetCol.AER] == st.mflag
    # sem precipitacao mensuravel (<= 1e-4) o tipo de precipitacao vira 1
    assert metm[1, 1, MetCol.PT] == 1 and metm[1, 1, MetCol.PT2] == 1

    # ncols > 32: le 1..32 mais 34 e 35, saltando a 33
    st, loc = _run_loc(_LAND, met=_met(
        "c35.met",
        cabecalho="3 35 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 2.0",
        corpo="".join(_linha_met(h, extra=" 15.0 7") for h in (12, 13, 14))))
    metm = loc["metm"]
    assert st.ncols == 35
    assert metm[1, 1, MetCol.HR] == q(12.0)
    # PONTO DE ATENCAO 15: neste ramo o laco de quantizacao vai de 1 a 34,
    # mas a leitura cobre 1..32 mais 34 e 35. Os dois conjuntos nao batem, e
    # as duas pontas sao observaveis:
    #   * a coluna 34 (VIS) e lida E quantizada;
    #   * a coluna 35 (AER) e lida e NAO quantizada -- fica com o valor cru,
    #     que e diferente do quantizado (ver ponto 16);
    #   * a coluna 33 (HI) nao e lida mas E quantizada: o sentinela -999.0
    #     sai de la como -999.0000000000001, ou seja a coluna 33 deixa de
    #     ser IGUAL a `mflag`. E por isso que todo teste de dado faltante no
    #     codigo usa tolerancia -- a linha 1176, por exemplo, faz
    #     `aint(dabs(met(k,ip_zen)-mflag)*1d5)*1d-5 > eps` -- e nunca
    #     igualdade.
    assert metm[1, 1, MetCol.VIS] == q(15.0) != 15.0      # lida e quantizada
    assert metm[1, 1, MetCol.AER] == 7.0 != q(7.0)        # lida, nao quantizada
    assert metm[1, 1, MetCol.HI] == q(st.mflag)           # nao lida, quantizada
    assert metm[1, 1, MetCol.HI] != st.mflag
    assert abs(metm[1, 1, MetCol.HI] - st.mflag) < 1e-12  # mas passa por tolerancia

    # precipitacao de verdade: o tipo fica como veio do arquivo
    st, loc = _run_loc(_LAND, met=_met(
        "chuva.met",
        corpo="".join(_linha_met(h, prec="2.5") for h in (12, 13, 14))))
    assert loc["metm"][1, 1, MetCol.PREC] == 2.5
    assert loc["metm"][1, 1, MetCol.PT] == 0
    assert loc["metm"][1, 1, MetCol.PT2] == 1        # prec2 continua 0.0

    # O limite e `<=`, nao `<`: prec exatamente 1e-4 AINDA zera o tipo. E o
    # 1e-4 sobrevive a quantizacao, entao a comparacao e exata.
    assert q(1e-4) == 1e-4
    st, loc = _run_loc(_LAND, met=_met(
        "limite.met",
        corpo="".join(_linha_met(h, prec="0.0001") for h in (12, 13, 14))))
    assert loc["metm"][1, 1, MetCol.PREC] == 1e-4
    assert loc["metm"][1, 1, MetCol.PT] == 1

    # descarte da ultima linha por salto de ano (815-818), sempre olhando
    # `metm(1,...)` -- sitio 1 fixo, independente de met_count
    st, loc = _run_loc(_LAND, met=_met(
        "salto.met",
        corpo=_linha_met(12) + _linha_met(13) + _linha_met(14, ano=2015)))
    assert loc["icnti"] == 2 and st.iend == 2
    # salto de exatamente 1 ano NAO e descartado (`> 1`, nao `>= 1`)
    st, loc = _run_loc(_LAND, met=_met(
        "salto1.met",
        corpo=_linha_met(12) + _linha_met(13) + _linha_met(14, ano=2011)))
    assert loc["icnti"] == 3 and st.iend == 3

    # met inexistente: ramo de erro do original (divergencia D1)
    _para(_LAND.replace("@MET@", os.path.join(d, "nao_existe.met")),
          "Error opening met file")

    # ---- locais do STEP 1 que so agora ficaram verificaveis ----
    st, loc = _run_loc(_LAND)
    assert loc["output1"] == "saida.fst" and loc["output5"] == "saida.snw"
    assert loc["surfdepth"] == 5e-2
    assert loc["ttest"] == "n" and loc["mtest"] == "n"
    # Ponto de atencao 13 no `fasstusersoil`: sem `!` na linha e sem leitura
    # (lcount == 0), o strip da linha 629 deixa o nome VAZIO -- e a linha
    # 878, ja no STEP 3, o restaura para o padrao. A captura acontece depois
    # dela, entao o que se ve aqui e o nome ja consertado.
    assert loc["lcount"] == 0 and loc["fasstusersoil"] == "fasstusersoil.inp"
    # ramo multiponto: agora da para ver os dois nomes lidos
    _MULTI = ("@MET@ !\n0 1 0.5 1 !\n1.5 0.0 1 1 0 100 0 !\n"
              "complexo.vitd !\nmeusolo.inp !\n")
    # `tolera=(IndexError,)`: ver o ponto de atencao 5 mais abaixo -- com
    # `num_frq == 0` o laco de lacunas do STEP 3 indexa `metm` com -1 e
    # `_ix` levanta. A captura do fim do hot start ja aconteceu, entao
    # `loc` vem preenchido.
    st, loc = _run_loc(_MULTI, tolera=(IndexError,))
    assert loc["vitd_input"] == "complexo.vitd"
    assert loc["fasstusersoil"] == "meusolo.inp"

    # ---- STEP 3: unidades de arquivo e hot start (linhas 825-1084) ----

    # As 16 unidades do STEP 3 ficam em campos de `FasstState`, e os nomes
    # estao em `_ALLOCATABLE_FIELDS` para que `reset()` nao as zere entre
    # pontos de grade.
    from .state import _UNIT_FIELDS
    assert len(_UNIT_FIELDS) == 16
    from .state import _ALLOCATABLE_FIELDS
    assert _UNIT_FIELDS <= _ALLOCATABLE_FIELDS

    st, loc = _run_loc(_LAND)
    # sempre abertas: as 4 tabelas de propriedades, o solo do usuario, a
    # saida principal, o inferred.out e os dois binarios de hot start
    for campo in ("unit30", "unit31", "unit32", "unit33", "unit90",
                  "unit2", "unit10", "unit17", "unit18"):
        assert getattr(st, campo) is not None, campo
    # `nprint == 1` no cabecalho -> unidade 3 aberta; `flprint == 1` ->
    # unidade 24; `sprint == 0` -> unidade 55 fica fechada; `vegint == 0`
    # (nao > EPS) -> unidade 28 fica fechada
    assert st.unit3 is not None and st.unit24 is not None
    assert st.unit55 is None and st.unit28 is None
    # ponto unico: nem a unidade 5 nem os binarios da execucao anterior
    assert st.unit5 is None and st.unit11 is None and st.unit12 is None

    # os handles sobrevivem a um `reset()` (o laco multiponto chama
    # `zero_parameters` por ponto, muito depois destes `open`)
    handles = {c: getattr(st, c) for c in _UNIT_FIELDS}
    st.reset()
    assert all(getattr(st, c) is handles[c] for c in _UNIT_FIELDS)

    # `npos`: posicao do ultimo '.' do nome do met, e a base dos nomes
    # derivados. `met_file_name[:npos]` INCLUI o ponto.
    st, loc = _run_loc(_LAND)
    assert loc["met_file_name"][loc["npos"] - 1] == "."
    assert loc["met_file_name"][loc["npos"]:] == "met"
    assert loc["fnlen"] == loc["npos"] + 3        # atribuido e nunca lido
    # os dois binarios de hot start saem do nome do met sem a extensao
    assert loc["ro1"].endswith("_ro_soil.dat")
    assert loc["ro2"].endswith("_ro_time.dat")
    assert loc["ro1"][:-len("_ro_soil.dat")] == loc["met_file_name"][:loc["npos"] - 1]
    assert os.path.exists(loc["ro1"]) and os.path.exists(loc["ro2"])
    assert os.path.exists(os.path.join(d, "inferred.out"))
    assert os.path.exists(os.path.join(d, "saida.fst"))      # output1

    # Cabecalho de 8 linhas da saida de fluxos (948-959), escrito na hora.
    # O handle continua aberto (o driver so fecha no STEP 6), e Python
    # bufferiza -- dai o flush antes de ler.
    st.unit24.f.flush()
    with open(os.path.join(d, "saida.flx")) as f:
        flx = f.read()
    assert flx.startswith("Surface Flux File\nConventions used: \n")
    assert "+ => flux towards surface" in flx
    assert "sheat > 0 -> Tair > Tgr => surface warms, air cools" in flx
    assert "lheat < 0 -> qair < qgr => evaporation; surface cools, air warms" in flx
    assert len(flx.rstrip("\n").split("\n")) == 8

    # sem '.' no nome do met, `npos` volta a `mlen` (linha 894)
    st, loc = _run_loc(_LAND, met=_met("semponto"))
    assert loc["npos"] == loc["mlen"] == len(loc["met_file_name"])

    # O laco de busca do '.' comeca EM `mlen`, nao em `mlen - 1`: so um
    # nome cujo ULTIMO caractere e '.' distingue as duas coisas, e para
    # isso ele precisa ter um '.' anterior tambem. Com 'a.b.' o original
    # acha o ponto final (npos = 4); comecando um antes acharia o do
    # meio (npos = 2) e as saidas sairiam "a.fst" em vez de "a.b.fst".
    st, loc = _run_loc(_LAND, met=_met("a.b."))
    assert loc["npos"] == loc["mlen"] == len(loc["met_file_name"])
    # `output1` nao serve aqui: no ponto unico ele vem do proprio `.inp`.
    # `ro1` e derivado do nome do met e serve.
    assert loc["ro1"].endswith("a.b_ro_soil.dat")


    # o casamento de indice VITD: `kk` e o sitio casado, e `i` sai do laco
    # valendo kk+1 -- e o `i = i - 1` da linha 1063 o traz de volta a kk
    # (primeira metade do ponto de atencao 5)
    st, loc = _run_loc(_LAND)
    assert loc["tflag"] == 1
    assert loc["kk"] == 1 and loc["i"] == loc["kk"]
    assert loc["vitd_indexm"][loc["kk"]] == st.vitd_index == 1234
    # `fdi = freq_id` e depois `freq_id = fdi`: ida e volta sem efeito no
    # ponto unico. `freq_id` vale 1, posto por `zero_parameters`.
    assert loc["fdi"] == st.freq_id == 1

    # vitd_index que nao casa com nenhum sitio do met: para o programa
    _SEM_CASO = _CABECALHO.replace("1234                         ! vitd_index",
                                   "9999                         ! vitd_index")
    _para(_SEM_CASO + _TERRA + "0 ! vegl\n0 ! vegh\n" + _SOLO,
          "Something wrong with input meteorological file")

    # `nprint`/`flprint`/`sprint`/`vegint` controlam quais saidas abrem
    _CAB2 = _CABECALHO.replace("1.5 0.0 1 1 0 100 0", "1.5 0.5 0 0 1 100 0")
    st, loc = _run_loc(_CAB2 + _TERRA + "0 ! vegl\n0 ! vegh\n" + _SOLO)
    assert loc["nprint"] == 0 and st.unit3 is None
    assert loc["flprint"] == 0 and st.unit24 is None
    assert loc["sprint"] == 1 and st.unit55 is not None
    assert loc["vegint"] == 0.5 and st.unit28 is not None     # vegint > EPS

    # A guarda da vegetacao e `vegint > eps`, nao `>= eps`: com vegint
    # valendo EXATAMENTE eps a unidade 28 NAO abre. O `repr` de
    # `sys.float_info.epsilon` faz round-trip, e esta leitura nao passa
    # por `anint`, entao o valor chega inteiro pelo `.inp`.
    _CAB_EPS = _CABECALHO.replace("1.5 0.0 1 1 0 100 0",
                                  "1.5 %r 0 0 1 100 0" % EPS)
    st, loc = _run_loc(_CAB_EPS + _TERRA + "0 ! vegl\n0 ! vegh\n" + _SOLO)
    assert loc["vegint"] == EPS and st.unit28 is None


    # ---- multiponto: unidade 5 e os nomes de saida derivados do met ----
    # `complexo.vitd` tem uma linha so, entao `num_frq` sai 0 e o laco de
    # `read_complex` nao roda -- da para atravessar o STEP 3 inteiro.
    st, loc = _run_loc(_MULTI, tolera=(IndexError,))
    assert st.unit5 is not None
    assert loc["num_frq"] == 0
    base = loc["met_file_name"][:loc["npos"]]
    assert loc["output1"] == base + "fst" and loc["output2"] == base + "grd"
    assert loc["output3"] == base + "flx" and loc["output4"] == base + "vgt"
    assert loc["output5"] == base + "snw"
    assert os.path.exists(loc["output1"])
    # PONTO DE ATENCAO 5, segunda metade: no multiponto o `i` que chega ao
    # `i = i - 1` e o contador de LINHAS da unidade 5, nao um indice de
    # sitio. Com N linhas no arquivo complexo: o laco de contagem sai com
    # `i = N+1` (a ultima volta e a que ve o EOF), `i = i - 2` deixa `N-1`,
    # o laco de `read_complex` nao mexe em `i` quando nao roda, e o
    # `i = i - 1` da linha 1063 fecha em `N-2`. Com N = 1 da -1 -- um
    # indice NEGATIVO, que o Dia 8 vai usar para indexar `metm`.
    assert loc["i"] == -1 and loc["kk"] == 0
    assert loc["tflag"] == 0        # nenhum casamento, e aqui nao ha `stop`
    # ...e e com esse -1 que o laco de lacunas (1095-1099) indexa `metm`.
    # Em NumPy -1 pegaria o ULTIMO sitio em silencio; `_ix` levanta
    # (divergencia D4). Este e o unico lugar em que o demo pede
    # `tolera=(IndexError,)`, e e o que prova que D4 esta ativa.
    try:
        _exec(_MULTI, None, FasstState(), None)
    except IndexError as e:
        assert "ponto de atencao 5" in str(e)
    else:
        raise AssertionError("esperado IndexError no laco de lacunas")

    # com mais de uma linha no complexo, `read_complex` e chamada e levanta
    vitd3 = os.path.join(d, "tres.vitd")
    with open(vitd3, "w") as f:
        f.write("cabecalho\nregistro 1\nregistro 2\n")
    st, loc = _run_loc(_MULTI.replace("complexo.vitd", "tres.vitd"))
    assert loc == {}, "read_complex deveria ter levantado antes da captura"
    assert st.unit5 is not None      # mas as unidades ja estavam abertas

    # ---- hot start (infer_test == 1): unidades 11 e 12 ----
    _HOT = ("@MET@ !\n1 0 0.5 1 !\nvelho_ro_soil.dat ! ro1o\n"
            "1.5 0.0 1 1 0 100 0 !\n1234 !\n" + "saida.fst !\n" * 5 +
            "0.0 0.0 !\n0.01 !\n0 !\n")
    st, loc = _run_loc(_HOT)
    # `read_old_data` e stub e levanta, entao a captura nao acontece --
    # mas as unidades 11 e 12 ja foram abertas antes dela
    assert loc == {}
    assert st.unit11 is not None and st.unit12 is not None
    assert st.istart == 2
    # ponto de atencao 11: o nome da unidade 12 vem de cortar 9 caracteres
    # ("_soil.dat") do nome da 11, sem validacao
    assert st.unit11.path == "velho_ro_soil.dat"
    assert st.unit12.path == "velho_ro_time.dat"
    # e sao binarios ('rb'), nao texto
    assert "b" in st.unit11.f.mode and "b" in st.unit17.f.mode

    # ponto 11 na pratica: um nome que NAO termina em "_soil.dat" gera um
    # nome truncado em lugar nenhum. "abc.dat" (7 caracteres) -> [:-2] ->
    # "abc.d" sem os 9, ou seja len-9 = -2 -> "abcde"[:-2]
    nome_curto = "curto.dat"                                  # 9 caracteres
    with open(os.path.join(d, nome_curto), "wb") as f:
        f.write(b"")
    with open(os.path.join(d, "_time.dat"), "wb") as f:
        f.write(b"")
    st, loc = _run_loc(_HOT.replace("velho_ro_soil.dat", nome_curto))
    assert st.unit12.path == "_time.dat"     # len("curto.dat") - 9 == 0

    # `inferred.out` e aberto com `status='unknown'` e o comentario da linha
    # 984 diz "erase it first if it already exists": o conteudo anterior
    # tem de ir embora.
    with open(os.path.join(d, "inferred.out"), "w") as f:
        f.write("sobra de uma execucao anterior\n")
    st, loc = _run_loc(_LAND)
    st.unit10.f.flush()
    with open(os.path.join(d, "inferred.out")) as f:
        assert f.read() == ""

    # ---- o resto do bloco de hot start, com `read_old_data` postico ----
    # Com o stub real a execucao levanta na linha 1069 e nada depois dela e
    # alcancavel: a chamada em si, o laco de `t1o`, o `t1n` e a busca de
    # `mpos`. Trocar a funcao por uma que registra os argumentos e devolve
    # um `timeo` sintetico cobre os quatro.
    _rod_real = globals()["read_old_data"]

    def _hot(horas, met=None, hot=None):
        """Roda o caminho de hot start com `read_old_data` postico.

        `horas` sao as horas dos registros antigos; `moverlapr` sai do
        tamanho da lista (no original e `read_old_data` que o escreve, e
        ninguem o zera, entao da para semea-lo).
        """
        chamadas = []

        def _falso(d1i_, ii_, wstart_, sdensi_, phie_, timeo_):
            chamadas.append((d1i_, ii_, wstart_, sdensi_, phie_, timeo_))
            t = [[0.0] * 5 for _ in range(MOVERLAP + 1)]
            for kk_, hr_ in enumerate(horas, start=1):
                t[kk_][1] = 2010.0
                t[kk_][2] = 100.0
                t[kk_][3] = float(hr_)
                t[kk_][4] = 0.0
            return 123.0, 0.5, t

        semente, loc = FasstState(), {}
        semente.moverlapr = len(horas)
        globals()["read_old_data"] = _falso
        try:
            _exec(hot if hot else _HOT, met, semente, loc)
        finally:
            globals()["read_old_data"] = _rod_real
        return semente, loc, chamadas

    # a chamada da linha 1069, argumento por argumento
    st, loc, chamadas = _hot((10, 11, 13, 14))
    assert len(chamadas) == 1
    d1i_, ii_, wstart_, sdensi_, phie_, timeo_ = chamadas[0]
    assert d1i_ == 3                             # `d1i = 3`, linha 1066
    assert ii_ == st.dstart                      # o 2o argumento e `dstart`
    # `max(moverlap, istart)`: MOVERLAP e 15 e istart vale 2 com
    # infer_test == 1, entao o max e sempre o 15
    assert wstart_ == MOVERLAP == 15 and st.istart == 2
    assert len(timeo_) == MOVERLAP + 1           # timeo(moverlap,4)
    # os valores de retorno voltam para as locais do driver
    assert loc["sdensi"] == 123.0 and loc["phie"] == 0.5

    # `t1o` recebe as datas antigas; `t1n` e a data da 1a linha do met do
    # sitio casado (`metm(kk,1,1..4)`), que e 2010/100 12:00
    assert [loc["t1o"][k] for k in (1, 2, 3, 4)] == [
        met_date(2010.0, 100.0, h, 0.0) for h in (10.0, 11.0, 13.0, 14.0)]
    # `loc["_hot"]`, nao `loc`: o laco de lacunas logo abaixo reescreve
    # `t1n` (e `t1oo`) com o par de datas da ULTIMA volta dele.
    assert loc["_hot"]["t1n"] == met_date(
        *(loc["metm"][loc["kk"], 1, c] for c in (1, 2, 3, 4)))
    # a busca de `mpos`: o ultimo k com t1o[k] < t1n <= t1o[k+1].
    # 11h < 12h <= 13h -> k = 2.
    assert st.mpos == 2

    # O 2o argumento e `dstart`, nao `istart`: com timstep == 3 e dois
    # sitios eles passam a ser diferentes (dstart = min(14, icnti) = 3,
    # istart = 2), e e o dstart que chega la.
    met2 = _met("hot3.met",
                cabecalho="3 32 32.5 -114.0 50.0 1234 2 -7.0 3.0 -999.0 2.0",
                corpo="".join(_linha_met(h) for h in (12, 13, 14))
                      + "0 0 40.0 -80.0 100.0 5678 0 0.0 0.0 0.0\n"
                      + "Year DOY Hr Min ...\n  unidades\n"
                      + "".join(_linha_met(h, ano=2011) for h in (12, 13, 14)))
    st, loc, chamadas = _hot((10, 11, 13, 14), met=met2)
    assert st.dstart == 3 and st.istart == 2 and st.dstart != st.istart
    assert chamadas[0][1] == 3

    # Comparacoes exatas: com um registro antigo CAINDO em cima de `t1n`
    # (12h), o `<` estrito do lado esquerdo e o `>=` do direito importam.
    # Original: k=1 (10h < 12h <= 12h) -> mpos = 1.
    st, loc, _ = _hot((10, 12, 13, 14))
    assert st.mpos == 1

    # quando nenhum par cerca `t1n`, `mpos` cai no fallback `= moverlapr`
    st, loc, _ = _hot((18, 19, 20, 21))
    assert st.mpos == st.moverlapr == 4

    # `moverlapr` no maximo (= MOVERLAP): o laco vai so ate moverlapr-1
    # justamente para nao ler `t1o(moverlapr+1)`, que estouraria o array
    # declarado `t1o(moverlap)`.
    st, loc, _ = _hot(tuple(range(1, MOVERLAP + 1)))
    assert st.moverlapr == MOVERLAP == 15
    assert len(loc["t1o"]) == MOVERLAP + 1       # indices 0..15
    assert st.mpos == 11                         # 11h < 12h <= 12h

    # O caso que realmente exercita o limite do laco: `t1n` DEPOIS de todos
    # os registros antigos. Ai a primeira metade da condicao e verdadeira
    # ate k = moverlapr, e so o limite `moverlapr - 1` evita a leitura de
    # `t1o[moverlapr + 1]`. Com met comecando as 23h e 15 registros antigos
    # de 1h a 15h, nenhum par cerca `t1n` -> fallback.
    met23 = _met("hot23.met",
                 corpo="".join(_linha_met(h) for h in (23, 0, 1)))
    st, loc, _ = _hot(tuple(range(1, MOVERLAP + 1)), met=met23)
    assert loc["_hot"]["t1n"] == met_date(2010.0, 100.0, 23.0, 0.0)
    assert loc["t1o"][MOVERLAP] < loc["_hot"]["t1n"]   # NAO curto-circuita
    assert st.mpos == st.moverlapr == MOVERLAP

    # `freq_id = fdi` (linha 1067): no ponto unico e ida e volta, mas o
    # valor tem de chegar inteiro -- `zero_parameters` o deixa em 1.
    st, loc, _ = _hot((10, 11, 13, 14))
    assert loc["fdi"] == 1 and st.freq_id == 1

    # ---- STEP 3: lacunas de tempo, maxlines e alocacao (1086-1153) ----

    # Um met com as horas que a lista pedir, e o `.inp` de terra firme
    # correspondente. `timstep` sai do cabecalho global do met.
    def _lacunas(horas, timstep=1.0, ano=2010, doys=None, minutos=None):
        nome = "gap%d.met" % len(_feitos)
        _feitos.append(nome)
        if doys is None:
            doys = [100] * len(horas)
        if minutos is None:
            minutos = [0] * len(horas)
        cab = ("%d 32 32.5 -114.0 50.0 1234 1 -7.0 %g -999.0 2.0"
               % (len(horas), timstep))
        anos = ano if isinstance(ano, (list, tuple)) else [ano] * len(horas)
        corpo = "".join(_linha_met(h, ano=a, doy=dd, minuto=mm)
                        for h, a, dd, mm in zip(horas, anos, doys, minutos))
        met = _met(nome, cabecalho=cab, corpo=corpo)
        return _run_loc(_LAND, met=met)

    _feitos = []

    # sem lacuna: tres horas consecutivas -> mtsteps = 0
    st, loc = _lacunas((12, 13, 14))
    assert loc["tstps"] == 24 and loc["rtstps"] == 24.0   # 24 / timstep
    assert loc["mtsteps"] == 0
    assert st.iend == 3 and st.maxlines == 3               # iend + mtsteps
    # `daylim` da ultima volta: todas as voltas ficaram no mesmo ano
    assert loc["daylim"] == 1.042e0

    # uma hora faltando (12, 14, 15). O produto da primeira volta e
    # 1.999999974 -- logo ABAIXO de 2 por causa do ponto 16 -- e o
    # `int()` o trunca para 1, que por sorte e o numero certo de
    # registros faltantes.
    st, loc = _lacunas((12, 14, 15))
    assert loc["mtsteps"] == 1 and st.maxlines == st.iend + 1 == 4

    # ...mas o acerto E sorte: entre 20h e 22h o erro do ponto 16 muda de
    # sinal (produto = 2.000000852) e a MESMA lacuna de um registro conta
    # 2. Ou seja `mtsteps` depende de QUAIS horas a lacuna cai entre.
    st, loc = _lacunas((20, 22, 23))
    assert loc["mtsteps"] == 2 and st.maxlines == 5

    # lacuna de tres horas: produto 3.999999945 -> 3
    st, loc = _lacunas((12, 16, 17))
    assert loc["mtsteps"] == 3 and st.maxlines == 6

    # PONTO DE ATENCAO 9: `timstep <= 0` deixa `rtstps` no zero da
    # declaracao (o `if` da 1092 nao dispara), `tstps` sai 0 e todo
    # produto sai 0 -- a lacuna de tres horas NAO e detectada.
    st, loc = _lacunas((12, 16, 17), timstep=0.0)
    assert loc["rtstps"] == 0.0 and loc["tstps"] == 0
    assert loc["mtsteps"] == 0 and st.maxlines == st.iend
    st, loc = _lacunas((12, 16, 17), timstep=-1.0)
    assert loc["tstps"] == 0 and loc["mtsteps"] == 0

    # virada de ano: `daylim` sobe de 1.042 para 25.07 (linha 1106).
    # Entre 2010/365 23h e 2011/001 00h o produto e 25.000000099 --
    # porque `met_date` divide o dia-do-ano por 366, nao por 365 -- e
    # fica logo ABAIXO de 25.07, entao a virada normal nao conta lacuna.
    st, loc = _lacunas((23, 0, 1), ano=(2010, 2011, 2011),
                       doys=(365, 1, 1))
    assert loc["daylim"] == 1.042e0        # ultima volta: ja no mesmo ano
    assert loc["mtsteps"] == 0
    # uma hora faltando NA virada passa de 25.07 e conta 25 -- o vies de
    # +24 embutido no 366 vai inteiro para `maxlines`.
    st, loc = _lacunas((22, 0, 1), ano=(2010, 2011, 2011),
                       doys=(365, 1, 1))
    assert loc["mtsteps"] == 25 and st.maxlines == st.iend + 25

    # PONTO DE ATENCAO 19 (novo): o ramo bissexto (fator 367) e
    # INALCANCAVEL para qualquer ano real. `mody` e calculado sobre o ano
    # JA QUANTIZADO pelo ponto 16, e `anint(2012*1d15)*1d-15` devolve
    # 2012.0000000000002274 -> `mody` = 2.27e-13, muito acima de
    # `eps` = 2.22e-16. Logo 2012 cai no ramo 366, igual a 2010.
    st, loc = _lacunas((20, 22, 23), ano=2012)
    assert loc["mody"] != 0.0 and abs(loc["mody"]) > EPS
    assert loc["mtsteps"] == 2                 # identico ao caso de 2010
    # Os unicos anos que passam o teste sao os multiplos de 4 pequenos o
    # suficiente para `ano*1d15` caber exato num double (< 2**53): 0, 4,
    # 8 e 16 -- note que 12 NAO passa (1.78e-15 > eps). Com ano 4 o ramo
    # 367 fica alcancavel, e entre 10h e 12h ele MUDA o resultado: o
    # produto e 2.000000851 com 367 e 1.994551258 com 366.
    st, loc = _lacunas((10, 12, 13), ano=4)
    assert loc["mody"] == 0.0
    assert loc["mtsteps"] == 2                 # com 366 daria 1
    # `modyo` (ponto de atencao 17) e calculado e nunca lido: aqui o ano
    # anterior e o mesmo, entao ele acaba valendo o mesmo que `mody`.
    assert loc["modyo"] == loc["mody"]

    # `maxlines` no hot start (1122): ganha `moverlapr + 1`. O `_hot` usa
    # `read_old_data` postico, que e quem escreve `moverlapr`.
    st, loc, _ = _hot((10, 11, 13, 14))
    assert st.moverlapr == 4
    assert st.infer_test == 1
    assert st.maxlines == st.iend + loc["mtsteps"] + st.moverlapr + 1

    # ---- a alocacao (1126-1148) ----
    # As formas tem de ser exatamente as que o slicing de
    # `zero_mxl_params` assume -- e a prova e preencher tudo com valor
    # nao-nulo, chamar `zero_mxl_params` e conferir os zeros.
    st, loc = _lacunas((12, 14, 15))
    n = st.maxlines
    assert n == 4
    for nome in _RUN_REAL_1D:
        a = getattr(st, nome)
        assert a is not None and a.shape == (n + 1,), nome
        assert a.dtype == np.float64, nome
    assert st.slushy.shape == (n + 1,) and st.slushy.dtype == np.int64
    assert st.sstate.shape == (n + 1,) and st.sstate.dtype.kind == "U"
    assert st.met.shape == (n + 1, MAXCOL + 1)
    assert st.dmet1.shape == (n + 1, 14)
    assert st.airt.shape == (3, n + 1)
    assert st.canopy_temp.shape == (NCLAYERS + 1, n + 1)
    # os 40 nomes de `_ALLOCATABLE_FIELDS` que nao sao unidades: nenhum
    # pode ter ficado em None
    from .state import _ALLOCATABLE_FIELDS, _UNIT_FIELDS
    for nome in _ALLOCATABLE_FIELDS - _UNIT_FIELDS:
        assert getattr(st, nome) is not None, nome

    # ordem obrigatoria: alocar ANTES de `zero_mxl_params`, que itera
    # sobre estas arrays. Preenche com lixo e confere o zeramento.
    for nome in _RUN_REAL_1D:
        getattr(st, nome)[:] = 7.5
    st.met[:, :] = 7.5
    st.dmet1[:, :] = 7.5
    st.airt[:, :] = 7.5
    st.canopy_temp[:, :] = 7.5
    st.sstate[:] = "x"
    zero_mxl_params(st)
    for nome in _RUN_REAL_1D:
        assert not getattr(st, nome)[1:n + 1].any(), nome
    assert not st.met[1:n + 1, 1:MAXCOL + 1].any()
    assert not st.airt[1:3, 1:n + 1].any()
    assert not st.canopy_temp[1:NCLAYERS + 1, 1:n + 1].any()
    assert list(st.sstate[1:n + 1]) == [""] * n
    # ponto 2 de `zero_mxl_params`: as colunas 11..13 de `dmet1` NAO sao
    # zeradas -- preservado, e so observavel porque a alocacao as criou
    assert not st.dmet1[1:n + 1, 1:11].any()
    assert (st.dmet1[1:n + 1, 11:14] == 7.5).all()
    # e `slushy` nao e tocado (ponto 9 de `zero_mxl_params`)
    assert (st.slushy == 0).all()

    # ---- multiponto COM registros: `read_complex` postico ----
    # Todo fixture multiponto acima tem `num_frq == 0`, que e justamente
    # o caso degenerado do ponto de atencao 5. Com um `read_complex`
    # postico o laco de registros VITD roda de verdade, `i` termina no
    # sitio casado e o laco de lacunas indexa `metm` com indice valido.
    _rc_real = globals()["read_complex"]

    def _multi(nlinhas, met=None):
        """Roda o multiponto com N linhas no `.vitd` e read_complex falso.

        `num_frq` sai `N - 1` (ver o ponto de atencao 5 acima), entao o
        laco chama o postico `N - 1` vezes.
        """
        import inspect
        chamadas = []

        def _falso(state_, wtype_, io_, lcount_, nm0_, wvel_, wdepth_,
                   mtest_, sname_):
            # `bind` contra a assinatura REAL: sem isto o postico
            # esconderia um erro de aridade na chamada do driver.
            inspect.signature(_rc_real).bind(
                state_, wtype_, io_, lcount_, nm0_, wvel_, wdepth_,
                mtest_, sname_)
            chamadas.append(state_.freq_id)
            # o real le UM registro VITD e escreve estes dois campos
            state_.vitd_index = 1234          # = pindex do met do fixture
            state_.freq_id = 10 * len(chamadas)
            return (wtype_, 0, lcount_, nm0_, wvel_, wdepth_, mtest_,
                    sname_)

        nome = "multi%d.vitd" % nlinhas
        with open(os.path.join(d, nome), "w") as f:
            f.write("cabecalho do arquivo complexo\n")
            for kk_ in range(nlinhas - 1):
                f.write("registro %d\n" % (kk_ + 1))
        st_, loc_ = FasstState(), {}
        globals()["read_complex"] = _falso
        try:
            _exec(_MULTI.replace("complexo.vitd", nome), met, st_, loc_)
        finally:
            globals()["read_complex"] = _rc_real
        return st_, loc_, chamadas

    # 3 linhas -> num_frq = 2 -> duas chamadas
    st, loc, chamadas = _multi(3)
    assert loc["num_frq"] == 2 and len(chamadas) == 2
    # `fdi = freq_id` so na PRIMEIRA volta (`if j == 1`, linha 1047): o
    # postico poe 10 na 1a chamada e 20 na 2a, e `fdi` fica com o 10.
    assert loc["fdi"] == 10 and st.freq_id == 20
    # o casamento aconteceu, e agora `i` SAI valendo o sitio casado --
    # ao contrario do caso `num_frq == 0`, em que ele vale -1
    assert loc["tflag"] == 1 and loc["kk"] == 1 and loc["i"] == 1
    # e por isso o laco de lacunas rodou sobre o sitio de verdade
    assert loc["mtsteps"] == 0 and st.maxlines == st.iend == 3
    # a unidade 5 voltou ao topo depois do laco (`rewind(5)`, linha 1053)
    assert st.unit5 is not None

    # uma lacuna no met, agora pelo caminho multiponto
    cab = "3 32 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 2.0"
    met_gap = _met("multigap.met", cabecalho=cab,
                   corpo="".join(_linha_met(h) for h in (12, 14, 15)))
    st, loc, chamadas = _multi(3, met=met_gap)
    assert loc["mtsteps"] == 1 and st.maxlines == 4

    # 2 linhas -> num_frq = 1 -> uma chamada, e `fdi` e o unico freq_id
    st, loc, chamadas = _multi(2)
    assert loc["num_frq"] == 1 and len(chamadas) == 1
    assert loc["fdi"] == 10 and st.freq_id == 10

    # ---- quatro cenarios que a 1a rodada de mutacao do Dia 8 pediu ----

    # (1) `tstps = int(rtstps)`, nao `nint`/`round`. Com os timsteps
    #     usuais (1, 3) o quociente 24/timstep e inteiro e os dois dao o
    #     mesmo. Com `timstep = 9` da 2.666...: `int` deixa 2, `round`
    #     deixaria 3 -- e com tstps = 3 a lacuna de 9 horas passaria a
    #     ser detectada (mtsteps = 2).
    st, loc = _lacunas((12, 21, 30), timstep=9.0)
    assert loc["rtstps"] == 24e0 / 9e0 and loc["tstps"] == 2
    assert loc["mtsteps"] == 0 and st.maxlines == st.iend

    # (2) `mody` sai do ano do registro NOVO (k), nao do anterior (k-1).
    #     So e observavel quando os dois anos caem em ramos diferentes, e
    #     pelo ponto 19 isso exige anos pequenos: 4 -> 5. O ano 4 da
    #     mody = 0 (ramo 367) e o 5 da mody = 1 (ramo 366).
    st, loc = _lacunas((12, 0, 1), ano=(4, 5, 5), doys=(365, 1, 1))
    assert loc["mody"] == 1.0                 # ultima volta, ano 5
    assert loc["mtsteps"] == 59                # com mody de k-1 daria 60

    # (3) O fator da GUARDA e o mesmo do `int()` -- a expressao aparece
    #     duplicada no original, e as duas copias tem de usar o MESMO
    #     fator. Com 1h02.4min entre dois registros o produto da 1.0400
    #     com 366 (abaixo de daylim = 1.042, nao conta lacuna) e 1.0428
    #     com 367 (acima, contaria). A coluna de minutos e lida como
    #     real, entao 2.4 passa.
    st, loc = _lacunas((12, 13, 14), minutos=(0, 2.4, 2.4))
    assert loc["mtsteps"] == 0           # com a guarda em 367 daria 1
    #     e o espelho, no ramo bissexto: ano 4 com 1h02.6min da 1.04333
    #     com 367 (acima de daylim, conta 1) e 1.04050 com 366 (abaixo).
    st, loc = _lacunas((12, 13, 14), ano=4, minutos=(0, 2.6, 2.6))
    assert loc["mody"] == 0.0 and loc["mtsteps"] == 1

    # (4) `maxlines` do hot start soma os TRES termos. Com o met padrao
    #     `mtsteps` e 0 e o termo desaparece; com uma lacuna ele aparece.
    cab_hot = "3 32 32.5 -114.0 50.0 1234 1 -7.0 1.0 -999.0 2.0"
    met_hot_gap = _met("hotgap.met", cabecalho=cab_hot,
                       corpo="".join(_linha_met(h) for h in (12, 14, 15)))
    st, loc, _ = _hot((10, 11, 13, 14), met=met_hot_gap)
    assert loc["mtsteps"] == 1 and st.moverlapr == 4
    assert st.maxlines == st.iend + 1 + 4 + 1

    print("demo OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
    else:
        main(sys.argv)
