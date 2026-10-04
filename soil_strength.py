"""
soil_strength.py -- traducao de soil_strength.F90

Fonte Fortran: soil_strength.F90 (763 linhas). Uma unica sub-rotina
publica (`soil_strength`), sem nenhuma sub-rotina chamada ("no
subroutines called" no cabecalho do .F90). Estima a resistencia do
solo em cada no do perfil -- indice de cone (CI), indice de cone
remoldado (RCI) e CBR (California Bearing Ratio) -- a partir da
umidade, do tipo de solo (USCS) e, quando disponivel, da cabeca de
pressao, e agrega esses valores em duas faixas de profundidade
(camada superficial, ate 6 polegadas, e camada logo abaixo, de 6 a 12
polegadas) para produzir tres saidas: a resistencia media da camada
superficial (ss_strengtho), da camada abaixo (ll_strengtho) e o CBR
medio da camada superficial (scbro).

Assinatura original:
    subroutine soil_strength(nstr_flag,ss_strengtho,ll_strengtho,scbro)
    integer(ip),intent(in):: nstr_flag(maxn)
    real(dp),intent(out):: ss_strengtho,ll_strengtho,scbro

Saidas calculadas mas nao devolvidas pela sub-rotina
------------------------------------------------------
O cabecalho do .F90 ainda tem uma linha comentada que revela uma
assinatura mais antiga, com mais saidas:
    !real(dp),intent(out):: sci,srci,ss_strength,ll_strength,scbr
Na assinatura atual, sci/srci/ss_strength/ll_strength/scbr (o metodo
"Steve Grant, George Mason", baseado em indices de cone geotecnicos
classicos) sao variaveis puramente locais -- sao calculadas (mais de
100 linhas de formulas de mecanica dos solos) mas nunca saem da
sub-rotina. O mesmo vale para scio, srcio, olscbr e olscbrsg. Das ~10
metricas agregadas calculadas no final da sub-rotina, so tres
realmente saem: ss_strengtho, ll_strengtho, scbro (o metodo "SMSPII").

Alem disso, `ols_flag` -- que controla tanto a escrita do arquivo de
diagnostico OLS_soil_strength.dat quanto o calculo de cbr_r1/cbr_r2/
cbr_r3 (curvas CBR-vs-umidade de Rosa) e do k-NN -- e uma constante
local zerada no topo da sub-rotina e nunca reatribuida em nenhum outro
lugar do arquivo. Ou seja, `if(ols_flag==1)` nunca e verdadeiro: nem o
arquivo e aberto/escrito, nem cbr_r1/cbr_r2/cbr_r3/knn_cbr sao
calculados de fato.

Tudo isso foi traduzido fielmente mesmo assim (inclusive o arquivo de
saida OLS_soil_strength.dat e o k-NN comentado no proprio Fortran) --
ver "Pontos de atencao".

Arquitetura
-------------
`soil_strength` recebe `state: FasstState` como argumento explicito.
`EPS` e `PI` vem de `fasst.constants`; `mflag` e campo real de
`FasstState`. `dense` esta listada no cabecalho do .F90 ("uses the
functions: dense") mas, assim como `soilhumid` em soil_moisture.F90,
nunca e de fato chamada em nenhum lugar do corpo desta sub-rotina --
por isso nao e importada aqui.

Unico parametro: nstr_flag (intent(in), array) -- nunca escrito nesta
sub-rotina, lido diretamente sem copia. As tres saidas (ss_strengtho,
ll_strengtho, scbro) sao devolvidas em tupla; nao ha nada intent(inout)
para copiar.

As tabelas de coeficientes (DATA do Fortran: smsp, scourse, cbr_min,
coefci, coefrci, ols, PL, R1, R2, R3) foram transcritas com indexacao
1-based identica ao Fortran, coluna por coluna -- cada linha do DATA
original do Fortran corresponde exatamente a uma coluna da tabela
(smsp(6,16), coefci(4,14), coefrci(4,14), ols(2,20), R1/R2/R3(3,20)),
entao cada tabela em Python e uma lista de colunas (1-based, indice 0
sem uso), e cada coluna e por sua vez uma lista 1-based dos valores
daquela coluna -- TABELA[coluna][linha] == tabela(linha,coluna) do
Fortran, sem nenhuma transposicao adicional necessaria (diferente das
matrizes de coeficientes de insol() em module_radiation, que exigiam
compensar o preenchimento column-major do Fortran).

Duas convencoes de indice por tipo de solo coexistem neste arquivo,
copiadas exatamente do Fortran: R1/R2/R3/ols/PL/cbr_min/scourse sao
indexadas diretamente por ntype(i) (1=GW,2=GP,...,18=EV,19=UK,20=UF,
e cbr_min/scourse continuam ate 30 para cobrir CO/AS/RO/agua/ar/neve);
coefci/coefrci sao indexadas por ntype(i)-4 (comecam em SW=5); smsp e
indexada por ntype(i)-2 (comeca em GM=3). Nenhuma dessas foi
"corrigida" para uma convencao unica -- os offsets sao exatamente os
do original.

Pontos de atencao
-------------------
* Provavel erro de digitacao do Fortran original, preservado sem
  correcao: no ramo SMSPII com gelo (linhas 405-406 do .F90), `strci`
  usa produto (`strci*(1-ice) + ice_str*ice`), mas `rci` usa
  exponenciacao (`rci**(1-ice) + ice_str*ice`); as linhas equivalentes
  de `cisg` e `rcisg` (483 e 539) usam produto. `rci` so alimenta
  `tot_rcio`/`srcio`, que a sub-rotina nao devolve, entao o efeito nas
  tres saidas reais e nulo.
* `ols_flag` e uma constante local (sempre 0) -- todo o bloco de
  escrita do arquivo OLS_soil_strength.dat, o calculo de cbr_r1/
  cbr_r2/cbr_r3 (curvas de Rosa) e o metodo k-NN (ja comentado no
  proprio Fortran) sao inatingiveis na pratica. Traduzidos fielmente
  mesmo assim.
* `plt` nao e reatribuida no ramo `nstr_flag(i)==1` (dentro do bloco
  morto de `ols_flag==1`) -- carrega o valor do ultimo no processado
  no ramo `nstr_flag(i)==0 ou ==3` (ou 0.0, se nenhum no ainda tiver
  passado por ali). Se `ols_flag` fosse ligado e o primeiro no com
  `scourse!=0/3` e `ntype>2` e `ntype!=15` do perfil tivesse
  `nstr_flag==1`, `plt` ainda estaria em 0.0 e `t1 = (mc/plt)/...`
  faria uma divisao por zero. Em Fortran, isso silenciosamente vira
  Infinity (sem travar o programa) e a checagem seguinte
  (`abs(t1) < 50`) descarta o resultado; em Python, `mc/plt` com
  `plt==0.0` levantaria `ZeroDivisionError`. Preservado sem guarda de
  divisao por zero, por fidelidade estrita a formula; como o bloco
  inteiro e inatingivel hoje (ols_flag=0), isso nao afeta a execucao
  atual -- mas quebraria se alguem reativasse ols_flag sem tratar esse
  caso. Ver Pendencias.
* `E/sigma0` e `E/sigma0t`, nos blocos de CI e RCI (unico trecho de
  codigo ativo -- fora do bloco morto de `ols_flag` -- com uma divisao
  cujo denominador depende de dado de entrada): a condicao
  `dabs(a)>eps and E/sigma0>0.0` calcula `E/sigma0` mesmo quando
  `sigma0` pode ser 0.0 para alguma combinacao de tipo de solo/umidade.
  Protegido aqui com `_fdiv`, que devolve +-inf/nan em vez de levantar
  `ZeroDivisionError` -- mesma convencao de `sp_humid.py`/`snow.py`/
  `module_lowveg.py`/`module_canopy.py`/`module_radiation.py`.
* `close(25)` no final do Fortran fecha uma unidade de arquivo que
  nunca e aberta em nenhum lugar ativo deste arquivo -- so aparece um
  `open(unit=25,...)` comentado, junto do restante do metodo k-NN.
  Traduzido como state.unit25.close() por simetria com
  state.unit54.close().
* Das ~10 metricas finais calculadas (Steve Grant: ss_strength, sci,
  srci, scbr, olscbrsg; SMSPII: ss_strengtho, scio, srcio, scbro,
  olscbr; mais ll_strength/ll_strengtho), so ss_strengtho, ll_strengtho
  e scbro saem da sub-rotina -- ver acima. Todas as demais foram
  mantidas como variaveis locais (calculadas, nunca retornadas), fieis
  a assinatura atual do Fortran.
"""

import math

import numpy as np

from .constants import EPS, PI, SPFLAG
from .state import FasstState

# ============================================================
# Tabelas de coeficientes (DATA do Fortran, linhas 60-226)
# ============================================================

CONVERT = 0.0142230  # convert psi to cm
M_FRICTION = 1.0  # friction factor (nome Fortran: m)
SGBETA = 165.0 * 3.141592653589793 / 180.0  # angle of the cone penetrometer
ICE_STR = 5.0 * 14.2233433  # ice strength (psi)

# scourse(30): classificacao grosso(1)/fino(2)/pavimento-rocha-neve(3)/ar-agua(0)
SCOURSE = [0,
           1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 1, 0, 3,
           3, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# cbr_min(30): CBR minimo por tipo de solo/material (1..30, ver docstring)
CBR_MIN = [0.0,
           60.0, 35.0, 30.0, 15.0, 20.0, 15.0, 20.0, 10.0,   # GW,GP,GM,GC,SW,SP,SM,SC
           5.0, 5.0, 1.0, 2.0, 1.0, 0.0, 15.0, 5.0,           # ML,CL,OL,CH,MH,OH,Pt,MC
           20.0, 0.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0,       # CM,EV,XX,CO,AS,XX,XX,XX
           100.0, 0.0, 0.0, 0.0, 0.0, 1.0]                    # RO,WA,AI,XX,XX,SN

# smsp(6,16): CIa,CIb,RCIa,RCIb,CImax,RCImax -- indexada por ntype(i)-2 (GM=col.1)
_smsp_cols = [
    [8.749, -1.1949, 12.542, -2.955, 750.0, 750.0],    # GM
    [9.056, -1.3566, 12.542, -2.955, 750.0, 750.0],    # GC
    [3.987, 0.815, 3.987, 0.815, 750.0, 750.0],         # SW
    [3.987, 0.815, 3.987, 0.815, 750.0, 750.0],         # SP
    [8.749, -1.1949, 12.542, -2.955, 750.0, 750.0],    # SM
    [9.056, -1.3566, 12.542, -2.955, 750.0, 750.0],    # SC
    [10.225, -1.565, 11.936, -2.407, 300.0, 300.0],     # ML
    [10.998, -1.848, 15.506, -3.530, 300.0, 300.0],     # CL
    [10.977, -1.754, 17.399, -3.584, 300.0, 300.0],     # OL
    [13.641, -2.417, 13.686, -2.705, 300.0, 300.0],     # CH
    [12.321, -2.044, 23.641, -5.191, 300.0, 300.0],     # MH
    [13.046, -2.172, 12.189, -1.942, 300.0, 300.0],     # OH
    [0.0, 0.0, 0.0, 0.0, 15.0, 15.0],                   # Pt
    [9.056, -1.3566, 12.542, -2.955, 750.0, 750.0],    # MC (SMSC)
    [9.454, -1.385, 14.236, -3.137, 300.0, 300.0],      # CM (CLML)
    [3.987, 0.815, 3.987, 0.815, 750.0, 750.0],         # EV
]
SMSP = [None] + [[None] + col for col in _smsp_cols]

# coefci(4,14): a,b,c,lambda -- indexada por ntype(i)-4 (SW=col.1)
_coefci_cols = [
    [0.0, 0.0, 0.0, 0.0],                    # SW
    [0.0, 0.0, 0.0, 0.0],                    # SP
    [90.2617, 57.9645, -248.0, 1.1812],      # SM
    [13.862, 175.1, 3080.4, 0.4803],         # SC
    [-5.9508, 74.1596, -100.0, 0.4877],      # ML
    [42.7166, 2.935, 11.8003, 1.0988],       # CL
    [118.1, 2.4527, -506.3, 0.9892],         # OL
    [33.9469, 2.0355, -129.5, 1.3644],       # CH
    [75.6996, 1.5039, -206.9, 1.0707],       # MH
    [87.8898, 0.0117, -210.2, 1.2245],       # OH
    [0.0, 0.0, 0.0, 0.0],                    # PT
    [13.862, 175.1, 3080.4, 0.4803],         # MC - SMSC
    [-4.3898, 75.6072, -500.0, 0.4880],      # CM - CLML
    [0.0, 0.0, 0.0, 0.0],                    # EV
]
COEFCI = [None] + [[None] + col for col in _coefci_cols]

# coefrci(4,14): a,b,c,lambda -- indexada por ntype(i)-4 (SW=col.1)
_coefrci_cols = [
    [0.0, 0.0, 0.0, 0.0],                       # SW
    [0.0, 0.0, 0.0, 0.0],                       # SP
    [3.47e-3, 2.5e-4, -0.0216, 2.8795],          # SM
    [59.9976, 1.9953, -54.8093, 0.6354],         # SC
    [-2.3143, 42.0819, -100.0, 0.3557],          # ML
    [0.0, 0.0, 0.0, 0.0],                       # CL
    [0.0, 0.0, 0.0, 0.0],                       # OL
    [13.8991, 105.0, -17.6119, 0.4388],          # CH
    [0.0, 0.0, 0.0, 0.0],                       # MH
    [0.0, 0.0, 0.0, 0.0],                       # OH
    [0.0, 0.0, 0.0, 0.0],                       # PT
    [59.9976, 1.9953, -54.809, 0.6354],          # MC - SMSC
    [-8.8264, 44.0451, -100.0, 0.3553],          # CM - CLML
    [0.0, 0.0, 0.0, 0.0],                       # EV
]
COEFRCI = [None] + [[None] + col for col in _coefrci_cols]

# ols(2,20): CBR = a*CI**b -- indexada diretamente por ntype(i) (GW=col.1..UF=col.20)
_ols_cols = [
    [1.1392, 0.4896], [1.1392, 0.4896], [1.1392, 0.4896], [1.1392, 0.4896],  # GW,GP,GM,GC
    [1.1392, 0.4896], [1.1392, 0.4896], [1.1392, 0.4896], [1.1392, 0.4896],  # SW,SP,SM,SC
    [0.1111, 0.7390], [0.1266, 0.6986], [0.1305, 0.6776], [0.1264, 0.6979],  # ML,CL,OL,CH
    [0.0820, 0.7174], [0.1305, 0.6776], [0.2985, 0.5358],                    # MH,OH,Pt
    [1.1392, 0.4896], [0.1281, 0.6984], [1.1392, 0.4896],                    # MC,CM,EV
    [1.1392, 0.4896], [0.1305, 0.6776],                                     # UK,UF
]
OLS = [None] + [[None] + col for col in _ols_cols]

# PL(20): limite de plasticidade -- indexada diretamente por ntype(i)
PL = [None,
      999.0, 999.0, 27.0, 14.0, 999.0, 999.0, 27.5, 16.0, 27.0,   # GW,GP,GM,GC,SW,SP,SM,SC,ML
      19.5, 42.0, 24.0, 42.5, 71.0, 166.0, 18.5, 17.0, 999.0,      # CL,OL,CH,MH,OH,Pt,MC,CM,EV
      27.5, 24.0]                                                  # UK,UF

# R1(3,20): CBR = a + b*MC**c -- indexada diretamente por ntype(i)
_r1_cols = [
    [999.0, 999.0, 999.0],       # GW
    [7.7814, -0.04242, 1.649],    # GP
    [999.0, 999.0, 999.0],       # GM
    [999.0, 999.0, 999.0],       # GC
    [999.0, 999.0, 999.0],       # SW
    [999.0, 999.0, 999.0],       # SP
    [7.7814, -0.042428, 1.649],   # SM
    [999.0, 999.0, 999.0],       # SC
    [-1.815, 5480.0, -2.336],     # ML
    [-5.671, 135.0, -0.9031],     # CL
    [999.0, 999.0, 999.0],       # OL
    [-4.52, 77.7, -0.66],         # CH
    [-3.184, 733.0, -1.287],      # MH
    [999.0, 999.0, 999.0],       # OH
    [999.0, 999.0, 999.0],       # Pt
    [7.7814, -0.042428, 1.649],   # MC
    [-6.09, 139.08, -0.8953],     # CM
    [999.0, 999.0, 999.0],       # EV
    [-119.625, 130.8733, -0.01994],  # UK
    [0.5062, 405.96e9, -1.591],    # UF
]
R1 = [None] + [[None] + col for col in _r1_cols]

# R2(3,20): CBR = a + b*exp(-(MC/PL)/c) -- indexada diretamente por ntype(i)
_r2_cols = [
    [999.0, 999.0, 999.0],     # GW
    [999.0, 999.0, 999.0],     # GP
    [999.0, 999.0, 999.0],     # GM
    [999.0, 999.0, 999.0],     # GC
    [999.0, 999.0, 999.0],     # SW
    [999.0, 999.0, 999.0],     # SP
    [-0.418, 3771.0, 0.1129],   # SM
    [999.0, 999.0, 999.0],     # SC
    [-0.169, 165.5, 0.2117],    # ML
    [-0.5509, 52.32, 0.3828],   # CL
    [999.0, 999.0, 999.0],     # OL
    [-0.5957, 26.97, 0.6256],   # CH
    [-2.356, 24.86, 0.5914],    # MH
    [999.0, 999.0, 999.0],     # OH
    [999.0, 999.0, 999.0],     # Pt
    [-0.418, 3771.0, 0.1129],   # MC
    [0.0, 72.12, 0.3076],       # CM
    [999.0, 999.0, 999.0],     # EV
    [0.02801, 32.46, 0.426],    # UK
    [-0.03409, 31.94, 0.5042],  # UF
]
R2 = [None] + [[None] + col for col in _r2_cols]

# R3(3,20): CBR = 10**(a + b*exp(-MC/c)) -- indexada diretamente por ntype(i)
_r3_cols = [
    [999.0, 999.0, 999.0],                    # GW
    [999.0, 999.0, 999.0],                    # GP
    [999.0, 999.0, 999.0],                    # GM
    [999.0, 999.0, 999.0],                    # GC
    [999.0, 999.0, 999.0],                    # SW
    [999.0, 999.0, 999.0],                    # SP
    [-1008791.0, 1008792.22, 17866100.0],      # SM
    [999.0, 999.0, 999.0],                    # SC
    [-1.273, 8.091, 12.22],                    # ML
    [-3.257, 5.495, 47.16],                    # CL
    [999.0, 999.0, 999.0],                    # OL
    [-0.8406, 3.108, 34.82],                   # CH
    [-1.955, 5.031, 50.8],                     # MH
    [999.0, 999.0, 999.0],                    # OH
    [999.0, 999.0, 999.0],                    # Pt
    [-1008791.0, 1008792.22, 17866100.0],      # MC
    [-2.453, 5.037, 33.53],                    # CM
    [999.0, 999.0, 999.0],                    # EV
    [-1.779, 2.628, 101.9],                    # UK
    [-0.9324, 2.033, 49.99],                   # UF
]
R3 = [None] + [[None] + col for col in _r3_cols]


def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _fdiv(x, y):
    """Divisao IEEE-754 pura: ``x / 0.0`` -> +-inf / nan, sem levantar.

    Mesma funcao de `sp_humid.py`/`snow.py`/`module_lowveg.py`/
    `module_canopy.py`/`module_radiation.py`.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


def soil_strength(state: FasstState, nstr_flag):
    """
    Traducao fiel de `subroutine soil_strength(nstr_flag,ss_strengtho,
    ll_strengtho,scbro)` -- ver docstring do modulo para as saidas
    calculadas mas nunca devolvidas pela sub-rotina (dado ols_flag=0
    fixo e a assinatura reduzida).
    """
    nnodes = state.nnodes

    # ---------------- zero-out variables (linhas 236-299) ----------------
    ols_flag = 0
    scount = 0
    lcount = 0
    strci = 0.0
    rci = 0.0
    mc = 0.0
    tot_ci = 0.0
    tot_rci = 0.0
    sstrengtho = 0.0
    tot_strengtho = 0.0
    tot_cio = 0.0
    tot_rcio = 0.0
    cbr = 0.0
    lstrengtho = 0.0
    tot_cbro = 0.0
    ss_strengtho = 0.0
    ll_strengtho = 0.0
    E = 0.0
    D = 0.0
    cisg = 0.0
    rcisg = 0.0
    sigma0 = 0.0
    lambda_ = 0.0
    a = 0.0
    b = 0.0
    c = 0.0
    smt = 0.0
    pht = 0.0
    w = 0.0
    sigma0t = 0.0
    cisgt = 0.0
    rcisgt = 0.0
    cisgmax = 0.0
    rcisgmax = 0.0
    cbrsg = 0.0
    sstrength = 0.0
    ss_strength = 0.0
    tot_strength = 0.0
    sci = 0.0
    srci = 0.0
    lstrength = 0.0
    ll_strength = 0.0
    scbr = 0.0
    tot_cbr = 0.0
    ols_cbr = 0.0
    tot_olscbr = 0.0
    olscbr = 0.0
    ols_cbrsg = 0.0
    tot_olscbrsg = 0.0
    olscbrsg = 0.0
    plt = 0.0
    cbr_r1 = 0.0
    cbr_r2 = 0.0
    cbr_r3 = 0.0
    t1 = 0.0
    t2 = 0.0

    scio = 0.0
    srcio = 0.0
    scbro = 0.0

    # st(11,nnodes): valores por no, para o arquivo OLS (so usado se ols_flag==1)
    st = [None] + [_new_array(nnodes) for _ in range(11)]

    # ---------------- open the output file; write headers (linhas 305-358) ----------------
    # NUNCA executa com ols_flag fixo em 0 -- ver docstring do modulo.
    if state.iw == state.istart and ols_flag == 1:
        state.unit54 = open('OLS_soil_strength.dat', 'w')

        state.unit54.write("USCS = USCS soil type\n")
        state.unit54.write("Z(m) = depth in meters\n")
        state.unit54.write("THETA = actaul volumetric water content\n")
        state.unit54.write("WC = 100 * THETA/porosity (vol/vol)\n")
        state.unit54.write("MC = 100 * THETA/dry density (mass/mass)\n")
        state.unit54.write("SMSP_CI = cone index based on SMSPII curves\n")
        state.unit54.write("SG_CI = cone index based on Steve Grant curves\n")
        state.unit54.write("CBR_CI = Debrah CBR vs CI curves using SMSPII CI values\n")
        state.unit54.write("CBR_CI_SG = Debrah CBR vs CI curves using Steve Grant CI values\n")
        state.unit54.write("CBR_MC_R1 = Rosa CBR vs MC curve 1\n")
        state.unit54.write("CBR_MC_R2 = Rosa CBR vs MC curve 2\n")
        state.unit54.write("CBR_MC_R3 = Rosa CBR vs MC curve 3\n")
        state.unit54.write("CBR_KNN = k-nearest neighbor CBR\n")
        state.unit54.write("CBR_JRAC = JRAC CBR using SMSPII CI curves\n")
        state.unit54.write("CBR_JRACSG = JRAC CBR using Steve Grant CI curves\n")
        state.unit54.write("Dry Density (g/cm^3)\n")
        state.unit54.write(
            "NOTE: No SMSPII CI vs THETA curves exist for gravels or peat. "
            "Therefore,\n      a default value of 150.0 is used for gravels, "
            "15.0 for unfrozen peat\n      and 300.0 for frozen peat.\n"
        )
        state.unit54.write(" \n")
        state.unit54.write(
            "{:10d} {:6d} {:3d} {:10.6f} {:11.6f} {:11.6f} {:8d} {:8d} "
            "{:6.2f} {:5.2f} {:8.2f}\n".format(
                state.freq_id, (state.iend - state.istart + 1), 14, state.lat, state.mlong,
                state.elev, state.vitd_index, state.met_count, state.timeoffset,
                state.timstep, state.mflag,
            )
        )
        state.unit54.write(
            " USCS     Z(m)      THETA      WC      MC      "
            " SMSP_CI        SG_CI       CBR_CI    CBR_CI_SG    "
            "CBR_MC_R1    CBR_MC_R2    CBR_MC_R3      CBR_KNN      "
            "CBR_JRAC    CBR_JRACSG     DENSITY\n"
        )
        # bloco do k-NN (abertura da unidade 25, arvore k-d) comentado no
        # proprio Fortran original -- nao traduzido, ver docstring do modulo.

    # ---------------- laco principal por no (linhas 360-694) ----------------
    for i in range(1, nnodes + 1):
        if SCOURSE[state.ntype[i]] == 3:  # pavements, rock, snow
            strci = 300.0
            rci = 300.0
            cbr = CBR_MIN[state.ntype[i]]
            ols_cbr = cbr
            ols_cbrsg = cbr
        elif SCOURSE[state.ntype[i]] == 0:  # air, water
            strci = 0.0
            rci = 0.0
            cbr = CBR_MIN[state.ntype[i]]
            ols_cbr = cbr
            ols_cbrsg = cbr
        elif state.ntype[i] <= 2:  # GP, GW
            strci = 150.0
            rci = 150.0
            cbr = CBR_MIN[state.ntype[i]]
            ols_cbr = OLS[state.ntype[i]][1] * strci ** OLS[state.ntype[i]][2]
            ols_cbrsg = ols_cbr
        elif state.ntype[i] == 15:  # peat
            strci = SMSP[state.ntype[i] - 2][5]
            rci = SMSP[state.ntype[i] - 2][6]
            if state.ice[i] > 0.0:  # frozen peat
                strci = 300.0
                rci = 300.0
            cbr = CBR_MIN[state.ntype[i]]
            ols_cbr = cbr
            ols_cbrsg = cbr
        else:  # all other soils
            # SMSPII method modified for frozen soil
            mc = 100.0 * state.soil_moist[i] / state.nsoilp[i][1]  # % by weight
            if mc <= 0.0:
                mc = 1e-3  # saftey net!!!!!!!!

            t1 = SMSP[state.ntype[i] - 2][1] + SMSP[state.ntype[i] - 2][2] * math.log(mc)
            if abs(t1) < 50.0:
                strci = math.exp(t1)
            if state.ntype[i] <= 4:
                strci = strci * 1.1
            if strci > SMSP[state.ntype[i] - 2][5]:
                strci = SMSP[state.ntype[i] - 2][5]

            t1 = SMSP[state.ntype[i] - 2][3] + SMSP[state.ntype[i] - 2][4] * math.log(mc)
            if abs(t1) < 50.0:
                rci = math.exp(t1)
            if state.ntype[i] <= 4:
                rci = rci * 1.1
            if rci > SMSP[state.ntype[i] - 2][6]:
                rci = SMSP[state.ntype[i] - 2][6]

            if state.ice[i] > EPS:
                strci = min(300.0, strci * (1.0 - state.ice[i]) + ICE_STR * state.ice[i])
                rci = min(300.0, rci ** (1.0 - state.ice[i]) + ICE_STR * state.ice[i])

            cbr = 2e-5 * strci * strci + 6e-3 * strci + 0.129
            if cbr > 100.0:
                cbr = 100.0
            cbr = (cbr + CBR_MIN[state.ntype[i]]) * 0.5

            if nstr_flag[i] == 0:
                ols_cbr = OLS[state.ntype[i]][1] * (strci ** OLS[state.ntype[i]][2])
            elif nstr_flag[i] == 1:  # unknown, coarse-grained
                ols_cbr = OLS[19][1] * (strci ** OLS[19][2])
            else:  # unknown, fine-grained
                ols_cbr = OLS[20][1] * (strci ** OLS[20][2])

            # Steve Grant, George Mason Method
            # need another head and soil moisture value for this method
            smt = 1.01 * state.soil_moist[i]
            dct = math.cos(PI - SGBETA) / math.sin(PI - SGBETA)

            if smt >= state.nsoilp[i][9]:  # saturated
                pht = 0.0
            elif abs(smt - state.nsoilp[i][8]) < 1e-3 or smt <= state.nsoilp[i][15]:
                pht = state.pheadmin[i]  # totally unsaturated (cm)
            else:  # general conditions
                w = (smt - state.nsoilp[i][8]) / (state.nsoilp[i][9] - state.nsoilp[i][8])  # unitless
                if w > EPS and abs(w ** (-1.0 / state.nsoilp[i][12]) - 1.0) > EPS:
                    pht = (w ** (-1.0 / state.nsoilp[i][12]) - 1.0) ** (1.0 / state.nsoilp[i][11])
                pht = -(1.0 / state.nsoilp[i][10]) * pht  # cm

            # CI
            if state.ntype[i] > 4:
                lambda_ = COEFCI[state.ntype[i] - 4][4]
                E = -(math.sin(lambda_) * math.cos(SGBETA)
                      + M_FRICTION * math.cos(lambda_) * math.sin(SGBETA)) / (
                    math.cos(lambda_) - math.cos(SGBETA)
                )
                D = (math.sin(lambda_) + M_FRICTION * math.sin(SGBETA)) / (
                    math.cos(lambda_) - math.cos(SGBETA)
                )

                a = COEFCI[state.ntype[i] - 4][1]
                if state.soil_moist[i] < state.nsoilp[i][9]:
                    c = 0.0
                    b = COEFCI[state.ntype[i] - 4][2]
                else:
                    c = COEFCI[state.ntype[i] - 4][3]
                    b = 0.0

                sigma0 = a + b * (abs(state.phead[i]) * CONVERT) + c * (
                    state.soil_moist[i] - state.nsoilp[i][9]
                )
                sigma0t = a + b * (abs(pht) * CONVERT) + c * (smt - state.nsoilp[i][9])

                if abs(a) > EPS and _fdiv(E, sigma0) > 0.0:
                    cisg = (sigma0 / math.sqrt(3.0)) * (
                        3.0 * PI - 2.0 * SGBETA + math.asin(M_FRICTION)
                        + math.log(E / (math.sqrt(3.0) * sigma0))
                        + D * 0.5 + M_FRICTION * dct - math.sqrt(1.0 - M_FRICTION * M_FRICTION)
                    )
                    cisg = max(0.0, min(999.0, cisg))
                else:
                    cisg = 0.0

                if abs(a) > EPS and _fdiv(E, sigma0t) > 0.0:
                    cisgt = (sigma0t / math.sqrt(3.0)) * (
                        3.0 * PI - 2.0 * SGBETA + math.asin(M_FRICTION)
                        + math.log(E / (math.sqrt(3.0) * sigma0t))
                        + D * 0.5 + M_FRICTION * dct - math.sqrt(1.0 - M_FRICTION * M_FRICTION)
                    )
                    cisgt = max(0.0, min(999.0, cisgt))
                else:
                    cisgt = 0.0

                # cisgmax = SMSP[state.ntype[i]-2][5]  # temporario -- comentado no Fortran
                # if cisgt-cisg>0 and cisg<0.8*cisgmax: cisg = max(0,min(999,0.8*cisgmax))  # comentado

                if state.ice[i] > EPS:
                    cisg = min(300.0, cisg * (1.0 - state.ice[i]) + ICE_STR * state.ice[i])

                cbrsg = 2e-5 * cisg * cisg + 6e-3 * cisg + 0.129
                if cbrsg > 100.0:
                    cbrsg = 100.0
                cbrsg = (cbrsg + CBR_MIN[state.ntype[i]]) * 0.5

                if nstr_flag[i] == 0:
                    if cisg > EPS:
                        ols_cbrsg = OLS[state.ntype[i]][1] * cisg ** OLS[state.ntype[i]][2]
                else:
                    if cisg > EPS:
                        ols_cbr = OLS[19][1] * cisg ** OLS[19][2]

                # RCI
                lambda_ = COEFRCI[state.ntype[i] - 4][4]
                E = -(math.sin(lambda_) * math.cos(SGBETA)
                      + M_FRICTION * math.cos(lambda_) * math.sin(SGBETA)) / (
                    math.cos(lambda_) - math.cos(SGBETA)
                )
                D = (math.sin(lambda_) + M_FRICTION * math.sin(SGBETA)) / (
                    math.cos(lambda_) - math.cos(SGBETA)
                )

                a = COEFRCI[state.ntype[i] - 4][1]
                b = COEFRCI[state.ntype[i] - 4][2]
                if state.soil_moist[i] < state.nsoilp[i][9]:
                    c = 0.0
                else:
                    c = COEFRCI[state.ntype[i] - 4][3]

                sigma0 = a + b * (-state.phead[i] * CONVERT) + c * (
                    state.soil_moist[i] - state.nsoilp[i][9]
                )
                sigma0t = a + b * (-pht * CONVERT) + c * (smt - state.nsoilp[i][9])

                if abs(a) > EPS and _fdiv(E, sigma0) > 0.0:
                    rcisg = (sigma0 / math.sqrt(3.0)) * (
                        3.0 * PI - 2.0 * SGBETA + math.asin(M_FRICTION)
                        + math.log(E / (math.sqrt(3.0) * sigma0))
                        + D * 0.5 + M_FRICTION * dct - math.sqrt(1.0 - M_FRICTION * M_FRICTION)
                    )
                    rcisg = max(0.0, min(999.0, rcisg))
                else:
                    rcisg = 0.0

                if abs(a) > EPS and _fdiv(E, sigma0t) > 0.0:
                    rcisgt = (sigma0t / math.sqrt(3.0)) * (
                        3.0 * PI - 2.0 * SGBETA + math.asin(M_FRICTION)
                        + math.log(E / (math.sqrt(3.0) * sigma0t))
                        + D * 0.5 + M_FRICTION * dct - math.sqrt(1.0 - M_FRICTION * M_FRICTION)
                    )
                    rcisgt = max(0.0, min(999.0, rcisgt))
                else:
                    rcisgt = 0.0

                rcisgmax = SMSP[state.ntype[i] - 2][6]  # temporary
                if rcisgt - rcisg > 0.0 and rcisg < 0.8 * rcisgmax:
                    rcisg = max(0.0, min(999.0, 0.8 * rcisgmax))

                if state.ice[i] > EPS:
                    rcisg = min(300.0, rcisg * (1.0 - state.ice[i]) + ICE_STR * state.ice[i])
            # fim do bloco Steve Grant (ntype(i) > 4)

        cisg = max(0.0, cisg)
        strci = max(0.0, strci)
        rcisg = max(0.0, rcisg)
        rci = max(0.0, rci)
        cbr = max(0.0, cbr)

        # determine strength based on fine/coarse and depth
        if SCOURSE[state.ntype[i]] == 1:  # coarse
            sstrength = cisg
            sstrengtho = strci
        elif SCOURSE[state.ntype[i]] == 2:  # fine
            sstrength = rcisg
            sstrengtho = rci
        else:  # peat, pavements, air, snow
            sstrength = strci
            sstrengtho = strci

        if state.elev - state.nz[i] <= 0.1524:  # 6"
            scount += 1
            tot_strength += sstrength
            tot_ci += cisg
            tot_rci += rcisg
            tot_cbr += cbrsg
            tot_strengtho += sstrengtho
            tot_cio += strci
            tot_rcio += rci
            tot_cbro += cbr
            tot_olscbr += ols_cbr
            tot_olscbrsg += ols_cbrsg
        elif 0.1524 < state.elev - state.nz[i] <= 0.3048:  # 6" - 12"
            lstrength += sstrength
            lstrengtho += sstrengtho
            lcount += 1
        elif state.elev - state.nz[i] > 0.3048 and lcount == 0:
            lstrength += sstrength
            lstrengtho += sstrengtho
            lcount += 1

        # ---- bloco OLS (curvas de Rosa CBR-vs-umidade + k-NN) ----
        # NUNCA executa com ols_flag fixo em 0 -- ver docstring do modulo.
        if ols_flag == 1:
            cbr_r1 = SPFLAG
            cbr_r2 = SPFLAG
            cbr_r3 = SPFLAG

            if nstr_flag[i] == 0 or nstr_flag[i] == 3:
                if abs(R1[state.ntype[i]][1] - SPFLAG) > EPS:
                    if mc > EPS:
                        cbr_r1 = R1[state.ntype[i]][1] + R1[state.ntype[i]][2] * (
                            mc ** R1[state.ntype[i]][3]
                        )
                else:
                    if mc > EPS:
                        cbr_r1 = R1[19][1] + R1[19][2] * (mc ** R1[19][3])

                if abs(state.nsoilp[i][22] - SPFLAG) > EPS:
                    plt = state.nsoilp[i][22]
                else:
                    plt = PL[state.ntype[i]]

                t1 = (mc / plt) / R2[state.ntype[i]][3]
                t2 = mc / R3[state.ntype[i]][3]
                if abs(plt - SPFLAG) > EPS:
                    if abs(t1) < 50.0:
                        cbr_r2 = R2[state.ntype[i]][1] + R2[state.ntype[i]][2] * math.exp(-t1)

                    if abs(t2) < 50.0:
                        cbr_r3 = 10.0 ** (R3[state.ntype[i]][1] + R3[state.ntype[i]][2] * math.exp(-t2))
            elif nstr_flag[i] == 1:  # unknown, coarse-grained
                if mc > EPS:
                    cbr_r1 = R1[19][1] + R1[19][2] * (mc ** R1[19][3])

                # `plt` nao e redefinida neste ramo -- carrega o valor do
                # ultimo no processado no ramo nstr_flag==0/3 (ou 0.0). Se
                # for 0.0 aqui, a divisao abaixo levanta ZeroDivisionError
                # em Python (Fortran produziria Infinity silenciosamente) --
                # ver docstring do modulo, "Pontos de atencao". Bloco
                # inatingivel hoje (ols_flag=0), preservado por fidelidade.
                t1 = (mc / plt) / R2[19][3]
                t2 = (mc / PL[19]) / R2[19][3]
                if abs(state.nsoilp[i][22] - SPFLAG) > EPS:
                    if abs(t1) < 50.0:
                        cbr_r2 = R2[19][1] + R2[19][2] * math.exp(-t1)
                else:
                    if abs(t1) < 50.0:
                        cbr_r2 = R2[19][1] + R2[19][2] * math.exp(-t2)

                t1 = mc / R3[19][3]
                if abs(t1) < 50.0:
                    cbr_r3 = 10.0 ** (R3[19][1] + R3[19][2] * math.exp(-t1))
            elif nstr_flag[i] == 2:  # unknown, fine-grained
                if mc > EPS:
                    cbr_r1 = R1[20][1] + R1[20][2] * (mc ** R1[20][3])

                t1 = (mc / plt) / R2[20][3]
                t2 = (mc / PL[20]) / R2[20][3]
                if abs(state.nsoilp[i][22] - SPFLAG) > EPS:
                    if abs(t1) < 50.0:
                        cbr_r2 = R2[20][1] + R2[20][2] * math.exp(-t1)
                else:
                    if abs(t2) < 50.0:
                        cbr_r2 = R2[20][1] + R2[20][2] * math.exp(-t2)

                t1 = mc / R3[20][3]
                if abs(t1) < 50.0:
                    cbr_r3 = 10.0 ** (R3[20][1] + R3[20][2] * math.exp(-t1))

            # k-nn method: inteiramente comentado no proprio Fortran
            # original (abertura de arquivo, arvore k-d, busca dos k
            # vizinhos mais proximos) -- nao traduzido; apenas o valor
            # sentinela que sobra dessa desativacao:
            knn_cbr = SPFLAG

            # write OLS info to temporary array to reverse output
            st[1][i] = mc
            st[2][i] = strci
            st[3][i] = cisg
            st[4][i] = ols_cbr
            st[5][i] = ols_cbrsg
            st[6][i] = cbr_r1
            st[7][i] = cbr_r2
            st[8][i] = cbr_r3
            st[9][i] = knn_cbr
            st[10][i] = cbr
            st[11][i] = cbrsg
    # fim do laco principal por no

    # ---------------- output the data (linhas 696-711) ----------------
    # NUNCA executa com ols_flag fixo em 0 -- ver docstring do modulo.
    if ols_flag == 1:
        for j in range(1, nnodes + 1):
            i = nnodes - j + 1
            state.unit54.write(
                "  {:2s}  {:10.4f} {:8.4f} {:8.2f}{:8.2f} "
                "{:12.2f} {:12.2f} {:12.2f} {:12.2f} {:12.2f} "
                "{:12.2f} {:12.2f} {:12.2f} {:12.2f} {:12.2f} "
                "{:12.3f}\n".format(
                    state.node_type[i], state.nz[i], state.soil_moist[i],
                    100.0 * state.soil_moist[i] / state.nsoilp[i][2],
                    st[1][i], st[2][i], st[3][i], st[4][i], st[5][i],
                    st[6][i], st[7][i], st[8][i], st[9][i], st[10][i],
                    st[11][i], state.nsoilp[i][1],
                )
            )

        state.unit54.write(" \n")

    # ---------------- Steve Grant, George Mason method (linhas 713-733) ----------------
    # calculado por fidelidade, mas NUNCA retornado -- ver docstring do modulo.
    if scount != 0:
        f1 = 1.0 / float(scount)
        ss_strength = tot_strength * f1
        sci = tot_ci * f1
        srci = tot_rci * f1
        scbr = tot_cbr * f1
        olscbrsg = tot_olscbrsg * f1
    else:
        ss_strength = state.mflag
        sci = state.mflag
        srci = state.mflag
        scbr = state.mflag
        olscbrsg = state.mflag

    if lcount != 0:
        ll_strength = lstrength / float(lcount)
    else:
        ll_strength = state.mflag

    # ---------------- SMSPII method (linhas 735-754) -- esta e a que sai da sub-rotina ----------------
    if scount != 0:
        ss_strengtho = tot_strengtho / float(scount)
        scio = tot_cio / float(scount)
        srcio = tot_rcio / float(scount)
        scbro = tot_cbro / float(scount)
        olscbr = tot_olscbr / float(scount)
    else:
        ss_strengtho = state.mflag
        scio = state.mflag
        srcio = state.mflag
        scbro = state.mflag
        olscbr = state.mflag

    if lcount != 0:
        ll_strengtho = lstrengtho / float(lcount)
    else:
        ll_strengtho = state.mflag

    # ---------------- close output files (linhas 756-760) ----------------
    if state.iw == state.iend and state.single_multi_flag == 0:
        # state.unit25 nunca e de fato aberta nesta versao (ver "Pontos de
        # atencao") -- fechada aqui apenas por simetria com o Fortran.
        if hasattr(state, "unit25"):
            state.unit25.close()
        if hasattr(state, "unit54"):
            state.unit54.close()

    return ss_strengtho, ll_strengtho, scbro
