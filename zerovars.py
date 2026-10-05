"""
zerovars.py -- traducao de module_zerovars.F90

Fonte Fortran: module_zerovars.F90 (295 linhas). Duas sub-rotinas publicas,
ambas "zeradoras" de estado de modulo de `fasst_global`:

  * `zero_parameters` (linhas 5-218) -- zera/inicializa o estado de TAMANHO
    FIXO (perfil de solo, neve, vegetacao, variaveis de calculo) e devolve
    os 11 `intent(out)` do driver (flags, parametros de agua, `sname`).
  * `zero_mxl_params` (linhas 221-293) -- zera os arrays ALOCAVEIS
    dimensionados por `maxlines` (um elemento por linha do arquivo met).

Assinaturas originais:
    subroutine zero_parameters(tflag,nt0,nm0,stflag,lcount,wtype,wvel,&
                               wdepth,ttest,mtest,sname)
    integer(ip),intent(out):: tflag,nt0,nm0,stflag,lcount,wtype
    real(dp),intent(out):: wvel,wdepth
    character(len=1),intent(out):: ttest,mtest
    character(len=4),intent(out):: sname(maxl)

    subroutine zero_mxl_params            (sem argumentos)

Arquitetura
-------------
Mesma convencao de `icethick` / `soil_moisture`: o estado de modulo e o
`state: FasstState` explicito, mutado no lugar; os escalares `intent(out)`
sao DEVOLVIDOS. `zero_parameters(state)` devolve a `ZeroParamsOut`
(NamedTuple) com os 11 argumentos de saida; `zero_mxl_params(state)` nao
devolve nada.

Isto NAO e o mesmo que `FasstState.reset()` (state.py): `reset()` e uma
limpeza GERAL ("instancia nova", mais sentinelas), ao passo que aqui cada
escrita do Fortran e reproduzida 1-para-1, incluindo o que o Fortran NAO
zera (pontos 1 e 2) -- e o que NAO zera e justamente o que importa para
quem quer saber o que vaza entre pontos no modo multi-ponto.

Fan-out ZERO: "no subroutines called" (linha 12) -- so atribuicoes e
lacos. `spflag`, `maxl`, `maxn`, `maxp`, `extran`, `nclayers`, `maxcol`
sao `parameter` de fasst_global e vieram para `fasst.constants`.

Fan-in: `fasst_driver.F90` -- `zero_parameters` nas linhas 269 (ponto unico,
logo DEPOIS da leitura de `gwl`/`mstflag` na linha 266 -- ver ponto 3) e
1289 (multi-ponto, antes de `read_complex`); `zero_mxl_params` nas linhas
1167 (ponto unico) e 1292 (multi-ponto), ambas depois do `allocate`
(linhas 1126-1148). O driver ainda nao foi transcrito, entao nenhum call
site Python existe.

Pontos que NAO devem ser "limpos" (numeracao = doc de transcricao)
--------------------------------------------------------------------
 1. `zero_parameters` NAO e uma limpeza total: so toca os campos que
    lista. Ficam de fora, p.ex., `gwl`, `mstflag`, `mflag`, `timstep`,
    `iw`, `maxlines`, `lat`, `elev`, `laif`, `dzveg`, `storcl/storcs`,
    `delzsi`, `stll/stls`, `ilail` e todos os alocaveis (inventario
    completo e medido no doc de transcricao). Preservado: NAO zerar mais
    do que o Fortran zera.
 2. LACUNAS de dimensao: `sn_stat(15)` so tem os elementos 1..13 zerados
    (linhas 74-76) e `dmet1(maxlines,13)` so tem as colunas 1..10 zeradas
    (linhas 282-284). `sn_stat(14:15)` e `dmet1(:,11:13)` conservam o valor
    ANTERIOR (= vazam entre pontos). `FasstState.reset()` divergiu de
    proposito (zera `sn_stat` inteiro); aqui a lacuna e PRESERVADA.
 3. `zero_parameters` NAO toca `gwl` nem `mstflag`, que o driver le
    imediatamente antes (fasst_driver.F90:266). Se alguem "completar" a
    zeragem, esses dois valores lidos do arquivo de entrada seriam
    apagados.
 4. Escritas duplas, ordem importa: `iswe = 0d0` e depois `iswe = spflag`
    (linhas 72-73); `rho_fac(i) = 0d0` e depois `= 1d0` (194-195);
    `albf = 1d0 - fola` LE `fola` ja zerado (linha 92), logo `albf == 1.0`
    e nao "albedo da folhagem". Reproduzidas na mesma ordem.
 5. Sentinelas nao-nulos: `freq_id = 1`; `iswe, ftemp, isigfl, iepf, ifola,
    ihfol = spflag` (999.0); `isoilp/soilp(maxl,maxp) = spflag`;
    `ifoliage_type, ilai, iclump, irho, itau, ialp, ieps, idzveg
    (nclayers) = spflag`; `wtype = int(spflag)` = 999, `wdepth = wvel =
    spflag`. Nao sao zeros: quem tratar "zerado" como "== 0" erra.
 6. Strings de comprimento fixo: Fortran zera com brancos (`'  '`, `' '`,
    `'    '`); aqui viram `""` (convencao de state.py -- comparar sempre com
    `.strip()`). `ttest`/`mtest = 'n'` sao devolvidos como `"n"`.
 7. Arrays 1-based com slot fantasma 0 (convencao de state.py): todo laco
    `do i=1,n` vira fatia `[1:n+1]` e o slot 0 NUNCA e escrito. `avect`
    (0:4, nclayers) e zerado inteiro (i=0..4, j=1..nclayers cobre o array
    todo; ver o ponto 8 do doc de fasst_global sobre o layout de `avect`).
 8. `zero_mxl_params` limita os lacos por `state.maxlines` (campo de
    estado), NAO pelo shape dos arrays: array alocado maior que `maxlines`
    conserva as linhas excedentes. `airt` e `(2, maxlines)` e
    `canopy_temp` e `(nclayers, maxlines)` (fasst_driver.F90:1128-1129),
    ou seja, em Python `airt[1:3, 1:maxlines+1]`.
 9. `zero_mxl_params` NAO zera `slushy` (alocavel inteiro), apesar de ser
    alocado junto com os demais (fasst_driver.F90:1136). Em Fortran o
    conteudo e indefinido apos o `allocate`; em Python o campo nem e
    tocado. Preservado; quem aloca `slushy` decide o valor inicial.
10. Alocaveis ainda nao alocados (`None`): `zero_mxl_params` levanta
    `TypeError` (em Fortran seria comportamento indefinido). Chame-a
    depois de alocar, como o driver faz.

Ver o doc de transcricao para o restante (PENDENCIAS, VALIDACAO).
"""

from typing import NamedTuple

import numpy as np

from .constants import EXTRAN, MAXCOL, MAXL, MAXN, MAXP, NCLAYERS, SPFLAG
from .state import FasstState, _chars


class ZeroParamsOut(NamedTuple):
    """Os 11 argumentos `intent(out)` de `zero_parameters`."""

    tflag: int
    nt0: int
    nm0: int
    stflag: int
    lcount: int
    wtype: int
    wvel: float
    wdepth: float
    ttest: str
    mtest: str
    sname: np.ndarray  # character(len=4) sname(maxl), slot fantasma 0


# module_zerovars.F90:223-229 -- `use fasst_global,only:` dos arrays
# alocaveis 1-D de float64 zerados no laco `do i=1,maxlines`.
_MXL_REAL_1D = (
    "sdens", "surfci", "surfrci", "surfcbr", "surfice", "surfmoist",
    "surficep", "surfmoistp", "surfemis", "surfemisf", "surfemisc", "surfd",
    "ft", "tt", "frthick", "twthick",
    "cheat", "cheat1", "sdown", "sup", "irdown", "irup", "pheat1", "lheat",
    "evap_heat", "sheat", "melt", "tmelt", "lhes", "lheatf", "cevap",
    "ponding", "tot_moist", "tot_thick",
)


def zero_parameters(state: FasstState) -> ZeroParamsOut:
    """Zera o estado de tamanho fixo de `state` (so o que o Fortran zera).

    Devolve os 11 `intent(out)` (flags, parametros de agua, `sname`).
    """
    n1 = MAXN + 1
    nx1 = MAXN + EXTRAN + 1

    state.freq_id = 1                                                  # complex id

    # ---------------- various flags (linhas 28-40) ----------------
    tflag = 0
    nt0 = 0
    nm0 = 0
    stflag = 0
    lcount = 0
    state.ecount = 0
    state.rough = 0.0
    ttest = "n"
    mtest = "n"
    sname = _chars(MAXL, 4)                                            # '    ' -> "" (ponto 6)

    # ---------------- initial soil conditions (linhas 42-55) ----------------
    state.nlayers = 0
    state.nnodes = 0
    state.ntot = 0
    state.ntemp = 0
    state.icase = 0
    state.sgralbedo = 0.0
    state.sgremis = 0.0
    state.slope_fasst = 0.0
    state.aspect = 0.0
    state.albedo_fasst = 0.0
    state.emis = 0.0
    state.isurfoldg = 0.0
    state.isurfoldf = 0.0

    # ---------------- snow & ice parameters (linhas 57-76) ----------------
    state.sn_istat[1] = 0
    state.sn_istat[2] = 0
    state.hsaccum = 0.0
    state.hi = 0.0
    state.dsnow = 0.0
    state.refreeze = 0.0
    state.newsd = 0.0
    state.atopf = 0.0
    state.km[1:EXTRAN + 1] = 0.0
    state.sphm[1:EXTRAN + 1] = 0.0
    state.hm = 0.0
    state.refreezei = 0.0
    state.iswe = 0.0
    state.iswe = SPFLAG                                                # ponto 4: escrita dupla
    state.sn_stat[1:14] = 0.0                                          # ponto 2: so 1..13 de 15

    # ---------------- vegetation parameters (linhas 78-113) ----------------
    state.isigfl = SPFLAG                                              # initial low vegetation density
    state.iepf = SPFLAG                                                # initial low veg emissivity
    state.ifola = SPFLAG                                               # initial low veg absorptivity
    state.ihfol = SPFLAG                                               # initial low veg height (cm)
    state.vegl_type = 0
    state.vegh_type = 0
    state.veg_flagl = 0                                                # low vegetation flag (0 = none; 1 = yes)
    state.veg_flagh = 0                                                # high vegetation flag (0 = no trees; 1 = yes)
    state.izh = 0.0                                                    # initial canopy height (m)
    state.isigfh = 0.0                                                 # initial canopy density
    state.hfol = 0.0                                                   # vegetation height (m)
    state.fola = 0.0                                                   # vegetation absorptivity
    state.epf = 0.0                                                    # vegetation emissivity
    state.albf = 1.0 - state.fola                                      # ponto 4: fola ja zerado -> 1.0
    state.lail = 0.0                                                   # vegetation leaf area index (m^2/m^2)
    state.sigfl = 0.0
    state.state = 0.0
    state.ftemp = SPFLAG
    state.z0l = 0.0
    state.zd = 0.0
    state.chnf = 0.0
    state.sqrt_chnf = 0.0
    state.trmlm = 0.0
    state.trmhm = 0.0
    state.rsl = 0.0
    state.storll = 0.0
    state.storls = 0.0
    state.sigfh = 0.0
    state.uaf = 0.0
    state.hfol_tot = 0.0

    state.water_flag = 0
    wtype = int(SPFLAG)                                                # water type: 0 = rivers, 1 = lakes
    wdepth = SPFLAG                                                    # water depth (m)
    wvel = SPFLAG                                                      # water velocity (m/s)

    state.avect[0:5, :] = 0.0                                          # i=0..4, j=1..nclayers: array todo (ponto 7)

    # ---------------- calculation variables (linhas 121-131) ----------------
    state.step = 0
    state.stepi = 0
    state.error_code = 0
    state.error_type = 0
    state.vitd_index = 0
    state.deltat_fasst = 0.0
    state.deltati = 0.0
    state.hpond = 0.0
    state.toptemp = 0.0
    state.ptemp = 0.0

    # ---------------- do i=1,maxn (linhas 133-153) ----------------
    state.ntype[1:n1] = 0
    state.icourse[1:n1] = 0
    state.zti[1:n1] = 0.0
    state.tm[1:n1] = 0.0
    state.zm[1:n1] = 0.0
    state.sm[1:n1] = 0.0
    state.soil_moist[1:n1] = 0.0
    state.nz[1:n1] = 0.0
    state.nzi[1:n1] = 0.0
    state.pheadmin[1:n1] = 0.0
    state.frl[1:n1] = 0.0
    state.frh[1:n1] = 0.0
    state.sinkr[1:n1] = 0.0
    state.delzs[1:n1] = 0.0
    state.bftm[1:n1] = 0.0
    state.nsoilp[1:n1, 1:MAXP + 1] = 0.0

    # ---------------- do i=1,maxn+extran (linhas 155-178) ----------------
    state.grthcond[1:nx1] = 0.0
    state.grspheat[1:nx1] = 0.0
    state.stt[1:nx1] = 0.0
    state.ice[1:nx1] = 0.0
    state.phead[1:nx1] = 0.0
    state.wvc[1:nx1] = 0.0
    state.khl[1:nx1] = 0.0
    state.khu[1:nx1] = 0.0
    state.vin[1:nx1] = 0.0
    state.flowu[1:nx1] = 0.0
    state.flowl[1:nx1] = 0.0
    state.fv1[1:nx1] = 0.0
    state.source[1:nx1] = 0.0
    state.sink[1:nx1] = 0.0
    state.too[1:nx1] = 0.0
    state.smoo[1:nx1] = 0.0
    state.woo[1:nx1] = 0.0
    state.ioo[1:nx1] = 0.0
    state.phoo[1:nx1] = 0.0
    state.dsmdh[1:nx1] = 0.0
    state.node_type[1:nx1] = ""                                        # ' ' (ponto 6)
    state.nclass[1:nx1] = ""

    # ---------------- do i=1,nclayers (linhas 180-189) ----------------
    c1 = NCLAYERS + 1
    state.ifoliage_type[1:c1] = SPFLAG
    state.ilai[1:c1] = SPFLAG
    state.iclump[1:c1] = SPFLAG
    state.irho[1:c1] = SPFLAG
    state.itau[1:c1] = SPFLAG
    state.ialp[1:c1] = SPFLAG
    state.ieps[1:c1] = SPFLAG
    state.idzveg[1:c1] = SPFLAG

    # ---------------- do i=1,maxl (linhas 191-202) ----------------
    l1 = MAXL + 1
    state.soiltype[1:l1] = 0
    state.lthick[1:l1] = 0.0
    state.rho_fac[1:l1] = 0.0
    state.rho_fac[1:l1] = 1.0                                          # ponto 4: escrita dupla
    state.stype[1:l1] = ""                                             # '  '
    state.sclass[1:l1] = ""                                            # ' '
    state.isoilp[1:l1, 1:MAXP + 1] = SPFLAG
    state.soilp[1:l1, 1:MAXP + 1] = SPFLAG

    # ---------------- stor, veg_prp, rk (linhas 204-216) ----------------
    state.stor[1:2 * NCLAYERS + 1] = 0.0
    state.veg_prp[1:19, 1:18] = 0.0
    state.rk[1:19, 1:n1] = 0.0

    return ZeroParamsOut(tflag, nt0, nm0, stflag, lcount, wtype, wvel,
                         wdepth, ttest, mtest, sname)


def zero_mxl_params(state: FasstState) -> None:
    """Zera os arrays alocaveis dimensionados por `state.maxlines`.

    Os alocaveis precisam estar alocados (ponto 10). Colunas 11..13 de
    `dmet1` NAO sao zeradas (ponto 2) e `slushy` nao e tocado (ponto 9).
    """
    n1 = state.maxlines + 1                                            # ponto 8: limite = maxlines, nao o shape

    for name in _MXL_REAL_1D:
        getattr(state, name)[1:n1] = 0.0

    state.airt[1:3, 1:n1] = 0.0                                        # airt(1,i), airt(2,i)
    state.sstate[1:n1] = ""                                            # ' '
    state.met[1:n1, 1:MAXCOL + 1] = 0.0
    state.dmet1[1:n1, 1:11] = 0.0                                      # ponto 2: so colunas 1..10 de 13
    state.canopy_temp[1:NCLAYERS + 1, 1:n1] = 0.0
