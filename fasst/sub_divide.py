"""
sub_divide.py -- traducao de sub_divide.F90

Fonte Fortran: sub_divide.F90 (163 linhas). Uma unica sub-rotina
publica (`sub_divide`), sem nenhuma sub-rotina chamada ("no
subroutines called" no cabecalho do .F90). Interpola as variaveis
meteorologicas do passo de tempo atual para uma subdivisao interna
("setorizacao") do intervalo entre dois registros meteorologicos
consecutivos -- usada quando o passo de integracao do FASST e menor
que o passo do arquivo meteorologico, precisando de valores
intermediarios entre duas leituras. Tambem aplica a correcao de
radiacao solar/infravermelha por inclinacao do terreno (Jordan, 1991),
quando o sitio tem uma inclinacao (`slope_fasst`) diferente de zero.

Assinatura original:
    subroutine sub_divide(ii,delvar,delvar1,dmet)
    integer(ip),intent(in):: ii
    real(dp),intent(in):: delvar(maxcol),delvar1(12)
    real(dp),intent(out):: dmet(13)

Arquitetura
-------------
`sub_divide` recebe `state: FasstState` como argumento explicito. `EPS`
e `PI` vem de `fasst.constants`; `mflag`, `slope_fasst`, `albedo_fasst`,
`sloper`, `aspect`, `met`, `dmet1`, `mstflag`, `iw` sao campos de
`FasstState`. `ip_az`/`ip_zen`/`ip_tsoil`/`ip_irup` (indices de coluna
do meteorologico, no Fortran) sao `MetCol.AZ`/`MetCol.ZEN`/
`MetCol.TSOIL`/`MetCol.IRUP` -- usados tanto para indexar `met`/`dmet1`
quanto `delvar`, que tem a mesma largura de coluna que `met`
(`delvar(maxcol)`).

`dmet` e o unico parametro `intent(out)`; `delvar`/`delvar1` sao
`intent(in)` puros. A funcao devolve `dmet` diretamente (nao em tupla
de 1 elemento), mesma convencao usada para outras funcoes deste
projeto com uma unica saida (por exemplo, `emisatm` em
module_radiation.py).

Indexacao 1-based identica ao Fortran (posicao 0 sem uso) em `delvar`,
`delvar1` e `dmet`.

Pontos de atencao
-------------------
* Provavel bug do Fortran original, preservado sem correcao: no ramo
  em que `solzen` (em radianos) excede 85 graus, o denominador da
  correcao de radiacao direta usa `dcos(85d0)` -- 85.0 SEM a conversao
  `*f2` (graus para radianos) que todo o resto da sub-rotina aplica
  consistentemente (`90d0*f2`, `85d0*f2` nas comparacoes de limiar
  logo acima, e o proprio `dcos(solzen)` do ramo irmao, onde `solzen`
  ja foi convertido). `dcos(85d0)` calcula o cosseno de 85 RADIANOS
  (um angulo que deu varias voltas completas), nao de 85 graus --
  quase certamente deveria ser `dcos(85d0*f2)`. Traduzido literalmente
  (`math.cos(85.0)`, nao `math.cos(85.0*f2)`).
* `difcoeff`, `sdifcorr1`, `sdifcorr2` sao zeradas no topo da
  sub-rotina mas nunca atribuidas de novo em nenhum lugar do codigo
  ativo -- a formula que as usaria (uma correcao alternativa para a
  radiacao difusa, citada sem referencia) esta inteira comentada no
  Fortran original. `dmet(10)` (radiacao difusa) e recalculada por um
  caminho mais simples (`dmet(1) - dmet(9)`, total menos direta) logo
  depois. Codigo morto inofensivo, preservado por fidelidade.
* A correcao de radiacao de onda longa por inclinacao (Sellers et al.
  1995, citada no comentario) tambem esta inteira comentada -- `dmet(2)`
  (radiacao IV incidente) nunca e ajustada por inclinacao do terreno,
  apesar do modelo ter sido citado para isso.
* `dmet(6)`/`dmet(7)` vem sempre de `delvar1(6)`/`delvar1(7)`, nos dois
  ramos (`iw>=2` e `iw<2`) -- no ramo `iw<2`, uma alternativa comentada
  (`dmet1(iw,6)`/`dmet1(iw,7)`) mostra que essa fonte foi considerada e
  descartada.
"""

import math

from .constants import EPS, PI
from .indices import MetCol
from .state import FasstState


def _new_array(n, fill=0.0):
    """1-based com folga no indice 0 -- ver fasst/state.py, `_reals`/`_ints`."""
    return [fill] * (n + 1)


def _aint(x, p=0):
    """Fortran ``AINT(x*10**p)*10**-p`` -- trunca em direção a zero."""
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


def sub_divide(state: FasstState, ii, delvar, delvar1):
    # ---------------- zero-out variables (linhas 21-34) ----------------
    solzen = 0.0
    saz = 0.0
    sdir = 0.0
    sdif = 0.0
    delaz = 0.0
    difcoeff = 0.0  # nunca reatribuida em codigo ativo -- ver docstring
    sdircorr = 0.0
    sdifcorr1 = 0.0  # nunca reatribuida em codigo ativo -- ver docstring
    sdifcorr2 = 0.0  # nunca reatribuida em codigo ativo -- ver docstring

    dmet = _new_array(13)

    # ---------------- sectorize the meteorological parameters (linhas 37-96) ----------------
    ltest = float(ii)  # -1
    if ltest < 0.0:
        ltest = 0.0

    if state.iw >= 2:
        solzen = state.met[state.iw - 1][MetCol.ZEN] + ltest * delvar[MetCol.ZEN]
        saz = state.met[state.iw - 1][MetCol.AZ] + ltest * delvar[MetCol.AZ]

        dmet[1] = state.dmet1[state.iw - 1][1] + ltest * delvar1[1]  # total incoming solar
        dmet[2] = state.dmet1[state.iw - 1][2] + ltest * delvar1[2]  # incoming ir
        dmet[3] = state.dmet1[state.iw - 1][3] + ltest * delvar1[3]  # wind speed
        dmet[4] = state.dmet1[state.iw - 1][4] + ltest * delvar1[4]  # air temperature (K)
        dmet[5] = state.dmet1[state.iw - 1][5] + ltest * delvar1[5]  # relative humidity

        dmet[6] = delvar1[6]
        dmet[7] = delvar1[7]

        if _aint(abs(delvar1[8] - state.mflag) * 1e5) * 1e-5 > EPS:  # reflected solar
            dmet[8] = state.dmet1[state.iw - 1][8] + ltest * delvar1[8]
        else:
            dmet[8] = state.mflag

        dmet[9] = state.dmet1[state.iw - 1][9] + ltest * delvar1[9]  # direct solar
        dmet[10] = state.dmet1[state.iw - 1][10] + ltest * delvar1[10]  # diffuse solar
        dmet[11] = state.dmet1[state.iw - 1][11] + ltest * delvar1[11]  # air pressure

        if _aint(abs(delvar1[12] - state.mflag) * 1e5) * 1e-5 > EPS:  # upwelling ir
            dmet[12] = state.met[state.iw - 1][MetCol.IRUP] + ltest * delvar1[12]
        else:
            dmet[12] = state.mflag

        if (state.mstflag == 1
                and _aint(abs(delvar[MetCol.TSOIL] - state.mflag) * 1e5) * 1e-5 > EPS):
            # measured soil surface temp
            dmet[13] = state.met[state.iw - 1][MetCol.TSOIL] + ltest * delvar[MetCol.TSOIL]
        else:
            dmet[13] = state.mflag

    else:
        solzen = state.met[state.iw][MetCol.ZEN]
        saz = state.met[state.iw][MetCol.AZ]

        dmet[1] = state.dmet1[state.iw][1]  # total incoming solar
        dmet[2] = state.dmet1[state.iw][2]  # incoming IR
        dmet[3] = state.dmet1[state.iw][3]  # wind speed
        dmet[4] = state.dmet1[state.iw][4]  # air temperature
        dmet[5] = state.dmet1[state.iw][5]  # relative humidity

        # dmet[6] = state.dmet1[state.iw][6]  -- comentado no Fortran, ver docstring
        # dmet[7] = state.dmet1[state.iw][7]  -- comentado no Fortran, ver docstring
        dmet[6] = delvar1[6]
        dmet[7] = delvar1[7]

        dmet[8] = state.dmet1[state.iw][8]  # reflected solar
        dmet[9] = state.dmet1[state.iw][9]  # direct incoming solar
        dmet[10] = state.dmet1[state.iw][10]  # diffuse incoming solar
        dmet[11] = state.dmet1[state.iw][11]  # air pressure
        dmet[12] = state.dmet1[state.iw][12]  # reflected & emitted IR
        dmet[13] = state.dmet1[state.iw][13]  # measured soil surface temp

    if dmet[1] <= 1e-2:
        dmet[1] = 0.0
    dmet[2] = max(0.0, dmet[2])

    sdir = max(0.0, dmet[9])
    sdif = max(0.0, dmet[10])

    # ---------------- CALCULATE THE SHORT and LONG WAVE RADIATION TERM (linhas 104-150) ----------------
    # correct for a surface slope based on Jordan, R., CRREL Rep. 91-16, p.25.
    # solar zenith angle[solzen(degrees)]
    # solar azimuthal angle; + = clockwise [saz(degrees from north)]
    # aspect: 180 = true South; 360 = true North
    if state.slope_fasst > EPS:
        f2 = PI / 180.0
        solzen = solzen * f2

        if solzen >= 90.0 * f2:
            dmet[1] = 0.0
        elif solzen < 90.0 * f2 and dmet[1] > 0.0:
            if state.sloper > EPS and state.sloper < 90.0 * f2:
                delaz = abs(state.aspect - saz) * f2
                # short wave
                if solzen <= 85.0 * f2:
                    sdircorr = (
                        (math.sin(solzen) * math.sin(state.sloper) * math.cos(delaz))
                        / math.cos(solzen)
                    )
                else:
                    # dcos(85d0) sem conversao para radianos -- ver "Pontos de
                    # atencao" no docstring do modulo; preservado literalmente.
                    sdircorr = (
                        (math.sin(solzen) * math.sin(state.sloper) * math.cos(delaz))
                        / math.cos(85.0)
                    )

                sdir = max(0.0, sdir * (math.cos(state.sloper) + sdircorr))
                dmet[9] = sdir

                # difcoeff = 1d0 + 0.5d0*sin(solzen) + 2d0*sin(2d0*solzen)
                # sdifcorr1 = (pi - delaz)*(difcoeff + dcos(sloper))
                # sdifcorr2 = delaz*(1d0 + difcoeff*dcos(sloper))
                # sdif = dmax1(0d0,sdif*((sdifcorr1+sdifcorr2)/(pi*(1d0+difcoeff))))
                # -- comentado no Fortran, ver docstring do modulo.

                pre_slope = dmet[1]
                dmet[1] = (
                    sdir + sdif
                    + dmet[1] * state.albedo_fasst * (1.0 - math.cos(state.sloper)) * 0.5
                )
                if dmet[1] < 1e-2:
                    dmet[1] = 0.0

                dmet[10] = dmet[1] - dmet[9]

                if (_aint(abs(dmet[8] - state.mflag) * 1e5) * 1e-5 > EPS
                        and pre_slope > EPS):
                    dmet[8] = dmet[8] * dmet[1] / pre_slope

                # long wave (Sellers et al. (1995), JGR 100(D12), p.25607-25629)
                # dmet(2) = dmax1(0d0,dmet(2)*((pi - sloper)/pi)/dcos(sloper))
                # -- comentado no Fortran, ver docstring do modulo.

    if abs(dmet[1]) <= EPS:  # no up-solar if none coming in
        dmet[8] = 0.0
    elif (dmet[8] >= dmet[1]
            and _aint(abs(dmet[8] - state.mflag) * 1e5) * 1e-5 > EPS):
        # can't have more going out than coming in
        dmet[8] = state.albedo_fasst * dmet[1]

    for i in range(1, 14):
        dmet[i] = _anint(dmet[i], 20)

    return dmet
