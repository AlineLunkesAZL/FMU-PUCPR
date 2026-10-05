"""
us_soil_tools.py -- traducao de US_soil_tools.F90

Fonte Fortran: US_soil_tools.F90 (413 linhas, 2 sub-rotinas):
`get_user_soil_params` e `upr_case` (apendice da primeira, chamada por
ela e tambem por initsurface.F90).

`upr_case` converte as letras minusculas ASCII (a-z) de uma cadeia em
maiusculas. `get_user_soil_params` procura, num arquivo de tipos de solo
definidos pelo usuario (unidade 90), o solo cujo nome de 4 caracteres
coincide com `name`, le seus 22 parametros (um por linha, formato
list-directed), valida cada faixa, completa os parametros que podem ser
derivados dos demais e, para solos que nao sejam materiais especiais,
estima a umidade minima/maxima e os parametros de van Genuchten com a
funcao de pedotransferencia de Vereecken et al. (1989). Se o solo nao e
encontrado, devolve um padrao (areia siltosa, USCS 7).

Assinaturas originais:
    subroutine get_user_soil_params(name,usclass,ustype,rhod,poros,alb,em,
        qtz,kdry,sathdrc,minwc,maxwc,vgbpres,vgexp,spheat,vgm,orgf,psand,
        psilt,pclay,pcarbon,pl,p200)
    subroutine upr_case(d1,al)

Arquitetura
-------------
Nenhuma das duas sub-rotinas usa campos de estado: `get_user_soil_params`
importa de `fasst_global` apenas as constantes `eps` e `spflag`
(`EPS`/`SPFLAG` de `fasst.constants`). Por isso as duas sao funcoes
independentes, sem `state: FasstState`.

O arquivo da unidade 90 e recebido explicitamente (`unit90`, objeto de
arquivo ja aberto, com `readline()`), em vez de lido de um numero de
unidade global. A posicao de leitura e mantida pelo proprio objeto:
quem chama (initsurface.py) faz `unit90.seek(0)` depois de cada chamada,
equivalente ao `rewind(90)` do Fortran.

`get_user_soil_params` devolve uma tupla na mesma ordem da lista de
argumentos do Fortran:
    (name, usclass, ustype, rhod, poros, alb, em, qtz, kdry, sathdrc,
     minwc, maxwc, vgbpres, vgexp, spheat, vgm, orgf, psand, psilt,
     pclay, pcarbon, pl, p200)
`name` (`intent(inout)`) volta em maiusculas, com 4 caracteres.
`usclass` e 'USCS' ou 'USDA'. `upr_case` devolve a cadeia convertida.

Formato do arquivo (reconstruido a partir do codigo de leitura): cada
solo ocupa 1 linha de nome (4 primeiros caracteres), 22 linhas de dados
(ustype, rhod, rhoid, thetas, poros, alb, em, qtz, orgf, kdry, spheat,
sathdrc, minwc, maxwc, vgbpres, vgexp, psand, psilt, pclay, pcarbon, pl,
p200) e 1 linha em branco (`nslines = 23` linhas puladas por solo que
nao coincide). Cada leitura de dado le o primeiro valor da proxima linha
nao vazia (leitura list-directed: separadores sao espaco e virgula, e o
restante da linha e ignorado); o expoente pode ser escrito com `d`/`D`.

Pontos de atencao
-------------------
* Provavel inconsistencia de unidades no Fortran original, preservada
  sem correcao: a funcao de pedotransferencia de Vereecken et al. (1989)
  (linhas 297-330 do .F90) recebe `psand`, `pclay` e `pcarbon`, que o
  proprio codigo valida na faixa [0, 1] (fracoes) nas linhas 222-248,
  mas os coeficientes (por exemplo `2.5d-2*psand`, `-2.3d-2*pclay`,
  `1.5d-4*psand*psand`) correspondem, na forma usual dessa funcao, a
  teores em percentual (0-100). Com fracoes, os termos de areia, argila
  e carbono ficam cerca de 100 vezes menores do que o esperado, e `vgbpres`
  e `vgexp` dependem quase so da densidade (`rhod`). Confirmar contra a
  publicacao antes de qualquer correcao.
* O arquivo e percorrido uma unica vez, sequencialmente: a funcao nao
  rebobina a unidade 90 (quem chama faz `seek(0)` depois). Se a unidade
  nao estiver no inicio, solos anteriores a posicao atual nao sao
  encontrados.
* Quando a leitura do nome chega ao fim do arquivo, o Fortran deixa
  `sname` com o valor anterior e prossegue pela comparacao; como esse
  valor ja foi comparado e rejeitado, o resultado e o mesmo que sair do
  laco. A traducao reproduz esse caminho.
* A mensagem do caso "solo nao encontrado" diz "defaulting to sand", mas
  o padrao atribuido e `ustype = 7` (areia siltosa, "silty sand", de
  Graves, MA). Discrepancia de texto do proprio Fortran.
* Uma leitura que falha (valor nao numerico), com valor nulo (virgula
  inicial ou `/`) ou que encontra o fim do arquivo deixa a variavel com
  o valor que ela tinha, sem interromper -- a leitura Fortran usa
  `iostat=`, que suprime o aborto. Esse valor e `SPFLAG` para os
  parametros devolvidos, mas 0.0 para `rhoid` e `thetas` (variaveis
  locais zeradas no topo), que entao disparam a verificacao de faixa. As
  verificacoes seguem rodando sobre o valor retido.
* Caminho de erro: qualquer parametro fora da faixa liga `ptest`; ao
  final, `ptest == 1` imprime a mensagem e interrompe a execucao
  (`stop`), levantando `SystemExit` aqui. `sathdrc`, `minwc`, `maxwc`,
  `vgbpres` e `vgexp` (este so no ramo `vgexp <= 0`) so contam como
  erro se o solo nao for material especial (`tflag`, ustype 20, 21, 25
  ou 26: concreto, asfalto, rocha, agua).
* `rhoid/thetas` e `rhod/rhoid` (linhas 268 e 271) podem dividir por
  zero quando `thetas` ou `rhoid` sao nulos (valores ja sinalizados como
  fora da faixa por `ptest`); o Fortran produz Infinity e segue ate o
  `stop`. `_fdiv` reproduz isso em vez de levantar `ZeroDivisionError`
  antes da mensagem de erro.
* Ordem das derivacoes: `maxwc = poros` (linha 285) roda antes da
  pedotransferencia. Por isso a estimativa de `maxwc` pela PTF
  (`0.81 - 0.283*rhod + 1e-3*pclay`) so e usada quando `poros` tambem nao
  foi informado nem derivado; com `poros` conhecido, `maxwc` ja vale
  `poros` e a PTF nao o altera.
* `ustype` fora de `[-99, 99]` e tratado como codigo USDA (`usclass =
  'USDA'`): valores positivos passam por `map_usda_soiltype_to_uscs`
  (`fasst.functions`), que levanta `ValueError` fora da faixa da tabela
  (divergencia documentada daquela funcao em relacao ao Fortran).
"""

import math
import re

import numpy as np

from .constants import EPS, SPFLAG
from .functions import map_usda_soiltype_to_uscs


def _fdiv(x, y):
    """Divisao IEEE-754 pura: ``x / 0.0`` -> +-inf / nan, sem levantar."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.float64(x) / np.float64(y))


def upr_case(d1, al):
    """Fortran ``upr_case(d1, al)``: converte a-z em A-Z nos primeiros ``d1``
    caracteres de ``al`` (a conversao e ASCII: codigos 97-122 menos 32)."""
    head = "".join(chr(ord(c) - 32) if "a" <= c <= "z" else c for c in al[:d1])
    return head + al[d1:]


def _isflag(x):
    return abs(x - SPFLAG) <= EPS


def _notflag(x):
    return abs(x - SPFLAG) > EPS


def _read_token(unit90):
    """Equivalente a ``read(90,*,iostat=io) x``: devolve ``(token, io)``.

    Pula linhas vazias (a leitura list-directed continua no registro
    seguinte ate encontrar um valor). ``token`` e ``None`` no fim do arquivo
    (``io == -1``) ou quando o primeiro item e um valor nulo (virgula
    inicial ou ``/``), caso em que a variavel fica inalterada.
    """
    while True:
        line = unit90.readline()
        if line == "":
            return None, -1
        text = line.strip("\r\n").strip()
        if text == "":
            continue
        if text[0] in ",/":
            return None, 0
        return re.split(r"[\s,/]+", text, maxsplit=1)[0], 0


def _read_real(unit90, current):
    """``read(90,*,iostat=io) x`` para variavel real; mantem ``current`` se a
    leitura nao entregar um valor valido."""
    tok, _io = _read_token(unit90)
    if tok is None:
        return current
    try:
        return float(tok.replace("d", "e").replace("D", "e"))
    except ValueError:
        return current


def _read_int(unit90, current):
    """``read(90,*,iostat=io) k`` para variavel inteira."""
    tok, _io = _read_token(unit90)
    if tok is None:
        return current
    try:
        return int(tok)
    except ValueError:
        return current


def get_user_soil_params(unit90, name):
    # ---------------- zero-out variables (linhas 38-73) ----------------
    ptest = 0
    tflag = 0
    nslines = 23  # number of data lines + blank
    io = 0
    rhoid = 0.0
    thetas = 0.0
    lnn = 0.0
    lnalpha = 0.0
    sname = "    "

    ustype = 0
    rhod = poros = alb = em = qtz = kdry = sathdrc = SPFLAG
    minwc = maxwc = vgbpres = vgexp = spheat = vgm = SPFLAG
    orgf = psand = psilt = pclay = pcarbon = pl = p200 = SPFLAG
    usclass = "  "

    # Make sure the name is in caps
    name = name[:4].ljust(4)
    d1 = len(name.rstrip())
    name = upr_case(d1, name)
    print("User soil name  ", name)

    while io != -1:
        line = unit90.readline()
        if line == "":
            io = -1  # fim do arquivo: sname mantem o valor anterior
        else:
            io = 0
            sname = line.strip("\r\n")[:4].ljust(4)

        d1 = len(sname.rstrip())
        sname = upr_case(d1, sname)

        if name == sname:
            ustype = _read_int(unit90, ustype)
            if abs(ustype) < 100:  # USCS classification
                usclass = "USCS"
            else:
                usclass = "USDA"
                if ustype > 0:
                    lis = ustype - 100
                    ustype = map_usda_soiltype_to_uscs(lis)
                else:
                    ustype = -(abs(ustype) - 100)

            if ustype in (20, 21, 25, 26):
                tflag = 1

            rhod = _read_real(unit90, rhod)  # bulk dry density (g/cm^3)
            if rhod <= 0.0 and _notflag(rhod):
                print(name, " bulk dry density out of range")
                ptest = 1

            rhoid = _read_real(unit90, rhoid)  # intrinsic dry density (g/cm^3)
            if rhoid <= 0.0 and _notflag(rhoid):
                print(name, " intrinsic dry density out of range")
                ptest = 1

            thetas = _read_real(unit90, thetas)  # solids fraction (0.0 - 1.0)
            if (thetas <= 0.0 or thetas > 1.0) and _notflag(thetas):
                print(name, " solids fraction out of range")
                ptest = 1

            poros = _read_real(unit90, poros)  # porosity (0.0 - 1.0)
            if (poros <= 0.0 or poros > 1.0) and _notflag(poros):
                print(name, " porosity out of range")
                ptest = 1

            alb = _read_real(unit90, alb)  # albedo (0.0 - 1.0)
            if (alb <= 0.0 or alb > 1.0) and _notflag(alb):
                print(name, " albedo out of range")
                ptest = 1

            em = _read_real(unit90, em)  # emissivity (0.0 - 1.0)
            if (em <= 0.0 or em > 1.0) and _notflag(em):
                print(name, " emissivity out of range")
                ptest = 1

            qtz = _read_real(unit90, qtz)  # quartz content (0.0 - 1.0)
            if (qtz < 0.0 or qtz > 1.0) and _notflag(qtz):
                print(name, " quartz content out of range")
                ptest = 1

            orgf = _read_real(unit90, orgf)  # organic fraction of solids (0.0 - 1.0)
            if (orgf < 0.0 or orgf > 1.0) and _notflag(orgf):
                print(name, " organic fraction out of range")
                ptest = 1

            kdry = _read_real(unit90, kdry)  # dry thermal conductivity (W/m*K)
            if kdry <= 0.0 and _notflag(kdry):
                print(name, " dry thermal conductivity out of range")
                ptest = 1

            spheat = _read_real(unit90, spheat)  # specific heat (J/kg*K)
            if (spheat <= 0.0 or spheat > 2200.0) and _notflag(spheat):
                print(name, " specific heat out of range")
                ptest = 1

            sathdrc = _read_real(unit90, sathdrc)  # saturated hydraulic conductivity (cm/s)
            if sathdrc <= 0.0 and _notflag(sathdrc):
                if tflag != 1:
                    print(name, " saturated hydraulic conductivity out of range")
                    ptest = 1

            minwc = _read_real(unit90, minwc)  # residual (minimum) water content (vol/vol)
            if (minwc <= 0.0 or minwc > poros) and _notflag(minwc):
                if tflag != 1:
                    print(name, " min water content out of range")
                    ptest = 1

            maxwc = _read_real(unit90, maxwc)  # maximum water content (vol/vol)
            if (maxwc <= 0.0 or maxwc > poros) and _notflag(maxwc):
                if tflag != 1:
                    print(name, " max water content out of range")
                    ptest = 1

            vgbpres = _read_real(unit90, vgbpres)  # van Genuchten bubbling pressure head (cm)
            if vgbpres <= 0.0 and _notflag(vgbpres):
                if tflag != 1:
                    print(name, " van Genuchten pressure head out of range")
                    ptest = 1
            if _notflag(vgbpres) and vgbpres > EPS:
                vgbpres = 1.0 / vgbpres  # 1/cm

            vgexp = _read_real(unit90, vgexp)  # van Genuchten exponent, n
            if vgexp <= 0.0 and _notflag(vgexp):
                if tflag != 1:
                    print(name, " van Genuchten exponent out of range")
                    ptest = 1
            elif abs(vgexp - 1.0) <= EPS:
                print(name, " van Genuchten exponent can not equal 1")  # causes a divide by 0 error
                ptest = 1
            if _notflag(vgexp) and vgexp > EPS:
                vgm = 1.0 - 1.0 / vgexp  # unitless

            psand = _read_real(unit90, psand)  # percent sand
            if (psand < 0.0 or psand > 1.0) and _notflag(psand):
                print(name, " percent sand out of range")
                ptest = 1

            psilt = _read_real(unit90, psilt)  # percent silt
            if (psilt < 0.0 or psilt > 1.0) and _notflag(psilt):
                print(name, " percent silt out of range")
                ptest = 1

            pclay = _read_real(unit90, pclay)  # percent clay
            if (pclay < 0.0 or pclay > 1.0) and _notflag(pclay):
                print(name, " percent clay out of range")
                ptest = 1

            pcarbon = _read_real(unit90, pcarbon)  # carbon content
            if (pcarbon < 0.0 or pcarbon > 1.0) and _notflag(pcarbon):
                print(name, " percent carbon out of range")
                ptest = 1

            pl = _read_real(unit90, pl)  # plastic limit
            if (pl < 0.0 or pl > 2e2) and _notflag(pl):
                print(name, " plastic limit out of range")
                ptest = 1

            p200 = _read_real(unit90, p200)  # percent fines passing #200 sieve
            if (p200 < 0.0 or p200 > 1e2) and _notflag(p200):
                print(name, " percent passing #200 sieve out of range")
                ptest = 1

            # See if any parameters can be calculated from others
            if _isflag(rhod) and (_notflag(rhoid) and _notflag(thetas)):
                rhod = rhoid * thetas

            if _isflag(rhoid) and (_notflag(rhod) and _notflag(thetas)):
                rhoid = _fdiv(rhod, thetas)

            if _isflag(thetas) and (_notflag(rhod) and _notflag(rhoid)):
                thetas = _fdiv(rhod, rhoid)

            if _isflag(poros) and _notflag(thetas):
                poros = 1.0 - thetas
                if poros <= 0.0 or poros > 1.0:
                    print(name, " porosity out of range")
                    ptest = 1

            if _isflag(rhod) and (_notflag(rhoid) and _notflag(poros)):
                rhod = rhoid * (1.0 - poros)

            if _notflag(poros) and _isflag(maxwc):
                maxwc = poros

            if _isflag(orgf) and _notflag(pcarbon):
                orgf = 2.0 * pcarbon  # SSSA recommendation
                if orgf > 1.0:
                    orgf = 1.0

            if _isflag(qtz) and _notflag(psand):
                qtz = psand  # Peters-Lidard et al., J. Atmos. Sci., V55 (1998)

            # calculate max water content, min water content, van Genuchten
            # parameters using PTF of Vereecken et al. (1989) "Estimating the
            # soil moisture retention characteristic from texture, bulk density
            # and carbon content", Soil Science, V.148, No.6, p.389-403
            # (ver "Pontos de atencao" sobre as unidades de psand/pclay/pcarbon)
            if tflag != 1 and ustype != 30:
                if _notflag(psand) and _notflag(pclay) and _notflag(pcarbon):
                    if _isflag(psilt):
                        psilt = 1.0 - (psand + pclay)

                    if _isflag(minwc):
                        minwc = 1.5e-2 + 5e-3 * pclay + 1.4e-2 * pcarbon

                    if _isflag(maxwc):
                        maxwc = 0.81 - 0.283 * rhod + 1e-3 * pclay

                    if maxwc > poros:
                        maxwc = poros

                    if _isflag(vgbpres):
                        lnalpha = (-2.486 + 2.5e-2 * psand - 0.351 * pcarbon
                                   - 2.617 * rhod - 2.3e-2 * pclay)
                        if abs(lnalpha) < 5e1:
                            vgbpres = math.exp(lnalpha)  # 1/cm

                    if _isflag(vgexp):
                        lnn = 5.3e-2 - 9e-3 * psand - 1.3e-2 * pclay + 1.5e-4 * psand * psand
                        if abs(lnn) < 5e1 and abs(lnn) > EPS:
                            vgexp = math.exp(lnn)
                            vgm = 1.0 - 1.0 / vgexp

            if ptest == 1:
                print(" Soil Input Parameters out of range, STOPPING")
                raise SystemExit("get_user_soil_params: Soil Input Parameters out of range")
            return (name, usclass, ustype, rhod, poros, alb, em, qtz, kdry, sathdrc,
                    minwc, maxwc, vgbpres, vgexp, spheat, vgm, orgf, psand, psilt,
                    pclay, pcarbon, pl, p200)
        else:
            for _ir in range(nslines):  # skip to the next soil type
                unit90.readline()

    # Unable to find the user defined soil type; silty-sand will be used as a default
    print(" !! Can not find your soil, defaulting to sand !! ")  # Graves silty-sand, MA
    usclass = "USCS"
    ustype = 7  # silty sand
    rhod = 1.49  # Guyman et al. (1993), p.55
    rhoid = 2.73  # Guyman et al. (1993), p.55  # noqa: F841 (so usada para thetas)
    thetas = rhod / rhoid
    poros = 1.0 - thetas
    alb = 0.35  # Sullivan et al. (1997), p.38
    em = 0.92  # Sullivan et al. (1997), p.51
    qtz = 0.8  # Tarnawski (1997), p.96
    orgf = 0.0
    kdry = 0.831  # Jordan (2000)   (0.281 comentado no Fortran)
    spheat = 830.0
    sathdrc = 0.000533  # Guyman et al. (1993), p.55
    minwc = 0.001
    maxwc = poros
    vgbpres = 1.0 / 23.5474  # Jordan (2000)
    vgexp = 1.50  # Jordan (2000)
    vgm = 1.0 - 1.0 / vgexp  # unitless
    psand = SPFLAG  # Loamy Sand, Cospy et al. (1984), p.683
    psilt = SPFLAG
    pclay = SPFLAG
    pcarbon = SPFLAG  # Vereecken et al. (1989), p.391
    pl = SPFLAG
    p200 = SPFLAG
    return (name, usclass, ustype, rhod, poros, alb, em, qtz, kdry, sathdrc,
            minwc, maxwc, vgbpres, vgexp, spheat, vgm, orgf, psand, psilt,
            pclay, pcarbon, pl, p200)
