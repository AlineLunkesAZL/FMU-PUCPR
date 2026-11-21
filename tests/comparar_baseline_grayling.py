"""Teste de comparacao: FASST via adaptador da FMU x baseline do Grayling (issue #2).

Roda `FasstAdaptador.passo` hora a hora com `tests/inputs/grayling/grayling_1992.met`
(casos A e B), grava a saida com os nomes de coluna dos CSVs de `tests/reference/` e compara
com eles (tolerancia absoluta 0,01, a do `comparar_com_fortran.py`).

    python -B tests/comparar_baseline_grayling.py              # o teste
    python -B tests/comparar_baseline_grayling.py --controle   # baseline x baseline (deve passar)

Sai com codigo 1 se algum caso diverge, apontando o primeiro instante e a primeira coluna.
Com o modelo provisorio (`ModeloSubstituto`) o teste FALHA de proposito; passa a valer de
verdade quando o adaptador chamar a `FasstPasso` (issue #3).

Convencoes:
- passo k usa a linha k do .met (meteorologia do FIM do passo, como o FASST); nao ha
  deslocamento de serie, que so existe nos testes com fmpy.
- contorno="fasst_original": fundo com 0,15 W/m2 + lencol freatico. T_laje nao e usado
  nesse contorno; vai um valor finito qualquer.
- saidas gravadas em tests/saida/ (fora do git).
"""
import csv
import math
import os
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RAIZ)

from fasst_adaptador import FasstAdaptador  # noqa: E402

TOL = 0.01
DT = 3600.0
T_LAJE = 295.15  # [K] ignorado com contorno="fasst_original"
MET = os.path.join(RAIZ, "tests", "inputs", "grayling", "grayling_1992.met")
REF = os.path.join(RAIZ, "tests", "reference")
SAIDA = os.path.join(RAIZ, "tests", "saida")
TEMPO = ["time_s", "Year", "JD", "Hr", "M"]

# saida do adaptador -> coluna da baseline, por arquivo. PROVISORIO ate os Anexos A e B.
# Sinal do flx: positivo = para a superficie, logo LHeat = -LE (LE > 0 = evaporacao).
MAPA = {
    "fst": {"Grtemp": lambda s: s["T_sub"], "Grmoistv": lambda s: s["theta"]},
    "flx": {"LHeat": lambda s: -s["LE"]},
    "vgt": {"LVegTmp": lambda s: s["T_folha"]},
}
# parametros provisorios do ModeloSubstituto; a FasstPasso le os do fassti.inp de cada caso
CASOS = {
    "A": dict(params=dict(LAI=0.0, h_sub=2.0, theta0=0.12, T0=291.45), arquivos=["fst", "flx"]),
    "B": dict(params=dict(LAI=2.0, h_sub=2.0, theta0=0.12, T0=291.45), arquivos=["fst", "flx", "vgt"]),
}


def ler_csv(caminho):
    with open(caminho, newline="") as f:
        linhas = list(csv.reader(f))
    return linhas[0], linhas[1:]


def ler_met(caminho):
    """Linha 1 = cabecalho numerico, linha 2 = nomes, linha 3 = unidades, depois dados."""
    with open(caminho) as f:
        f.readline()
        nomes = f.readline().split()
        f.readline()
        return [dict(zip(nomes, map(float, l.split()))) for l in f if l.strip()]


def entradas(m):
    """Linha do .met -> argumentos do adaptador em SI."""
    for k in ("tmp", "rh", "ws", "prec", "tsol", "ir"):
        if m[k] == 9999.0:
            raise ValueError(f"{k} ausente (9999) em JD {m['doy']:.0f} {m['hr']:.0f}h")
    return dict(T_ar=m["tmp"] + 273.15, UR=m["rh"] / 100.0, G_sol=m["tsol"], L_ceu=m["ir"],
                vento=m["ws"], chuva=m["prec"] / DT, T_laje=T_LAJE)  # chuva: mm/passo -> kg/m2s


def rodar_adaptador(caso, tempos):
    """Devolve {arquivo: [linhas]} com as colunas de tempo + as colunas mapeadas."""
    met = ler_met(MET)
    a = FasstAdaptador(**CASOS[caso]["params"], contorno="fasst_original")
    saida = {arq: [] for arq in CASOS[caso]["arquivos"]}
    for k, t in enumerate(tempos):
        m = met[k]
        if (int(m["doy"]), int(m["hr"])) != (int(t[2]), int(t[3])):
            raise AssertionError(f"met e baseline desalinhados na linha {k}: {m['doy']} {m['hr']} x {t[2]} {t[3]}")
        s = a.passo(DT, **entradas(m))
        for arq in saida:
            saida[arq].append(list(t) + [f(s) for f in MAPA[arq].values()])
    os.makedirs(SAIDA, exist_ok=True)
    for arq, linhas in saida.items():
        with open(os.path.join(SAIDA, f"grayling-{caso}-{arq}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(TEMPO + list(MAPA[arq]))
            w.writerows(linhas)
    return {arq: (TEMPO + list(MAPA[arq]), linhas) for arq, linhas in saida.items()}


def primeira_divergencia(ref, obtido):
    """(linha, coluna, esperado, obtido) da 1a divergencia > TOL, por instante e ordem de coluna."""
    cab_r, lin_r = ref
    cab_o, lin_o = obtido
    cols = [c for c in cab_r if c in cab_o and c not in TEMPO]
    ir, io = [cab_r.index(c) for c in cols], [cab_o.index(c) for c in cols]
    if len(lin_r) != len(lin_o):
        return (min(len(lin_r), len(lin_o)), "numero de linhas", len(lin_r), len(lin_o))
    for k, (r, o) in enumerate(zip(lin_r, lin_o)):
        for c, i, j in zip(cols, ir, io):
            try:
                esp, obt = float(r[i]), float(o[j])
            except ValueError:  # colunas de texto (FTU, USCS): igualdade exata
                if r[i] != o[j]:
                    return (k, c, r[i], o[j])
                continue
            if not abs(esp - obt) <= TOL:  # NaN/inf contam como divergencia
                return (k, c, esp, obt)
    return None


def main(controle=False):
    falhou = False
    for caso in CASOS:
        refs = {arq: ler_csv(os.path.join(REF, f"2026-11-12-grayling-{caso}-{arq}.csv"))
                for arq in CASOS[caso]["arquivos"]}
        tempos = [l[:5] for l in refs["fst"][1]]
        obtidos = refs if controle else rodar_adaptador(caso, tempos)
        divs = {arq: primeira_divergencia(refs[arq], obtidos[arq]) for arq in refs}
        divs = {a: d for a, d in divs.items() if d}
        if not divs:
            print(f"caso {caso}: IGUAL A BASELINE ({', '.join(refs)}; tol {TOL})")
            continue
        falhou = True
        arq, (k, col, esp, obt) = min(divs.items(), key=lambda x: x[1][0])
        t = refs[arq][1][k] if k < len(refs[arq][1]) else ["?"] * 5
        print(f"caso {caso}: DIVERGE em {arq}.{col}, linha {k} (time_s {t[0]}, dia {t[2]} {t[3]}h): "
              f"esperado {esp}, obtido {obt}, dif {abs(esp - obt) if isinstance(esp, float) else '-'}")
        for a, (k2, c2, _, _) in divs.items():
            if a != arq:
                print(f"         tambem diverge em {a}.{c2} a partir da linha {k2}")
    return 1 if falhou else 0


if __name__ == "__main__":
    sys.exit(main(controle="--controle" in sys.argv))
