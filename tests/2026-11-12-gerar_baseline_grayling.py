"""Gera a baseline do FASST Python nos casos Grayling A (solo nu) e B (grama).

    python -B tests/2026-11-12-gerar_baseline_grayling.py <pacote_ref> <tabelas> <saida>

  <pacote_ref>  pasta referencia_fortran_grayling (grayling_1992.met, FASST_default*.inp,
                A_solo_nu/ e B_grama/ com fassti.inp e ref_*.csv)
  <tabelas>     pasta com MODIS_IGBP_vegprop.inp, UMD_LDAS_vegprop.inp e fasstusersoil.inp
                (o driver Python exige esses arquivos; o pacote nao os traz)
  <saida>       onde gravar os CSVs

Cada caso roda num diretorio temporario proprio, num processo novo (python -B, sem
__pycache__). As saidas saida.fst/.flx/.vgt/.grd viram CSV com os MESMOS nomes de coluna
dos ref_*.csv (compativel com comparar_com_fortran.py), mais uma coluna `time_s` (segundos
desde o primeiro registro do met). Os valores sao copiados como o FASST os imprimiu, sem
reformatar. So entra o periodo comparado com o Fortran: dias 259 a 285 (648 h, sem neve).
"""
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # tem `fasst/`
CASOS = {"A": "A_solo_nu", "B": "B_grama"}
SAIDAS = ("fst", "flx", "vgt", "grd")
JD_INI, JD_FIM = 259, 285
TABELAS = ("MODIS_IGBP_vegprop.inp", "UMD_LDAS_vegprop.inp", "fasstusersoil.inp")


def rodar(pacote, tabelas, caso, tmp):
    for f in ("grayling_1992.met", "FASST_defaultlayerprop.inp", "FASST_defaultvegprop.inp"):
        shutil.copy(os.path.join(pacote, f), tmp)
    shutil.copy(os.path.join(pacote, CASOS[caso], "fassti.inp"), tmp)
    for f in TABELAS:
        shutil.copy(os.path.join(tabelas, f), tmp)
    env = dict(os.environ, PYTHONPATH=RAIZ, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-B", "-c",
                        "from fasst.fasst_driver import main; main(['fasst_driver', 'fassti.inp'])"],
                       cwd=tmp, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"caso {caso}: FASST falhou (exit {r.returncode})\n{r.stdout}\n{r.stderr}")


def para_csv(pacote, caso, tmp, saida, nome):
    with open(os.path.join(pacote, CASOS[caso], f"ref_{nome}.csv"), newline="") as f:
        cab = next(csv.reader(f))
    linhas, t0 = [], None
    for ln in open(os.path.join(tmp, f"saida.{nome}")):
        t = ln.split()
        if not (t and t[0].isdigit() and len(t[0]) == 4):
            continue                                  # cabecalhos
        if len(t) != len(cab):
            sys.exit(f"saida.{nome}: {len(t)} valores numa linha, esperado {len(cab)}:\n{ln}")
        s = int(t[1]) * 86400 + int(t[2]) * 3600 + int(t[3]) * 60
        t0 = s if t0 is None else t0
        if JD_INI <= int(t[1]) <= JD_FIM:
            linhas.append([str(s - t0)] + t)
    destino = os.path.join(saida, f"2026-11-12-grayling-{caso}-{nome}.csv")
    with open(destino, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["time_s"] + cab)
        w.writerows(linhas)
    return destino, len(linhas), len(cab) + 1


def main():
    pacote, tabelas, saida = sys.argv[1:4]
    os.makedirs(saida, exist_ok=True)
    for caso in CASOS:
        tmp = tempfile.mkdtemp(prefix=f"fasst_{caso}_")
        rodar(pacote, tabelas, caso, tmp)
        for nome in SAIDAS:
            p, n, c = para_csv(pacote, caso, tmp, saida, nome)
            sha = hashlib.sha256(open(p, "rb").read()).hexdigest()
            print(f"{os.path.basename(p)}  {n} linhas x {c} colunas  sha256={sha}")
        shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
