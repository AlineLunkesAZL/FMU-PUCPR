# Baseline de referência: FASST Python, caso Grayling (MI, 1992)

Gerado pelo Python, validado contra o Fortran compilado.

Esta é a "resposta certa" para os testes de regressão da FMU. Os CSVs são a saída do FASST
em Python nos dois casos do pacote `referencia_fortran_grayling` (enviado pela orientadora).
Comparados com os `ref_*.csv` do FASST Fortran, coincidem em todas as colunas, com diferença
máxima 0,0 na precisão que o Fortran grava (648 linhas por arquivo e 7 776 no perfil).

## Ambiente

| item | valor |
|---|---|
| commit | `b448b449a6ee7049884e89ca886c0371e94b4701` |
| `git status --porcelain` antes de criar `tests/` | `?? docs/doc-transcricao-fasst_driver_V1.docx` (não afeta a execução) |
| Python | 3.14.10 |
| NumPy | 2.5.2 |
| sistema | Windows-10-10.0.19045-SP0 |
| pacotes | `2026-11-12-requirements-baseline.txt` (`pip freeze`) |

## Comando exato

Executado da raiz do repositório:

```
python -B tests/2026-11-12-gerar_baseline_grayling.py "C:/Users/gelen/Área de Trabalho/referencia_fortran_grayling" exemplo tests/reference
```

O script copia as entradas de cada caso para um diretório temporário e roda
`fasst.fasst_driver.main(['fasst_driver', 'fassti.inp'])` num processo novo. Depois
converte `saida.fst/.flx/.vgt/.grd` em CSV. Cada execução completa (A e B) leva cerca de 8 min.

## Arquivos de entrada

| arquivo | origem | sha256 |
|---|---|---|
| `grayling_1992.met` | pacote da orientadora | `31feeacd50a1356628e08e1890a86ead847b2c1da2b71791fbf03958c8d5f10c` |
| `FASST_defaultlayerprop.inp` | pacote | `0fee8545c1d53863ca849a998667dc3f71644a08d790e18291c5998531195d88` |
| `FASST_defaultvegprop.inp` | pacote | `e2e9f1b32c28eeb4faedd60675f084b6a8fe7712a1110322d16747366e3a34e5` |
| `A_solo_nu/fassti.inp` | pacote | `7d82b4fbdc04f87ca2c9cefae5ee2679a9763e86f90b5e8ba017c33b8ff077be` |
| `B_grama/fassti.inp` | pacote | `1b52bee75a219367ba33e646f74f0c750d519044a0672d767130b63383ad0477` |
| `MODIS_IGBP_vegprop.inp` | `exemplo/` deste repo | `d3331b5a207f24434ffdee5cd08d9964c3dd645373007edf7f35e1ba0312eb5d` |
| `UMD_LDAS_vegprop.inp` | `exemplo/` deste repo | `8aa4bd0c2293558c8e97b04d677b785b6238e07c5dc4218b6e2b368c2d1a2002` |
| `fasstusersoil.inp` | `exemplo/` deste repo (vazio) | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |

Os três últimos não vêm no pacote, mas o driver Python exige que existam. No caso A, o
resultado é idêntico com essas tabelas e com arquivos vazios. O pacote e `exemplo/` não estão
versionados no git; por isso os sha256 acima.

## Saídas

Período: dias 259 a 285 de 1992 (15/09 a 11/10), 648 horas, sem neve. É o mesmo período da
referência Fortran. A partir do dia 286 cai neve e o Python e o Fortran divergem.

Colunas: `time_s` (segundos desde o primeiro registro do met, 1992-259 00:00), seguida dos mesmos
nomes de coluna dos `ref_*.csv`. Assim os arquivos servem direto no `comparar_com_fortran.py`.
Os valores são os números impressos pelo FASST, sem reformatar.

| arquivo | linhas (sem cabeçalho) | colunas | sha256 |
|---|---|---|---|
| `2026-11-12-grayling-A-fst.csv` | 648 | 49 | `b70a1a4b670374261764e3a5e376ee9fd2750f7fb635ea0896940647026993c8` |
| `2026-11-12-grayling-A-flx.csv` | 648 | 18 | `fe8304f338b655d6f854c705a975451e913e4b265e8cb6b2a8ea041a86d4db1b` |
| `2026-11-12-grayling-A-vgt.csv` | 648 | 22 | `4e9d20eb06d61a2450ce3805140b0802a5861f39a3f23dfb1c2d9ccf9266b065` |
| `2026-11-12-grayling-A-grd.csv` | 7776 | 14 | `1c3c1e99154d4e97f99cee5d321abaad48ea2ff43cd45591aeaaa545284be4cb` |
| `2026-11-12-grayling-B-fst.csv` | 648 | 49 | `2b03c395f89f992d9c4fba4c3b1e11e37aa5d8bded0481fbd492e72874af0e5e` |
| `2026-11-12-grayling-B-flx.csv` | 648 | 18 | `c0daeb0a7a13495a5531bbb1a019a1abb6d091b35cf33759e7fa1898e1e0566f` |
| `2026-11-12-grayling-B-vgt.csv` | 648 | 22 | `545a132a246bbe22eb120c635ce05fdaaec8da0d61716f4c88622fd69d11a8de` |
| `2026-11-12-grayling-B-grd.csv` | 7776 | 14 | `5ff01211660b059727f5f02eedd5cc841c93621ed54da683f040be0a97393d5c` |

`fst` = superfície, `flx` = fluxos de energia, `vgt` = vegetação, `grd` = perfil do solo
(12 nós por hora, chave `N`).

## Determinismo

O comando foi rodado duas vezes, em diretórios de saída diferentes. Os 8 CSVs saíram idênticos
byte a byte (mesmo sha256).

## Limite de precisão

Os CSVs têm a precisão impressa pelo FASST: 2 a 5 casas decimais, conforme a coluna. Servem
para regressão com tolerância da ordem de 0,01 (a do `comparar_com_fortran.py`). Não servem
para conferir diferenças de ~1e-10.
