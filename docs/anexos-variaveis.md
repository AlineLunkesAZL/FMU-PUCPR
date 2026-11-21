# Anexos A e B: variáveis da FMU e correspondência com o FASST

Os nomes da FMU são os de `greenroof_fmu.py` (`ENTRADAS`, `PARAMETROS`, `SAIDAS`) e os
argumentos de `fasst_adaptador.py`. As colunas do FASST são as dos `ref_*.csv` do Fortran e
dos CSVs de `tests/reference/`. É isso que `tests/comparar_baseline_grayling.py` usa.

**Estado: provisório.** As variáveis abaixo são as do modelo provisório (`ModeloSubstituto`).
Quando o adaptador passar a chamar a `FasstPasso` (#3), esta lista será revista. Os parâmetros
reais do FASST (solo, van Genuchten, vegetação, número de camadas) entram aqui e no
`greenroof_fmu.py`.

## Anexo A: entradas e parâmetros

### Entradas (a cada passo; no acoplamento vêm do DOMUS)

| FMU | Unidade | FASST (`.met`) | Unidade no `.met` | Conversão |
|---|---|---|---|---|
| `T_ar` | K | `tmp` (col. 6) | °C | `tmp + 273,15` |
| `UR` | 0 a 1 | `rh` (col. 7) | % | `rh / 100` |
| `G_sol` | W/m² | `tsol` (col. 23) | W/m² | nenhuma |
| `L_ceu` | W/m² | `ir` (col. 27) | W/m² | nenhuma |
| `vento` | m/s | `ws` (col. 8) | m/s | nenhuma |
| `chuva` | kg/(m²·s) | `prec` (col. 10) | mm por passo | `prec / dt` (o adaptador volta com `chuva · dt`) |
| `T_laje` | K | não existe | — | só com `contorno="laje"`. Na regressão (`"fasst_original"`) o fundo é fluxo fixo de 0,15 W/m² + lençol freático, e `T_laje` é ignorada |

O passo k usa a linha k do `.met`: a meteorologia do **fim** do passo, como no FASST.

### Parâmetros (fixos na inicialização)

| FMU | Unidade | FASST | Caso A (solo nu) | Caso B (grama) |
|---|---|---|---|---|
| `LAI` | — | `vegl_type` no `fassti.inp` → colunas 6–7 (LAI máx./mín.) de `FASST_defaultvegprop.inp` | 0 (`vegl_type = 0`) | tipo 2 (grama curta): máx. 2,0, mín. 0,5 |
| `h_sub` | m | espessura da camada no `fassti.inp` (linha `6 2.0`: USCS SP, 2 m) | 2,0 | 2,0 |
| `theta0` | m³/m³ | umidades iniciais do `fassti.inp` | `n`: o FASST calcula | `n`: o FASST calcula |
| `T0` | K | temperaturas iniciais do `fassti.inp` | `n`: o FASST calcula | `n`: o FASST calcula |

Como os dois casos usam `n` (sem valores iniciais), o teste passa ao modelo provisório
`theta0 = 0,12` e `T0 = 291,45 K` (a temperatura do ar no primeiro registro). A `FasstPasso`
deve calcular esses valores como o Fortran.

## Anexo B: saídas

| FMU | Unidade | FASST (coluna da baseline) | Unidade | Correspondência | Casos |
|---|---|---|---|---|---|
| `T_sub` | K | `fst.Grtemp` | K | igual | A, B |
| `theta` | m³/m³ | `fst.Grmoistv` | m³/m³ | igual | A, B |
| `LE` | W/m² | `flx.LHeat` | W/m² | `LHeat = −LE` (no `flx`, positivo = para a superfície; `LE > 0` = evaporação) | A, B |
| `T_folha` | K | `vgt.LVegTmp` | K | igual | B (em A vale 9999, sem vegetação) |
| `q_laje` | W/m² | não existe | — | sem comparação: no FASST original o fundo é fixo em 0,15 W/m² | — |
| `escoamento` | kg/(m²·s) | não existe nos CSVs | — | sem comparação | — |

O perfil do solo (`grd`: `Grtemp_K`, `Water`, `Ice` por nó, 12 nós por hora) não é
comparado, porque o adaptador não devolve o perfil. Se a `FasstPasso` passar a expor o perfil,
ele entra aqui e no `MAPA` do teste.

Fonte única no código: `MAPA` e `CASOS` em `tests/comparar_baseline_grayling.py`. Se esta
tabela mudar, aquele dicionário muda junto.
