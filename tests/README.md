# Testes

## Comparação contra a baseline do Grayling (issue #2)

Da raiz do repositório:

```
python -B tests/comparar_baseline_grayling.py
```

O script roda `FasstAdaptador.passo` (`fasst_adaptador.py`) hora a hora com
`inputs/grayling/grayling_1992.met`, nos casos A (solo nu) e B (grama), durante as 648 h da
baseline. A saída vai para `tests/saida/grayling-<caso>-<fst|flx|vgt>.csv`, com os nomes de
coluna de `reference/`, e é comparada com a baseline com tolerância absoluta de 0,01. Se algum
caso diverge, o script sai com código 1 e mostra o primeiro instante e a primeira coluna que
divergem.

`--controle` compara a baseline com ela mesma, em todas as colunas, e tem que passar. Serve
para conferir o comparador.

O GitHub Actions roda os dois comandos a cada push e PR
(`.github/workflows/comparacao-baseline.yml`).

**Estado atual.** O adaptador ainda usa o `ModeloSubstituto`, então o teste falha de
propósito logo na linha 0. Ele passa a valer quando o adaptador chamar a `FasstPasso`
(issue #3).

Colunas comparadas hoje (lista completa em `docs/anexos-variaveis.md`, Anexos A e B):

| saída do adaptador | coluna da baseline | casos |
|---|---|---|
| `T_sub` [K] | `fst.Grtemp` | A, B |
| `theta` [m³/m³] | `fst.Grmoistv` | A, B |
| `-LE` [W/m²] | `flx.LHeat` (positivo = para a superfície) | A, B |
| `T_folha` [K] | `vgt.LVegTmp` | B (em A vale 9999, sem vegetação) |

`q_laje` e `escoamento` não têm coluna na baseline, e `grd` (perfil por nó) não é comparado.

Entradas (`.met` → adaptador): `tmp` °C → `T_ar` K; `rh` % → `UR` (0 a 1); `tsol` → `G_sol`;
`ir` → `L_ceu`; `ws` → `vento`; `prec` mm/passo ÷ 3600 → `chuva` kg/(m²·s). O passo k usa
a linha k do `.met`, que é a meteorologia do fim do passo, sem deslocar a série.
