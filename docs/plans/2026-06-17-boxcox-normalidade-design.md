# Box-Cox: diagnóstico de normalidade e transformação da resposta

## Motivação

O `llm()` ajusta o modelo sem checar o pressuposto de normalidade/estabilização
de variância. A ideia é, ao construir o modelo, rodar um Box-Cox com a mesma
fórmula e avisar quando λ = 1 não estiver no intervalo de máxima verossimilhança
(IC 95%), sugerindo a transformação normalizadora — e oferecer uma função para
aplicá-la.

## Decisões (brainstorming)

- **Integrado ao `llm()`**, não uma função `boxcox()` separada (o nome `boxcox`
  fica reservado para `MASS::boxcox`). Quando o pressuposto é violado, o `llm()`
  emite `warning()`.
- **Cálculo via `MASS::boxcox`** (não na mão).
- **`transformar_resposta()` sem `lambda`/`tipo` usa a transformação padrão
  sugerida** (o λ "redondo" mais próximo dentro do IC).

## Diagnóstico (dentro de `llm()`)

Novo parâmetro `checar_boxcox = TRUE`. Quando há `y` positivo:

1. `MASS::boxcox(formula, data = dados, lambda = seq(-2, 2, 0.01), plotit = FALSE)`.
2. `lambda_hat` = argmax do perfil.
3. IC 95% = `{λ : ℓ(λ) ≥ max(ℓ) − qchisq(0.95, 1)/2}`; `ic = range`.
4. `contem_1` = `1 ∈ ic`.
5. `lambda_sug`/`tipo_sug` = λ redondo (`-1, -0.5, 0, 0.5, 1, 2`) mais próximo de
   `lambda_hat` **dentro** do IC (se nenhum cair no IC, o mais próximo no geral).

Resultado guardado em `m$boxcox` (lista com `ok`, `motivo`, `lambda_hat`, `ic`,
`contem_1`, `lambda_sug`, `tipo_sug`). Se `contem_1 == FALSE`, `warning()`:

> Pressuposto Box-Cox: 1 ∉ IC95% (λ̂=…, IC=[…, …]). Sugestão: `log` (λ=0).
> Use transformar_resposta(m) para aplicar.

Sem `y`, com `y ≤ 0`, ou falha do MASS → guarda o motivo em `m$boxcox`, sem
warning. `print.llm` mostra uma linha do diagnóstico.

## Aplicação — `transformar_resposta(m, lambda = NULL, tipo = NULL)`

Precedência: `tipo` → `lambda` → `m$boxcox$tipo_sug` (padrão sugerida).

| `tipo`          | transformação | λ equiv. |
|-----------------|---------------|----------|
| `identidade`    | y             | 1        |
| `log`           | log(y)        | 0        |
| `raiz`          | √y            | 0.5      |
| `inversa_raiz`  | 1/√y          | −0.5     |
| `inversa`       | 1/y           | −1       |
| `quadrado`      | y²            | 2        |

- `tipo`: forma **literal** (log(y), √y, 1/y…).
- `lambda`: forma **Box-Cox** `(y^λ − 1)/λ` (log(y) se λ = 0).

Reajusta chamando `llm()` sobre `dados` com a coluna resposta substituída (mesmo
nome de variável). O novo objeto guarda `m$transformacao` (`tipo`, `lambda`,
`rotulo`) e o check do Box-Cox roda de novo no `y` transformado.

Erros: sem `y`; `y ≤ 0` para transformações que exigem positividade
(log/raiz/inversa/inversa_raiz/lambda); `tipo` desconhecido; padrão sem
diagnóstico disponível.

## Arquivos

- `R/boxcox.R` — `.diagnostico_boxcox()`, `transformar_resposta()`, helpers.
- `R/llm.R` — parâmetro `checar_boxcox`, chamada do diagnóstico, warning, print.
- `tests/testthat/test-boxcox.R`.
- `README` — seção de uso.
