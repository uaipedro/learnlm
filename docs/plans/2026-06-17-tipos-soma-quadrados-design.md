# learnlm — tipos de soma de quadrados (I, II, III)

Estende `anava()` com `tipo = c("I","II","III")` e `explicar = FALSE`. Sem Tipo IV.
Cada tipo é a mesma máquina de SQ de hipótese (`SQ(L) = (Lb)'(LGL')⁻¹(Lb)`,
df = posto de `L`), mudando qual `L` se testa por termo.

## Construção por tipo

- **Tipo I (sequencial):** `SQ(T) = R(T | termos anteriores)`, via reduções `R(·)`.
  Depende da ordem da fórmula. É o que `anava()` já faz.
- **Tipo II (marginalidade):** `SQ(T) = R(T | termos que NÃO contêm T)`. Via
  `.reducao_modelo`. Independe da ordem. `.contem(T, o)` = todos os fatores de T
  estão em o.
- **Tipo III (parcial):** `L` montado de médias de célula equiponderadas. Para
  cada fator: matriz de contraste soma-zero `(l-1)×l` se o fator está no termo,
  ou linha de média `(1/l)·1ᵀ` se não está. `C_T` = produto de Kronecker dessas
  matrizes (na ordem das células de `expand.grid`); `L = C_T · M`, com `M` as
  funções estimáveis de média de célula (`.lambda_celula`). df = ∏(l_f−1) sobre
  os fatores do termo.

## Resíduo / Total

Iguais para os três tipos (vêm do ajuste do modelo completo). Em desbalanceado
as SQ dos termos II/III **não somam** a SQ do modelo — anoto isso na impressão.

## explicar = TRUE

Tabela + um bloco por termo (estilo `passos()`): tipo, hipótese em palavras,
conjunto de ajuste. Tipo III mostra "médias marginais equiponderadas iguais" e
se `L` é estimável.

## Células vazias / estimabilidade

`L` do Tipo III usa médias de célula que podem não ser estimáveis (fronteira do
Tipo IV). Reuso `eh_estimavel`/`.validar_L`: se não estimável, aviso claro e a
linha sai `NA`. Tipo II por reduções ainda calcula, com ressalva.

## Restrições

- Tipo III exige termos só de fatores (sem covariável numérica) — guardo com erro
  claro. Tipo I/II valem com numérico.
- `tipo=` (numerador) e `erros=` (denominador, split-plot) são eixos ortogonais e
  combinam.

## Validação

- Invariante interno: balanceado / sem interação ⇒ `I == II == III` (teste).
- Tipo II/III contra `car::Anova(lm, type=2/3)` num desbalanceado clássico.
- Tipo III contra médias marginais equiponderadas calculadas à mão.

## Arquivos

`anava.R` (tipo=, explicar=, `.contem`, `.termos_tipo1/2/3`, `.l_tipo3`,
`.contr_centro`), `tests/testthat/test-tipos-sq.R`, README, `man/*`,
`DESCRIPTION` → 0.1.3.
