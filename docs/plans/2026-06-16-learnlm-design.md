# learnlm — desenho do pacote

Pacote R didático para **modelos lineares de posto incompleto**, voltado a
graduação e pós de universidades federais. Material-base: curso GES122 (UFLA).

## Princípios

1. **Toda função interna é pública e standalone.** A mesma `inversa_condicional()`
   que `llm()` usa por dentro pode ser chamada numa matriz crua.
2. **A matriz X é a COMPLETA (super-parametrizada), nunca a restrita.** Resolve a
   dor de `model.matrix`, que já entrega X com casela de referência embutida.
   `matriz_modelo()` devolve X com uma coluna por nível (μ, cada bᵢ, cada τⱼ, cada
   interação) — posto incompleto de propósito.
3. **Cada cálculo carrega explicação e verificação.** Alimenta `passos(m)` e a
   exploração peça a peça (`explicar = TRUE` nas funções de álgebra).

## Decisões de design (brainstorming)

- Interação: **objeto inspecionável** — `llm()` devolve objeto S3 rico.
- Idioma: **tudo em português** (funções, argumentos, saídas).
- API em **snake_case** (decisão explícita do usuário; estilo moderno de pacotes R).
  Notação estatística preservada (`X`, `XtX`, `beta_hat`).
- Nome: pacote `learnlm`, função principal `llm()`.
- Stack: S3 puro; única dependência `MASS` (para `ginv` / Moore-Penrose).

## Escopo

### Fase 1 (v0.1 — esta entrega)
- `matriz_modelo()` — X completa a partir da fórmula, com descritores de coluna.
- `posto()`, `inversa_condicional()`, `inversa_minimos_quadrados()`,
  `inversa_moore_penrose()` (via `MASS::ginv`).
- `llm()` — monta X, X'X, X'y, posto, deficiência; objeto S3.
- `betas()` — solução do SEN sob `restricao = nenhuma | soma_zero |
  casela_referencia`, mostrando que Xβ é invariante.
- `eh_estimavel()`, `projetor()`.
- `passos()` — narra a dedução inteira; cada passo é recuperável.

### Fase 2 (depois)
Inferência: formas quadráticas, distribuição normal multivariada, IC e regiões.

### Fase 3 (depois)
ANAVA, DIC/DBC, fatorial, reduções de SQ (notação R(·)).

## Verificações que ancoram a corretude

- `A G A = A` para a inversa condicional; `AG` simétrica para a de mínimos quadrados.
- `P = X(X'X)⁻X'` simétrica e idempotente; invariante à escolha da g-inversa.
- `Xβ` idêntico para todas as restrições e g-inversas (lição central).
- Em `μ + τ` (one-way): contrastes estimáveis, parâmetros individuais não.
