# learnlm — estratos de erro (split-plot / parcelas subdivididas)

Extensão para testar fontes e contrastes contra um **quadrado médio de erro
designado**, em vez do resíduo global. Cobre split-plot e delineamentos com mais
de um estrato de erro, mantendo a maquinaria de posto incompleto (g-inversas,
reduções R(·)) — válida para desbalanceado e contrastes, ao contrário de `aov`.

## Motivação

A decomposição de SQ do pacote já está correta para split-plot (cada SQ bate com
os estratos de `aov(... + Error())`). O único gap é que `anava()` e a inferência
sempre dividem pelo `QME` global. Falta poder designar o erro.

## Interface

Erro especificado por **nome de termo do modelo** (string).

```r
anava(m, erros = c(bloco = "bloco:irrig", irrig = "bloco:irrig"))

teste_t("irrig1 - irrig2", m, erro = "bloco:irrig")
teste_F(c("irrig1 - irrig2"), m, erro = "bloco:irrig")
intervalo_confianca("irrig1 - irrig2", m, erro = "bloco:irrig")
regiao_confianca(L, m, erro = "bloco:irrig")
```

Regras da tabela `anava`:

- `erros` é named character (fonte -> termo de erro). Default `NULL` = atual.
- Fontes não listadas usam o `Residuo` global.
- Termos que aparecem como erro de outra fonte viram estrato de erro puro
  (`F`/`pvalor` = `NA`).
- Coluna `erro` mostra, por fonte, contra qual QM foi testada.
- Ordem da fórmula é mantida; o QM de cada termo sai da decomposição sequencial
  R(·) — logo a ordem importa (Tipo I), como já é no pacote.

Compatibilidade: `erros`/`erro` default `NULL` não altera o uso atual.

## Internos

- `.decomposicao(ml)` — núcleo extraído de `anava()`; devolve por termo
  `GL/SQ/QM` e o resíduo (`SQE/glE/QME`).
- `.qm_erro(ml, fonte)` — `NULL`/`"Residuo"` -> `(QME, glE)` global; senão
  localiza o termo e devolve `(QM, gl)`. Erro claro se o nome não existir.
- Inferência: troca só o par `(s², gl)`. `var_unit = λ'(X'X)⁻λ` inalterado;
  `ep = sqrt(QM_erro · var_unit)`, quantis com `gl_erro`.

## Ressalva (documentada)

Reproduz exatamente os testes por estratos em delineamentos **balanceados**. Para
split-plot **desbalanceado** o rigoroso é modelo misto (REML, `lme4`); o
`λ'(X'X)⁻λ` com QM de erro designado é aproximação didática.

## Testes

`test-estratos.R`: split-plot balanceado vs. oráculo calculado à mão (não `aov`);
`teste_t(erro=)` com SE clássico de parcela; erro com termo inexistente;
`erros = NULL` idêntico ao atual.
