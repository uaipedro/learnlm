.reducao_modelo <- function(rotulos, resposta, dados, intercepto = TRUE) {
  formula <- if (length(rotulos) == 0) {
    stats::reformulate("1", response = resposta)
  } else {
    stats::reformulate(rotulos, response = resposta)
  }
  X <- matriz_modelo(formula, dados, intercepto)
  if (ncol(X) == 0) return(list(R = 0, posto = 0))
  y <- .resposta(formula, dados)
  XtX <- t(X) %*% X
  Xty <- as.numeric(t(X) %*% y)
  b0 <- as.numeric(inversa_moore_penrose(XtX) %*% Xty)
  list(R = sum(b0 * Xty), posto = posto(X))
}

#' Reducao de soma de quadrados R(.)
#'
#' Calcula `R(mu, termos) = b0' X'y` para o modelo com os termos indicados. Com
#' `ajustado_por`, devolve a reducao parcial `R(termos | mu, ajustado_por) =
#' R(mu, ajustado_por, termos) - R(mu, ajustado_por)`, na notacao do curso.
#'
#' @param objeto Objeto `llm`.
#' @param termos Vetor de nomes de termos (p.ex. `"trat"` ou `c("a", "a:b")`).
#'   `NULL` (padrao) usa todos os termos do modelo.
#' @param ajustado_por Vetor de termos pelos quais ajustar (entram antes).
#' @return Lista com `R` (valor da reducao) e `gl` (graus de liberdade).
#' @export
#' @examples
#' dados <- data.frame(
#'   bloco = factor(rep(1:3, times = 4)),
#'   trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
#'   y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
#' )
#' m <- llm(y ~ bloco + trat, dados)
#' reducao(m, "trat", ajustado_por = "bloco")  # R(tau | mu, b)
reducao <- function(objeto, termos = NULL, ajustado_por = NULL) {
  ml <- .como_llm(objeto)
  resp <- all.vars(ml$formula)[1]
  if (is.null(termos)) termos <- attr(stats::terms(ml$formula), "term.labels")

  if (is.null(ajustado_por)) {
    r <- .reducao_modelo(termos, resp, ml$dados, ml$intercepto)
    return(list(R = r$R, gl = r$posto))
  }
  completo <- .reducao_modelo(union(ajustado_por, termos), resp, ml$dados, ml$intercepto)
  base <- .reducao_modelo(ajustado_por, resp, ml$dados, ml$intercepto)
  list(R = completo$R - base$R, gl = completo$posto - base$posto)
}

.decomposicao <- function(ml) {
  resp <- all.vars(ml$formula)[1]
  rotulos <- attr(stats::terms(ml$formula), "term.labels")
  if (length(rotulos) == 0) stop("O modelo nao tem termos para a ANAVA.", call. = FALSE)

  Rmu <- .reducao_modelo(character(0), resp, ml$dados, ml$intercepto)
  acum <- Rmu
  SQ <- numeric(length(rotulos))
  GL <- numeric(length(rotulos))
  for (k in seq_along(rotulos)) {
    atual <- .reducao_modelo(rotulos[seq_len(k)], resp, ml$dados, ml$intercepto)
    SQ[k] <- atual$R - acum$R
    GL[k] <- atual$posto - acum$posto
    acum <- atual
  }
  SQE <- sum(ml$y^2) - acum$R
  glE <- ml$n - ml$posto
  list(rotulos = rotulos, SQ = SQ, GL = GL, QM = SQ / GL,
       SQE = SQE, glE = glE, QME = SQE / glE,
       SQTotal = sum(ml$y^2) - Rmu$R)
}

.qm_erro <- function(ml, fonte) {
  if (is.null(fonte) || identical(fonte, "Residuo")) {
    v <- qme(ml)
    return(list(QM = v$QME, gl = v$gl, nome = "Residuo"))
  }
  dec <- .decomposicao(ml)
  i <- match(fonte, dec$rotulos)
  if (is.na(i)) {
    stop(sprintf("Termo de erro '%s' nao existe no modelo. Termos: %s.",
                 fonte, paste(dec$rotulos, collapse = ", ")), call. = FALSE)
  }
  list(QM = dec$QM[i], gl = dec$GL[i], nome = fonte)
}

#' Tabela da analise de variancia (ANAVA)
#'
#' Monta a ANAVA sequencial (cada termo ajustado pelos anteriores, na ordem da
#' formula), com as somas de quadrados obtidas por reducoes `R(.)`. Equivale a
#' decomposicao do tipo I e reproduz `anova(lm(...))`. Serve para DIC, DBC e
#' fatoriais. Com `erros`, cada fonte e testada contra um estrato de erro
#' designado (split-plot / parcelas subdivididas).
#'
#' @param objeto Objeto `llm`.
#' @param erros Vetor nomeado opcional `fonte = termo_de_erro` (p.ex.
#'   `c(irrig = "bloco:irrig")`). Fontes nao listadas usam o `Residuo` global.
#'   Termos usados como erro viram estrato de erro puro (sem `F`).
#' @return `data.frame` de classe `anava_llm` com colunas `GL`, `SQ`, `QM`, `F`,
#'   `pvalor`, uma linha por termo mais `Residuo` e `Total`. Com `erros`, ganha a
#'   coluna `erro` (contra qual QM cada fonte foi testada).
#' @details Os estratos de erro reproduzem exatamente os testes classicos em
#'   delineamentos *balanceados*. Para split-plot *desbalanceado* o rigoroso e
#'   modelo misto (REML, `lme4`); aqui e uma aproximacao didatica.
#' @export
#' @examples
#' dados <- data.frame(
#'   bloco = factor(rep(1:3, times = 4)),
#'   trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
#'   y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
#' )
#' anava(llm(y ~ bloco + trat, dados))
anava <- function(objeto, erros = NULL) {
  ml <- .como_llm(objeto)
  if (is.null(ml$y)) stop("O modelo nao tem resposta y.", call. = FALSE)
  dec <- .decomposicao(ml)
  rotulos <- dec$rotulos

  if (!is.null(erros)) {
    bad_src <- setdiff(names(erros), rotulos)
    if (length(bad_src)) {
      stop(sprintf("Fonte(s) inexistente(s) em 'erros': %s.",
                   paste(bad_src, collapse = ", ")), call. = FALSE)
    }
    bad_err <- setdiff(unname(erros), c(rotulos, "Residuo"))
    if (length(bad_err)) {
      stop(sprintf("Termo(s) de erro inexistente(s): %s.",
                   paste(bad_err, collapse = ", ")), call. = FALSE)
    }
  }
  termos_erro <- intersect(unname(erros), rotulos)

  Fobs <- pvalor <- numeric(length(rotulos))
  col_erro <- character(length(rotulos))
  for (k in seq_along(rotulos)) {
    if (rotulos[k] %in% termos_erro) {
      Fobs[k] <- NA
      pvalor[k] <- NA
      col_erro[k] <- "(erro)"
    } else {
      en <- if (!is.null(erros) && rotulos[k] %in% names(erros)) erros[[rotulos[k]]] else "Residuo"
      qe <- .qm_erro(ml, en)
      Fobs[k] <- dec$QM[k] / qe$QM
      pvalor[k] <- stats::pf(Fobs[k], dec$GL[k], qe$gl, lower.tail = FALSE)
      col_erro[k] <- qe$nome
    }
  }

  tab <- data.frame(
    GL = c(dec$GL, dec$glE, ml$n - 1),
    SQ = c(dec$SQ, dec$SQE, dec$SQTotal),
    QM = c(dec$QM, dec$QME, NA),
    F = c(Fobs, NA, NA),
    pvalor = c(pvalor, NA, NA),
    row.names = c(rotulos, "Residuo", "Total")
  )
  if (!is.null(erros)) tab$erro <- c(col_erro, "", "")
  class(tab) <- c("anava_llm", "data.frame")
  tab
}

.imprimir_tab_anava <- function(x, digits) {
  y <- as.data.frame(x)
  y$SQ <- round(y$SQ, digits)
  y$QM <- round(y$QM, digits)
  y$F <- round(y$F, digits)
  y$pvalor <- format.pval(y$pvalor, digits = digits, na.form = "")
  y$F[is.na(x$F)] <- ""
  y$QM[is.na(x$QM)] <- ""
  print(y)
}

#' @export
print.anava_llm <- function(x, digits = 4, ...) {
  cat("Analise de Variancia (ANAVA sequencial, SQ por reducoes)\n\n")
  .imprimir_tab_anava(x, digits)
  invisible(x)
}
