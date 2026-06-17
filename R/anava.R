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

.contem <- function(termo, outro) {
  ft <- strsplit(termo, ":", fixed = TRUE)[[1]]
  fo <- strsplit(outro, ":", fixed = TRUE)[[1]]
  all(ft %in% fo)
}

.fatores_modelo <- function(ml) {
  rot <- attr(stats::terms(ml$formula), "term.labels")
  unique(unlist(strsplit(rot, ":", fixed = TRUE)))
}

.contr_centro <- function(l) {
  C <- matrix(0, l - 1, l)
  for (i in seq_len(l - 1)) {
    C[i, i] <- 1
    C[i, i + 1] <- -1
  }
  C
}

# L do Tipo III: contrastes das medias de celula equiponderadas.
.l_tipo3 <- function(ml, termo) {
  fs <- .fatores_modelo(ml)
  if (any(vapply(fs, function(f) is.numeric(ml$dados[[f]]), logical(1)))) {
    stop("Tipo III exige termos apenas de fatores (sem covariavel numerica). ",
         "Use tipo = 'I' ou 'II'.", call. = FALSE)
  }
  niveis <- lapply(fs, function(f) levels(as.factor(ml$dados[[f]])))
  names(niveis) <- fs
  fT <- strsplit(termo, ":", fixed = TRUE)[[1]]

  grade <- expand.grid(niveis, stringsAsFactors = FALSE)
  M <- t(apply(grade, 1, function(linha) {
    alvo <- stats::setNames(as.character(unlist(linha)), fs)
    .lambda_celula(ml, alvo)
  }))

  mats <- lapply(fs, function(f) {
    l <- length(niveis[[f]])
    if (f %in% fT) .contr_centro(l) else matrix(1 / l, 1, l)
  })
  C <- Reduce(kronecker, rev(mats))
  C %*% M
}

.termos_sq <- function(ml, tipo) {
  resp <- all.vars(ml$formula)[1]
  rotulos <- attr(stats::terms(ml$formula), "term.labels")
  if (length(rotulos) == 0) stop("O modelo nao tem termos para a ANAVA.", call. = FALSE)
  SQ <- GL <- numeric(length(rotulos))
  hip <- character(length(rotulos))
  estim <- rep(TRUE, length(rotulos))

  if (tipo == "I") {
    acum <- .reducao_modelo(character(0), resp, ml$dados, ml$intercepto)
    for (k in seq_along(rotulos)) {
      atual <- .reducao_modelo(rotulos[seq_len(k)], resp, ml$dados, ml$intercepto)
      SQ[k] <- atual$R - acum$R
      GL[k] <- atual$posto - acum$posto
      ant <- if (k == 1) "mu" else paste(c("mu", rotulos[seq_len(k - 1)]), collapse = ", ")
      hip[k] <- sprintf("R(%s | %s)  [ajustado pelos anteriores]", rotulos[k], ant)
      acum <- atual
    }
  } else if (tipo == "II") {
    for (k in seq_along(rotulos)) {
      tm <- rotulos[k]
      U <- rotulos[vapply(rotulos, function(o) o != tm && !.contem(tm, o), logical(1))]
      base <- .reducao_modelo(U, resp, ml$dados, ml$intercepto)
      comT <- .reducao_modelo(union(U, tm), resp, ml$dados, ml$intercepto)
      SQ[k] <- comT$R - base$R
      GL[k] <- comT$posto - base$posto
      hip[k] <- sprintf("R(%s | %s)  [ajustado por quem NAO contem %s]",
                        tm, paste(c("mu", U), collapse = ", "), tm)
    }
  } else {
    nao_estim <- character(0)
    for (k in seq_along(rotulos)) {
      tm <- rotulos[k]
      L <- .l_tipo3(ml, tm)
      ok <- all(vapply(seq_len(nrow(L)),
                       function(i) isTRUE(eh_estimavel(L[i, ], ml)), logical(1)))
      estim[k] <- ok
      GL[k] <- nrow(L)
      outros <- setdiff(rotulos, tm)
      hip[k] <- sprintf("medias marginais equiponderadas de %s iguais  [ajustado por %s]",
                        tm, if (length(outros)) paste(outros, collapse = ", ") else "nada")
      if (ok) {
        SQ[k] <- .sq_hipotese(ml, L)$SQ
      } else {
        SQ[k] <- NA
        nao_estim <- c(nao_estim, tm)
      }
    }
    if (length(nao_estim) > 0) {
      warning(sprintf("Tipo III: funcoes nao estimaveis (celula vazia?) em: %s; SQ = NA.",
                      paste(nao_estim, collapse = ", ")), call. = FALSE)
    }
  }
  list(rotulos = rotulos, SQ = SQ, GL = GL, QM = SQ / GL, hip = hip, estim = estim)
}

#' Tabela da analise de variancia (ANAVA)
#'
#' Monta a tabela de ANAVA com as somas de quadrados de um dos tres tipos
#' (notacao SAS), todas pela mesma maquina de SQ de hipotese do pacote:
#' \describe{
#'   \item{`"I"`}{sequencial: `R(termo | anteriores)`; depende da ordem; reproduz
#'     `anova(lm(...))`.}
#'   \item{`"II"`}{marginalidade: `R(termo | termos que NAO o contem)`; nao
#'     depende da ordem.}
#'   \item{`"III"`}{parcial: medias marginais equiponderadas (contra todos os
#'     outros termos, inclusive interacoes); exige termos so de fatores.}
#' }
#' Com `erros`, cada fonte e testada contra um estrato de erro designado
#' (split-plot). `tipo` (numerador) e `erros` (denominador) sao independentes.
#'
#' @param objeto Objeto `llm`.
#' @param tipo Tipo de soma de quadrados: `"I"` (padrao), `"II"` ou `"III"`.
#' @param erros Vetor nomeado opcional `fonte = termo_de_erro` (p.ex.
#'   `c(irrig = "bloco:irrig")`). Fontes nao listadas usam o `Residuo` global.
#'   Termos usados como erro viram estrato de erro puro (sem `F`).
#' @param explicar Se `TRUE`, imprime a tabela e a hipotese testada por fonte, e
#'   devolve a tabela invisivelmente.
#' @return `data.frame` de classe `anava_llm` com colunas `GL`, `SQ`, `QM`, `F`,
#'   `pvalor`, uma linha por termo mais `Residuo` e `Total`. Com `erros`, ganha a
#'   coluna `erro` (contra qual QM cada fonte foi testada).
#' @details No Tipo I as SQ dos termos somam a SQ do modelo; nos Tipos II/III, em
#'   dados desbalanceados, nao somam (sao ajustes parciais). Em dados balanceados
#'   sem interacao os tres tipos coincidem. Os estratos de erro reproduzem os
#'   testes classicos em delineamentos *balanceados*; para split-plot
#'   *desbalanceado* o rigoroso e modelo misto (REML, `lme4`).
#' @export
#' @examples
#' dados <- data.frame(
#'   bloco = factor(rep(1:3, times = 4)),
#'   trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
#'   y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
#' )
#' m <- llm(y ~ bloco + trat, dados)
#' anava(m)
#' anava(m, tipo = "III", explicar = TRUE)
anava <- function(objeto, tipo = c("I", "II", "III"), erros = NULL,
                  explicar = FALSE) {
  ml <- .como_llm(objeto)
  if (is.null(ml$y)) stop("O modelo nao tem resposta y.", call. = FALSE)
  tipo <- match.arg(tipo)
  dec <- .decomposicao(ml)
  tt <- .termos_sq(ml, tipo)
  rotulos <- tt$rotulos

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
      Fobs[k] <- tt$QM[k] / qe$QM
      pvalor[k] <- stats::pf(Fobs[k], tt$GL[k], qe$gl, lower.tail = FALSE)
      col_erro[k] <- qe$nome
    }
  }

  tab <- data.frame(
    GL = c(tt$GL, dec$glE, ml$n - 1),
    SQ = c(tt$SQ, dec$SQE, dec$SQTotal),
    QM = c(tt$QM, dec$QME, NA),
    F = c(Fobs, NA, NA),
    pvalor = c(pvalor, NA, NA),
    row.names = c(rotulos, "Residuo", "Total")
  )
  if (!is.null(erros)) tab$erro <- c(col_erro, "", "")
  attr(tab, "tipo") <- tipo
  attr(tab, "hipoteses") <- stats::setNames(tt$hip, rotulos)
  class(tab) <- c("anava_llm", "data.frame")

  if (explicar) {
    print(tab)
    .explicar_anava(tab)
    return(invisible(tab))
  }
  tab
}

.explicar_anava <- function(tab) {
  hip <- attr(tab, "hipoteses")
  if (is.null(hip)) return(invisible(NULL))
  cat(sprintf("\nHipotese de cada fonte (Tipo %s):\n", attr(tab, "tipo")))
  for (tm in names(hip)) {
    cat(sprintf("  %-16s %s\n", tm, hip[[tm]]))
  }
  invisible(NULL)
}

.imprimir_tab_anava <- function(x, digits) {
  y <- as.data.frame(x)
  y$SQ <- round(y$SQ, digits)
  y$QM <- round(y$QM, digits)
  y$F <- round(y$F, digits)
  y$pvalor <- format.pval(y$pvalor, digits = digits, na.form = "")
  y$F[is.na(x$F)] <- ""
  y$QM[is.na(x$QM)] <- ""
  y$SQ[is.na(x$SQ)] <- ""
  print(y)
}

#' @export
print.anava_llm <- function(x, digits = 4, ...) {
  tipo <- attr(x, "tipo")
  if (is.null(tipo)) {
    cat("Analise de Variancia (ANAVA sequencial, SQ por reducoes)\n\n")
  } else {
    cat(sprintf("Analise de Variancia (SQ Tipo %s)\n\n", tipo))
  }
  .imprimir_tab_anava(x, digits)
  invisible(x)
}
