.lambda_celula <- function(ml, niveis_alvo) {
  lambda <- numeric(ml$p)
  for (j in seq_len(ml$p)) {
    d <- ml$descr[[j]]
    if (d$tipo == "intercepto") {
      lambda[j] <- 1
      next
    }
    if (d$tipo == "numerico") next
    fats <- names(d$niveis)
    casa <- all(fats %in% names(niveis_alvo)) &&
      all(vapply(fats, function(f) {
        identical(unname(d$niveis[f]), unname(niveis_alvo[f]))
      }, logical(1)))
    if (casa) lambda[j] <- 1
  }
  lambda
}

.sq_hipotese <- function(ml, L) {
  .validar_L(ml, L, "contrastes do desdobramento")
  b <- .solucao(ml)
  Lb <- as.numeric(L %*% b)
  W <- L %*% ml$G %*% t(L)
  list(SQ = as.numeric(t(Lb) %*% solve(W) %*% Lb), gl = nrow(L))
}

#' Desdobramento da interacao
#'
#' Estuda os efeitos de um fator dentro de cada nivel de outro, como manda o
#' material quando a interacao e significativa ("avaliar os efeitos de um fator
#' dentro de cada nivel do outro fator"). Para cada nivel de `dentro_de`, monta
#' os contrastes entre as medias de celula do `fator` e testa por forma
#' quadratica. Vale a identidade `soma_j SQ(fator | dentro_de = j) = SQ(fator) +
#' SQ(fator:dentro_de)`.
#'
#' @param objeto Objeto `llm` de um modelo com interacao entre os dois fatores.
#' @param fator Nome do fator cujos efeitos se estuda.
#' @param dentro_de Nome do fator cujos niveis sao fixados.
#' @return `data.frame` de classe `desdobramento_llm` com uma linha por nivel de
#'   `dentro_de`, mais `Residuo`.
#' @export
#' @examples
#' fat <- data.frame(
#'   tipo   = factor(rep(c("1", "2", "3"), each = 4)),
#'   metodo = factor(rep(rep(c("1", "2"), each = 2), times = 3)),
#'   y = c(39.02, 38.79, 38.96, 39.01, 35.74, 35.41,
#'         35.58, 35.52, 37.02, 36.00, 35.70, 36.04)
#' )
#' m <- llm(y ~ tipo + metodo + tipo:metodo, fat)
#' desdobramento(m, "tipo", dentro_de = "metodo")
desdobramento <- function(objeto, fator, dentro_de) {
  ml <- .como_llm(objeto)
  if (is.null(ml$y)) stop("O modelo nao tem resposta y.", call. = FALSE)
  niveis_f <- levels(as.factor(ml$dados[[fator]]))
  niveis_d <- levels(as.factor(ml$dados[[dentro_de]]))
  v <- qme(ml)

  linhas <- list()
  for (d in niveis_d) {
    celulas <- lapply(niveis_f, function(i) {
      alvo <- c(stats::setNames(i, fator), stats::setNames(d, dentro_de))
      .lambda_celula(ml, alvo)
    })
    L <- do.call(rbind, lapply(seq_len(length(celulas) - 1), function(k) {
      celulas[[k]] - celulas[[k + 1]]
    }))
    sq <- .sq_hipotese(ml, L)
    QM <- sq$SQ / sq$gl
    Fobs <- QM / v$S2
    linhas[[length(linhas) + 1L]] <- data.frame(
      GL = sq$gl, SQ = sq$SQ, QM = QM, F = Fobs,
      pvalor = stats::pf(Fobs, sq$gl, v$gl, lower.tail = FALSE),
      row.names = sprintf("%s / %s=%s", fator, dentro_de, d)
    )
  }
  res <- data.frame(GL = v$gl, SQ = v$SQE, QM = v$S2, F = NA, pvalor = NA,
                    row.names = "Residuo")
  tab <- rbind(do.call(rbind, linhas), res)
  attr(tab, "info") <- list(fator = fator, dentro_de = dentro_de)
  class(tab) <- c("desdobramento_llm", "anava_llm", "data.frame")
  tab
}

#' @export
print.desdobramento_llm <- function(x, digits = 4, ...) {
  info <- attr(x, "info")
  cat(sprintf("Desdobramento: %s dentro de cada nivel de %s\n\n",
              info$fator, info$dentro_de))
  .imprimir_tab_anava(x, digits)
  invisible(x)
}
