.indicadora_fator <- function(x, nome) {
  x <- as.factor(x)
  niveis <- levels(x)
  M <- vapply(niveis, function(l) as.integer(x == l), integer(length(x)))
  M <- matrix(M, nrow = length(x))
  colnames(M) <- paste0(nome, niveis)
  M
}

.colunas_termo <- function(rotulo, dados) {
  vars <- strsplit(rotulo, ":", fixed = TRUE)[[1]]
  infos <- lapply(vars, function(v) {
    x <- dados[[v]]
    if (is.numeric(x)) {
      list(v = v, tipo = "numerico", niveis = NA_character_,
           ind = matrix(x, ncol = 1, dimnames = list(NULL, v)))
    } else {
      x <- as.factor(x)
      list(v = v, tipo = "fator", niveis = levels(x),
           ind = .indicadora_fator(x, v))
    }
  })

  grade <- expand.grid(lapply(infos, function(i) seq_len(ncol(i$ind))))
  n <- nrow(dados)
  M <- matrix(1, nrow = n, ncol = nrow(grade))
  nomes <- character(nrow(grade))
  descr <- vector("list", nrow(grade))

  for (k in seq_len(nrow(grade))) {
    col <- rep(1, n)
    niveis_k <- character(0)
    rotulos_k <- character(0)
    for (j in seq_along(infos)) {
      cj <- grade[k, j]
      col <- col * infos[[j]]$ind[, cj]
      if (infos[[j]]$tipo == "fator") {
        niv <- infos[[j]]$niveis[cj]
        niveis_k[infos[[j]]$v] <- niv
        rotulos_k <- c(rotulos_k, paste0(infos[[j]]$v, niv))
      } else {
        rotulos_k <- c(rotulos_k, infos[[j]]$v)
      }
    }
    M[, k] <- col
    nomes[k] <- paste(rotulos_k, collapse = ":")
    descr[[k]] <- list(
      termo = rotulo,
      niveis = niveis_k,
      tipo = if (length(infos) > 1) "interacao" else infos[[1]]$tipo
    )
  }
  colnames(M) <- nomes
  list(M = M, descr = descr)
}

#' Matriz de delineamento completa (super-parametrizada)
#'
#' Constroi a matriz X com **uma coluna por nivel** de cada fator (mais
#' interacoes), sem casela de referencia e sem restricoes. Diferente de
#' `model.matrix`, que devolve a matriz ja restrita. Em modelos de posto
#' incompleto, esta e a X que aparece no material teorico.
#'
#' @param formula Formula do modelo, p.ex. `y ~ bloco + trat + bloco:trat`.
#' @param dados `data.frame` com as variaveis da formula.
#' @param intercepto Se `TRUE` (padrao), inclui a coluna `mu` (intercepto).
#' @return Matriz X de classe `matriz_modelo`, com o atributo `descr`
#'   (descritores de cada coluna: termo, niveis, tipo).
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 2)),
#'   y = c(10, 12, 20, 22, 30, 28)
#' )
#' matriz_modelo(y ~ trat, dados)
matriz_modelo <- function(formula, dados, intercepto = TRUE) {
  tt <- stats::terms(formula)
  rotulos <- attr(tt, "term.labels")
  n <- nrow(dados)

  if (length(rotulos) == 0) {
    X <- matrix(numeric(0), nrow = n, ncol = 0)
    descr <- list()
  } else {
    partes <- lapply(rotulos, .colunas_termo, dados = dados)
    X <- do.call(cbind, lapply(partes, `[[`, "M"))
    descr <- do.call(c, lapply(partes, `[[`, "descr"))
  }

  if (intercepto) {
    X <- cbind(mu = rep(1, n), X)
    descr <- c(
      list(list(termo = "(intercepto)", niveis = character(0), tipo = "intercepto")),
      descr
    )
  }
  rownames(X) <- NULL
  attr(X, "descr") <- descr
  attr(X, "formula") <- formula
  class(X) <- c("matriz_modelo", "matrix", "array")
  X
}

#' @export
print.matriz_modelo <- function(x, ...) {
  r <- posto(x)
  p <- ncol(x)
  cat(sprintf("Matriz de delineamento %d x %d\n", nrow(x), p))
  cat(sprintf("posto = %d, %s (deficiencia = %d)\n", r,
              if (r < p) "posto INCOMPLETO" else "posto completo", p - r))
  y <- unclass(x)
  attr(y, "descr") <- NULL
  attr(y, "formula") <- NULL
  print(y, ...)
  invisible(x)
}

.resposta <- function(formula, dados) {
  tt <- stats::terms(formula)
  if (attr(tt, "response") == 0) return(NULL)
  yname <- all.vars(formula)[1]
  as.numeric(dados[[yname]])
}

.niveis_referencia <- function(formula, dados) {
  vars <- unique(unlist(strsplit(attr(stats::terms(formula), "term.labels"),
                                 ":", fixed = TRUE)))
  ref <- list()
  for (v in vars) {
    x <- dados[[v]]
    if (!is.numeric(x)) ref[[v]] <- levels(as.factor(x))[1]
  }
  ref
}
