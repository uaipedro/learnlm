#' Ajusta e dichava um modelo linear passo a passo
#'
#' Monta a matriz de delineamento completa, o sistema de equacoes normais
#' (`X'X b = X'y`), calcula posto e deficiencia e registra a deducao em
#' [passos()]. O objeto resultante expoe cada peca (`m$X`, `m$XtX`, `m$Xty`, ...)
#' e serve de entrada para [betas()], [projetor()] e [eh_estimavel()].
#'
#' @param formula Formula do modelo, p.ex. `y ~ bloco + trat`.
#' @param dados `data.frame` com as variaveis.
#' @param intercepto Se `TRUE` (padrao), inclui o intercepto `mu`.
#' @return Objeto de classe `llm`.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 2)),
#'   y = c(10, 12, 20, 22, 30, 28)
#' )
#' m <- llm(y ~ trat, dados)
#' m$X
#' betas(m, restricao = "soma_zero")
llm <- function(formula, dados, intercepto = TRUE) {
  X <- matriz_modelo(formula, dados, intercepto)
  y <- .resposta(formula, dados)
  XtX <- t(X) %*% X
  Xty <- if (!is.null(y)) as.numeric(t(X) %*% y) else NULL
  r <- posto(X)
  p <- ncol(X)

  obj <- list(
    formula = formula,
    dados = dados,
    X = X,
    y = y,
    XtX = XtX,
    Xty = Xty,
    posto = r,
    p = p,
    n = nrow(X),
    deficiencia = p - r,
    descr = attr(X, "descr"),
    niveis_ref = .niveis_referencia(formula, dados)
  )
  class(obj) <- "llm"
  obj$passos <- .construir_passos(obj)
  obj
}

#' @export
print.llm <- function(x, ...) {
  cat("Modelo linear (learnlm)\n")
  cat("Formula: "); print(x$formula)
  cat(sprintf("Observacoes: %d   Parametros: %d   posto(X): %d\n",
              x$n, x$p, x$posto))
  if (x$deficiencia > 0) {
    cat(sprintf("Posto INCOMPLETO: deficiencia = %d.\n", x$deficiencia))
    cat("=> Os parametros individuais nao sao estimaveis; ha infinitas solucoes b.\n")
  } else {
    cat("Posto completo: solucao unica.\n")
  }
  cat("\nUse passos(m) para ver a deducao, ou m$X, m$XtX, m$Xty para as pecas.\n")
  cat("betas(m, restricao = ...), projetor(m), eh_estimavel(lambda, m).\n")
  invisible(x)
}

.como_llm <- function(objeto) {
  if (inherits(objeto, "llm")) return(objeto)
  stop("Esperava um objeto 'llm' (criado por llm()).", call. = FALSE)
}

.pegar_X <- function(objeto) {
  if (inherits(objeto, "llm")) return(objeto$X)
  as.matrix(objeto)
}
