#' Distribuicao normal multivariada
#'
#' Cria um objeto representando `N_p(mu, Sigma)`, com metodos para a densidade,
#' transformacoes lineares (`A Y + b ~ N(A mu + b, A Sigma A')`) e marginais.
#'
#' @param mu Vetor de medias.
#' @param Sigma Matriz de covariancias (padrao identidade).
#' @return Objeto de classe `normal_mv`.
#' @export
#' @examples
#' Y <- normal_multivariada(c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2))
#' densidade(Y, c(0, 0))
#' transformar(Y, matrix(c(1, 1), 1))  # soma das componentes
normal_multivariada <- function(mu, Sigma = diag(length(mu))) {
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  if (nrow(Sigma) != length(mu)) {
    stop("Sigma deve ser quadrada de ordem length(mu).", call. = FALSE)
  }
  structure(list(mu = mu, Sigma = Sigma, dim = length(mu)), class = "normal_mv")
}

#' @export
print.normal_mv <- function(x, ...) {
  cat(sprintf("Normal multivariada N_%d(mu, Sigma)\n", x$dim))
  cat("mu:", round(x$mu, 4), "\n")
  cat("Sigma:\n")
  print(round(x$Sigma, 4))
  invisible(x)
}

#' Densidade da normal multivariada
#'
#' @param objeto Objeto `normal_mv`.
#' @param y Ponto (vetor) onde avaliar a densidade.
#' @return Valor da densidade em `y`.
#' @export
densidade <- function(objeto, y) {
  if (!inherits(objeto, "normal_mv")) {
    stop("Esperava um objeto 'normal_mv'.", call. = FALSE)
  }
  y <- as.numeric(y)
  d <- objeto$dim
  dif <- y - objeto$mu
  quad <- as.numeric(t(dif) %*% solve(objeto$Sigma) %*% dif)
  (2 * pi)^(-d / 2) * det(objeto$Sigma)^(-1 / 2) * exp(-quad / 2)
}

#' Transformacao linear de uma normal multivariada
#'
#' Devolve a distribuicao de `A Y + b`, que e `N(A mu + b, A Sigma A')`.
#'
#' @param objeto Objeto `normal_mv`.
#' @param A Matriz da transformacao.
#' @param b Vetor de deslocamento (padrao zeros).
#' @return Novo objeto `normal_mv`.
#' @export
transformar <- function(objeto, A, b = NULL) {
  if (!inherits(objeto, "normal_mv")) {
    stop("Esperava um objeto 'normal_mv'.", call. = FALSE)
  }
  A <- rbind(A)
  if (is.null(b)) b <- numeric(nrow(A))
  normal_multivariada(as.numeric(A %*% objeto$mu) + b, A %*% objeto$Sigma %*% t(A))
}

#' Marginal de uma normal multivariada
#'
#' @param objeto Objeto `normal_mv`.
#' @param indices Indices das componentes a manter.
#' @return Novo objeto `normal_mv`.
#' @export
marginal <- function(objeto, indices) {
  if (!inherits(objeto, "normal_mv")) {
    stop("Esperava um objeto 'normal_mv'.", call. = FALSE)
  }
  normal_multivariada(objeto$mu[indices], objeto$Sigma[indices, indices, drop = FALSE])
}

#' Distribuicao de uma funcao estimavel
#'
#' Para `lambda'beta` estimavel, `lambda'b ~ N(lambda'beta, sigma^2 lambda'(X'X)^-
#' lambda)`. Devolve a normal correspondente, com a media estimada e a variancia
#' estimada (usando `s^2 = QME`).
#'
#' @param lambda Vetor de coeficientes ou string de contraste.
#' @param objeto Objeto `llm`.
#' @return Objeto `normal_mv` unidimensional.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' distribuicao_estimavel("trat1 - trat3", llm(y ~ trat, dados))
distribuicao_estimavel <- function(lambda, objeto) {
  ml <- .como_llm(objeto)
  lambda <- .como_lambda(lambda, ml)
  if (!eh_estimavel(lambda, ml)) {
    warning("lambda nao e estimavel.", call. = FALSE)
  }
  est <- .estimativa_lambda(ml, lambda)
  v <- qme(ml)
  normal_multivariada(est$valor, matrix(v$S2 * est$var_unit, 1, 1))
}
