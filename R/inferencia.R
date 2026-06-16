#' Estimativa da variancia residual (QME)
#'
#' Calcula a soma de quadrados do residuo `SQE = y'(I - P) y`, seus graus de
#' liberdade `n - posto(X)` e o estimador `s^2 = QME = SQE / (n - posto(X))`,
#' que e nao-viesado para `sigma^2`.
#'
#' @param objeto Objeto `llm`.
#' @return Lista com `SQE`, `gl`, `S2` (= `QME`).
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' qme(llm(y ~ trat, dados))
qme <- function(objeto) {
  ml <- .como_llm(objeto)
  if (is.null(ml$y)) stop("O modelo nao tem resposta y.", call. = FALSE)
  M <- projetor(ml, residuos = TRUE)
  SQE <- as.numeric(t(ml$y) %*% M %*% ml$y)
  gl <- ml$n - ml$posto
  list(SQE = SQE, gl = gl, S2 = SQE / gl, QME = SQE / gl)
}

.estimativa_lambda <- function(ml, lambda) {
  b <- .solucao(ml)
  list(
    valor = sum(lambda * b),
    var_unit = as.numeric(t(lambda) %*% ml$G %*% lambda)
  )
}

.montar_L <- function(ml, L) {
  if (is.character(L)) {
    return(do.call(rbind, lapply(L, function(e) as.numeric(contraste(ml, e)))))
  }
  rbind(L)
}

.validar_L <- function(ml, L, contexto = "L") {
  nao_estim <- which(!vapply(seq_len(nrow(L)),
                             function(i) isTRUE(eh_estimavel(L[i, ], ml)), logical(1)))
  if (length(nao_estim) > 0) {
    warning(sprintf("Linha(s) %s de %s nao sao estimaveis; o resultado nao tem interpretacao valida.",
                    paste(nao_estim, collapse = ", "), contexto), call. = FALSE)
  }
  if (posto(L) < nrow(L)) {
    stop(sprintf("As linhas de %s nao sao linearmente independentes (posto %d < %d linhas).",
                 contexto, posto(L), nrow(L)), call. = FALSE)
  }
  length(nao_estim) == 0
}

.resultado_inferencia <- function(vec, tipo, expr = NULL, estimavel = NULL) {
  attr(vec, "tipo") <- tipo
  attr(vec, "expr") <- expr
  attr(vec, "estimavel") <- estimavel
  class(vec) <- c("inferencia_llm", "numeric")
  vec
}

#' @export
print.inferencia_llm <- function(x, digits = 4, ...) {
  titulo <- switch(attr(x, "tipo"),
    ic = "Intervalo de confianca",
    t  = "Teste t  (H0: lambda'beta = valor)",
    F  = "Teste F (Wald)  (H0: L beta = m)")
  cat(titulo, "\n")
  if (!is.null(attr(x, "expr"))) cat("Funcao:", paste(attr(x, "expr"), collapse = "; "), "\n")
  if (!is.null(attr(x, "estimavel"))) cat("Estimavel:", isTRUE(attr(x, "estimavel")), "\n")
  vals <- unclass(x)
  attributes(vals) <- list(names = names(vals))
  print(round(vals, digits))
  invisible(x)
}

#' Intervalo de confianca para uma funcao estimavel
#'
#' Para `lambda'beta` estimavel, devolve `lambda'b +/- t_{alpha/2; n-r} *
#' sqrt(s^2 * lambda'(X'X)^- lambda)`. Avisa se `lambda` nao for estimavel.
#'
#' @param lambda Vetor de coeficientes (comprimento = colunas de X).
#' @param objeto Objeto `llm`.
#' @param nivel Nivel de confianca (padrao 0.95).
#' @return Vetor nomeado com `estimativa`, `ep`, `gl`, `inferior`, `superior`.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' m <- llm(y ~ trat, dados)
#' intervalo_confianca(c(0, 1, 0, -1), m)
intervalo_confianca <- function(lambda, objeto, nivel = 0.95) {
  ml <- .como_llm(objeto)
  expr <- if (is.character(lambda)) lambda else NULL
  lambda <- .como_lambda(lambda, ml)
  estimavel <- isTRUE(eh_estimavel(lambda, ml))
  if (!estimavel) {
    warning("lambda nao e estimavel; o IC nao tem interpretacao valida.",
            call. = FALSE)
  }
  est <- .estimativa_lambda(ml, lambda)
  v <- qme(ml)
  ep <- sqrt(v$S2 * est$var_unit)
  tcrit <- stats::qt(1 - (1 - nivel) / 2, v$gl)
  .resultado_inferencia(
    c(estimativa = est$valor, ep = ep, gl = v$gl,
      inferior = est$valor - tcrit * ep,
      superior = est$valor + tcrit * ep),
    tipo = "ic", expr = expr, estimavel = estimavel)
}

#' Teste t para H0: lambda'beta = valor
#'
#' @param lambda Vetor de coeficientes (comprimento = colunas de X).
#' @param objeto Objeto `llm`.
#' @param valor Valor sob a hipotese nula (padrao 0).
#' @return Vetor nomeado com `estimativa`, `ep`, `t`, `gl`, `pvalor`.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' teste_t("trat1 - trat3", llm(y ~ trat, dados))
teste_t <- function(lambda, objeto, valor = 0) {
  ml <- .como_llm(objeto)
  expr <- if (is.character(lambda)) lambda else NULL
  lambda <- .como_lambda(lambda, ml)
  estimavel <- isTRUE(eh_estimavel(lambda, ml))
  if (!estimavel) {
    warning("lambda nao e estimavel; o teste nao tem interpretacao valida.",
            call. = FALSE)
  }
  est <- .estimativa_lambda(ml, lambda)
  v <- qme(ml)
  ep <- sqrt(v$S2 * est$var_unit)
  t <- (est$valor - valor) / ep
  .resultado_inferencia(
    c(estimativa = est$valor, ep = ep, t = t, gl = v$gl,
      pvalor = 2 * stats::pt(-abs(t), v$gl)),
    tipo = "t", expr = expr, estimavel = estimavel)
}

#' Teste F (Wald) para H0: L beta = m
#'
#' Testa simultaneamente `s` funcoes estimaveis (linhas de `L`). A estatistica e
#' `F = (Lb - m)' [L (X'X)^- L']^{-1} (Lb - m) / (s * s^2) ~ F(s, n-r)`.
#'
#' @param L Matriz `s x p` (cada linha uma funcao estimavel), ou um vetor de
#'   expressoes de contraste (p.ex. `c("trat1 - trat2", "trat2 - trat3")`).
#' @param objeto Objeto `llm`.
#' @param m Vetor `s` sob H0 (padrao zeros).
#' @return Vetor nomeado com `F`, `gl_num`, `gl_den`, `pvalor`.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' m <- llm(y ~ trat, dados)
#' teste_F(c("trat1 - trat2", "trat2 - trat3"), m)   # H0: sem efeito de trat
teste_F <- function(L, objeto, m = NULL) {
  ml <- .como_llm(objeto)
  expr <- if (is.character(L)) L else NULL
  L <- .montar_L(ml, L)
  if (is.null(m)) m <- numeric(nrow(L))
  estimavel <- .validar_L(ml, L, "L")
  b <- .solucao(ml)
  Lb <- as.numeric(L %*% b)
  W <- L %*% ml$G %*% t(L)
  d <- Lb - m
  Q <- as.numeric(t(d) %*% solve(W) %*% d)
  v <- qme(ml)
  s <- nrow(L)
  Fobs <- Q / (s * v$S2)
  .resultado_inferencia(
    c(F = Fobs, gl_num = s, gl_den = v$gl,
      pvalor = stats::pf(Fobs, s, v$gl, lower.tail = FALSE)),
    tipo = "F", expr = expr, estimavel = estimavel)
}

#' Regiao de confianca para um conjunto de funcoes estimaveis
#'
#' Descreve o elipsoide `(Lb - Lbeta)' [L (X'X)^- L']^{-1} (Lb - Lbeta) <=
#' s * s^2 * F_{alpha; s, n-r}`.
#'
#' @param L Matriz `s x p` (cada linha uma funcao estimavel), ou um vetor de
#'   expressoes de contraste.
#' @param objeto Objeto `llm`.
#' @param nivel Nivel de confianca (padrao 0.95).
#' @return Lista com `centro` (Lb), `W` = `L (X'X)^- L'`, `constante` (lado
#'   direito da desigualdade) e `dentro()`, funcao que testa se um ponto pertence
#'   a regiao.
#' @export
regiao_confianca <- function(L, objeto, nivel = 0.95) {
  ml <- .como_llm(objeto)
  L <- .montar_L(ml, L)
  .validar_L(ml, L, "L")
  b <- .solucao(ml)
  centro <- as.numeric(L %*% b)
  W <- L %*% ml$G %*% t(L)
  v <- qme(ml)
  s <- nrow(L)
  Fcrit <- stats::qf(nivel, s, v$gl)
  constante <- s * v$S2 * Fcrit
  Winv <- solve(W)
  reg <- list(
    centro = centro,
    W = W,
    constante = constante,
    nivel = nivel,
    gl = c(s, v$gl),
    dentro = function(ponto) {
      d <- as.numeric(ponto) - centro
      as.numeric(t(d) %*% Winv %*% d) <= constante
    }
  )
  if (s == 2) reg$pontos <- .pontos_elipse(centro, Winv, constante)
  class(reg) <- "regiao_llm"
  reg
}

#' Forma quadratica y'Ay
#'
#' @param y Vetor.
#' @param A Matriz simetrica.
#' @return Escalar `y'Ay`.
#' @export
#' @examples
#' forma_quadratica(c(1, 2, 3), diag(3))   # = 14
forma_quadratica <- function(y, A) {
  y <- as.numeric(y)
  as.numeric(t(y) %*% A %*% y)
}

#' Esperanca de uma forma quadratica
#'
#' `E[y'Ay] = tr(A Sigma) + mu' A mu`.
#'
#' @param A Matriz simetrica.
#' @param mu Vetor de medias.
#' @param Sigma Matriz de covariancias (padrao identidade).
#' @return Escalar `E[y'Ay]`.
#' @export
#' @examples
#' esperanca_forma_quadratica(diag(3), mu = c(1, 2, 3))   # tr(A) + mu'mu
esperanca_forma_quadratica <- function(A, mu, Sigma = diag(length(mu))) {
  mu <- as.numeric(mu)
  sum(diag(A %*% Sigma)) + as.numeric(t(mu) %*% A %*% mu)
}

#' Condicao de distribuicao qui-quadrado de uma forma quadratica
#'
#' `y'Ay / sigma^2` segue qui-quadrado se e somente se `A Sigma` for idempotente.
#' Os graus de liberdade sao `posto(A)`.
#'
#' @details Pressupoe `y ~ N(mu, sigma^2 Sigma)`. A funcao checa apenas a condicao
#'   distribucional (`A Sigma` idempotente); a qui-quadrado e nao-central salvo se
#'   `A mu = 0`.
#'
#' @param A Matriz simetrica.
#' @param Sigma Matriz de covariancias (padrao identidade).
#' @param tol Tolerancia numerica.
#' @return Lista com `qui_quadrado` (logico) e `gl` = `posto(A)`.
#' @export
#' @examples
#' P <- projetor(llm(y ~ trat, data.frame(
#'   trat = factor(rep(c("A", "B"), each = 3)), y = c(1, 2, 3, 4, 5, 6))))
#' forma_quadratica_qui2(P)   # projetor: idempotente => qui-quadrado
forma_quadratica_qui2 <- function(A, Sigma = diag(nrow(A)), tol = 1e-8) {
  AS <- A %*% Sigma
  idem <- max(abs(AS %*% AS - AS)) < tol
  list(qui_quadrado = idem, gl = posto(A))
}
