.pontos_elipse <- function(centro, Winv, k, n = 200) {
  e <- eigen(Winv, symmetric = TRUE)
  V <- e$vectors
  d <- e$values
  t <- seq(0, 2 * pi, length.out = n)
  w <- rbind(sqrt(k) * cos(t), sqrt(k) * sin(t))
  pts <- V %*% diag(1 / sqrt(d)) %*% w
  data.frame(x = centro[1] + pts[1, ], y = centro[2] + pts[2, ])
}

#' @export
print.regiao_llm <- function(x, ...) {
  cat(sprintf("Regiao de confianca (%.0f%%) para %d funcoes estimaveis\n",
              100 * x$nivel, x$gl[1]))
  cat("Centro (Lb):", round(x$centro, 4), "\n")
  cat(sprintf("Constante (s * s^2 * F_{%g; %d, %d}): %.4f\n",
              x$nivel, x$gl[1], x$gl[2], x$constante))
  if (!is.null(x$pontos)) {
    cat("Use plot() para ver a elipse, ou $pontos para os limites.\n")
  }
  invisible(x)
}

#' Grafico da regiao de confianca (2 funcoes)
#'
#' Desenha a elipse de confianca quando ha exatamente duas funcoes estimaveis.
#'
#' @param x Objeto `regiao_llm` (de [regiao_confianca()]).
#' @param rotulos Rotulos dos eixos (padrao `c("psi1", "psi2")`).
#' @param ... Passado a `plot`.
#' @return Invisivelmente, o objeto.
#' @export
plot.regiao_llm <- function(x, rotulos = c("psi1", "psi2"), ...) {
  if (is.null(x$pontos)) {
    stop("O grafico so esta disponivel para regioes de 2 funcoes.", call. = FALSE)
  }
  graphics::plot(x$pontos$x, x$pontos$y, type = "l",
                 xlab = rotulos[1], ylab = rotulos[2],
                 main = sprintf("Regiao de confianca %.0f%%", 100 * x$nivel), ...)
  graphics::points(x$centro[1], x$centro[2], pch = 19)
  graphics::abline(h = x$centro[2], v = x$centro[1], lty = 3, col = "grey60")
  invisible(x)
}
