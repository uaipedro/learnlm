#' Posto de uma matriz
#'
#' @param A Matriz numerica.
#' @param tol Tolerancia para o posto.
#' @return Inteiro com o posto de `A`.
#' @export
#' @examples
#' posto(matrix(c(1, 2, 2, 4), 2))
posto <- function(A, tol = 1e-7) {
  qr(as.matrix(A), tol = tol)$rank
}

#' Inversa generalizada condicional
#'
#' Constroi uma g-inversa (g1) tal que `A %*% G %*% A == A`, pelo metodo da
#' submatriz nao-singular: escolhe r (= `posto(A)`) linhas e colunas independentes,
#' inverte o bloco e posiciona o resultado. Nao e unica.
#'
#' @param A Matriz numerica.
#' @param explicar Se `TRUE`, imprime a construcao e a verificacao e devolve a
#'   g-inversa invisivelmente.
#' @param tol Tolerancia para o posto.
#' @return Matriz `G` (g-inversa condicional de `A`).
#' @export
#' @examples
#' A <- matrix(c(1, 1, 1, 0, 1, 0, 1, 0, 1), 3)
#' G <- inversa_condicional(A)
#' all.equal(A %*% G %*% A, A)
inversa_condicional <- function(A, explicar = FALSE, tol = 1e-7) {
  A <- as.matrix(A)
  r <- qr(A, tol = tol)$rank
  cols_ind <- qr(A, tol = tol)$pivot[seq_len(r)]
  linhas_ind <- qr(t(A), tol = tol)$pivot[seq_len(r)]
  M <- A[linhas_ind, cols_ind, drop = FALSE]
  G <- matrix(0, ncol(A), nrow(A))
  G[cols_ind, linhas_ind] <- solve(M)

  if (explicar) {
    cat("== Inversa generalizada condicional ==\n")
    cat(sprintf("posto(A) = %d (de %d x %d)\n", r, nrow(A), ncol(A)))
    cat("Linhas independentes escolhidas:", linhas_ind, "\n")
    cat("Colunas independentes escolhidas:", cols_ind, "\n")
    cat("Bloco nao-singular M = A[linhas, colunas] invertido e reposicionado.\n")
    cat("Verificacao A G A = A:",
        isTRUE(all.equal(A %*% G %*% A, A, check.attributes = FALSE)), "\n")
    return(invisible(G))
  }
  G
}

#' Inversa generalizada de minimos quadrados
#'
#' G-inversa `G` tal que `A %*% G %*% A == A` e `A %*% G` e simetrica. Construida
#' como `(A'A)^- A'`, usando a g-inversa condicional de `A'A`.
#'
#' @param A Matriz numerica.
#' @param explicar Se `TRUE`, imprime a verificacao e devolve `G` invisivelmente.
#' @return Matriz `G` (g-inversa de minimos quadrados de `A`).
#' @export
#' @examples
#' A <- matrix(c(1, 1, 1, 0, 1, 0, 1, 0, 1), 3)
#' inversa_minimos_quadrados(A, explicar = TRUE)
inversa_minimos_quadrados <- function(A, explicar = FALSE) {
  A <- as.matrix(A)
  AtA <- t(A) %*% A
  G <- inversa_condicional(AtA) %*% t(A)

  if (explicar) {
    AG <- A %*% G
    cat("== Inversa generalizada de minimos quadrados ==\n")
    cat("G = (A'A)^- A'\n")
    cat("Verificacao A G A = A:",
        isTRUE(all.equal(A %*% G %*% A, A, check.attributes = FALSE)), "\n")
    cat("Verificacao A G simetrica:",
        isTRUE(all.equal(AG, t(AG), check.attributes = FALSE)), "\n")
    return(invisible(G))
  }
  G
}

#' Inversa de Moore-Penrose
#'
#' G-inversa unica que satisfaz as quatro condicoes de Penrose, obtida via
#' decomposicao em valores singulares (`MASS::ginv`).
#'
#' @param A Matriz numerica.
#' @param explicar Se `TRUE`, imprime as verificacoes e devolve `G` invisivelmente.
#' @return Matriz `G` (inversa de Moore-Penrose de `A`).
#' @export
#' @examples
#' A <- matrix(c(1, 1, 1, 0, 1, 0, 1, 0, 1), 3)
#' inversa_moore_penrose(A, explicar = TRUE)
inversa_moore_penrose <- function(A, explicar = FALSE) {
  A <- as.matrix(A)
  G <- MASS::ginv(A)

  if (explicar) {
    cat("== Inversa de Moore-Penrose (via SVD) ==\n")
    cat("Condicoes de Penrose:\n")
    cat("  1. A G A = A :",
        isTRUE(all.equal(A %*% G %*% A, A, check.attributes = FALSE)), "\n")
    cat("  2. G A G = G :",
        isTRUE(all.equal(G %*% A %*% G, G, check.attributes = FALSE)), "\n")
    AG <- A %*% G
    GA <- G %*% A
    cat("  3. (A G)' = A G :",
        isTRUE(all.equal(AG, t(AG), check.attributes = FALSE)), "\n")
    cat("  4. (G A)' = G A :",
        isTRUE(all.equal(GA, t(GA), check.attributes = FALSE)), "\n")
    return(invisible(G))
  }
  G
}
