#' Solucoes do sistema de equacoes normais
#'
#' Resolve `X'X b = X'y` de tres modos didaticos: sem restricao (usando uma
#' g-inversa de `X'X`), ou impondo uma restricao (soma zero ou casela de
#' referencia). Em todos os casos os valores ajustados `Xb` sao identicos — o que
#' muda e a solucao `b`. Essa invariancia e a licao central de modelos de posto
#' incompleto.
#'
#' @param objeto Objeto `llm`.
#' @param restricao Um de `"nenhuma"`, `"soma_zero"`, `"casela_referencia"`.
#' @param inversa Quando `restricao = "nenhuma"`, qual g-inversa de `X'X` usar:
#'   `"condicional"`, `"moore_penrose"` ou `"minimos_quadrados"`.
#' @return Vetor `b` nomeado pelas colunas de X. O atributo `ajustados` traz `Xb`.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 2)),
#'   y = c(10, 12, 20, 22, 30, 28)
#' )
#' m <- llm(y ~ trat, dados)
#' betas(m, restricao = "nenhuma")
#' betas(m, restricao = "soma_zero")
#' betas(m, restricao = "casela_referencia")
betas <- function(objeto,
                  restricao = c("nenhuma", "soma_zero", "casela_referencia"),
                  inversa = c("condicional", "moore_penrose", "minimos_quadrados")) {
  ml <- .como_llm(objeto)
  if (is.null(ml$y)) stop("O modelo nao tem resposta y.", call. = FALSE)
  restricao <- match.arg(restricao)

  if (restricao == "nenhuma") {
    inversa <- match.arg(inversa)
    G <- switch(inversa,
      condicional = ml$G,
      moore_penrose = inversa_moore_penrose(ml$XtX),
      minimos_quadrados = inversa_minimos_quadrados(ml$XtX)
    )
    beta <- as.numeric(G %*% ml$Xty)
  } else {
    H <- .matriz_restricao(ml, restricao)
    beta <- .betas_restrito(ml$X, ml$y, H)
  }

  names(beta) <- colnames(ml$X)
  attr(beta, "ajustados") <- as.numeric(ml$X %*% beta)
  attr(beta, "restricao") <- restricao
  beta
}

.betas_restrito <- function(X, y, H) {
  XtX <- t(X) %*% X
  Xty <- as.numeric(t(X) %*% y)
  p <- ncol(X)
  q <- nrow(H)
  if (posto(rbind(X, H)) < p) {
    stop(paste("As restricoes nao identificam a solucao: posto([X; H]) < numero de",
               "parametros. Faltam restricoes (ou elas sao dependentes das colunas de X)."),
         call. = FALSE)
  }
  K <- rbind(cbind(XtX, t(H)), cbind(H, matrix(0, q, q)))
  rhs <- c(Xty, numeric(q))
  sol <- solve(K, rhs)
  sol[seq_len(p)]
}

.matriz_restricao <- function(ml, restricao) {
  H <- switch(restricao,
    soma_zero = .restricao_soma_zero(ml),
    casela_referencia = .restricao_casela(ml)
  )
  .linhas_independentes(H)
}

.linhas_independentes <- function(H) {
  if (is.null(H) || nrow(H) == 0) return(H)
  q <- qr(t(H))
  H[q$pivot[seq_len(q$rank)], , drop = FALSE]
}

.restricao_soma_zero <- function(ml) {
  descr <- ml$descr
  p <- ml$p
  linhas <- list()
  termos <- unique(vapply(descr, `[[`, "", "termo"))
  for (tm in termos) {
    idx <- which(vapply(descr, function(d) d$termo == tm, logical(1)))
    d1 <- descr[[idx[1]]]
    if (d1$tipo %in% c("intercepto", "numerico")) next
    fatores <- names(d1$niveis)
    for (f in fatores) {
      outros <- setdiff(fatores, f)
      chaves <- vapply(idx, function(i) {
        paste(descr[[i]]$niveis[outros], collapse = "|")
      }, "")
      for (ch in unique(chaves)) {
        linha <- numeric(p)
        linha[idx[chaves == ch]] <- 1
        linhas[[length(linhas) + 1L]] <- linha
      }
    }
  }
  do.call(rbind, linhas)
}

.restricao_casela <- function(ml) {
  descr <- ml$descr
  p <- ml$p
  ref <- ml$niveis_ref
  linhas <- list()
  for (i in seq_len(p)) {
    d <- descr[[i]]
    if (d$tipo %in% c("intercepto", "numerico")) next
    fatores <- names(d$niveis)
    eh_ref <- any(vapply(fatores, function(f) {
      identical(unname(d$niveis[f]), ref[[f]])
    }, logical(1)))
    if (eh_ref) {
      linha <- numeric(p)
      linha[i] <- 1
      linhas[[length(linhas) + 1L]] <- linha
    }
  }
  do.call(rbind, linhas)
}

#' Projetor ortogonal no espaco coluna de X
#'
#' Calcula `P = X (X'X)^- X'`, a matriz que projeta `y` no espaco gerado pelas
#' colunas de X. `P` e simetrica, idempotente e invariante a escolha da
#' g-inversa. Os valores ajustados sao `Py` e os residuos `(I - P) y`.
#'
#' @param objeto Objeto `llm` ou matriz X.
#' @param residuos Se `TRUE`, devolve `M = I - P` (projetor dos residuos).
#' @param explicar Se `TRUE`, imprime as verificacoes e devolve a matriz
#'   invisivelmente.
#' @return Matriz `P` (ou `M`).
#' @export
projetor <- function(objeto, residuos = FALSE, explicar = FALSE) {
  if (inherits(objeto, "llm")) {
    P <- objeto$P
  } else {
    X <- .pegar_X(objeto)
    P <- X %*% inversa_condicional(t(X) %*% X) %*% t(X)
  }
  saida <- if (residuos) diag(nrow(P)) - P else P

  if (explicar) {
    cat(sprintf("== Projetor %s ==\n", if (residuos) "dos residuos M = I - P" else "P = X(X'X)^- X'"))
    cat("Simetrica:", isTRUE(all.equal(P, t(P), check.attributes = FALSE)), "\n")
    cat("Idempotente (P^2 = P):",
        isTRUE(all.equal(P %*% P, P, check.attributes = FALSE)), "\n")
    cat("posto(P) =", posto(P), "= posto(X)\n")
    return(invisible(saida))
  }
  saida
}

#' Testa se uma funcao linear dos parametros e estimavel
#'
#' `lambda' b` e estimavel se e somente se `lambda` pertence ao espaco linha de X,
#' isto e, `H lambda = lambda` com `H = (X'X)^- (X'X)`. Em modelos de posto
#' incompleto, parametros individuais nao sao estimaveis, mas contrastes sim.
#'
#' @param lambda Vetor de coeficientes (comprimento = numero de colunas de X).
#' @param objeto Objeto `llm` ou matriz X.
#' @param tol Tolerancia numerica.
#' @return `TRUE`/`FALSE`. O atributo `residuo` traz `max|H lambda - lambda|`.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 2)),
#'   y = c(10, 12, 20, 22, 30, 28)
#' )
#' m <- llm(y ~ trat, dados)
#' eh_estimavel(c(0, 1, 0, 0), m)   # tau_A sozinho: NAO estimavel
#' eh_estimavel(c(0, 1, -1, 0), m)  # contraste tau_A - tau_B: estimavel
eh_estimavel <- function(lambda, objeto, tol = 1e-8) {
  X <- .pegar_X(objeto)
  XtX <- if (inherits(objeto, "llm")) objeto$XtX else t(X) %*% X
  H <- XtX %*% .g_inversa(objeto)
  lambda <- if (is.character(lambda)) .como_lambda(lambda, objeto) else as.numeric(lambda)
  if (length(lambda) != ncol(X)) {
    stop(sprintf("lambda deve ter comprimento %d (colunas de X).", ncol(X)),
         call. = FALSE)
  }
  residuo <- max(abs(as.numeric(H %*% lambda) - lambda))
  resultado <- residuo < tol
  attr(resultado, "residuo") <- residuo
  resultado
}
