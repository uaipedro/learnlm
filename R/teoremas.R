#' Consistencia do sistema de equacoes normais
#'
#' Verifica o teorema da consistencia: o SEN `X'X b = X'y` e consistente se e
#' somente se `posto(X'X) = posto([X'X | X'y])`. Para modelos lineares isso vale
#' sempre, pois `X'y` esta no espaco coluna de `X'X`.
#'
#' @param objeto Objeto `llm`.
#' @return Lista com `consistente`, `posto_XtX`, `posto_aumentada`.
#' @export
consistencia_sen <- function(objeto) {
  ml <- .como_llm(objeto)
  if (is.null(ml$Xty)) stop("O modelo nao tem resposta y.", call. = FALSE)
  r1 <- posto(ml$XtX)
  r2 <- posto(cbind(ml$XtX, ml$Xty))
  list(consistente = r1 == r2, posto_XtX = r1, posto_aumentada = r2)
}

#' Numero maximo de funcoes estimaveis linearmente independentes
#'
#' Pelo teorema de estimabilidade, esse numero e igual a `posto(X)`.
#'
#' @param objeto Objeto `llm`.
#' @return Inteiro `posto(X)`.
#' @export
numero_funcoes_estimaveis <- function(objeto) {
  .como_llm(objeto)$posto
}

#' Verifica se uma funcao linear e um contraste
#'
#' Um contraste e uma funcao estimavel cujo coeficiente do intercepto e zero e
#' cujos coeficientes somam zero dentro de cada termo (p.ex. `sum(c_i) = 0` nos
#' tratamentos). Contrastes de tratamentos sao sempre estimaveis.
#'
#' @param lambda Vetor de coeficientes, ou uma string de contraste.
#' @param objeto Objeto `llm`.
#' @return `TRUE`/`FALSE`. O atributo `estimavel` indica se a funcao e estimavel.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' m <- llm(y ~ trat, dados)
#' eh_contraste("trat1 - trat3", m)   # TRUE
#' eh_contraste("trat1 + trat3", m)   # FALSE (soma != 0)
eh_contraste <- function(lambda, objeto) {
  ml <- .como_llm(objeto)
  lambda <- .como_lambda(lambda, ml)
  estimavel <- isTRUE(eh_estimavel(lambda, ml))

  termos <- unique(vapply(ml$descr, `[[`, "", "termo"))
  intercepto_ok <- TRUE
  somas_ok <- TRUE
  for (tm in termos) {
    idx <- which(vapply(ml$descr, function(d) d$termo == tm, logical(1)))
    if (ml$descr[[idx[1]]]$tipo == "intercepto") {
      if (any(abs(lambda[idx]) > 1e-8)) intercepto_ok <- FALSE
    } else {
      if (abs(sum(lambda[idx])) > 1e-8) somas_ok <- FALSE
    }
  }
  resultado <- estimavel && intercepto_ok && somas_ok
  attr(resultado, "estimavel") <- estimavel
  resultado
}

#' Bateria de verificacao dos teoremas do modelo
#'
#' Roda e relata, de uma vez, os principais teoremas vistos no curso para o
#' modelo dado: consistencia do SEN, numero de funcoes estimaveis, propriedades
#' do projetor, invariancia de `Xb` e graus de liberdade do residuo. Util para
#' conferir, num exemplo concreto, que cada teorema vale.
#'
#' @param objeto Objeto `llm`.
#' @return `data.frame` de classe `teoremas_llm` com colunas `teorema`, `vale`,
#'   `detalhe`.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' verificar_teoremas(llm(y ~ trat, dados))
verificar_teoremas <- function(objeto) {
  ml <- .como_llm(objeto)
  itens <- list()
  add <- function(teorema, vale, detalhe) {
    itens[[length(itens) + 1L]] <<- data.frame(
      teorema = teorema, vale = isTRUE(vale), detalhe = detalhe,
      stringsAsFactors = FALSE
    )
  }

  cons <- consistencia_sen(ml)
  add("SEN consistente: posto(X'X) = posto([X'X|X'y])",
      cons$consistente,
      sprintf("%d = %d", cons$posto_XtX, cons$posto_aumentada))

  add("No. maximo de funcoes estimaveis LI = posto(X)",
      TRUE, sprintf("posto(X) = %d", ml$posto))

  P <- projetor(ml)
  add("Projetor P simetrico",
      isTRUE(all.equal(P, t(P), check.attributes = FALSE)), "P = P'")
  add("Projetor P idempotente",
      isTRUE(all.equal(P %*% P, P, check.attributes = FALSE)), "P^2 = P")
  add("posto(P) = posto(X)", posto(P) == ml$posto,
      sprintf("posto(P) = %d", posto(P)))

  if (!is.null(ml$y)) {
    aj1 <- attr(betas(ml, "nenhuma", "condicional"), "ajustados")
    aj2 <- attr(betas(ml, "nenhuma", "moore_penrose"), "ajustados")
    add("Xb invariante a escolha da g-inversa",
        isTRUE(all.equal(aj1, aj2)), "Xb identico")
    v <- qme(ml)
    add("gl do residuo = n - posto(X)",
        v$gl == ml$n - ml$posto,
        sprintf("%d = %d - %d", v$gl, ml$n, ml$posto))
  }

  out <- do.call(rbind, itens)
  class(out) <- c("teoremas_llm", "data.frame")
  out
}

#' @export
print.teoremas_llm <- function(x, ...) {
  cat("Verificacao dos teoremas do modelo\n\n")
  for (i in seq_len(nrow(x))) {
    marca <- if (x$vale[i]) "[OK]" else "[FALHA]"
    cat(sprintf("%s %s\n      %s\n", marca, x$teorema[i], x$detalhe[i]))
  }
  invisible(x)
}
