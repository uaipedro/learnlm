.parse_combinacao <- function(expr) {
  expr <- gsub("\\s+", "", expr)
  # limitacao conhecida: coeficientes em notacao cientifica (2e-3) quebram aqui,
  # pois o '-' do expoente vira separador de termo. Use 0.002 no lugar.
  expr <- gsub("-", "+-", expr)
  partes <- strsplit(expr, "+", fixed = TRUE)[[1]]
  partes[partes != ""]
}

.coef_token <- function(parte) {
  mm <- regmatches(parte, regexec("^(-?)([0-9]*\\.?[0-9]*)\\*?(.+)$", parte))[[1]]
  sinal <- mm[2]
  num <- mm[3]
  token <- mm[4]
  if (num == "") num <- "1"
  list(coef = as.numeric(paste0(sinal, num)), token = token)
}

.split_indices <- function(resto, k) {
  if (grepl("[^A-Za-z0-9]", resto)) {
    return(strsplit(resto, "[^A-Za-z0-9]+")[[1]])
  }
  if (k == 1) return(resto)
  strsplit(resto, "")[[1]]
}

.nivel_por_indice <- function(ml, f, idx) {
  niveis <- levels(as.factor(ml$dados[[f]]))
  if (idx %in% niveis) return(idx)
  pos <- suppressWarnings(as.integer(idx))
  if (is.na(pos) || pos < 1 || pos > length(niveis)) {
    stop(sprintf("Indice '%s' invalido para o fator '%s'.", idx, f), call. = FALSE)
  }
  niveis[pos]
}

.resolver_token <- function(token, ml, simbolos) {
  cols <- colnames(ml$X)
  if (token %in% cols) return(match(token, cols))

  sym <- sub("^([A-Za-z]+).*$", "\\1", token)
  resto <- sub("^[A-Za-z]+", "", token)
  termo <- if (!is.null(simbolos) && sym %in% names(simbolos)) simbolos[[sym]] else sym
  fatores <- strsplit(termo, ":", fixed = TRUE)[[1]]
  idxs <- .split_indices(resto, length(fatores))
  if (length(idxs) != length(fatores)) {
    stop(sprintf("Nao consegui mapear '%s' para o termo '%s' (use simbolos= ou separadores).",
                 token, termo), call. = FALSE)
  }
  alvo <- mapply(function(f, idx) .nivel_por_indice(ml, f, idx), fatores, idxs)

  for (j in seq_len(ml$p)) {
    d <- ml$descr[[j]]
    if (!identical(sort(names(d$niveis)), sort(fatores))) next
    casa <- all(vapply(fatores, function(f) {
      identical(unname(d$niveis[f]), unname(alvo[f]))
    }, logical(1)))
    if (casa) return(j)
  }
  stop(sprintf("Coluna para '%s' nao encontrada (celula inexistente?).", token),
       call. = FALSE)
}

#' Constroi um contraste/funcao linear a partir de uma string
#'
#' Traduz uma expressao legivel como `"trat1 - trat3"` ou `"g12 - g21"` no vetor
#' `lambda` correspondente as colunas de X. Tokens podem ser:
#' nomes exatos de coluna (`"tratA"`, `"tipo1:metodo2"`); ou um simbolo seguido
#' de indices de nivel (`"trat1"` = primeiro nivel de `trat`). Use `simbolos`
#' para apelidar termos (p.ex. `c(t = "trat", g = "tipo:metodo")`), de modo que
#' `"t1"` e `"g12"` funcionem. Coeficientes sao aceitos: `"2*t1 - 0.5*t3"`.
#'
#' @param objeto Objeto `llm`.
#' @param expr String com a combinacao linear.
#' @param simbolos Vetor nomeado opcional, mapeando simbolo -> termo.
#' @return Vetor `lambda` nomeado, de classe `contraste`, com atributos `expr` e
#'   `estimavel`. Pode ser passado a [eh_estimavel()], [intervalo_confianca()],
#'   [teste_t()] e [teste_F()].
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
#' )
#' m <- llm(y ~ trat, dados)
#' contraste(m, "trat1 - trat3")
#' contraste(m, "t1 - t3", simbolos = c(t = "trat"))
contraste <- function(objeto, expr, simbolos = NULL) {
  ml <- .como_llm(objeto)
  partes <- .parse_combinacao(expr)
  lambda <- numeric(ml$p)
  names(lambda) <- colnames(ml$X)
  for (parte in partes) {
    ct <- .coef_token(parte)
    j <- .resolver_token(ct$token, ml, simbolos)
    lambda[j] <- lambda[j] + ct$coef
  }
  attr(lambda, "expr") <- expr
  attr(lambda, "estimavel") <- isTRUE(eh_estimavel(lambda, ml))
  class(lambda) <- c("contraste", "numeric")
  lambda
}

#' @export
print.contraste <- function(x, ...) {
  cat("Contraste:", attr(x, "expr"), "\n")
  v <- unclass(x)
  nz <- v[abs(v) > 1e-12]
  cat("lambda (coeficientes nao-nulos):\n")
  print(nz)
  cat("Estimavel:", isTRUE(attr(x, "estimavel")), "\n")
  invisible(x)
}

.como_lambda <- function(x, objeto) {
  if (is.character(x)) {
    return(as.numeric(contraste(objeto, x)))
  }
  as.numeric(x)
}
