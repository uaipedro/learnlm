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
#' @param checar_boxcox Se `TRUE` (padrao), roda o diagnostico de Box-Cox sobre a
#'   resposta e avisa quando `lambda = 1` fica fora do IC95% (veja `m$boxcox` e
#'   [transformar_resposta()]).
#' @param fix_boxcox Se `NULL` (padrao), apenas avisa. Se `"tipo"` ou `"lambda"`
#'   e o pressuposto for violado, ja devolve o modelo corrigido — `"tipo"` aplica
#'   a transformacao padrao sugerida e `"lambda"` aplica a forma de Box-Cox com o
#'   `lambda_hat` estimado (via [transformar_resposta()]).
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
llm <- function(formula, dados, intercepto = TRUE, checar_boxcox = TRUE,
                fix_boxcox = NULL) {
  if (!is.null(fix_boxcox) && !fix_boxcox %in% c("tipo", "lambda")) {
    stop("fix_boxcox deve ser NULL, 'tipo' ou 'lambda'.", call. = FALSE)
  }
  X <- matriz_modelo(formula, dados, intercepto)
  y <- .resposta(formula, dados)
  XtX <- t(X) %*% X
  Xty <- if (!is.null(y)) as.numeric(t(X) %*% y) else NULL
  r <- posto(X)
  p <- ncol(X)
  G <- inversa_condicional(XtX)

  obj <- list(
    formula = formula,
    dados = dados,
    intercepto = intercepto,
    X = X,
    y = y,
    XtX = XtX,
    Xty = Xty,
    G = G,
    P = X %*% G %*% t(X),
    posto = r,
    p = p,
    n = nrow(X),
    deficiencia = p - r,
    descr = attr(X, "descr"),
    niveis_ref = .niveis_referencia(formula, dados)
  )
  class(obj) <- "llm"
  obj$passos <- .construir_passos(obj)

  if (checar_boxcox || !is.null(fix_boxcox)) {
    obj$boxcox <- .diagnostico_boxcox(formula, dados)
    if (isTRUE(obj$boxcox$ok) && !obj$boxcox$contem_1) {
      if (is.null(fix_boxcox)) {
        warning(.aviso_boxcox(obj$boxcox), call. = FALSE)
      } else {
        corrigido <- suppressWarnings(
          if (fix_boxcox == "lambda") {
            transformar_resposta(obj, lambda = obj$boxcox$lambda_hat)
          } else {
            transformar_resposta(obj)
          }
        )
        status <- if (isTRUE(corrigido$boxcox$contem_1)) {
          "1 agora esta no IC95%"
        } else {
          sprintf("1 ainda fora do IC95%% (lambda_hat = %.2f)",
                  corrigido$boxcox$lambda_hat)
        }
        message(sprintf("Box-Cox: resposta corrigida (%s); %s.",
                        corrigido$transformacao$rotulo, status))
        return(corrigido)
      }
    }
  }
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
  if (!is.null(x$transformacao)) {
    cat(sprintf("Resposta transformada: %s\n", x$transformacao$rotulo))
  }
  bc <- x$boxcox
  if (!is.null(bc)) {
    if (isTRUE(bc$ok)) {
      cat(sprintf("Box-Cox: lambda_hat = %.2f, IC%.0f%% = [%.2f, %.2f]; 1 %s no IC",
                  bc$lambda_hat, 100 * bc$nivel, bc$ic[1], bc$ic[2],
                  if (bc$contem_1) "ESTA" else "NAO esta"))
      if (!bc$contem_1) {
        cat(sprintf(" -> sugestao: %s (lambda = %g)", bc$tipo_sug, bc$lambda_sug))
      }
      cat("\n")
    } else {
      cat(sprintf("Box-Cox: nao avaliado (%s)\n", bc$motivo))
    }
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

.g_inversa <- function(objeto) {
  if (inherits(objeto, "llm")) return(objeto$G)
  X <- as.matrix(objeto)
  inversa_condicional(t(X) %*% X)
}

.solucao <- function(ml) as.numeric(ml$G %*% ml$Xty)
