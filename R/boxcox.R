.redondos_boxcox <- c(
  identidade = 1, log = 0, raiz = 0.5,
  inversa_raiz = -0.5, inversa = -1, quadrado = 2
)

.diagnostico_boxcox <- function(formula, dados, nivel = 0.95) {
  y <- .resposta(formula, dados)
  if (is.null(y)) {
    return(list(ok = FALSE, motivo = "sem resposta y"))
  }
  if (any(y <= 0)) {
    return(list(ok = FALSE, motivo = "resposta nao-positiva (Box-Cox exige y > 0)"))
  }

  bc <- tryCatch(
    MASS::boxcox(formula, data = dados, lambda = seq(-2, 2, 0.01), plotit = FALSE),
    error = function(e) NULL
  )
  if (is.null(bc)) {
    return(list(ok = FALSE, motivo = "MASS::boxcox falhou para esta formula"))
  }

  lambda_hat <- bc$x[which.max(bc$y)]
  limite <- max(bc$y) - stats::qchisq(nivel, 1) / 2
  dentro <- bc$x[bc$y >= limite]
  ic <- range(dentro)
  contem_1 <- ic[1] <= 1 && 1 <= ic[2]

  no_ic <- .redondos_boxcox[.redondos_boxcox >= ic[1] & .redondos_boxcox <= ic[2]]
  candidatos <- if (length(no_ic) > 0) no_ic else .redondos_boxcox
  i_sug <- which.min(abs(candidatos - lambda_hat))
  tipo_sug <- names(candidatos)[i_sug]
  lambda_sug <- unname(candidatos[i_sug])

  list(
    ok = TRUE, motivo = NULL, nivel = nivel,
    lambda_hat = lambda_hat, ic = ic, contem_1 = contem_1,
    tipo_sug = tipo_sug, lambda_sug = lambda_sug
  )
}

.aviso_boxcox <- function(bc) {
  sprintf(
    paste0(
      "Pressuposto Box-Cox: 1 nao esta no IC%.0f%% (lambda_hat = %.2f, ",
      "IC = [%.2f, %.2f]). Sugestao: '%s' (lambda = %g). ",
      "Use transformar_resposta(m) para aplicar."
    ),
    100 * bc$nivel, bc$lambda_hat, bc$ic[1], bc$ic[2], bc$tipo_sug, bc$lambda_sug
  )
}

.transforma_y <- function(y, tipo = NULL, lambda = NULL) {
  exige_positivo <- function() {
    if (any(y <= 0)) {
      stop("Esta transformacao exige y > 0.", call. = FALSE)
    }
  }
  if (!is.null(tipo)) {
    if (!tipo %in% names(.redondos_boxcox)) {
      stop(sprintf("tipo '%s' desconhecido. Use um de: %s.", tipo,
                   paste(names(.redondos_boxcox), collapse = ", ")),
           call. = FALSE)
    }
    switch(tipo,
      identidade   = list(valor = y,           rotulo = "y"),
      log          = { exige_positivo(); list(valor = log(y),       rotulo = "log(y)") },
      raiz         = { exige_positivo(); list(valor = sqrt(y),      rotulo = "sqrt(y)") },
      inversa_raiz = { exige_positivo(); list(valor = 1 / sqrt(y),  rotulo = "1/sqrt(y)") },
      inversa      = { exige_positivo(); list(valor = 1 / y,        rotulo = "1/y") },
      quadrado     = list(valor = y^2,         rotulo = "y^2")
    )
  } else {
    exige_positivo()
    if (lambda == 0) {
      list(valor = log(y), rotulo = "log(y)")
    } else {
      list(valor = (y^lambda - 1) / lambda,
           rotulo = sprintf("(y^%g - 1)/%g", lambda, lambda))
    }
  }
}

#' Transforma a resposta e reajusta o modelo
#'
#' Aplica uma transformacao normalizadora a `y` e reajusta o `llm`. Sem `lambda`
#' nem `tipo`, usa a transformacao padrao sugerida pelo diagnostico de Box-Cox
#' guardado em `m$boxcox`.
#'
#' @param m Objeto `llm`.
#' @param lambda Expoente de Box-Cox; aplica `(y^lambda - 1)/lambda` (ou `log(y)`
#'   se `lambda = 0`).
#' @param tipo Transformacao literal: `"identidade"`, `"log"`, `"raiz"`,
#'   `"inversa_raiz"`, `"inversa"` ou `"quadrado"`.
#' @return Novo objeto `llm`, com `m$transformacao` registrando o que foi aplicado.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 4)),
#'   y = c(1, 2, 1, 3, 12, 18, 15, 20, 80, 95, 70, 110)
#' )
#' m <- llm(y ~ trat, dados)
#' m2 <- transformar_resposta(m, tipo = "log")
#' m2$transformacao
transformar_resposta <- function(m, lambda = NULL, tipo = NULL) {
  ml <- .como_llm(m)
  if (is.null(ml$y)) stop("O modelo nao tem resposta y.", call. = FALSE)

  if (is.null(tipo) && is.null(lambda)) {
    if (is.null(ml$boxcox) || !isTRUE(ml$boxcox$ok)) {
      stop(paste("Sem diagnostico de Box-Cox disponivel para sugerir a",
                 "transformacao. Informe lambda = ou tipo =."), call. = FALSE)
    }
    tipo <- ml$boxcox$tipo_sug
  }

  tr <- .transforma_y(ml$y, tipo = tipo, lambda = lambda)

  yname <- all.vars(ml$formula)[1]
  dados_novo <- ml$dados
  dados_novo[[yname]] <- tr$valor

  m_novo <- llm(ml$formula, dados_novo, intercepto = ml$intercepto,
                checar_boxcox = !is.null(ml$boxcox))
  m_novo$transformacao <- list(tipo = tipo, lambda = lambda, rotulo = tr$rotulo)
  m_novo
}
