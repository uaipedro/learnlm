.passo <- function(numero, titulo, explicacao, objeto = NULL, verificacao = NULL) {
  structure(
    list(numero = numero, titulo = titulo, explicacao = explicacao,
         objeto = objeto, verificacao = verificacao),
    class = "passo_llm"
  )
}

#' @export
print.passo_llm <- function(x, ...) {
  cat(sprintf("== Passo %d: %s ==\n", x$numero, x$titulo))
  cat(strwrap(x$explicacao, width = 76), sep = "\n")
  cat("\n")
  if (!is.null(x$objeto)) {
    cat("\n")
    print(x$objeto)
  }
  if (!is.null(x$verificacao)) {
    cat("\nVerificacao:", x$verificacao, "\n")
  }
  invisible(x)
}

.construir_passos <- function(ml) {
  ps <- list()

  ps[[1]] <- .passo(
    1, "Matriz de delineamento X (completa)",
    paste(
      "X tem uma coluna por nivel de cada fator (super-parametrizada), sem casela",
      "de referencia. Cada linha e uma observacao. Em modelos de posto incompleto",
      sprintf("o posto de X (%d) e menor que o numero de parametros (%d).",
              ml$posto, ml$p)
    ),
    objeto = ml$X
  )

  ps[[2]] <- .passo(
    2, "Sistema de equacoes normais  X'X b = X'y",
    paste(
      "Minimizar ||y - Xb||^2 leva a X'X b = X'y. A matriz X'X e quadrada",
      sprintf("(%d x %d) e simetrica.", ml$p, ml$p),
      if (ml$deficiencia > 0)
        "Como X tem posto incompleto, X'X e SINGULAR: nao tem inversa usual."
      else "X'X e nao-singular."
    ),
    objeto = ml$XtX
  )

  ps[[3]] <- .passo(
    3, "Posto e deficiencia",
    paste(
      sprintf("posto(X) = %d, parametros p = %d.", ml$posto, ml$p),
      if (ml$deficiencia > 0)
        sprintf(paste("Deficiencia = %d: existem infinitas solucoes b, e sao",
                      "necessarias g-inversas ou restricoes para obter uma."),
                ml$deficiencia)
      else "Posto completo: solucao unica."
    ),
    verificacao = sprintf("posto(X) = posto(X'X) = %d", posto(ml$XtX))
  )

  if (ml$deficiencia > 0 && !is.null(ml$y)) {
    bc <- .solucao(ml)
    bmp <- as.numeric(inversa_moore_penrose(ml$XtX) %*% ml$Xty)
    aj_c <- as.numeric(ml$X %*% bc)
    aj_mp <- as.numeric(ml$X %*% bmp)
    comp <- rbind(condicional = bc, moore_penrose = bmp)
    colnames(comp) <- colnames(ml$X)

    ps[[4]] <- .passo(
      4, "Solucoes b com g-inversas diferentes",
      paste(
        "Usando g-inversas diferentes de X'X obtem-se solucoes b DIFERENTES",
        "(b = (X'X)^- X'y). A condicional e a de Moore-Penrose dao vetores",
        "distintos — mas os valores ajustados Xb sao identicos."
      ),
      objeto = round(comp, 4),
      verificacao = sprintf("Xb identico para as duas g-inversas: %s",
                            isTRUE(all.equal(aj_c, aj_mp)))
    )

    rest <- list()
    for (r in c("soma_zero", "casela_referencia")) {
      rest[[r]] <- tryCatch(
        .betas_restrito(ml$X, ml$y, .matriz_restricao(ml, r)),
        error = function(e) NULL
      )
    }
    rest <- Filter(Negate(is.null), rest)

    if (length(rest) > 0) {
      comp_r <- do.call(rbind, rest)
      colnames(comp_r) <- colnames(ml$X)
      aj_r <- as.numeric(ml$X %*% rest[[1]])
      ps[[5]] <- .passo(
        5, "Solucoes b sob restricoes",
        paste(
          "Impor uma restricao (soma zero ou casela de referencia) seleciona UMA",
          "solucao. Cada restricao da um b diferente, mas — de novo — os valores",
          "ajustados Xb sao os mesmos da solucao por g-inversa."
        ),
        objeto = round(comp_r, 4),
        verificacao = sprintf("Xb (restricao) = Xb (g-inversa): %s",
                              isTRUE(all.equal(aj_r, aj_c)))
      )
    } else {
      # celulas vazias: as restricoes canonicas nao identificam a solucao
      ps[[5]] <- .passo(
        5, "Solucoes b sob restricoes (indisponivel)",
        paste(
          "As restricoes canonicas (soma zero / casela) NAO identificam a solucao",
          "neste modelo — tipico de delineamentos com celulas vazias. Use a solucao",
          "por g-inversa e trabalhe apenas com funcoes estimaveis (veja eh_estimavel)."
        )
      )
    }
  }

  ps
}

#' Narra a deducao do modelo passo a passo
#'
#' Imprime os passos registrados por [llm()]: matriz X, sistema de equacoes
#' normais, posto/deficiencia, solucoes por g-inversas e por restricoes.
#'
#' @param objeto Objeto `llm` ou `permutacao_llm` (de [teste_permutacao()]).
#' @param qual Opcional. Indice(s) numerico(s) dos passos, ou um texto para
#'   filtrar pelo titulo (p.ex. `"restric"`, `"nula"`).
#' @return Invisivelmente, a lista de passos.
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 2)),
#'   y = c(10, 12, 20, 22, 30, 28)
#' )
#' m <- llm(y ~ trat, dados)
#' passos(m)
#' passos(m, qual = 1)
passos <- function(objeto, qual = NULL) {
  if (inherits(objeto, "permutacao_llm")) return(.passos_permutacao(objeto, qual))
  ml <- .como_llm(objeto)
  ps <- ml$passos
  alvo <- ps
  if (!is.null(qual)) {
    if (is.numeric(qual)) {
      alvo <- ps[qual]
    } else {
      alvo <- Filter(function(p) grepl(qual, p$titulo, ignore.case = TRUE), ps)
    }
  }
  for (p in alvo) {
    print(p)
    cat("\n")
  }
  invisible(ml$passos)
}
