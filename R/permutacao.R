# --------------------------------------------------------------------------- #
# Teste de permutacao (Freedman-Lane) para contrastes L'beta = m              #
# --------------------------------------------------------------------------- #

# ---- helpers internos ---------------------------------------------------- #

.proj_hipotese <- function(ml, L) {
  # Subesp. da hipotese: P_H = X G L' (L G X'X G L')^{-1} L G X'
  U <- ml$X %*% ml$G %*% t(L)       # n x q
  UtU <- t(U) %*% U
  P_H <- U %*% solve(UtU) %*% t(U)
  list(P_H = P_H, P_red = ml$P - P_H, U = U)
}

.proj_modelo_aux <- function(rotulos, resp, dados, intercepto) {
  formula <- if (length(rotulos) == 0) {
    stats::reformulate("1", response = resp)
  } else {
    stats::reformulate(rotulos, response = resp)
  }
  X <- matriz_modelo(formula, dados, intercepto)
  G <- inversa_condicional(t(X) %*% X)
  list(P = X %*% G %*% t(X), posto = posto(X))
}

.proj_strato <- function(ml, termo) {
  # A_k = P_k - P_{k-1}: projector da SQ sequencial do termo
  resp <- all.vars(ml$formula)[1]
  rotulos <- attr(stats::terms(ml$formula), "term.labels")
  k <- match(termo, rotulos)
  if (is.na(k)) {
    stop(sprintf("Termo de erro '%s' nao existe no modelo. Termos: %s.",
                 termo, paste(rotulos, collapse = ", ")), call. = FALSE)
  }
  pk  <- .proj_modelo_aux(rotulos[seq_len(k)], resp, ml$dados, ml$intercepto)
  pkm1 <- if (k == 1L) {
    .proj_modelo_aux(character(0), resp, ml$dados, ml$intercepto)
  } else {
    .proj_modelo_aux(rotulos[seq_len(k - 1L)], resp, ml$dados, ml$intercepto)
  }
  list(A = pk$P - pkm1$P, gl = pk$posto - pkm1$posto)
}

.estrato_grupos <- function(ml, erro) {
  if (is.null(erro) || identical(erro, "Residuo")) return(rep(1L, ml$n))
  fatores <- strsplit(erro, ":", fixed = TRUE)[[1]]
  faltando <- setdiff(fatores, names(ml$dados))
  if (length(faltando) > 0) {
    stop(sprintf("Fator(es) nao encontrado(s) nos dados: %s.",
                 paste(faltando, collapse = ", ")), call. = FALSE)
  }
  as.integer(interaction(ml$dados[fatores], drop = TRUE))
}

.permutar_dentro <- function(n, grupos) {
  idx <- seq_len(n)
  for (g in unique(grupos)) {
    pos <- which(grupos == g)
    if (length(pos) > 1L) idx[pos] <- pos[sample(length(pos))]
  }
  idx
}

.todas_permutacoes <- function(n) {
  # Retorna matrix (n! x n) com todas as permutacoes de 1:n (linhas)
  if (n == 1L) return(matrix(1L, 1L, 1L))
  sub <- .todas_permutacoes(n - 1L)
  ns  <- nrow(sub)
  out <- matrix(0L, ns * n, n)
  for (i in seq_len(n)) {
    linhas <- ((i - 1L) * ns + 1L):(i * ns)
    outros <- seq_len(n)[-i]
    out[linhas, i] <- n
    out[linhas, outros] <- sub
  }
  out
}

.stat_perm <- function(y_perm, LGXt, m, W_sc, Winv, q, M_denom, gl_denom) {
  Lb  <- as.numeric(LGXt %*% y_perm)
  d   <- Lb - m
  My  <- as.numeric(M_denom %*% y_perm)
  QM  <- sum(My^2) / gl_denom
  if (q == 1L) {
    d / sqrt(QM * W_sc)
  } else {
    as.numeric(t(d) %*% Winv %*% d) / (q * QM)
  }
}

# ---- passos de narração (chamados por passos()) --------------------------- #

.passos_permutacao <- function(res, qual = NULL) {
  ps <- list()

  # Passo 1 — Hipotese
  ps[[1]] <- .passo(
    1, "Hipotese: L'beta = m",
    paste(
      "A hipotese H0: L'beta = m define um subesp. LINEAR do espaco de parametros.",
      "Diferente do teste-F classico (que assume normalidade), a significancia sera",
      "avaliada pela DISTRIBUICAO das permutacoes, sem pressupostos distribucionais."
    )
  )

  # Passo 2 — Modelo reduzido e projetor de hipotese
  ps[[2]] <- .passo(
    2, "Projetor de hipotese P_H e modelo reduzido",
    paste(
      "Montamos P_H = U (U'U)^{-1} U' com U = X G L' — que projeta sobre o",
      "subesp. que H0 'reivindica'. O modelo reduzido (sob H0) usa",
      "P_red = P - P_H. Os valores ajustados pelo reduzido sao yhat_red = P_red y.",
      "",
      "Insight: P_H y eh a parte de y 'explicada' pelo contraste; P_red y eh o",
      "que sobra sem o contraste. So P_H y 'some' quando H0 eh verdade."
    )
  )

  # Passo 3 — Residuos permutaveis
  ps[[3]] <- .passo(
    3, "Residuos permutaveis e_null = (I - P_red) y",
    paste(
      "e_null = y - yhat_null inclui os residuos do modelo cheio (M y) E a",
      "componente da hipotese (P_H y). Sob H0 (L'beta = m), essa segunda parte",
      "e apenas ruido sem estrutura — logo e_null eh intercambivel entre",
      sprintf("%d unidade(s)", if (res$esquema$n_estratos == 1L) res$n_perm else sum(res$esquema$tamanhos[res$esquema$tamanhos > 1L])),
      if (res$esquema$n_estratos > 1L) {
        sprintf(
          "dentro de cada um dos %d estratos ('%s').",
          res$esquema$n_estratos, res$esquema$erro
        )
      } else {
        "(permutacao livre — sem estratos)."
      },
      "",
      "Para um contraste par-a-par (ex.: tau_A - tau_B = 0), permutar e_null",
      "equivale exatamente a trocar os rotulos A e B — a intuicao geometrica",
      "do teste de permutacao classico."
    )
  )

  # Passo 4 — Distribuicao nula e p-valor
  esquema_txt <- if (res$exato) {
    sprintf("EXATA (todas as %d permutacoes)", res$n_perm)
  } else {
    sprintf("Monte Carlo (%d permutacoes aleatorias)", res$n_perm)
  }
  ps[[4]] <- .passo(
    4, "Distribuicao nula e p-valor",
    paste(
      sprintf("Distribuicao %s.", esquema_txt),
      sprintf("Estatistica observada: %s = %.4f.", res$estatistica, res$stat_obs),
      sprintf("p-valor = %.4f.", res$pvalor),
      "",
      "O p-valor eh a proporcao de permutacoes que produzem uma estatistica tao",
      if (res$estatistica == "t") "extrema (em valor absoluto) quanto a observada."
      else                        "grande (F >= F_obs) quanto a observada.",
      "",
      if (!is.null(res$esquema$erro)) {
        sprintf(
          "Nota: permutacao restrita ao estrato '%s'. Para contrastes entre parcelas",
          res$esquema$erro
        )
      } else NULL,
      if (!is.null(res$esquema$erro)) {
        "em split-plot, a unidade permutavel e a parcela inteira (nao implementado)."
      } else NULL
    )
  )

  ps <- Filter(Negate(is.null), ps)
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
  invisible(ps)
}

# ---- funcao exportada ----------------------------------------------------- #

#' Teste de permutacao (Freedman-Lane) para H0: L beta = m
#'
#' Testa uma ou mais funcoes estimaveis sem pressupor normalidade dos erros.
#' O esquema de Freedman-Lane permuta os residuos do modelo reduzido (sob H0),
#' preservando os valores ajustados e — quando `erro` e informado — restringindo
#' a permutacao ao interior de cada estrato.
#'
#' @param L Contraste ou conjunto de contrastes: string (p.ex. `"trat1 - trat2"`),
#'   vetor de strings, ou matriz `q x p`. Mesma interface de [teste_F()].
#' @param objeto Objeto `llm`.
#' @param m Valor(es) sob H0 (padrao zeros). Vetor de comprimento `q`.
#' @param erro Nome de um termo do modelo a usar como estrato de erro (p.ex.
#'   `"bloco:irrig"` em split-plot). `NULL` (padrao) permuta livremente e usa o
#'   QME global como denominador.
#' @param n_perm Numero de permutacoes Monte Carlo (ignorado se a enumeracao
#'   exata for viavel).
#' @param seed Semente para reproducibilidade do Monte Carlo. `NULL` nao fixa.
#' @return Objeto de classe `permutacao_llm`. Use [print()], [plot()] e [passos()].
#' @export
#' @examples
#' dados <- data.frame(
#'   trat = factor(rep(c("A", "B", "C"), each = 6)),
#'   y = c(7, 9, 8, 10, 8, 9, 41, 36, 39, 42, 40, 38, 18, 18, 18, 17, 19, 20)
#' )
#' m <- llm(y ~ trat, dados)
#' teste_permutacao("trat1 - trat3", m, seed = 42)
#' teste_permutacao(c("trat1 - trat2", "trat2 - trat3"), m, seed = 42)
teste_permutacao <- function(L, objeto, m = NULL, erro = NULL,
                             n_perm = 2000, seed = NULL) {
  ml <- .como_llm(objeto)
  if (is.null(ml$y)) stop("O modelo nao tem resposta y.", call. = FALSE)

  expr <- if (is.character(L)) L else NULL
  L    <- .montar_L(ml, L)
  if (is.null(m)) m <- numeric(nrow(L))
  m    <- as.numeric(m)
  estimavel <- .validar_L(ml, L, "L")
  q    <- nrow(L)

  # ---- Projetores ----
  ph    <- .proj_hipotese(ml, L)
  P_H   <- ph$P_H
  P_red <- ph$P_red
  U     <- ph$U

  # ---- Quantidades pre-computadas ----
  W    <- L %*% ml$G %*% t(L)
  Winv <- if (q == 1L) matrix(1 / as.numeric(W)) else solve(W)
  W_sc <- if (q == 1L) as.numeric(W) else NULL   # escalar para o caso t
  LGXt <- L %*% ml$G %*% t(ml$X)                # q x n

  # ---- Fit nulo e residuos permutaveis ----
  # yhat_null = P_red y + U Winv m  (projecao sob H0: L'beta = m)
  yhat_null <- as.numeric(P_red %*% ml$y) + as.numeric(U %*% Winv %*% m)
  e_null    <- ml$y - yhat_null

  # ---- Denominador (QM do estrato escolhido) ----
  if (is.null(erro) || identical(erro, "Residuo")) {
    M_denom  <- diag(ml$n) - ml$P
    gl_denom <- ml$n - ml$posto
  } else {
    .qm_erro(ml, erro)          # valida existencia do termo
    ps_strato <- .proj_strato(ml, erro)
    M_denom  <- ps_strato$A
    gl_denom <- ps_strato$gl
  }

  # ---- Estatistica observada ----
  stat_obs <- .stat_perm(ml$y, LGXt, m, W_sc, Winv, q, M_denom, gl_denom)

  # ---- Estratos de permutacao ----
  grupos <- .estrato_grupos(ml, erro)
  if (all(tabulate(grupos) <= 1L)) {
    warning("Todos os estratos tem tamanho 1: permutacao impossivel.", call. = FALSE)
  }

  # ---- Exato vs Monte Carlo ----
  limite_exato <- 50000L
  n_un <- length(unique(grupos))
  exato <- (n_un == 1L) && (factorial(ml$n) <= limite_exato)

  if (exato) {
    perms <- .todas_permutacoes(ml$n)
    nula  <- apply(perms, 1L, function(idx) {
      .stat_perm(yhat_null + e_null[idx], LGXt, m, W_sc, Winv, q, M_denom, gl_denom)
    })
    n_perm_ef <- nrow(perms)
    pvalor <- if (q == 1L) mean(abs(nula) >= abs(stat_obs)) else mean(nula >= stat_obs)
  } else {
    # Preserva estado do RNG
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else NULL
    if (!is.null(seed)) set.seed(seed)

    nula <- replicate(n_perm, {
      idx <- .permutar_dentro(ml$n, grupos)
      .stat_perm(yhat_null + e_null[idx], LGXt, m, W_sc, Winv, q, M_denom, gl_denom)
    })

    if (!is.null(seed)) {
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = globalenv()))
          rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }
    n_perm_ef <- n_perm
    pvalor <- if (q == 1L) {
      (1 + sum(abs(nula) >= abs(stat_obs))) / (1 + n_perm)
    } else {
      (1 + sum(nula >= stat_obs)) / (1 + n_perm)
    }
  }

  res <- list(
    estatistica = if (q == 1L) "t" else "F",
    stat_obs    = stat_obs,
    pvalor      = pvalor,
    n_perm      = n_perm_ef,
    exato       = exato,
    nula        = nula,
    gl          = if (q == 1L) c(gl = gl_denom)
                  else c(gl_num = q, gl_den = gl_denom),
    esquema     = list(erro = erro, n_estratos = n_un,
                       tamanhos = tabulate(grupos)),
    expr        = expr,
    estimavel   = estimavel,
    q           = q,
    m           = m,
    n           = ml$n
  )
  class(res) <- "permutacao_llm"
  res
}

#' @export
print.permutacao_llm <- function(x, digits = 4, ...) {
  tipo <- if (x$q == 1L) "Teste t de permutacao (Freedman-Lane)"
          else            "Teste F de permutacao (Freedman-Lane)"
  cat(tipo, "\n")
  if (!is.null(x$expr)) cat("Funcao:", paste(x$expr, collapse = "; "), "\n")
  if (!is.null(x$estimavel)) cat("Estimavel:", isTRUE(x$estimavel), "\n")
  cat(sprintf("Estatistica observada: %s = %.4f\n", x$estatistica, x$stat_obs))
  cat(sprintf("p-valor: %.4f  (%s, n = %d permutacoes)\n",
              x$pvalor,
              if (x$exato) "exato" else "Monte Carlo",
              x$n_perm))
  cat(sprintf("Esquema: %s, %d estrato(s)\n",
              if (is.null(x$esquema$erro)) "permutacao livre"
              else sprintf("dentro de '%s'", x$esquema$erro),
              x$esquema$n_estratos))
  cat("Use plot(res) para ver a distribuicao nula, passos(res) para a narrativa.\n")
  invisible(x)
}
