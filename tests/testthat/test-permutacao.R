dados_dic <- data.frame(
  trat = factor(rep(c("A", "B", "C"), each = 6)),
  y = c(7, 9, 8, 10, 8, 9, 41, 36, 39, 42, 40, 38, 18, 18, 18, 17, 19, 20)
)
m_dic <- suppressWarnings(llm(y ~ trat, dados_dic))

dados_sp <- data.frame(
  bloco  = factor(rep(1:4, each = 6)),
  irrig  = factor(rep(rep(c("Com", "Sem"), each = 3), times = 4)),
  variedade = factor(rep(c("V1", "V2", "V3"), times = 8)),
  y = c(28,32,24, 18,21,16, 30,35,25, 20,23,17,
        25,29,22, 15,19,14, 27,31,23, 17,21,15)
)
m_sp <- suppressWarnings(llm(y ~ bloco + irrig + bloco:irrig + variedade + irrig:variedade, dados_sp))

# ---- consistencia do projetor de hipotese --------------------------------- #

test_that("P_H satisfaz y'P_H y == SQ_H do Wald", {
  L   <- .montar_L(m_dic, "trat1 - trat3")
  ph  <- .proj_hipotese(m_dic, L)
  SQ_H_direto <- as.numeric(t(m_dic$y) %*% ph$P_H %*% m_dic$y)
  b   <- .solucao(m_dic)
  W   <- L %*% m_dic$G %*% t(L)
  Lb  <- as.numeric(L %*% b)
  SQ_H_wald <- as.numeric(Lb %*% (1 / W) %*% Lb)
  expect_equal(SQ_H_direto, SQ_H_wald, tolerance = 1e-8)
})

test_that("P_red e simetrico e idempotente", {
  L   <- .montar_L(m_dic, c("trat1 - trat2", "trat2 - trat3"))
  ph  <- .proj_hipotese(m_dic, L)
  P_r <- ph$P_red
  expect_true(isTRUE(all.equal(P_r, t(P_r), check.attributes = FALSE, tolerance = 1e-8)))
  expect_true(isTRUE(all.equal(P_r %*% P_r, P_r, check.attributes = FALSE, tolerance = 1e-8)))
})

# ---- p-valor ≈ parametrico sob dados normais ------------------------------ #

test_that("permutacao concorda com teste_F sob dados normais (balanceado)", {
  set.seed(20260618)
  dados_n <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 10)),
    y = c(rnorm(10, 50, 5), rnorm(10, 60, 5), rnorm(10, 70, 5))
  )
  m_n <- suppressWarnings(llm(y ~ trat, dados_n))
  res_p  <- teste_permutacao(c("trat1 - trat2", "trat2 - trat3"), m_n, n_perm = 2000, seed = 1)
  res_F  <- suppressWarnings(teste_F(c("trat1 - trat2", "trat2 - trat3"), m_n))
  # Ambos devem ser igualmente significativos (< 0.05)
  expect_lt(res_p$pvalor, 0.05)
  expect_lt(unname(res_F["pvalor"]), 0.05)
})

test_that("permutacao escalar concorda com teste_t (sinal e magnitude)", {
  res_p <- teste_permutacao("trat1 - trat3", m_dic, seed = 42)
  res_t <- suppressWarnings(teste_t("trat1 - trat3", m_dic))
  # direcao do efeito deve coincidir
  expect_equal(sign(res_p$stat_obs), sign(unname(res_t["t"])))
  # p-valor permutacao pequeno (efeito real)
  expect_lt(res_p$pvalor, 0.05)
})

# ---- reproducibilidade ---------------------------------------------------- #

test_that("mesmo seed produz mesmo p-valor", {
  p1 <- teste_permutacao("trat1 - trat2", m_dic, n_perm = 500, seed = 123)$pvalor
  p2 <- teste_permutacao("trat1 - trat2", m_dic, n_perm = 500, seed = 123)$pvalor
  expect_equal(p1, p2)
})

test_that("seeds diferentes podem dar p-valores diferentes", {
  p1 <- teste_permutacao("trat1 - trat2", m_dic, n_perm = 200, seed = 1)$nula
  p2 <- teste_permutacao("trat1 - trat2", m_dic, n_perm = 200, seed = 2)$nula
  expect_false(isTRUE(all.equal(p1, p2)))
})

# ---- enumeracao exata ----------------------------------------------------- #

test_that("enumeracao exata e ativada para n pequeno (estrato unico)", {
  dados_p <- data.frame(
    trat = factor(c("A", "A", "B", "B", "C", "C")),
    y = c(2.0, 2.5, 5.0, 6.0, 1.0, 1.5)
  )
  m_p <- suppressWarnings(llm(y ~ trat, dados_p))
  res <- suppressWarnings(teste_permutacao("trat1 - trat2", m_p))
  expect_true(res$exato)
  # n! = 720 <= 50000; n_perm deve ser 720
  expect_equal(res$n_perm, 720L)
  # p-valor deve ser racional com denominador 720
  expect_gte(res$pvalor, 0)
  expect_lte(res$pvalor, 1)
})

# ---- permutacao restrita ao estrato --------------------------------------- #

test_that("permutacao restrita nao mistura estratos (bloco:irrig)", {
  # Verifica que os grupos de permutacao sao correctos
  grupos <- .estrato_grupos(m_sp, "bloco:irrig")
  # deve haver 4 blocos x 2 irrigacoes = 8 estratos
  expect_equal(length(unique(grupos)), 8L)
  # dentro de cada estrato ha 3 variedades
  expect_true(all(tabulate(grupos) == 3L))

  # em cada permutacao, linhas de estratos distintos nao se misturam
  set.seed(99)
  idx <- .permutar_dentro(m_sp$n, grupos)
  for (g in unique(grupos)) {
    pos <- which(grupos == g)
    expect_true(all(idx[pos] %in% pos))
  }
})

test_that("teste com erro = bloco:irrig roda sem erro", {
  res <- suppressWarnings(
    teste_permutacao("variedade1 - variedade2", m_sp,
                     erro = "bloco:irrig", n_perm = 200, seed = 7)
  )
  expect_s3_class(res, "permutacao_llm")
  expect_equal(res$esquema$n_estratos, 8L)
  expect_gte(res$pvalor, 0)
  expect_lte(res$pvalor, 1)
})

# ---- interface e erros ---------------------------------------------------- #

test_that("L como vetor de strings e L como matriz dao o mesmo resultado", {
  r_str <- teste_permutacao(c("trat1 - trat2", "trat2 - trat3"),
                             m_dic, n_perm = 300, seed = 5)
  L <- do.call(rbind, lapply(c("trat1 - trat2", "trat2 - trat3"),
                              function(e) as.numeric(contraste(m_dic, e))))
  r_mat <- teste_permutacao(L, m_dic, n_perm = 300, seed = 5)
  expect_equal(r_str$stat_obs, r_mat$stat_obs, tolerance = 1e-10)
  expect_equal(r_str$pvalor,   r_mat$pvalor,   tolerance = 1e-10)
})

test_that("erro inexistente para", {
  expect_error(
    teste_permutacao("trat1 - trat2", m_dic, erro = "bloco"),
    "nao existe"
  )
})

test_that("modelo sem y para", {
  m_sem_y <- suppressWarnings(llm(~ trat, dados_dic))
  expect_error(teste_permutacao("trat1 - trat2", m_sem_y), "nao tem resposta")
})

test_that("L nao estimavel gera warning", {
  lambda_nao_estim <- numeric(m_dic$p)
  lambda_nao_estim[2] <- 1  # parametro individual: nao estimavel
  expect_warning(
    teste_permutacao(rbind(lambda_nao_estim), m_dic, n_perm = 50, seed = 1),
    "nao sao estimaveis"
  )
})

# ---- print e plot nao dao erro -------------------------------------------- #

test_that("print e plot funcionam sem erro", {
  res <- teste_permutacao("trat1 - trat3", m_dic, n_perm = 100, seed = 1)
  expect_output(print(res), "Freedman")
  expect_silent(plot(res))
})

test_that("passos() aceita permutacao_llm", {
  res <- teste_permutacao("trat1 - trat3", m_dic, n_perm = 100, seed = 1)
  expect_output(passos(res), "Freedman|hipotese|permut")
})
