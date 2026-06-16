criar_aula17 <- function() {
  llm(y ~ trat, data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
  ))
}

test_that("teste_F avisa quando uma linha de L nao e estimavel", {
  m <- criar_aula17()
  expect_warning(teste_F(rbind(c(1, 0, 0, 0)), m), "estimaveis")
})

test_that("teste_F barra L com linhas linearmente dependentes", {
  m <- criar_aula17()
  expect_error(teste_F(rbind(c(0, 1, -1, 0), c(0, 2, -2, 0)), m),
               "linearmente independentes")
})

test_that("teste_F aceita L como vetor de strings de contraste", {
  m <- criar_aula17()
  por_string <- teste_F(c("trat1 - trat2", "trat2 - trat3"), m)
  por_matriz <- teste_F(rbind(c(0, 1, -1, 0), c(0, 0, 1, -1)), m)
  expect_equal(unname(por_string["F"]), unname(por_matriz["F"]))
})

test_that("resultados de inferencia tem classe e print proprios", {
  m <- criar_aula17()
  expect_s3_class(teste_t("trat1 - trat3", m), "inferencia_llm")
  expect_s3_class(intervalo_confianca("trat1 - trat3", m), "inferencia_llm")
  expect_output(print(teste_t("trat1 - trat3", m)), "Estimavel: TRUE")
})

test_that("betas aceita a g-inversa de minimos quadrados com Xb invariante", {
  m <- criar_aula17()
  b_mq <- betas(m, "nenhuma", "minimos_quadrados")
  b_c  <- betas(m, "nenhuma", "condicional")
  expect_equal(attr(b_mq, "ajustados"), attr(b_c, "ajustados"))
})

test_that(".betas_restrito barra restricoes insuficientes", {
  m <- criar_aula17()
  H_contraste <- rbind(c(0, 1, -1, 0))
  expect_error(.betas_restrito(m$X, m$y, H_contraste), "identificam")
})

test_that("anava respeita intercepto = FALSE (Total nao corrigido)", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
  )
  tab <- anava(llm(y ~ trat, dados, intercepto = FALSE))
  expect_equal(tab["Total", "SQ"], sum(dados$y^2))
  expect_equal(tab["Total", "SQ"], tab["trat", "SQ"] + tab["Residuo", "SQ"])
})

test_that("g-inversa e projetor ficam cacheados no objeto llm", {
  m <- criar_aula17()
  expect_equal(m$XtX %*% m$G %*% m$XtX, m$XtX, ignore_attr = TRUE)
  expect_equal(m$P, m$X %*% m$G %*% t(m$X))
})
