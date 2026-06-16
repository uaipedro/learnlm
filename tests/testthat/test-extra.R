test_that("normal multivariada: densidade, transformacao e marginal", {
  Y <- normal_multivariada(c(0, 0))
  expect_equal(densidade(Y, c(0, 0)), 1 / (2 * pi))
  expect_equal(densidade(Y, c(1, 2)), dnorm(1) * dnorm(2))

  Yc <- normal_multivariada(c(1, 2), matrix(c(2, 0.5, 0.5, 1), 2))
  soma <- transformar(Yc, matrix(c(1, 1), 1))
  expect_equal(soma$mu, 3)
  expect_equal(as.numeric(soma$Sigma), 4)
  expect_equal(marginal(Yc, 1)$mu, 1)
})

test_that("distribuicao_estimavel concorda com o IC", {
  m <- llm(y ~ trat, data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
  ))
  dist <- distribuicao_estimavel("trat1 - trat3", m)
  ic <- intervalo_confianca("trat1 - trat3", m)
  expect_equal(unname(dist$mu), unname(ic["estimativa"]))
  expect_equal(sqrt(as.numeric(dist$Sigma)), unname(ic["ep"]))
})

test_that("regiao de confianca: pontos na borda e predicado dentro", {
  m <- llm(y ~ trat, data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
  ))
  reg <- regiao_confianca(rbind(c(0, 1, -1, 0), c(0, 1, 0, -1)), m)
  Winv <- solve(reg$W)
  borda <- apply(reg$pontos, 1, function(p) {
    d <- as.numeric(p) - reg$centro
    as.numeric(t(d) %*% Winv %*% d)
  })
  expect_equal(borda, rep(reg$constante, length(borda)))
  expect_true(reg$dentro(reg$centro))
  expect_false(reg$dentro(reg$centro + c(1e4, 1e4)))
})

test_that("desdobramento satisfaz a identidade SQ(A) + SQ(AxB)", {
  fat <- data.frame(
    tipo   = factor(rep(c("1", "2", "3"), each = 4)),
    metodo = factor(rep(rep(c("1", "2"), each = 2), times = 3)),
    y = c(39.02, 38.79, 38.96, 39.01, 35.74, 35.41,
          35.58, 35.52, 37.02, 36.00, 35.70, 36.04)
  )
  m <- llm(y ~ tipo + metodo + tipo:metodo, fat)
  tab <- anava(m)

  dt <- desdobramento(m, "tipo", dentro_de = "metodo")
  semres <- rownames(dt) != "Residuo"
  expect_equal(sum(dt$SQ[semres]),
               tab["tipo", "SQ"] + tab["tipo:metodo", "SQ"], tolerance = 1e-4)
  expect_equal(sum(dt$GL[semres]), tab["tipo", "GL"] + tab["tipo:metodo", "GL"])

  dm <- desdobramento(m, "metodo", dentro_de = "tipo")
  expect_equal(sum(dm$SQ[rownames(dm) != "Residuo"]),
               tab["metodo", "SQ"] + tab["tipo:metodo", "SQ"], tolerance = 1e-4)
})
