criar_aula17 <- function() {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
  )
  llm(y ~ trat, dados)
}

test_that("qme reproduz o oraculo da Aula17", {
  v <- qme(criar_aula17())
  expect_equal(v$SQE, 26.75)
  expect_equal(v$gl, 9)
  expect_equal(v$S2, 26.75 / 9)
})

test_that("intervalo de confianca bate com o oraculo da Aula17", {
  m <- criar_aula17()
  ic <- intervalo_confianca(c(0, 1, 0, -1), m)
  expect_equal(unname(ic["estimativa"]), -9.25)
  expect_equal(unname(ic["inferior"]), -12.008, tolerance = 1e-3)
  expect_equal(unname(ic["superior"]), -6.492, tolerance = 1e-3)

  ic1 <- intervalo_confianca(c(0, 0.5, -1, 0.5), m)
  expect_equal(unname(ic1["estimativa"]), -26.375)
})

test_that("teste_F com contrastes reproduz o F dos tratamentos", {
  m <- criar_aula17()
  L <- rbind(c(0, 1, -1, 0), c(0, 0, 1, -1))
  expect_equal(unname(teste_F(L, m)["F"]),
               anava(m)["trat", "F"], tolerance = 1e-8)
})

test_that("esperanca de forma quadratica e condicao qui-quadrado", {
  A <- crossprod(matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 1), 3)) / 3
  mu <- c(1, 2, 3)
  expect_equal(esperanca_forma_quadratica(A, mu),
               sum(diag(A)) + as.numeric(t(mu) %*% A %*% mu))

  M <- projetor(criar_aula17(), residuos = TRUE)
  cq <- forma_quadratica_qui2(M)
  expect_true(cq$qui_quadrado)
  expect_equal(cq$gl, 9)
})
