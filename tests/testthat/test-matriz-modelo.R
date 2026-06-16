test_that("matriz_modelo devolve a X completa, nao a restrita", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 2)),
    y = c(10, 12, 20, 22, 30, 28)
  )
  X <- matriz_modelo(y ~ trat, dados)
  expect_equal(colnames(X), c("mu", "tratA", "tratB", "tratC"))
  expect_equal(ncol(X), 4)
  expect_equal(posto(X), 3)
  expect_equal(ncol(model.matrix(y ~ trat, dados)), 3)
})

test_that("interacoes geram todas as celulas", {
  dados <- expand.grid(a = factor(c("a1", "a2")), b = factor(c("b1", "b2", "b3")))
  dados$y <- seq_len(nrow(dados))
  X <- matriz_modelo(y ~ a + b + a:b, dados)
  expect_equal(ncol(X), 1 + 2 + 3 + 6)
  expect_equal(posto(X), 6)
})
