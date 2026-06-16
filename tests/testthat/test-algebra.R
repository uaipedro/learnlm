test_that("g-inversas satisfazem suas propriedades", {
  A <- matrix(c(1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0), nrow = 4)

  Gc <- inversa_condicional(A)
  expect_equal(A %*% Gc %*% A, A, ignore_attr = TRUE)

  Gl <- inversa_minimos_quadrados(A)
  expect_equal(A %*% Gl %*% A, A, ignore_attr = TRUE)
  expect_equal(A %*% Gl, t(A %*% Gl), ignore_attr = TRUE)

  Gp <- inversa_moore_penrose(A)
  expect_equal(A %*% Gp %*% A, A, ignore_attr = TRUE)
  expect_equal(Gp %*% A %*% Gp, Gp, ignore_attr = TRUE)
})

test_that("posto identifica deficiencia", {
  expect_equal(posto(matrix(c(1, 2, 2, 4), 2)), 1)
  expect_equal(posto(diag(3)), 3)
})
