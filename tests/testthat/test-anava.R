test_that("ANAVA do DBC reproduz anova(lm)", {
  dados <- data.frame(
    bloco = factor(rep(1:3, times = 4)),
    trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
    y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
  )
  m <- llm(y ~ bloco + trat, dados)
  tab <- anava(m)
  ref <- anova(lm(y ~ bloco + trat, dados))

  expect_equal(tab["bloco", "SQ"], ref["bloco", "Sum Sq"])
  expect_equal(tab["trat", "SQ"], ref["trat", "Sum Sq"])
  expect_equal(tab["Residuo", "SQ"], ref["Residuals", "Sum Sq"])
  expect_equal(tab["trat", "F"], ref["trat", "F value"])
  expect_equal(c(tab["bloco", "GL"], tab["trat", "GL"], tab["Residuo", "GL"]),
               c(2, 3, 6))
})

test_that("ANAVA do fatorial reproduz o oraculo da Aula35-36", {
  fat <- data.frame(
    tipo   = factor(rep(c("1", "2", "3"), each = 4)),
    metodo = factor(rep(rep(c("1", "2"), each = 2), times = 3)),
    y = c(39.02, 38.79, 38.96, 39.01,
          35.74, 35.41, 35.58, 35.52,
          37.02, 36.00, 35.70, 36.04)
  )
  tab <- anava(llm(y ~ tipo + metodo + tipo:metodo, fat))
  expect_equal(tab["tipo", "SQ"], 25.9001, tolerance = 1e-3)
  expect_equal(tab["metodo", "SQ"], 0.1141, tolerance = 1e-3)
  expect_equal(tab["tipo:metodo", "SQ"], 0.3025, tolerance = 1e-3)
  expect_equal(tab["Residuo", "SQ"], 0.6620, tolerance = 1e-3)
  expect_equal(tab$GL[1:4], c(2, 1, 2, 6))
})

test_that("reducao parcial coincide com a SQ sequencial", {
  dados <- data.frame(
    bloco = factor(rep(1:3, times = 4)),
    trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
    y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
  )
  m <- llm(y ~ bloco + trat, dados)
  r <- reducao(m, "trat", ajustado_por = "bloco")
  expect_equal(r$R, anava(m)["trat", "SQ"])
  expect_equal(r$gl, 3)
})
