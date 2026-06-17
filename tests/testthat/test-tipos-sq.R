criar_desbalanceado <- function() {
  d <- data.frame(
    A = factor(c("a1","a1","a1","a1","a2","a2","a2","a3","a3","a3","a3","a3")),
    B = factor(c("b1","b2","b2","b1","b1","b2","b1","b1","b2","b2","b1","b2")),
    y = c(5, 7, 8, 6, 12, 15, 11, 20, 24, 23, 19, 25)
  )
  llm(y ~ A + B + A:B, d)
}

sq <- function(tab) {
  z <- as.data.frame(tab)
  stats::setNames(z$SQ, rownames(z))
}

test_that("tipo I (padrao) reproduz anova(lm) sequencial", {
  m <- criar_desbalanceado()
  alvo <- anova(lm(y ~ A + B + A:B, m$dados))
  s <- sq(anava(m))
  expect_equal(unname(s[c("A","B","A:B")]),
               unname(alvo[c("A","B","A:B"), "Sum Sq"]), tolerance = 1e-8)
})

test_that("tipo II e III batem com car::Anova", {
  skip_if_not_installed("car")
  m <- criar_desbalanceado()
  fit <- lm(y ~ A + B + A:B, m$dados,
            contrasts = list(A = contr.sum, B = contr.sum))

  s2 <- sq(anava(m, "II"))
  c2 <- car::Anova(fit, type = 2)
  tt <- intersect(names(s2), rownames(c2))
  expect_equal(unname(s2[tt]), unname(c2[tt, "Sum Sq"]), tolerance = 1e-6)

  s3 <- sq(anava(m, "III"))
  c3 <- car::Anova(fit, type = 3)
  tt3 <- intersect(names(s3), rownames(c3))
  expect_equal(unname(s3[tt3]), unname(c3[tt3, "Sum Sq"]), tolerance = 1e-6)
})

test_that("o termo de maior ordem tem a mesma SQ nos tres tipos", {
  m <- criar_desbalanceado()
  expect_equal(sq(anava(m, "I"))["A:B"], sq(anava(m, "II"))["A:B"])
  expect_equal(sq(anava(m, "II"))["A:B"], sq(anava(m, "III"))["A:B"])
})

test_that("em dados balanceados os tres tipos coincidem", {
  db <- expand.grid(A = factor(c("a1","a2","a3")), B = factor(c("b1","b2")), rep = 1:3)
  db$y <- c(11,13,14, 12,15,16, 13,14,15, 12,13,16, 14,15,17, 13,16,18)
  mb <- llm(y ~ A + B + A:B, db)
  expect_equal(sq(anava(mb, "I")), sq(anava(mb, "II")))
  expect_equal(sq(anava(mb, "II")), sq(anava(mb, "III")))
})

test_that("explicar anexa hipoteses e tipo", {
  m <- criar_desbalanceado()
  tab <- anava(m, "III")
  expect_equal(attr(tab, "tipo"), "III")
  expect_named(attr(tab, "hipoteses"), c("A", "B", "A:B"))
  expect_output(print(anava(m, "II", explicar = TRUE)), "Hipotese de cada fonte")
})

test_that("tipo III avisa e da NA com celula vazia", {
  m <- criar_desbalanceado()
  d2 <- m$dados[!(m$dados$A == "a3" & m$dados$B == "b2"), ]
  m2 <- llm(y ~ A + B + A:B, d2)
  expect_warning(tab <- anava(m2, "III"), "nao estimaveis")
  expect_true(is.na(sq(tab)["A"]))
})
