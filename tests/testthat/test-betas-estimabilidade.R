criar_oneway <- function() {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 2)),
    y = c(10, 12, 20, 22, 30, 28)
  )
  llm(y ~ trat, dados)
}

test_that("Xb e invariante a g-inversas e restricoes", {
  m <- criar_oneway()
  b_cond <- betas(m, "nenhuma", "condicional")
  b_mp   <- betas(m, "nenhuma", "moore_penrose")
  b_sz   <- betas(m, "soma_zero")
  b_cr   <- betas(m, "casela_referencia")

  expect_equal(attr(b_cond, "ajustados"), attr(b_mp, "ajustados"))
  expect_equal(attr(b_cond, "ajustados"), attr(b_sz, "ajustados"))
  expect_equal(attr(b_cond, "ajustados"), attr(b_cr, "ajustados"))

  expect_false(isTRUE(all.equal(as.numeric(b_sz), as.numeric(b_cr))))
})

test_that("restricoes impoem o que prometem", {
  m <- criar_oneway()
  b_sz <- betas(m, "soma_zero")
  b_cr <- betas(m, "casela_referencia")
  expect_lt(abs(sum(b_sz[c("tratA", "tratB", "tratC")])), 1e-8)
  expect_lt(abs(b_cr["tratA"]), 1e-8)
})

test_that("Xb reproduz as medias por tratamento", {
  m <- criar_oneway()
  medias <- tapply(m$dados$y, m$dados$trat, mean)
  alvo <- as.numeric(medias[m$dados$trat])
  expect_equal(attr(betas(m, "soma_zero"), "ajustados"), alvo, ignore_attr = TRUE)
})

test_that("estimabilidade distingue parametros de contrastes", {
  m <- criar_oneway()
  expect_false(eh_estimavel(c(1, 0, 0, 0), m))
  expect_false(eh_estimavel(c(0, 1, 0, 0), m))
  expect_true(eh_estimavel(c(1, 1, 0, 0), m))
  expect_true(eh_estimavel(c(0, 1, -1, 0), m))
})

test_that("projetor e simetrico e idempotente", {
  m <- criar_oneway()
  P <- projetor(m)
  expect_equal(P, t(P), ignore_attr = TRUE)
  expect_equal(P %*% P, P, ignore_attr = TRUE)
  expect_equal(posto(P), posto(m$X))
})
