test_that("diagnostico de Box-Cox detecta 1 dentro do IC em dados normais", {
  set.seed(20260617)
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 8)),
    y = c(rnorm(8, 50, 5), rnorm(8, 60, 5), rnorm(8, 70, 5))
  )
  m <- expect_silent(llm(y ~ trat, dados))
  expect_true(m$boxcox$ok)
  expect_true(m$boxcox$contem_1)
})

test_that("dados com variancia crescente disparam warning e sugerem transformacao", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 5)),
    y = c(1, 2, 1, 3, 2, 12, 18, 15, 20, 14, 80, 95, 70, 110, 90)
  )
  expect_warning(llm(y ~ trat, dados), "Box-Cox")
  m <- suppressWarnings(llm(y ~ trat, dados))
  expect_false(m$boxcox$contem_1)
  expect_true(m$boxcox$tipo_sug %in% names(learnlm:::.redondos_boxcox))
})

test_that("checar_boxcox = FALSE pula o diagnostico", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B"), each = 4)),
    y = c(1, 2, 1, 3, 80, 95, 70, 110)
  )
  m <- expect_silent(llm(y ~ trat, dados, checar_boxcox = FALSE))
  expect_null(m$boxcox)
})

test_that("resposta nao-positiva nao avalia Box-Cox e nao avisa", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B"), each = 3)),
    y = c(-1, 0, 2, 5, 6, 7)
  )
  m <- expect_silent(llm(y ~ trat, dados))
  expect_false(m$boxcox$ok)
  expect_match(m$boxcox$motivo, "positiva")
})

test_that("transformar_resposta com tipo aplica a forma literal", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(1, 2, 1, 3, 12, 18, 15, 20, 80, 95, 70, 110)
  )
  m <- suppressWarnings(llm(y ~ trat, dados))
  m2 <- suppressWarnings(transformar_resposta(m, tipo = "log"))
  expect_equal(m2$transformacao$rotulo, "log(y)")
  expect_equal(m2$y, log(dados$y))
})

test_that("transformar_resposta com lambda aplica a forma de Box-Cox", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B"), each = 4)),
    y = c(1, 2, 3, 4, 12, 18, 15, 20)
  )
  m <- suppressWarnings(llm(y ~ trat, dados))
  m2 <- suppressWarnings(transformar_resposta(m, lambda = 0.5))
  expect_equal(m2$y, (dados$y^0.5 - 1) / 0.5)
  m0 <- suppressWarnings(transformar_resposta(m, lambda = 0))
  expect_equal(m0$y, log(dados$y))
})

test_that("transformar_resposta sem argumentos usa a sugestao do diagnostico", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 5)),
    y = c(1, 2, 1, 3, 2, 12, 18, 15, 20, 14, 80, 95, 70, 110, 90)
  )
  m <- suppressWarnings(llm(y ~ trat, dados))
  m2 <- suppressWarnings(transformar_resposta(m))
  expect_equal(m2$transformacao$tipo, m$boxcox$tipo_sug)
})

test_that("fix_boxcox = 'tipo' ja devolve o modelo corrigido", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 5)),
    y = c(1, 2, 1, 3, 2, 12, 18, 15, 20, 14, 80, 95, 70, 110, 90)
  )
  base <- suppressWarnings(llm(y ~ trat, dados))
  expect_message(llm(y ~ trat, dados, fix_boxcox = "tipo"), "corrigida")
  m <- suppressMessages(llm(y ~ trat, dados, fix_boxcox = "tipo"))
  expect_equal(m$transformacao$tipo, base$boxcox$tipo_sug)
})

test_that("fix_boxcox = 'lambda' usa o lambda_hat estimado", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 5)),
    y = c(1, 2, 1, 3, 2, 12, 18, 15, 20, 14, 80, 95, 70, 110, 90)
  )
  base <- suppressWarnings(llm(y ~ trat, dados))
  m <- suppressMessages(llm(y ~ trat, dados, fix_boxcox = "lambda"))
  expect_equal(m$transformacao$lambda, base$boxcox$lambda_hat)
})

test_that("fix_boxcox nao transforma quando o pressuposto e atendido", {
  set.seed(20260617)
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 8)),
    y = c(rnorm(8, 50, 5), rnorm(8, 60, 5), rnorm(8, 70, 5))
  )
  m <- expect_silent(llm(y ~ trat, dados, fix_boxcox = "tipo"))
  expect_null(m$transformacao)
})

test_that("fix_boxcox invalido da erro", {
  dados <- data.frame(trat = factor(c("A", "B")), y = c(1, 2))
  expect_error(llm(y ~ trat, dados, fix_boxcox = "raiz"),
               "NULL, 'tipo' ou 'lambda'")
})

test_that("tipo desconhecido e transformacao em y nao-positivo dao erro", {
  dados <- data.frame(
    trat = factor(rep(c("A", "B"), each = 3)),
    y = c(1, 2, 3, 4, 5, 6)
  )
  m <- suppressWarnings(llm(y ~ trat, dados))
  expect_error(transformar_resposta(m, tipo = "raiz_cubica"), "desconhecido")

  dados2 <- data.frame(
    trat = factor(rep(c("A", "B"), each = 3)),
    y = c(-1, 0, 2, 5, 6, 7)
  )
  m2 <- suppressWarnings(llm(y ~ trat, dados2))
  expect_error(transformar_resposta(m2, tipo = "log"), "y > 0")
})
