criar_split <- function() {
  sp <- expand.grid(variedade = factor(c("v1", "v2", "v3")),
                    bloco = factor(1:4),
                    irrig = factor(c("seco", "molhado")))
  sp$y <- c(53.1, 55.2, 57.0, 54.0, 56.1, 58.2, 52.2, 54.0, 56.1,
            55.0, 57.1, 59.0, 59.2, 61.0, 63.1, 60.1, 62.0, 64.2,
            58.1, 60.0, 62.1, 61.0, 63.2, 65.0)
  llm(y ~ bloco + irrig + bloco:irrig + variedade + irrig:variedade, sp)
}

test_that("erros = NULL nao muda a ANAVA atual", {
  m <- criar_split()
  expect_equal(unclass(anava(m)), unclass(anava(m, erros = NULL)))
})

test_that("anava com estrato de erro: F da parcela = QM(irrig)/QM(bloco:irrig)", {
  m <- criar_split()
  tab <- anava(m, erros = c(bloco = "bloco:irrig", irrig = "bloco:irrig"))

  F_irrig_manual <- tab["irrig", "QM"] / tab["bloco:irrig", "QM"]
  expect_equal(tab["irrig", "F"], F_irrig_manual)

  # termo usado como erro fica sem F (estrato de erro puro)
  expect_true(is.na(tab["bloco:irrig", "F"]))
  expect_equal(tab["bloco:irrig", "erro"], "(erro)")
  expect_equal(tab["irrig", "erro"], "bloco:irrig")

  # subparcela continua contra o residuo global
  expect_equal(tab["variedade", "F"], tab["variedade", "QM"] / tab["Residuo", "QM"])
  expect_equal(tab["variedade", "erro"], "Residuo")
})

# contraste de parcela (diferenca das medias marginais de irrig): media de linhas
# de X -> garantidamente estimavel; seu erro correto e o estrato de parcela.
lambda_irrig <- function(m) {
  seco <- m$dados$irrig == "seco"
  colMeans(m$X[seco, ]) - colMeans(m$X[!seco, ])
}

test_that("teste_t/teste_F com erro designado usam o QM do estrato", {
  m <- criar_split()
  lam <- lambda_irrig(m)
  expect_true(eh_estimavel(lam, m))

  tab <- anava(m, erros = c(irrig = "bloco:irrig"))
  QM_a <- tab["bloco:irrig", "QM"]
  gl_a <- tab["bloco:irrig", "GL"]

  tt <- teste_t(lam, m, erro = "bloco:irrig")
  expect_equal(unname(tt["gl"]), gl_a)

  # ep^2 = QM_erro * lambda'(X'X)^- lambda  -> recupera QM_erro
  var_unit <- as.numeric(t(lam) %*% m$G %*% lam)
  expect_equal(unname(tt["ep"])^2 / var_unit, QM_a)

  # teste_F de 1 contraste: F = t^2, e gl_den = gl do estrato
  tf <- teste_F(rbind(lam), m, erro = "bloco:irrig")
  expect_equal(unname(tf["F"]), unname(tt["t"])^2)
  expect_equal(unname(tf["gl_den"]), gl_a)
})

test_that("erro com nome de termo inexistente da mensagem clara", {
  m <- criar_split()
  lam <- lambda_irrig(m)
  expect_error(teste_t(lam, m, erro = "nao_existe"), "nao existe")
  expect_error(anava(m, erros = c(irrig = "nao_existe")), "inexistente")
  expect_error(anava(m, erros = c(nao_existe = "bloco:irrig")), "inexistente")
})
