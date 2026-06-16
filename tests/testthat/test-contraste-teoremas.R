criar_oneway <- function() {
  dados <- data.frame(
    trat = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
  )
  llm(y ~ trat, dados)
}

test_that("contraste por string monta o lambda certo", {
  m <- criar_oneway()
  expect_equal(as.numeric(contraste(m, "trat1 - trat3")), c(0, 1, 0, -1))
  expect_equal(as.numeric(contraste(m, "t1 - t3", simbolos = c(t = "trat"))),
               c(0, 1, 0, -1))
  expect_equal(as.numeric(contraste(m, "tratA - tratC")), c(0, 1, 0, -1))
  expect_equal(as.numeric(contraste(m, "2*trat1 - trat2 - trat3")),
               c(0, 2, -1, -1))
})

test_that("contraste de interacao por nome e por simbolo coincidem", {
  fat <- data.frame(
    tipo   = factor(rep(c("1", "2", "3"), each = 4)),
    metodo = factor(rep(rep(c("1", "2"), each = 2), times = 3)),
    y = c(39.02, 38.79, 38.96, 39.01, 35.74, 35.41,
          35.58, 35.52, 37.02, 36.00, 35.70, 36.04)
  )
  m <- llm(y ~ tipo + metodo + tipo:metodo, fat)
  por_nome <- contraste(m, "tipo1:metodo1 - tipo2:metodo1")
  por_simb <- contraste(m, "g11 - g21", simbolos = c(g = "tipo:metodo"))
  expect_equal(as.numeric(por_nome), as.numeric(por_simb))
})

test_that("funcoes de inferencia aceitam string de contraste", {
  m <- criar_oneway()
  expect_equal(unname(intervalo_confianca("trat1 - trat3", m)["estimativa"]), -9.25)
  expect_equal(unname(teste_t("trat1 - trat3", m)["estimativa"]), -9.25)
  expect_true(eh_estimavel("trat1 - trat3", m))
})

test_that("eh_contraste distingue contrastes de funcoes quaisquer", {
  m <- criar_oneway()
  expect_true(eh_contraste("trat1 - trat3", m))
  expect_false(eh_contraste("trat1 + trat3", m))
})

test_that("teoremas do modelo valem nos exemplos", {
  m <- criar_oneway()
  cons <- consistencia_sen(m)
  expect_true(cons$consistente)
  expect_equal(numero_funcoes_estimaveis(m), 3)
  expect_true(all(verificar_teoremas(m)$vale))
})
