library(MASS)
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

ok <- function(msg, cond) {
  cat(if (isTRUE(cond)) "  OK   " else "  FALHA", "-", msg, "\n")
  if (!isTRUE(cond)) stop("Verificacao falhou: ", msg, call. = FALSE)
}
perto <- function(a, b, tol = 1e-6) all(abs(a - b) < tol)

dados17 <- data.frame(
  trat = factor(rep(c("A", "B", "C"), each = 4)),
  y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
)
m17 <- llm(y ~ trat, dados17)

cat("== Contrastes por string ==\n")
c1 <- contraste(m17, "trat1 - trat3")
ok("trat1 - trat3 -> c(0,1,0,-1)", perto(as.numeric(c1), c(0, 1, 0, -1)))
ok("marcado como estimavel", isTRUE(attr(c1, "estimavel")))

c1b <- contraste(m17, "t1 - t3", simbolos = c(t = "trat"))
ok("simbolo customizado da o mesmo contraste", perto(as.numeric(c1), as.numeric(c1b)))

c1c <- contraste(m17, "tratA - tratC")
ok("nome de coluna direto da o mesmo contraste", perto(as.numeric(c1), as.numeric(c1c)))

c2 <- contraste(m17, "2*trat1 - trat2 - trat3")
ok("coeficientes: 2*trat1 - trat2 - trat3", perto(as.numeric(c2), c(0, 2, -1, -1)))

ok("intervalo_confianca aceita string direto",
   perto(intervalo_confianca("trat1 - trat3", m17)["estimativa"], -9.25))
ok("teste_t aceita string direto",
   perto(teste_t("trat1 - trat3", m17)["estimativa"], -9.25))

cat("\n== Contrastes em fatorial (interacao) ==\n")
fat <- data.frame(
  tipo   = factor(rep(c("1", "2", "3"), each = 4)),
  metodo = factor(rep(rep(c("1", "2"), each = 2), times = 3)),
  y = c(39.02, 38.79, 38.96, 39.01,
        35.74, 35.41, 35.58, 35.52,
        37.02, 36.00, 35.70, 36.04)
)
mF <- llm(y ~ tipo + metodo + tipo:metodo, fat)
g_nome <- contraste(mF, "tipo1:metodo1 - tipo2:metodo1")
g_simb <- contraste(mF, "g11 - g21", simbolos = c(g = "tipo:metodo"))
ok("interacao por nome de coluna == por simbolo g",
   perto(as.numeric(g_nome), as.numeric(g_simb)))
ok("contraste de interacao bem posicionado",
   sum(abs(as.numeric(g_nome))) == 2 &&
     g_nome["tipo1:metodo1"] == 1 && g_nome["tipo2:metodo1"] == -1)

cat("\n== Checagem de teoremas ==\n")
ok("eh_contraste: trat1 - trat3 (soma zero)", eh_contraste("trat1 - trat3", m17))
ok("eh_contraste: trat1 + trat3 NAO (soma != 0)", !eh_contraste("trat1 + trat3", m17))
ok("eh_estimavel aceita string", eh_estimavel("trat1 - trat3", m17))

cons <- consistencia_sen(m17)
ok("SEN consistente", cons$consistente)
ok("posto(X'X) = posto aumentada = 3",
   cons$posto_XtX == 3 && cons$posto_aumentada == 3)
ok("numero de funcoes estimaveis = posto(X) = 3",
   numero_funcoes_estimaveis(m17) == 3)

teo <- verificar_teoremas(m17)
ok("todos os teoremas valem no exemplo", all(teo$vale))

teoD <- verificar_teoremas(llm(y ~ bloco + trat, data.frame(
  bloco = factor(rep(1:3, times = 4)),
  trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
  y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
)))
ok("teoremas valem tambem no DBC", all(teoD$vale))

cat("\nTODAS AS VERIFICACOES DE CONTRASTES E TEOREMAS PASSARAM.\n")
