library(MASS)
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

ok <- function(msg, cond) {
  cat(if (isTRUE(cond)) "  OK   " else "  FALHA", "-", msg, "\n")
  if (!isTRUE(cond)) stop("Verificacao falhou: ", msg, call. = FALSE)
}
perto <- function(a, b, tol = 1e-6) all(abs(a - b) < tol)

cat("== Normal multivariada ==\n")
Y <- normal_multivariada(c(0, 0))
ok("densidade N(0,I) em (0,0) = 1/(2pi)", perto(densidade(Y, c(0, 0)), 1 / (2 * pi)))
ok("densidade N(0,I) em (1,2) = dnorm(1)*dnorm(2)",
   perto(densidade(Y, c(1, 2)), dnorm(1) * dnorm(2)))

Yc <- normal_multivariada(c(1, 2), matrix(c(2, 0.5, 0.5, 1), 2))
soma <- transformar(Yc, matrix(c(1, 1), 1))
ok("transformar (soma): media = 3", perto(soma$mu, 3))
ok("transformar (soma): variancia = 1'Sigma1 = 4",
   perto(as.numeric(soma$Sigma), 2 + 1 + 2 * 0.5))
mg <- marginal(Yc, 1)
ok("marginal 1: N(1, 2)", perto(mg$mu, 1) && perto(as.numeric(mg$Sigma), 2))

cat("\n== Distribuicao de funcao estimavel ==\n")
dados17 <- data.frame(
  trat = factor(rep(c("A", "B", "C"), each = 4)),
  y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
)
m17 <- llm(y ~ trat, dados17)
dist <- distribuicao_estimavel("trat1 - trat3", m17)
ic <- intervalo_confianca("trat1 - trat3", m17)
ok("media = estimativa do IC", perto(dist$mu, ic["estimativa"]))
ok("desvio = ep do IC", perto(sqrt(as.numeric(dist$Sigma)), ic["ep"]))

cat("\n== Regiao de confianca (elipse) ==\n")
L <- rbind(c(0, 1, -1, 0), c(0, 1, 0, -1))
reg <- regiao_confianca(L, m17)
ok("regiao tem pontos da elipse", !is.null(reg$pontos))
W <- reg$W
Winv <- solve(W)
checa_borda <- apply(reg$pontos, 1, function(p) {
  d <- as.numeric(p) - reg$centro
  as.numeric(t(d) %*% Winv %*% d)
})
ok("todos os pontos estao na borda (forma quad = constante)",
   perto(checa_borda, reg$constante, 1e-6))
ok("centro esta dentro da regiao", reg$dentro(reg$centro))
ok("ponto muito distante esta fora", !reg$dentro(reg$centro + c(1000, 1000)))

cat("\n== Desdobramento da interacao ==\n")
fat <- data.frame(
  tipo   = factor(rep(c("1", "2", "3"), each = 4)),
  metodo = factor(rep(rep(c("1", "2"), each = 2), times = 3)),
  y = c(39.02, 38.79, 38.96, 39.01,
        35.74, 35.41, 35.58, 35.52,
        37.02, 36.00, 35.70, 36.04)
)
mF <- llm(y ~ tipo + metodo + tipo:metodo, fat)
tabF <- anava(mF)

dt <- desdobramento(mF, "tipo", dentro_de = "metodo")
sq_dt <- sum(dt$SQ[rownames(dt) != "Residuo"])
ok("identidade: soma SQ(tipo/metodo) = SQ(tipo) + SQ(tipo:metodo)",
   perto(sq_dt, tabF["tipo", "SQ"] + tabF["tipo:metodo", "SQ"], 1e-4))
gl_dt <- sum(dt$GL[rownames(dt) != "Residuo"])
ok("identidade dos GL: 2+2 = 4", gl_dt == tabF["tipo", "GL"] + tabF["tipo:metodo", "GL"])

dm <- desdobramento(mF, "metodo", dentro_de = "tipo")
sq_dm <- sum(dm$SQ[rownames(dm) != "Residuo"])
ok("identidade (outra direcao): soma SQ(metodo/tipo) = SQ(metodo) + SQ(tipo:metodo)",
   perto(sq_dm, tabF["metodo", "SQ"] + tabF["tipo:metodo", "SQ"], 1e-4))
ok("desdobramento usa o QME do modelo completo",
   perto(dt["Residuo", "QM"], qme(mF)$S2))

cat("\nTODAS AS VERIFICACOES EXTRAS PASSARAM.\n")
