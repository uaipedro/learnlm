library(MASS)
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

ok <- function(msg, cond) {
  cat(if (isTRUE(cond)) "  OK   " else "  FALHA", "-", msg, "\n")
  if (!isTRUE(cond)) stop("Verificacao falhou: ", msg, call. = FALSE)
}
perto <- function(a, b, tol = 1e-3) all(abs(a - b) < tol)

cat("== FASE 2: inferencia (oraculo Aula17) ==\n")
dados17 <- data.frame(
  trat = factor(rep(c("A", "B", "C"), each = 4)),
  y = c(7, 9, 8, 10, 41, 36, 39, 42, 18, 18, 18, 17)
)
m17 <- llm(y ~ trat, dados17)
v <- qme(m17)
ok("SQE = 26.75", perto(v$SQE, 26.75))
ok("gl = 9", v$gl == 9)
ok("S2 = 2.9722", perto(v$S2, 2.9722))

ic2 <- intervalo_confianca(c(0, 1, 0, -1), m17)
ok("IC tau1-tau3: estimativa = -9.25", perto(ic2["estimativa"], -9.25))
ok("IC tau1-tau3: ep = 1.219", perto(ic2["ep"], 1.219))
ok("IC tau1-tau3: limites [-12.008; -6.492]",
   perto(ic2["inferior"], -12.008, 0.01) && perto(ic2["superior"], -6.492, 0.01))

ic1 <- intervalo_confianca(c(0, 0.5, -1, 0.5), m17)
ok("IC (1/2)tau1 - tau2 + (1/2)tau3: estimativa = -26.375",
   perto(ic1["estimativa"], -26.375))
ok("IC limites [-28.763; -23.987]",
   perto(ic1["inferior"], -28.763, 0.01) && perto(ic1["superior"], -23.987, 0.01))

tt <- teste_t(c(0, 1, 0, -1), m17)
ok("teste t: t = estimativa/ep", perto(tt["t"], ic2["estimativa"] / ic2["ep"]))
ok("teste t: pvalor < 0.05", tt["pvalor"] < 0.05)

cat("\n== FASE 2: formas quadraticas ==\n")
set.seed(1)
A <- crossprod(matrix(rnorm(9), 3)) / 3
mu <- c(1, 2, 3)
ok("E[y'Ay] = tr(A Sigma) + mu'A mu",
   perto(esperanca_forma_quadratica(A, mu),
         sum(diag(A)) + as.numeric(t(mu) %*% A %*% mu)))
P <- projetor(m17)
M <- projetor(m17, residuos = TRUE)
ok("y'(I-P)y/sigma^2 ~ qui-quadrado (M idempotente)",
   forma_quadratica_qui2(M)$qui_quadrado)
ok("gl da forma residual = n - posto = 9", forma_quadratica_qui2(M)$gl == 9)
ok("P I (I-P) = 0 (independencia das formas)",
   perto(P %*% M, matrix(0, 12, 12)))

cat("\n== FASE 2: teste F == ANAVA (one-way) ==\n")
L <- rbind(c(0, 1, -1, 0), c(0, 0, 1, -1))
tF <- teste_F(L, m17)
tab17 <- anava(m17)
ok("teste_F(contrastes) reproduz F dos tratamentos",
   perto(tF["F"], tab17["trat", "F"]))

cat("\n== FASE 3: ANAVA DBC vs anova(lm) ==\n")
dadosDBC <- data.frame(
  bloco = factor(rep(1:3, times = 4)),
  trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
  y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
)
mD <- llm(y ~ bloco + trat, dadosDBC)
tabD <- anava(mD)
ref <- stats::anova(stats::lm(y ~ bloco + trat, dadosDBC))
ok("DBC: SQ(bloco) = anova(lm)", perto(tabD["bloco", "SQ"], ref["bloco", "Sum Sq"], 1e-6))
ok("DBC: SQ(trat) = anova(lm)", perto(tabD["trat", "SQ"], ref["trat", "Sum Sq"], 1e-6))
ok("DBC: SQ(residuo) = anova(lm)", perto(tabD["Residuo", "SQ"], ref["Residuals", "Sum Sq"], 1e-6))
ok("DBC: GL conferem",
   tabD["bloco", "GL"] == 2 && tabD["trat", "GL"] == 3 && tabD["Residuo", "GL"] == 6)
ok("DBC: F(trat) = anova(lm)", perto(tabD["trat", "F"], ref["trat", "F value"], 1e-5))
ok("DBC: linhas somam o Total corrigido",
   perto(tabD["bloco", "SQ"] + tabD["trat", "SQ"] + tabD["Residuo", "SQ"],
         tabD["Total", "SQ"]))

cat("\n== FASE 3: reducoes R(.) ==\n")
Rtab <- reducao(mD, "trat", ajustado_por = "bloco")
ok("R(tau | mu, bloco) = SQ(trat) sequencial", perto(Rtab$R, tabD["trat", "SQ"], 1e-6))
ok("R(tau | mu, bloco) tem gl = 3", Rtab$gl == 3)

cat("\n== FASE 3: fatorial (oraculo Aula35-36) ==\n")
fat <- data.frame(
  tipo   = factor(rep(c("1", "2", "3"), each = 4)),
  metodo = factor(rep(rep(c("1", "2"), each = 2), times = 3)),
  y = c(39.02, 38.79, 38.96, 39.01,
        35.74, 35.41, 35.58, 35.52,
        37.02, 36.00, 35.70, 36.04)
)
mF <- llm(y ~ tipo + metodo + tipo:metodo, fat)
tabF <- anava(mF)
ok("Fatorial: SQ(tipo) = 25.9001", perto(tabF["tipo", "SQ"], 25.9001, 1e-3))
ok("Fatorial: SQ(metodo) = 0.1141", perto(tabF["metodo", "SQ"], 0.1141, 1e-3))
ok("Fatorial: SQ(tipo:metodo) = 0.3025", perto(tabF["tipo:metodo", "SQ"], 0.3025, 1e-3))
ok("Fatorial: SQ(residuo) = 0.6620", perto(tabF["Residuo", "SQ"], 0.6620, 1e-3))
ok("Fatorial: GL = 2,1,2,6", all(tabF$GL[1:4] == c(2, 1, 2, 6)))
ok("Fatorial: F(tipo) = 117.38", perto(tabF["tipo", "F"], 117.381, 1e-2))
refF <- stats::anova(stats::lm(y ~ tipo + metodo + tipo:metodo, fat))
ok("Fatorial: bate com anova(lm) em todas as SQ",
   perto(tabF$SQ[1:4], c(refF[1:3, "Sum Sq"], refF["Residuals", "Sum Sq"]), 1e-6))

cat("\nTODAS AS VERIFICACOES DAS FASES 2 e 3 PASSARAM.\n")
