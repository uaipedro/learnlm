library(MASS)

arquivos <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(arquivos, source))

ok <- function(msg, cond) {
  cat(if (isTRUE(cond)) "  OK   " else "  FALHA", "-", msg, "\n")
  if (!isTRUE(cond)) stop("Verificacao falhou: ", msg, call. = FALSE)
}

cat("== Algebra das g-inversas ==\n")
set.seed(1)
A <- matrix(c(1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0), nrow = 4)
Gc <- inversa_condicional(A)
ok("condicional: A G A = A", all.equal(A %*% Gc %*% A, A, check.attributes = FALSE))
Gl <- inversa_minimos_quadrados(A)
ok("minimos quadrados: A G A = A", all.equal(A %*% Gl %*% A, A, check.attributes = FALSE))
ok("minimos quadrados: A G simetrica",
   all.equal(A %*% Gl, t(A %*% Gl), check.attributes = FALSE))
Gp <- inversa_moore_penrose(A)
ok("moore-penrose: A G A = A", all.equal(A %*% Gp %*% A, A, check.attributes = FALSE))
ok("moore-penrose: G A G = G", all.equal(Gp %*% A %*% Gp, Gp, check.attributes = FALSE))

cat("\n== Matriz X completa (one-way, I = 3) ==\n")
dados1 <- data.frame(
  trat = factor(rep(c("A", "B", "C"), each = 2)),
  y = c(10, 12, 20, 22, 30, 28)
)
X1 <- matriz_modelo(y ~ trat, dados1)
ok("X tem 4 colunas (mu, tratA, tratB, tratC)", ncol(X1) == 4)
ok("colunas nomeadas corretamente",
   identical(colnames(X1), c("mu", "tratA", "tratB", "tratC")))
ok("posto incompleto: posto(X) = 3 < p = 4", posto(X1) == 3)
ok("model.matrix daria 3 colunas (restrita) — X aqui da 4",
   ncol(stats::model.matrix(y ~ trat, dados1)) == 3 && ncol(X1) == 4)

cat("\n== llm e invariancia de Xb ==\n")
m1 <- llm(y ~ trat, dados1)
ok("deficiencia = 1", m1$deficiencia == 1)
b_cond <- betas(m1, "nenhuma", "condicional")
b_mp   <- betas(m1, "nenhuma", "moore_penrose")
b_sz   <- betas(m1, "soma_zero")
b_cr   <- betas(m1, "casela_referencia")
medias <- tapply(dados1$y, dados1$trat, mean)
aj_alvo <- as.numeric(medias[dados1$trat])
ok("Xb (condicional) = medias por tratamento",
   all.equal(attr(b_cond, "ajustados"), aj_alvo, check.attributes = FALSE))
ok("Xb invariante: condicional = moore-penrose",
   all.equal(attr(b_cond, "ajustados"), attr(b_mp, "ajustados")))
ok("Xb invariante: condicional = soma_zero",
   all.equal(attr(b_cond, "ajustados"), attr(b_sz, "ajustados")))
ok("Xb invariante: condicional = casela_referencia",
   all.equal(attr(b_cond, "ajustados"), attr(b_cr, "ajustados")))
ok("solucoes b sao DIFERENTES entre restricoes",
   !isTRUE(all.equal(as.numeric(b_sz), as.numeric(b_cr))))
ok("soma_zero: tratA+tratB+tratC = 0",
   abs(sum(b_sz[c("tratA", "tratB", "tratC")])) < 1e-8)
ok("casela_referencia: tratA = 0 (nivel de referencia)",
   abs(b_cr["tratA"]) < 1e-8)

cat("\n== Projetor ==\n")
P <- projetor(m1)
ok("P simetrica", all.equal(P, t(P), check.attributes = FALSE))
ok("P idempotente", all.equal(P %*% P, P, check.attributes = FALSE))
ok("posto(P) = posto(X)", posto(P) == posto(m1$X))
ok("Py = valores ajustados", all.equal(as.numeric(P %*% m1$y), aj_alvo, check.attributes = FALSE))

cat("\n== Estimabilidade ==\n")
ok("mu sozinho NAO estimavel", !eh_estimavel(c(1, 0, 0, 0), m1))
ok("tauA sozinho NAO estimavel", !eh_estimavel(c(0, 1, 0, 0), m1))
ok("mu + tauA estimavel", eh_estimavel(c(1, 1, 0, 0), m1))
ok("contraste tauA - tauB estimavel", eh_estimavel(c(0, 1, -1, 0), m1))

cat("\n== DBC (bloco + trat) ==\n")
dadosDBC <- data.frame(
  bloco = factor(rep(1:3, times = 4)),
  trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
  y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
)
mD <- llm(y ~ bloco + trat, dadosDBC)
ok("X completa tem 1 + 3 + 4 = 8 colunas", ncol(mD$X) == 8)
ok("posto(X) = I + J - 1 = 6", posto(mD$X) == 6)
ok("deficiencia = 2", mD$deficiencia == 2)
bD_sz <- betas(mD, "soma_zero")
bD_cr <- betas(mD, "casela_referencia")
ok("DBC: Xb invariante (soma_zero = casela_referencia)",
   all.equal(attr(bD_sz, "ajustados"), attr(bD_cr, "ajustados")))
ok("DBC contraste de tratamento estimavel",
   eh_estimavel(c(0, 0, 0, 0, 1, -1, 0, 0), mD))

cat("\n== Fatorial com interacao (a*b) ==\n")
dadosF <- expand.grid(a = factor(c("a1", "a2")), b = factor(c("b1", "b2", "b3")),
                      rep = 1:2)
dadosF$y <- c(5, 7, 6, 9, 8, 11, 6, 8, 7, 10, 9, 12)
mF <- llm(y ~ a + b + a:b, dadosF)
ok("X completa: 1 + 2 + 3 + 6 = 12 colunas", ncol(mF$X) == 12)
ok("posto(X) = numero de celulas = 6", posto(mF$X) == 6)
bF_sz <- betas(mF, "soma_zero")
bF_cr <- betas(mF, "casela_referencia")
ok("Fatorial: Xb invariante entre restricoes",
   all.equal(attr(bF_sz, "ajustados"), attr(bF_cr, "ajustados")))
medias_cel <- tapply(dadosF$y, list(dadosF$a, dadosF$b), mean)
ok("Xb = medias de celula",
   all.equal(sort(unique(round(attr(bF_sz, "ajustados"), 6))),
             sort(unique(round(as.numeric(medias_cel), 6)))))

cat("\nTODAS AS VERIFICACOES PASSARAM.\n")
