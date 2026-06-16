# learnlm

Pacote R **didático** para aprender modelos lineares, com foco em **modelos de
posto incompleto**. Cada etapa da construção de um modelo é um objeto que pode
ser calculado, verificado e explorado isoladamente — e também montado de uma vez
com `llm()`, que narra a dedução passo a passo.

Material-base: curso GES122 (Modelos Lineares II, UFLA).

## Instalação

Direto do GitHub (instala o `MASS`, única dependência, junto):

```r
# install.packages("remotes")
remotes::install_github("uaipedro/learnlm")
```

Ou, a partir de um clone local, na raiz do projeto:

```sh
R CMD INSTALL .
```

## A ideia central

`model.matrix()` devolve a matriz **já restrita** (com casela de referência
embutida). Para estudar modelos de posto incompleto precisamos da matriz **X
completa** (super-parametrizada), com uma coluna por nível. É isso que
`matriz_modelo()` faz:

```r
library(learnlm)

dados <- data.frame(
  bloco = factor(rep(1:3, times = 4)),
  trat  = factor(rep(c("A", "B", "C", "D"), each = 3)),
  y = c(42, 45, 41, 50, 53, 49, 47, 48, 46, 39, 41, 38)
)

matriz_modelo(y ~ bloco + trat, dados)   # X 12 x 8, posto 6 (incompleto)
model.matrix(y ~ bloco + trat, dados)    # já restrita, 6 colunas
```

## Fluxo guiado

```r
m <- llm(y ~ bloco + trat, dados)

m                       # resumo: posto, deficiência
passos(m)               # a dedução inteira, passo a passo
passos(m, qual = 1)     # só a matriz X
passos(m, qual = "restric")  # só o passo das restrições

m$X; m$XtX; m$Xty       # as peças, acessíveis direto
```

## Resolvendo o sistema de equações normais

```r
betas(m, restricao = "nenhuma", inversa = "condicional")
betas(m, restricao = "nenhuma", inversa = "moore_penrose")
betas(m, restricao = "soma_zero")
betas(m, restricao = "casela_referencia")
```

Cada chamada dá um vetor `b` **diferente**, mas os valores ajustados `Xb` são
**idênticos** — a invariância que está no coração da teoria. O `Xb` fica no
atributo `ajustados`.

## Contrastes por string

Escreva o contraste como você escreveria no quadro:

```r
contraste(m, "trat1 - trat3")                    # 1º nível − 3º nível de trat
contraste(m, "t1 - t3", simbolos = c(t = "trat"))   # com símbolo apelidado
contraste(m, "tratA - tratC")                    # por nome de coluna
contraste(m, "2*trat1 - trat2 - trat3")          # com coeficientes

# interações:
contraste(m, "tipo1:metodo1 - tipo2:metodo1")
contraste(m, "g11 - g21", simbolos = c(g = "tipo:metodo"))
```

O contraste já vem marcado como estimável ou não. E as funções de inferência
aceitam a string direto:

```r
intervalo_confianca("trat1 - trat3", m)
teste_t("trat1 - trat3", m)
```

## Inferência (Fase 2)

```r
qme(m)                                  # SQE, gl, s² = QME
intervalo_confianca(c(0, 1, 0, -1), m)  # IC para função estimável
teste_t(c(0, 1, 0, -1), m)              # H0: λ'β = 0
teste_F(L, m)                           # várias funções (Wald)

reg <- regiao_confianca(L, m)           # elipsoide (centro, W, constante)
reg$dentro(c(0, 0))                     # o ponto está na região?
plot(reg)                               # elipse (quando L tem 2 linhas)

forma_quadratica(y, A)                  # y'Ay
esperanca_forma_quadratica(A, mu)       # tr(AΣ) + μ'Aμ
forma_quadratica_qui2(A)                # condição de qui-quadrado (AΣ idempotente)
```

### Distribuição da função estimável e normal multivariada

```r
distribuicao_estimavel("trat1 - trat3", m)   # N(λ'b, s² λ'(X'X)⁻λ)

Y <- normal_multivariada(c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2))
densidade(Y, c(0, 0))                   # f(y)
transformar(Y, matrix(c(1, 1), 1))      # distribuição de A Y + b
marginal(Y, 1)                          # marginal das componentes indicadas
```

## ANAVA e reduções de SQ (Fase 3)

```r
anava(m)                                # tabela ANAVA (SQ por reduções R(·))
reducao(m, "trat", ajustado_por = "bloco")   # R(τ | μ, bloco)
```

`anava()` usa reduções sequenciais `R(·)` e reproduz `anova(lm(...))`. Serve para
DIC, DBC e fatoriais.

Quando há interação, `desdobramento()` estuda um fator dentro de cada nível do
outro (vale `Σ_j SQ(fator | dentro = j) = SQ(fator) + SQ(fator:dentro)`):

```r
mf <- llm(y ~ tipo + metodo + tipo:metodo, fat)
desdobramento(mf, "tipo", dentro_de = "metodo")
```

## Checagem de teoremas

```r
consistencia_sen(m)            # posto(X'X) = posto([X'X | X'y])
numero_funcoes_estimaveis(m)   # = posto(X)
eh_contraste("trat1 - trat3", m)
verificar_teoremas(m)          # bateria: consistência, projetor, invariância, gl
```

## Ferramentas chamáveis diretamente

Toda função interna é pública e funciona sobre matrizes cruas:

```r
inversa_condicional(A, explicar = TRUE)        # mostra a construção e verifica AGA = A
inversa_minimos_quadrados(A)
inversa_moore_penrose(A)                        # via MASS::ginv
posto(A)
projetor(m, explicar = TRUE)                    # P = X(X'X)^- X', verifica simetria/idempotência
eh_estimavel(c(0, 1, -1, 0), m)                 # contraste estimável?
```

## Escopo

- **Fase 1:** matriz X completa, SEN, g-inversas, soluções por restrições,
  estimabilidade, projetor.
- **Fase 2:** inferência — σ²/QME, IC e testes t/F para funções estimáveis,
  região de confiança, formas quadráticas.
- **Fase 3:** ANAVA e reduções de SQ (notação R(·)) — DIC, DBC, fatorial.
- **Transversal:** contrastes por string e checagem de teoremas.

Validado contra os oráculos numéricos do material (IC da Aula 17, fatorial das
Aulas 35-36) e contra `anova(lm())`. Veja
`docs/plans/2026-06-16-learnlm-design.md` para o desenho.
