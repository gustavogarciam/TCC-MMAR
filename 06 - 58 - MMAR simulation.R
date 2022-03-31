# -----------
# Pacotes
# -----------

library(ggplot2)
library(latex2exp)
library(Rcpp)
library(gridExtra)
library(somebm)
library(readr)
library(dplyr)
library(lubridate)

# ----------------
# Leitura de dados
# ----------------

### GBPBRL ###

GBPBRL <- read_csv("D:/dataTCC/4-GBPBRL/GBPBRL.csv", col_names = FALSE)
price <- GBPBRL$X3

dt <- 10:20000
#H <- 2.156^-1

# -----------
# Log-Retorno
# -----------
X <- log(price) - log(price[1]) # X(t) = ln(P(t)/P(0))

# --------------------------------------------
# Função de partição Sq e Função de escala tau
# --------------------------------------------

Rcpp::cppFunction('

double Sq(NumericVector X, double q, int dt) {
  long double total = 0;
  int n = X.size();
  int k = 0;
  while (k + dt <= n-1) {
    total += pow(abs(X[k+dt] - X[k]), q);
    k += dt;
  }
  total += pow(abs(X[n-1] - X[k]), q);
  return total;
}

')

tau <- function(qs, X, dt){
  tq <- NULL
  i <- 1
  for(q in qs){
    teste <- NULL
    j <- 1
    for(t in dt){
      teste[j] <- Sq(X, q, t)
      j <- j+1
    }
    reg <- lm(log(teste) ~ log(dt))
    tq[i] <- reg$coefficients[2]
    i <- i+1
    print(q) # comente para não ver a progressão
  }
  result <- data.frame(q = qs, tau_q = tq)
  return(result)
}

tausimples <- function(q, X, dt){
  fpart <- rep(0, length(dt))
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X, q, t)
    i <- i+1
  }
  reg <- lm(log(fpart) ~ log(dt))
  return(reg$coefficients[2])
}

# ---------------
# Estimação de H
# ---------------

H <- (uniroot(tausimples, interval = c(0, 5), X = X, dt = dt)$root)^-1

# -------------------------------------------------------------------------
# Ajuste do espectro multifractal f: estimação de alpha_0, lambda e sigma^2
# -------------------------------------------------------------------------

q <- seq(0, 5, 0.05)

tauq <- tau(q,X,dt)

reg_tau <- lm(tauq$tau_q ~ poly(tauq$q,2, raw=T))

# Coeficientes da regressão
b0 <- reg_tau$coefficients[1]
b1 <- reg_tau$coefficients[2]
b2 <- reg_tau$coefficients[3]

# Parâmetros
alpha0 <- b1
lambda <- alpha0/H
sigma2 <- 2*(lambda-1)/log(2)
s <- sqrt(sigma2)

# ---------------------------------------------------
# Simulação do Tempo Estocástico de Negociação: theta
# ---------------------------------------------------

set.seed(0)

R     <- 10 # Iterações, a simulação será para o intervalo de tempo [0,T] com T = 2^R = 2^10
range <- 1 # Comprimento dos intervalos da cascata da medida
mass <- 1 # Massa da medida

k     <- 1
range <- 1/2


mass  <- rep(mass, each = 2) # Subintervalos com mesma medida original
mass  <- mass*rlnorm(length(mass), -log(2)*lambda, log(2)*s) 
# Cada subintervalo multiplicado por v.a. M log-normal
# Estes parâmetros se devem aos autores terem definido V = -log_b M como N(lambda, sigma2) em vez de V = ln M
# Parâmetros tais que EM = 1/b = 1/2
# Cada subintervalo tende a ter a medida cortada pela metade tal que o total será, em média, constante

# Teste de validade, descomente para executar
# mean(rlnorm(10^6, -log(2)*lambda, log(2)*s)) - 0.5

for (k in 2:R){
  range <- range/2
  mass  <- rep(mass, each = 2)
  mass  <- mass*rlnorm(length(mass), -log(2)*lambda, log(2)*s)
  #print(sum(mass)) # teste para ver se a massa parece ser constante em média
}

f.r        <- 1/range
y          <- mass*f.r

# theta é a F.D.A. da medida obtida
theta <- cumsum(y)/sum(y)

# ---------------------------------------------------------------------------
# Simulação do MMAR: Compondo theta com o Movimento Browniano Fracionário B_H
# ---------------------------------------------------------------------------

# theta[1:2^10]
# Com 1024 pontos, os primeiros valores de theta tendem a ser na 4ª casa decimal.
# Caso a medida tenha resultado em valores baixos nas primeiras observações é conveniente interpolar em [0,1] com variações de 0.00001 = 10^-5.

Rcpp::cppFunction('

NumericVector composeThetaBh(NumericVector Bh, NumericVector theta) {
  int n = theta.size();
  NumericVector result(n);
  int k;
  double t;
  for(int i = 0; i < n; i++) {
    k = floor(pow(10, 5)*theta[i]);
    t = pow(10, 5)*theta[i] - k;
    result[i] = (1-t)*Bh[k] + t*Bh[k+1];
  }
  return result;
}

')

# O C++ lida melhor com índice 0. P.ex., se theta[1] = 2x10^-6, 10^5*theta[1] = 0,2 e faremos
# 0,8 Bh[0] + 0,2 Bh[1]. O R tem problemas qdo índice da array é 0. Vamos adicionar o componente B_H (0) = 0 adiante.

Bh <- c(0,somebm::fbm(hurst = H, 10^5))
X.sim <- composeThetaBh(Bh, theta)


# --- Gráfico p/ os dados --- #
grafBh <- data.frame(t = 1:100002/100000, fbm = Bh)
gBh <- ggplot(data=grafBh, aes(x=t, y=fbm)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  ylab(TeX(r'($B_H (t)$)'))+
  xlab(TeX(r'($t$)'))+
  ggtitle("Movimento Browniano Fracionário")

gBh

grafXsim <- data.frame(t = 1:1024, Xs = X.sim)
gXsim <- ggplot(data=grafXsim, aes(x=t, y=Xs)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  ylab(TeX(r'($X_{sim} (t)$)'))+
  xlab(TeX(r'($t \, (Ticks)$)'))+
  ggtitle("Simulação MMAR - Câmbio BRL/GBP")

gXsim

grid.arrange(gBh, gXsim, nrow=2)
