library(ggplot2)
library(latex2exp)
library(Rcpp)
library(gridExtra)
library(somebm)
library(readr)
library(dplyr)
library(lubridate)
library(e1071)
library(reshape2)
library(rugarch)

# Dados do USDDM ###########################################

H <- 1.905^-1
alpha0 <- 0.558692
lambda <- 1.064308
sigma2 <- 0.1855545
s <- sqrt(sigma2)

# Funções usadas ############################################

genTheta_8192 <- function(lambda, s){
  R     <- 13
  range <- 1
  mass <- 1
  
  k     <- 1
  range <- 1/2
  
  mass  <- rep(mass, each = 2)
  mass  <- mass*rlnorm(length(mass), -log(2)*lambda, log(2)*s) 
  
  for (k in 2:R){
    range <- range/2
    mass  <- rep(mass, each = 2)
    mass  <- mass*rlnorm(length(mass), -log(2)*lambda, log(2)*s)
  }
  
  theta <- cumsum(mass)/sum(mass)
  
  return(theta)
}

genBh <- function(H){
  return(c(0,somebm::fbm(H, 10^6)))
}

Rcpp::cppFunction('

NumericVector composeThetaBh(NumericVector theta, NumericVector Bh) {
  int n = theta.size();
  NumericVector result(n);
  int k;
  double t;
  for(int i = 0; i < n; i++) {
    k = floor(pow(10, 6)*theta[i]);
    t = pow(10, 6)*theta[i] - k;
    result[i] = (1-t)*Bh[k] + t*Bh[k+1];
  }
  return result;
}

')

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
    #print(q) # comente para não ver a progressão
  }
  result <- data.frame(q = qs, tau_q = tq)
  return(result)
}


# ----------------------------
# Simulação do MMAR
# ----------------------------

# Simulação

set.seed(0)

XM <- matrix(0, nrow = 6118, ncol = 10000)

for(i in 1:10000){
  XM[,i] <- composeThetaBh(genTheta_8192(lambda, s), genBh(H))[1:6118]
  print(i) # Aviso: este processo é muito lento.
}


# -------------------------------------------------------
# É multifractal? Teste: consistência na estimação de tau
# -------------------------------------------------------

# Parâmetros para a função de partição

dt <- 1:180

q <- c(0.5, 1, 2, 3, 5)

tau.obs <- c(-0.71599535, -0.45015292,  0.04493628,  0.49617300,  1.28476661)
tau.teo <- c(-0.72508420, -0.45902868,  0.04650156,  0.51659071,  1.35044577)

# Tabela 3, Calvet & Fisher 2002, pdf pg 19 ou artigo pg 399

zero <- rep(0, 10000)
estim.tau <- data.frame(q0.5 = zero, q1 = zero, q2 = zero, q3 = zero, q5 = zero)

for(i in 1:10000){
  estim.tau[i,] <- tau(q, XM[,i], dt)$tau_q
  print(i)
}

colMeans(estim.tau)
tau.obs
tau.teo

colMeans(estim.tau) - tau.obs
colMeans(estim.tau) - tau.teo

rm(XG)

for(i in 1:5){
  print(sum(estim.tau[,i] < tau.obs[i])/10000)
}

percentis <- c(.5, 1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.5)/100

for(i in 1:5){
  print(quantile(estim.tau[,i], probs = percentis))
}

for(i in 1:5){
  print(sd(estim.tau[,i]))
}
