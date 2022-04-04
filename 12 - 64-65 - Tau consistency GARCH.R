library(ggplot2)
library(latex2exp)
library(Rcpp)
library(gridExtra)
library(somebm)
library(readr)
library(dplyr)
library(lubridate)
library(reshape2)
library(rugarch)

# ----------------
# Leitura de dados
# ----------------

### USDDM ###

USDDM <- read_csv("USDDM.csv")

# Appendix 9.3: 6119 observations for weekdays from 1973-06-04 to 1996-12-31
# Except Christmas and New Year's Day

USDDM$weekday <- wday(USDDM$Date) # 7 - Saturday, 1 - Sunday

USDDM <- USDDM %>% filter(weekday != 7)
USDDM <- USDDM %>% filter(weekday != 1)

isChristmas <- month(USDDM$Date) == 12 & day(USDDM$Date) == 25
USDDM <- USDDM[!isChristmas,]
isNewYear <- month(USDDM$Date) == 1 & day(USDDM$Date) == 1
USDDM <- USDDM[!isNewYear,]

# 6119 observations, as desired
# Appendix 9.3, JPY/USD data: results not substantially different when
# i) replacing missing values with previous day's price or ii) eliminating missing values
# Here it is assumed the same applies to DM/USD data

for(i in 2:nrow(USDDM)){
  ifelse(is.na(USDDM[i,2]), USDDM[i,2] <- USDDM[i-1,2], 1)
}

price <- USDDM$USDDM

# --------------
# Funções usadas
# --------------

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
# Ajuste e simulação do GARCH
# ----------------------------

# Retornos

y <- diff(log(price), 1)


# Ajuste GARCH(1,1)

garch <- ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                    mean.model = list(armaOrder=c(0,0),include.mean = FALSE))
fit.garch = ugarchfit(data = y, spec = garch)

coef(fit.garch)

# Simulação

set.seed(0)

sim.garch <- ugarchsim(fit.garch, n.sim = 6118, m.sim = 10000)

G <- sim.garch@simulation[["seriesSim"]]

rm(sim.garch)

# -------------------------------------------------------
# É multifractal? Teste: consistência na estimação de tau
# -------------------------------------------------------

# Parâmetros para a função de partição

dt <- 1:180

q <- c(0.5, 1, 2, 3, 5)

tau.obs <- c(-0.71599535, -0.45015292,  0.04493628,  0.49617300,  1.28476661)
tau.teo <- c(-0.72508420, -0.45902868,  0.04650156,  0.51659071,  1.35044577)

# Tabela 3, Calvet & Fisher 2002, artigo pg 399

XG <- G
for(i in 1:10000){
  XG[,i] <- cumsum(G[,i])
}

rm(G)

zero <- rep(0, 10000)
estim.tau <- data.frame(q0.5 = zero, q1 = zero, q2 = zero, q3 = zero, q5 = zero)

for(i in 1:10000){
  estim.tau[i,] <- tau(q, XG[,i], dt)$tau_q
  print(i)
}

rm(XG)

colMeans(estim.tau)
tau.obs
tau.teo

colMeans(estim.tau) - tau.obs
colMeans(estim.tau) - tau.teo

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
