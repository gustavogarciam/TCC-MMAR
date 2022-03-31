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

dt <- 1:180
# H <- 1.905^-1
# alpha0 <- 0.558692
# lambda <- 1.064308
# sigma2 <- 0.1855545
# s <- sqrt(sigma2)

### USDBRL diário ###

USDBRL <- read_csv("USDBRL.csv")
price <- USDBRL$USDBRL

dt <- 1:180
# H <- 1.808^-1
# alpha0 <- 0.6296443
# lambda <- 1.138397
# sigma2 <- 0.399329
# s <- sqrt(sigma2)

### USDBRL intraday ###

USDBRL <- read_csv("USDBRLt.csv", col_names = FALSE)
price <- USDBRL$X3

dt <- 100:62000
# H <- 1.842^-1

### USDCHF diário ###

library(readr)
USDCHF <- read_csv("FRB_H10.csv")
price <- na.omit(USDCHF$USDCHF)

dt <- 1:150

### USDCHF intraday ###

setwd("~/USDCHF")
data.1 <- read.csv("DAT_MT_USDCHF_M1_2000.csv",  header = FALSE)
data.a <- read.csv("DAT_MT_USDCHF_M1_2001.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2002.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2003.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2004.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2005.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2006.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2007.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2008.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2009.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2010.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2011.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2012.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2013.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_2014.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_201501.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_201502.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_201503.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_201504.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_201505.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)
data.a <- read.csv("DAT_MT_USDCHF_M1_201506.csv",  header = FALSE)
data.1 <- rbind(data.1,data.a)

price <- data.1$V3

dt <- 100:180000
# H <- 2.301^-1
# alpha0 <- 0.5895289
# lambda <- 1.356506
# sigma2 <- 1.028659
# s <- sqrt(sigma2)

### EURBRL ###

EURBRL <- read_csv("EURBRL.csv", col_names = FALSE)
price <- EURBRL$X3

dt <- 10:20000
# H <- 2.404^-1

### GBPBRL ###

GBPBRL <- read_csv("GBPBRL.csv", col_names = FALSE)
price <- GBPBRL$X3

dt <- 10:20000
# H <- 2.156^-1

  
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
