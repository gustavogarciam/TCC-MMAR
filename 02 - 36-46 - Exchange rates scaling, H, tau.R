# Based on Fisher, Calvet, Mandelbrot (1997):
# Multifractality of Deutschemark/US Dollar Exchange Rates

library(Rcpp)
library(readr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(reshape2)
library(latex2exp)

# Partition Function

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

# Data input. Adapt as necessary.

data <- read.csv(file.choose())
price <- data[,2]
X <- log(price)-log(price[1])

##############################
# Partition function Scaling #
##############################

dt <- 1:180 # Adapt as necessary
q <- c(seq(1.5, 2.5, 0.25), 3:5)

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X, qi, t)
    i <- i+1
  }
  print(qi) # Comment this out if you don't want to check progression.
  dados <- cbind(dados, fpart)
}

# This division makes all lines have the same intercept
for(i in 2:(length(q)+1)){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q15", "q175", "q2", "q225", "q25", "q3", "q4", "q5")

dados_long <- melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("1,5", "1,75", "2", "2,25", "2,5", "3", "4", "5"), each=length(dt))

# Replicating Figures 6/7.

ggplot(data=dados_long, aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  theme_bw()+
  ylab(TeX(r'($\log_{10} S_q$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("Função de partição câmbio aaa/bbb")



###########################################
# Scaling function graph and H estimation #
###########################################

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

# Using wide range of q-s to plot the scaling function tau(q)
# Section 5.1, pg 18, says 99 values for q were used, between 0.01 and 30
# Higher resolution near q = 2.

q <- c(0, 0.01, 0.05,
       seq(0.1, 1.9, 0.2),
       seq(1.90, 2.25, 0.05),
       seq(2.5, 5, 0.5),
       6:10,
       seq(12.5, 20, 2.5))
tauq <- tau(q,X,dt)

# H estimation
H <- (uniroot(tausimples, interval = c(0, 5), X = X, dt = dt)$root)^-1
tausimples(H^-1, X, dt)

# Replicating Figure 10.
ggplot(data=tauq,
       aes(x=q, y=tau_q)) +
  geom_abline(intercept = -1, slope = 0.5, linetype = "dotted")+ # simple Brownian motion scaling
  geom_path()+
  geom_point()+
  theme_bw()+
  ylab(TeX(r'($\tau (q)$)'))+
  xlab(TeX(r'($q$)'))+
  ggtitle("Função de escala - câmbio aaa/bbb")
