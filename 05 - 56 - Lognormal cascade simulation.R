# -- -- -- -- -- ---- -- - - - - -- - - - --- -- - - -- - ---- - - -- - - - - - - -- --
# R Code 
#
# Title:      Multifractality in Asset Returns: Theory and Evidence
#             A. Fisher; L. Calvet (2002) 
# Date:       24-February-2022
# Author:     Raul Matsushita
# Adaptation: Gustavo Garcia
# - --- - - -- - -- - -- - - - - - --- -- - - - -- - - - - -- - -- - - - - - - - - -- --

# ------------------------------------------------------------------------------
# --- call required packages ---------------------------------------------------
# --- and functions ------------------------------------------------------------
# ------------------------------------------------------------------------------

library(ggplot2)
library(latex2exp)
library(Rcpp)
library(gridExtra)

# ------------------------------------------------------------------------------
# --- set working folder -------------------------------------------------------
# ------------------------------------------------------------------------------
# change the path accordingly #

#getwd()
#setwd("")

# ------------------------------------------------------------------------------
# The Log-Normal Measure on [0,1]
# ------------------------------------------------------------------------------

lambda.gbpbrl <- 1.189
sigma.gbpbrl <- sqrt(0.545)

set.seed(0)

R     <- 10 # Total iterations
range <- 1 # delta_t length
mass <- 1

k     <- 1
range <- 1/2


mass  <- rep(mass, each = 2) # Creating subintervals with same original measure
mass  <- mass*rlnorm(length(mass), -lambda.gbpbrl*log(2),  sigma.gbpbrl*log(2)) # Each subinterval will be multiplied by M, a log-normal r.v.
# Parameters from GBPBRL spectrum such that EM = 1/2.
# Each subinterval is likely to have its mass halved compared to previous interval, the total will be constant on average

for (k in 2:R){
  range <- range/2
  mass  <- rep(mass, each = 2)
  mass  <- mass*rlnorm(length(mass), -lambda.gbpbrl*log(2),  sigma.gbpbrl*log(2))
  print(sum(mass)) # sanity check, uncomment to check that mass is almost-constant on average
}

f.r        <- 1/range
y          <- mass*f.r

# --------------------------------------------------------------------------

# Cumulative distribution function

theta <- cumsum(y)/sum(y)
t <- (1:2^R)/(2^R)
plotdf <- data.frame(t = t, theta=theta, y=y)

gmed <- ggplot(data=plotdf,
               aes(x=t, y=theta)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  ylab(TeX(r'($\theta (t)$)'))+
  xlab(TeX(r'($t$)'))+
  ggtitle("Tempo de transação estocástico - cascata Log-Normal")

gmed

gk10 <- ggplot(data=plotdf,
               aes(x=t, y=y)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  ylab("Densidade")+
  xlab(TeX(r'($t$)'))+
  ggtitle("Medida Log-Normal")

gk10

grid.arrange(gk10, gmed, nrow=2)

# --------------------------------------
# Checking if simulation is multifractal
# --------------------------------------

# ---------------------------
# Scaling function
# ---------------------------

q <- seq(0.5, 1.5, 0.25)
dt <- 1:100
dados <- data.frame(Delta_t=dt)

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

X.theta <- log(theta) - log(theta[1]) # In the MMAR, X(t) = log(P(t)/P(0))

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.theta, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}

# Estimates of tau(q)
FuncaoEscala <- data.frame(q = q, tauObs = 0)
for(i in 1:length(q)){
  reg <- lm(log(dados[,i+1])~log(dados[,1]))
  FuncaoEscala$tauObs[i] <- reg$coefficients[2]
}
FuncaoEscala

# Graph
# This division makes all lines have the same intercept
for(i in 2:(length(q)+1)){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q05", "q075", "q1", "q125", "q150")


library(ggplot2)
library(reshape2)
library(latex2exp)
dados_long <- melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("0,5", "0,75", "1", "1,25", "1,5"), each=length(dt))

gSq <- ggplot(data=dados_long,
              aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  theme_bw()+
  ylab(TeX(r'($\log_{10} S_q$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$  (lags) )'))+
  ggtitle("Função de partição - Tempo de transação estocástico")

gSq