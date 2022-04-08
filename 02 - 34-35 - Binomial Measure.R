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
# The Binomial Measure on [0,1] (stochastic case)
# ------------------------------------------------------------------------------

range <- 1 # delta_t length
mass <- 1
m.L   <- 0.6 # Mass shuffling proportion
m.R   <- 1 - m.L

# Let move.L = 1 if mass is to be transfered by m.L to the left and = 0 otherwise.
# For example, with move.L = [1,0] and m.L = 0.6 a mass [1, 1] will be spread as [0.6, 0.4, 0.4, 0.6].

binomShuffle <- function(move.L){
  n <- length(move.L)
  move.R <- 1 - move.L
  shuffle <- rep(0, 2*n)
  for(i in 1:n){ # Left is odd
    shuffle[2*i-1] <- m.L*move.L[i] + m.R*move.R[i]
  }
  for(i in 1:n){ # Right is even
    shuffle[2*i] <- m.R*move.L[i] + m.L*move.R[i]
  }
  return(shuffle)
}

# Test
binomShuffle(1)
binomShuffle(0)
binomShuffle(c(1,0))

k     <- 1
range <- 1/2

set.seed(0)
p <- 0.5
(mass  <- rep(mass, each = 2)*binomShuffle(rbinom(1, 1, p)) ) 
# Mass shared across intervals with 'range' length.
# Modify rbinom's 'p' parameter for unequal mass spliting probability.

R <- 4 # Number of iterations

for (k in 2:R){
  range <- range/2
  len.m <- length(mass)
  mass <- rep(mass, each = 2)*binomShuffle(rbinom(len.m, 1, p))
  # print(sum(mass)) # sanity check, uncomment to see that mass is conserved
}

f.r        <- 1/range
y          <- mass*f.r


# Graph until 4th iteration
R <- 4
range <- 1
mass <- 1

k     <- 1
range <- 1/2

set.seed(0)
p <- 0.5
(mass  <- rep(mass, each = 2)*binomShuffle(rbinom(1, 1, p)) ) 
# Mass shared across intervals with 'range' length.
# Modify rbinom's 'p' parameter for unequal mass spliting probability.

for (k in 2:R){
  range <- range/2
  len.m <- length(mass)
  mass <- rep(mass, each = 2)*binomShuffle(rbinom(len.m, 1, p))
  # print(sum(mass)) # sanity check, uncomment to see that mass is conserved
}

f.r        <- 1/range
y          <- mass*f.r

theta <- cumsum(y)/sum(y)
plot(theta)
t <- 1:16/16
plotdf <- data.frame(t = t, theta=theta, y=y)

gk4 <- ggplot(data=plotdf,
              aes(x=t, y=y)) +
  geom_step()+
  #geom_point()+
  theme_bw()+
  ylab("Densidade")+
  xlab(TeX(r'($t$)'))+
  ggtitle("Medida binomial - 4 iterações")


k     <- 1
range <- 1/2

set.seed(0)
p <- 0.5
(mass  <- rep(mass, each = 2)*binomShuffle(rbinom(1, 1, p)) ) 
# Mass shared across intervals with 'range' length.
# Modify rbinom's 'p' parameter for unequal mass spliting probability.

R <- 4 # Number of iterations

for (k in 2:R){
  range <- range/2
  len.m <- length(mass)
  mass <- rep(mass, each = 2)*binomShuffle(rbinom(len.m, 1, p))
  # print(sum(mass)) # sanity check, uncomment to see that mass is conserved
}

f.r        <- 1/range
y          <- mass*f.r


# Graph until 10th iteration
R <- 10
range <- 1
mass <- 1

k     <- 1
range <- 1/2

set.seed(0)
p <- 0.5
(mass  <- rep(mass, each = 2)*binomShuffle(rbinom(1, 1, p)) ) 
# Mass shared across intervals with 'range' length.
# Modify rbinom's 'p' parameter for unequal mass spliting probability.

for (k in 2:R){
  range <- range/2
  len.m <- length(mass)
  mass <- rep(mass, each = 2)*binomShuffle(rbinom(len.m, 1, p))
  # print(sum(mass)) # sanity check, uncomment to see that mass is conserved
}

f.r        <- 1/range
y          <- mass*f.r


theta <- cumsum(y)/sum(y)
plot(theta)
t <- 1:1024/1024
plotdf <- data.frame(t = t, theta=theta, y=y)

gmed <- ggplot(data=plotdf,
       aes(x=t, y=theta)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  ylab(TeX(r'($\theta (t)$)'))+
  xlab(TeX(r'($t$)'))+
  ggtitle("Tempo de transação estocástico - cascata binomial - 10 iterações")

gk10 <- ggplot(data=plotdf,
       aes(x=t, y=y)) +
  geom_step()+
  #geom_point()+
  theme_bw()+
  ylab("Densidade")+
  xlab(TeX(r'($t$)'))+
  ggtitle("Medida binomial - 10 iterações")



grid.arrange(gk4, gk10, gmed, nrow=3)



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

##################################
# tau(q) function for different q.

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

q <- c(0, 0.01, 0.05,
       seq(0.1, 1.9, 0.2),
       seq(1.91, 2.20, 0.01),
       seq(2.3, 5, 0.1),
       6:10,
       seq(12.5, 20, 2.5))
tauq <- tau(q,X.theta,dt)

# Replicating Figure 10
regfit <- data.frame(q = q, tauadj = predict(reg_tau))

gTau <- ggplot(data=tauq,
       aes(x=q, y=tau_q)) +
  geom_abline(intercept = -1, slope = 0.5, linetype = "dotted")+
  geom_path()+
  geom_point()+
  #ylim(-2, 10)+
  theme_bw()+
  ylab(TeX(r'($\tau_{\theta} (q)$)'))+
  xlab(TeX(r'($q$)'))+
  ggtitle("Função de escala - Tempo de transação estocástico")#+
#geom_line(color='red',data = regfit, aes(x=q, y=tauadj))

grid.arrange(gSq, gTau, nrow=2)
