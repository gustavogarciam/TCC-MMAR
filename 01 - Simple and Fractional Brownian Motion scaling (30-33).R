# Fisher, Calvet, Mandelbrot (1997):
# Multifractality of Deutschemark/US Dollar Exchange Rates

# Reproducing their Figure 2, scaling of Brownian Motion (Wiener process)
# H = 1/2, so ln Sq slope is expected to be (1/2)*q

# Section 2.3, page 6, says figures generated using T = 10^5 intervals
# Maximum Delta_t = 10^3 lags, 8 different q

q <- c(seq(1.5, 2.5, 0.25), 3:5)

dt <-1:1000

# How this works:
# Sq partitions [0,T] in N delta-t intervals, using continuous time.
# We have discrete data points, where X(0) = X[1] and X(T) = X[length(bm)]
# We want: length(bm)-1 = N*Delta_t, where N is a positive integer
# For delta_t = 1 day, we want Sq to take the non-overlapping intervals
# [1,2], [2,3], ... [len-1,len], total of N intervals.
# For delta_t = 2 days we want [1,3], [3,5], ... [len-2, len] and so on.
# The last interval will be X[length(bm)] - X[length(bm)-dt]

library(Rcpp)

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

# Generating a Brownian Motion
library(somebm)
set.seed(0)
x <- bm(n=10^5)

# For each q, storing data to plot the standardized partition function
dados <- data.frame(Delta_t=dt)
for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(x, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}

# Estimates of tau(q)
FuncaoEscala <- data.frame(q = q, tauObs = 0)
for(i in 1:8){
  reg <- lm(log(dados[,i+1])~log(dados[,1]))
  FuncaoEscala$tauObs[i] <- reg$coefficients[2]
}
FuncaoEscala$tauExpected <- H*q - 1
FuncaoEscala # Table 1

# Graph
# This division makes all lines have the same intercept
for(i in 2:9){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q15", "q175", "q2", "q225", "q25",
                  "q3", "q4", "q5")


library(ggplot2)
library(reshape2)
library(latex2exp)
dados_long <- melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("1,5", "1,75", "2", "2,25", "2,5", "3", "4", "5"), each=length(dt))

tau <- q/2 - 1

# Figure 2
ggplot(data=dados_long,
       aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  theme_bw()+
  geom_abline(intercept = 0, slope = tau[1], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[2], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[3], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[4], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[5], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[6], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[7], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[8], linetype = "dashed")+
  ylab(TeX(r'($\log_{10} (S_q)$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("Função de partição do movimento Browniano simples")

########################

####### Figure 4, for H = 0.4 and H = 0.8

####### Generating a Fractional Brownian Motion with H = 0.4
library(somebm)
set.seed(0)
H <- 0.4
x <- fbm(hurst = H, n=10^5)

# For each q, storing data to plot the standardized partition function
dados <- data.frame(Delta_t=dt)
for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(x, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}

# Estimates of tau(q)
FuncaoEscala <- data.frame(q = q, tauObs = 0)
for(i in 1:8){
  reg <- lm(log(dados[,i+1])~log(dados[,1]))
  FuncaoEscala$tauObs[i] <- reg$coefficients[2]
}
FuncaoEscala$tauExpected <- H*q - 1
FuncaoEscala # Table 2

# Graph
# This division makes all lines have the same intercept
for(i in 2:9){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q15", "q175", "q2", "q225", "q25",
                  "q3", "q4", "q5")


dados_long <- melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("1,5", "1,75", "2", "2,25", "2,5", "3", "4", "5"), each=length(dt))

tau <- H*q - 1

ggplot(data=dados_long,
       aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  theme_bw()+
  geom_abline(intercept = 0, slope = tau[1], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[2], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[3], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[4], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[5], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[6], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[7], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[8], linetype = "dashed")+
  ylab(TeX(r'($\log_{10} (S_q)$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("Função de partição do movimento Browniano fracionário, H = 0,4")

####### Generating a Fractional Brownian Motion with H = 0.8
set.seed(0)
H <- 0.8
x <- fbm(hurst = H, n=10^5)

# For each q, storing data to plot the standardized partition function
dados <- data.frame(Delta_t=dt)
for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(x, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}

# Estimates of tau(q)
FuncaoEscala <- data.frame(q = q, tauObs = 0)
for(i in 1:8){
  reg <- lm(log(dados[,i+1])~log(dados[,1]))
  FuncaoEscala$tauObs[i] <- reg$coefficients[2]
}
FuncaoEscala$tauExpected <- H*q - 1
FuncaoEscala # Table 3

# Graph
# This division makes all lines have the same intercept
for(i in 2:9){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q15", "q175", "q2", "q225", "q25",
                  "q3", "q4", "q5")

dados_long <- melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("1,5", "1,75", "2", "2,25", "2,5", "3", "4", "5"), each=length(dt))

tau <- H*q - 1

ggplot(data=dados_long,
       aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  theme_bw()+
  geom_abline(intercept = 0, slope = tau[1], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[2], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[3], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[4], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[5], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[6], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[7], linetype = "dashed")+
  geom_abline(intercept = 0, slope = tau[8], linetype = "dashed")+
  ylab(TeX(r'($\log_{10} (S_q)$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("Função de partição do movimento Browniano fracionário, H = 0,8")