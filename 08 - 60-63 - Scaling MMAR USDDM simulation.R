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

# Dados do USDDM ###########################################

H <- 1.905^-1
alpha0 <- 0.558692
lambda <- 1.064308
sigma2 <- 0.1855545
s <- sqrt(sigma2)

# Funções usadas ############################################

genTheta_100k <- function(lambda, s){
  R     <- 17
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
  return(c(0,somebm::fbm(H, 10^7)))
}

Rcpp::cppFunction('

NumericVector composeThetaBh(NumericVector theta, NumericVector Bh) {
  int n = theta.size();
  NumericVector result(n);
  int k;
  double t;
  for(int i = 0; i < n; i++) {
    k = floor(pow(10, 7)*theta[i]);
    t = pow(10, 7)*theta[i] - k;
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
    print(q) # comente para não ver a progressão
  }
  result <- data.frame(q = qs, tau_q = tq)
  return(result)
}

# Estimação ########################################################

set.seed(0)

X.sim <- composeThetaBh(genTheta_100k(lambda, s), genBh(H))[1:100000]

# É multifractal? Teste: comportamento linear de escala (log-log)

dt <- 1:1000

q <- c(0.5, 1, 2, 3, 5)

tau.mb <- q/2 - 1
tau.obs <- c(-0.71599535, -0.45015292,  0.04493628,  0.49617300,  1.28476661)
# tau.teo <- c(-0.72508420, -0.45902868,  0.04650156,  0.51659071,  1.35044577)

# Figura 9, Calvet & Fisher 2002, artigo pg 398
# Linhas pontilhadas: escala para o Mov Browniano Simples.
# Linhas tracejadas: escala com o MMAR e valor previsto com os dados originais (6119 dias, janela 1:180)


# ===========
# Gráfico 1 =
# ===========

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.sim, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}

# This division makes all lines have the same intercept
for(i in 2:(length(q)+1)){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q.5", "q1", "q2", "q3", "q5")

dados_long <- reshape2::melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("0,5", "1", "2", "3", "5"), each=length(dt))

g.mmar.1 <- ggplot(data=dados_long,
                   aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  geom_abline(intercept = 0, slope = tau.mb[1], linetype = "dotted")+ # pontilhadas: escala mov browniano simples
  geom_abline(intercept = 0, slope = tau.mb[2], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[3], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[4], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[5], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.obs[1], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[2], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[3], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[4], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[5], linetype = "longdash")+
  ylab(TeX(r'($\log_{10} S_q$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("MMAR - Simulação 1")

# ===========
# Gráfico 2 =
# ===========

X.sim <- composeThetaBh(genTheta_100k(lambda, s), genBh(H))[1:100000]

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.sim, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}


# This division makes all lines have the same intercept
for(i in 2:(length(q)+1)){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q.5", "q1", "q2", "q3", "q5")

dados_long <- reshape2::melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("0,5", "1", "2", "3", "5"), each=length(dt))

g.mmar.2 <- ggplot(data=dados_long,
                   aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  geom_abline(intercept = 0, slope = tau.mb[1], linetype = "dotted")+ # pontilhadas: escala mov browniano simples
  geom_abline(intercept = 0, slope = tau.mb[2], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[3], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[4], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[5], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.obs[1], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[2], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[3], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[4], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[5], linetype = "longdash")+
  ylab(TeX(r'($\log_{10} S_q$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("MMAR - Simulação 2")

# ===========
# Gráfico 3 #
# ===========

X.sim <- composeThetaBh(genTheta_100k(lambda, s), genBh(H))[1:100000]

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.sim, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}

# This division makes all lines have the same intercept
for(i in 2:(length(q)+1)){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q.5", "q1", "q2", "q3", "q5")

dados_long <- reshape2::melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("0,5", "1", "2", "3", "5"), each=length(dt))

g.mmar.3 <- ggplot(data=dados_long,
                   aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  geom_abline(intercept = 0, slope = tau.mb[1], linetype = "dotted")+ # pontilhadas: escala mov browniano simples
  geom_abline(intercept = 0, slope = tau.mb[2], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[3], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[4], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[5], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.obs[1], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[2], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[3], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[4], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[5], linetype = "longdash")+
  ylab(TeX(r'($\log_{10} S_q$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("MMAR - Simulação 3")


# ===========
# Gráfico 4 =
# ===========

X.sim <- composeThetaBh(genTheta_100k(lambda, s), genBh(H))[1:100000]

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.sim, qi, t)
    i <- i+1
  }
  dados <- cbind(dados, fpart)
}


# This division makes all lines have the same intercept
for(i in 2:(length(q)+1)){
  dados[,i] <- dados[,i]/dados[1,i]
}

names(dados) <- c("Delta_t", "q.5", "q1", "q2", "q3", "q5")

dados_long <- reshape2::melt(dados, id="Delta_t")
names(dados_long)[2] <- "q"
dados_long$q <- rep(c("0,5", "1", "2", "3", "5"), each=length(dt))


g.mmar.4 <- ggplot(data=dados_long,
                   aes(x=log10(Delta_t), y=log10(value), colour=q)) +
  geom_line()+
  #geom_point()+
  theme_bw()+
  geom_abline(intercept = 0, slope = tau.mb[1], linetype = "dotted")+ # pontilhadas: escala mov browniano simples
  geom_abline(intercept = 0, slope = tau.mb[2], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[3], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[4], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.mb[5], linetype = "dotted")+
  geom_abline(intercept = 0, slope = tau.obs[1], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[2], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[3], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[4], linetype = "longdash")+
  geom_abline(intercept = 0, slope = tau.obs[5], linetype = "longdash")+
  ylab(TeX(r'($\log_{10} S_q$)'))+
  xlab(TeX(r'($\log_{10} \Delta t$)'))+
  ggtitle("MMAR - Simulação 4")

# =================
# Resultado final =
# =================

grid.arrange(g.mmar.1, g.mmar.2, g.mmar.3, g.mmar.4, nrow = 2)
