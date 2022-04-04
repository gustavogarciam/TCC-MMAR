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
    print(q) # comente para não ver a progressão
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

sim.garch <- ugarchsim(fit.garch, n.sim = 100000, m.sim = 4)

G <- sim.garch@simulation[["seriesSim"]]

X.garch1 <- cumsum(G[,1])
X.garch2 <- cumsum(G[,2])
X.garch3 <- cumsum(G[,3])
X.garch4 <- cumsum(G[,4])

# ---------------------------------------------------------------
# É multifractal? Teste: comportamento linear de escala (log-log)
# ---------------------------------------------------------------

# Parâmetros para a função de partição

dt <- 1:1000

q <- c(0.5, 1, 2, 3, 5)

tau.mb <- q/2 - 1
tau.obs <- c(-0.71599535, -0.45015292,  0.04493628,  0.49617300,  1.28476661)
#tau.teo <- c(-0.72508420, -0.45902868,  0.04650156,  0.51659071,  1.35044577)

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
    fpart[i] <- Sq(X.garch1, qi, t)
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

g.g1 <- ggplot(data=dados_long,
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
  ggtitle("GARCH - Simulação 1")

# ===========
# Gráfico 2 =
# ===========

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.garch2, qi, t)
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

g.g2 <- ggplot(data=dados_long,
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
  ggtitle("GARCH - Simulação 2")

# ===========
# Gráfico 3 #
# ===========

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.garch3, qi, t)
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

g.g3 <- ggplot(data=dados_long,
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
  ggtitle("GARCH - Simulação 3")


# ===========
# Gráfico 4 =
# ===========

dados <- data.frame(Delta_t=dt)

for(qi in q){
  fpart <- NULL
  i <- 1
  for(t in dt){
    fpart[i] <- Sq(X.garch4, qi, t)
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


g.g4 <- ggplot(data=dados_long,
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
  ggtitle("GARCH - Simulação 4")

# =================
# Resultado final =
# =================

grid.arrange(g.g1, g.g2, g.g3, g.g4, nrow = 2)
