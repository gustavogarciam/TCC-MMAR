# Simulação de diversas trajetórias possíveis no MMAR

USDBRL <- readr::read_csv("USDBRL.csv")
priceB <- USDBRL$USDBRL[1:(5698-1024)]
priceF <- USDBRL$USDBRL[(5698-1023):5698]
X <- log(priceB) - log(priceB[1]) # X(t) = ln(P(t)/P(0))
dt <- 1:180

H <- 1.808^-1
alpha0 <- 0.6296443
lambda <- 1.138397
sigma2 <- 0.399329
s <- sqrt(sigma2)

genTheta <- function(lambda, s){
  R     <- 10
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
  return(c(0,somebm::fbm(H, 10^5)))
}

Rcpp::cppFunction('

NumericVector composeThetaBh(NumericVector theta, NumericVector Bh) {
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

escala <- sd(diff(X, lag = 1024))

nuvemSimulacoes <- matrix(0, nrow = 1024, ncol = 101)
nuvemSimulacoes[,1] <- 1:1024

set.seed(0)

for(simul in 1:100){
  nuvemSimulacoes[,simul+1] <- composeThetaBh(genTheta(lambda, s), genBh(H))
  print(simul) # comente para não ver a progressão
}

nuvemSimulacoes <- nuvemSimulacoes[,-1]
precoSimul <- priceB[4674]*exp(nuvemSimulacoes)
precoSimulreescala <- priceB[4674]*exp(escala*nuvemSimulacoes)



######

data <- as.Date(USDBRL$Date, "%d/%m/%Y")

plotsimuls <- data.frame(d = data)
plotsimuls$dolarAntes <- c(priceB, rep(NA, 1024))
plotsimuls$dolarDepois <- c(rep(NA, 4674), priceF)

for(sim in 1:100){
  plotsimuls[,sim+3] <- c(rep(NA, 4674), precoSimul[,sim])
}


plotsimulsescala <- data.frame(d = data)
plotsimulsescala$dolarAntes <- c(priceB, rep(NA, 1024))
plotsimulsescala$dolarDepois <- c(rep(NA, 4674), priceF)

for(sim in 1:100){
  plotsimulsescala[,sim+3] <- c(rep(NA, 4674), precoSimulreescala[,sim])
}


library(ggplot2)
library(reshape2)
library(gghighlight)

meltdf <- melt(plotsimuls,id="d")
g <- ggplot(meltdf)+
  geom_line(aes(d, value, colour = variable))+
  gghighlight((variable == "dolarDepois")  | (variable == "dolarAntes"), use_direct_label = FALSE)+
  theme_bw()+
  xlab("Data")+
  ylab("Cotação BRL/USD")+
  theme(legend.position="none")+
  ggtitle("Simulação MMAR sem ajuste de desvio padrão - BRL/USD")



meltdf <- melt(plotsimulsescala,id="d")
gesc <- ggplot(meltdf)+
  geom_line(aes(d, value, colour = variable))+
  gghighlight((variable == "dolarDepois")  | (variable == "dolarAntes"), use_direct_label = FALSE)+
  theme_bw()+
  xlab("Data")+
  ylab("Cotação BRL/USD")+
  theme(legend.position="none")+
  ggtitle("Simulação MMAR com ajuste de desvio padrão - BRL/USD")


library(gridExtra)
grid.arrange(g, gesc, nrow = 2)