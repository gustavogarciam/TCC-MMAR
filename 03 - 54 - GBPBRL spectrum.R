library(ggplot2)
library(reshape2)
library(latex2exp)

data.gbp <- read.csv("GBPBRL.csv",  header = FALSE)

price <- data.gbp$V3
X.gbp <- log(price) - log(price[1])


# Partition Function
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

# Scale Function
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
    print(q) # comment if you don't want to see progression
  }
  result <- data.frame(q = qs, tau_q = tq)
  return(result)
}

q <- seq(0, 5, 0.05)
dt <- 10:20000
tauq <- tau(q,X.gbp,dt)

(H <- 2.156^-1) # dt <- 10:20000

reg_tau <- lm(tauq$tau_q ~ poly(tauq$q,2, raw=T))
regfit <- data.frame(q = q, tauadj = predict(reg_tau))

b0 <- reg_tau$coefficients[1]
b1 <- reg_tau$coefficients[2]
b2 <- reg_tau$coefficients[3]

a <- seq(0, 1, 0.01)
fa <- -b0 + ((b1-a)^2)/(4*b2)
espectro <- data.frame(alpha = a, f = fa)

g1 <- ggplot(data=tauq,
       aes(x=q, y=tau_q)) +
  #geom_path()+
  geom_point()+
  theme_bw()+
  ylab(TeX(r'($\tau (q)$)'))+
  xlab(TeX(r'($q$)'))+
  ggtitle("Função de escala - câmbio BRL/GBP")+
  geom_line(color='red',data = regfit, aes(x=q, y=tauadj))
g1

g2 <- ggplot(data=espectro,
             aes(x=alpha, y=f)) +
  ylim(-0.001,1)+
  geom_line(color = 'red')+
  #geom_point()+
  theme_bw()+
  ylab(TeX(r'($f(\alpha)$)'))+
  xlab(TeX(r'($\alpha$)'))+
  ggtitle("Espectro Multifractal - câmbio BRL/GBP")

g2

library(gridExtra)
grid.arrange(g1, g2, nrow=2)


(alpha0 <- b1)
(lambda <- b1/H)
(sigma2 <- 2*(lambda-1)/log(2))

# sanity check
b0 - (-1)
b2 - (-H*(b1-H))

# sanity check: EM = 1/b?
set.seed(0)
M <- rlnorm(10^6, -log(2)*lambda, log(2)*(sigma2^0.5))
mean(M) - 1/2 # correto
