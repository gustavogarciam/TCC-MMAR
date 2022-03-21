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

X <- 1:10
# Simple test.
# When dt = 4, partition should be [1,5], [5,9] and [9,10].
# With q = 2, result should be (5-1)^2 + (9-5)^2 + (10-9)^2 = 33
Sq(X, 2, 4)