#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix yUpdate(NumericMatrix X, NumericMatrix Y, NumericMatrix weights,
                      double K, double alpha, double beta){
  
  NumericMatrix delta (Y.nrow(), Y.ncol());
  int m;
  int n;
  int p;
  
  // loop through points y on the circle
  for(int j = 0; j < Y.nrow(); j++){
    if(j == 0){
      m = Y.nrow() - 1;
      n = j;
      p = j + 1;
    } else if (j == Y.nrow() - 1) {
      m = j - 1;
      n = j;
      p = 0;
    } else {
      m = j - 1;
      n = j;
      p = j + 1;
    }
    // loop through nodes in the data set and add to update matrix based on weight
    for(int i = 0; i < X.nrow(); i++){
      delta(j,_) = delta(j,_) + weights(i,j) * (X(i,_) - Y(j,_));
    }
    delta(j,_) = alpha * delta(j,_);
    delta(j,_) = delta(j,_) + (beta * K * (Y(m,_) - (2 * Y(n,_)) + Y(p,_)));
  }
  
  return(delta);
}

// [[Rcpp::export]]
double getDistance(NumericVector x, NumericVector y){
  double d = 0;
  for(int i = 0; i < x.size(); i++){
    d += std::pow(x[i]- y[i], 2.0);
  }
  d = std::sqrt(d);
  return(d);
}

// [[Rcpp::export]]
NumericMatrix weightUpdate(NumericMatrix X, NumericMatrix Y, double K){
  
  NumericMatrix deltas (X.nrow(), Y.nrow());
  NumericVector delta (Y.nrow());
  double phi;
  double distance;
  double total;
  
  for(int i = 0; i < X.nrow(); i++){
    total = 0;
    for(int j = 0; j < Y.nrow(); j++){
      distance = getDistance(X(i,_), Y(j,_));
      phi = std::exp(-std::pow(distance, 2.0)/(2 * std::pow(K, 2.0)));
      total += phi;
      delta[j] = phi;
    }
    delta = delta / total;
    deltas(i,_) = delta;
  }
  
  return(deltas);
}
