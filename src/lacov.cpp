// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// local autocovariance with kernel matrix, sparse
//
// [[Rcpp::export]]
arma::mat lacov_kmat(const arma::mat & x,
                     arma::sp_mat k) {
  const int n = x.n_rows;
  const int p = x.n_cols;
  arma::mat L(p,p);
  L.zeros();
  
  arma::sp_mat::iterator it     = k.begin();
  arma::sp_mat::iterator it_end = k.end();
  
  for(; it != it_end; ++it)
  {
    L = L + (*it) * (x.row(it.row()).t()*x.row(it.col()));
  }  
  
  L = L / sqrt(n * accu(k % k));
  
  return L;
}

// ball kernel local autocovariance
// [[Rcpp::export]]
arma::mat lacov_ball(const arma::mat & coords,
                     const arma::vec & time,
                     const arma::mat & x,
                     const double & h,
                     const double & tau) {
  const int n = x.n_rows;
  const int p = x.n_cols;
  arma::mat L(p,p);
  L.zeros();
  double distance, lag, F;
  F = 0;
  
  if (tau == 0) {
    L = x.t() * x;  
    F = n;
  }
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      lag = fabs(time(i) - time(j));
      if(distance <= h && lag == tau){
        L += x.row(i).t() * x.row(j);
        L += x.row(j).t() * x.row(i);
        F += 2;
      }
    }
  }
  
  L = L / sqrt(n * F);
  
  return L;
}


// ring kernel local autocovariance
// [[Rcpp::export]]
arma::mat lacov_ring(const arma::mat & coords,
                     const arma::vec & time,
                     const arma::mat & x,
                     const double & h1,
                     const double & h2,
                     const double & tau) {
  const int n = x.n_rows;
  const int p = x.n_cols;
  arma::mat L(p,p);
  L.zeros();
  double distance, lag, F;
  F = 0;
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      lag = fabs(time(i) - time(j));
      if(h1 < distance && distance <= h2 && lag == tau){
        L += x.row(i).t() * x.row(j);
        L += x.row(j).t() * x.row(i);
        F += 2;
      }
    }
  }
  
  L = L / sqrt(n * F);
  
  return L;
}


// Gauss kernel local autocovariance
// [[Rcpp::export]]
arma::mat lacov_gauss(const arma::mat & coords,
                      const arma::vec & time,
                      const arma::mat & x,
                      const double & h,
                      const double & tau) {
  const int n = x.n_rows;
  const int p = x.n_cols;
  arma::mat L(p,p);
  L.zeros();
  double distance, lag, F, f;
  F = 0;
  
  if (tau == 0) {
    L = x.t() * x;  
    F = n;
  }
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      lag = fabs(time(i) - time(j));
      if(lag == tau){
        f = exp(- 0.5 * pow(distance, 2) / h);
        L += f * x.row(i).t() * x.row(j);
        L += f * x.row(j).t() * x.row(i);
        F += 2 * pow(f, 2);
      }
    }
  }
  
  L = L / sqrt(n * F);
  
  return L;
}
