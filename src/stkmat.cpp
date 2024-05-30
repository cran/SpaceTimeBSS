// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// computes the kernel matrix for a ring kernel
//
// [[Rcpp::export]]
arma::sp_mat stkmat_ring(const arma::mat & coords,
                         const arma::vec & time,
                         const double & h1,
                         const double & h2,
                         const double & tau) {
  const int n = coords.n_rows;
  double distance, lag;
  arma::sp_mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      lag = fabs(time(i) - time(j));
      if(h1 < distance && distance <= h2 && lag == tau){
        k(i,j) = k(j,i) = 1;
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for a ball kernel
//
// [[Rcpp::export]]
arma::sp_mat stkmat_ball(const arma::mat & coords,
                         const arma::vec & time,
                         const double & h,
                         const double & tau) {
  const int n = coords.n_rows;
  double distance, lag;
  arma::sp_mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      lag = fabs(time(i) - time(j));
      if(distance <= h && lag == tau){
        k(i,j) = k(j,i) = 1;
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for a Gauss kernel
//
// [[Rcpp::export]]
arma::sp_mat stkmat_gauss(const arma::mat & coords,
                          const arma::vec & time,
                          const double & h,
                          const double & tau) {
  const int n = coords.n_rows;
  double distance, lag;
  arma::sp_mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      lag = fabs(time(i) - time(j));
      if(lag == tau){
        k(i,j) = k(j,i) = exp(- 0.5 * pow(distance, 2) / h);
      }
    }
  }
  
  return k;
}
