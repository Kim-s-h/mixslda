#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

// Code adapted from Glynn et al. (2019)
// https://github.com/G-Lynn/DLTM


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat theta_step(double sig1, double sig2, arma::vec b1, arma::vec b2, arma::vec resp, arma::vec mem, arma::mat theta_p, arma::mat theta1, arma::mat theta0, arma::mat X_dk ){

  int K = theta0.n_cols;
  int D = theta0.n_rows;

  arma::vec U = arma::randu(D);

  for(int d = 0; d<D; d++){

        int d_class = mem(d);

        arma::mat theta0d = theta0.row(d);
        arma::mat theta1d = theta1.row(d);

        arma::vec term1_n;
        arma::vec term1_d;
        arma::vec term2_n;
        arma::vec term2_d;
        arma::vec term3_n;
        arma::vec term3_d;
        arma::vec term4_n;
        arma::vec term4_d;

        arma::vec g_theta0(K);
        arma::vec g_theta1(K);

        arma::vec crit;

        for(int k = 0; k<K; k++){
            g_theta0(k) = tgamma(theta0d(k));
            g_theta1(k) = tgamma(theta1d(k));
        }

        if(d_class==1){
            term1_n = (X_dk.row(d) + theta_p.row(d) - 1)*arma::trans(arma::log(theta1d));
            term2_n = 0.5*(1/sig1)*(-2*resp(d)*(theta1d*b1) + (theta1d*b1)*(theta1d*b1));
            term3_n = log(prod(g_theta1));
            term4_n = (theta1d-1)*arma::trans(arma::log(theta0d));

            term1_d = (X_dk.row(d) + theta_p.row(d) - 1)*arma::trans(arma::log(theta0d));
            term2_d = 0.5*(1/sig1)*(-2*resp(d)*(theta0d*b1) + (theta0d*b1)*(theta0d*b1));
            term3_d = log(prod(g_theta0));
            term4_d = (theta0d-1)*arma::trans(arma::log(theta1d));

            crit = term1_n - term2_n - term3_n + term4_n - term1_d + term2_d + term3_d - term4_d;
        } else if (d_class==2){
            term1_n = (X_dk.row(d) + theta_p.row(d) - 1)*arma::trans(arma::log(theta1d));
            term2_n = 0.5*(1/sig2)*(-2*resp(d)*(theta1d*b2) + (theta1d*b2)*(theta1d*b2));
            term3_n = log(prod(g_theta1));
            term4_n = (theta1d-1)*arma::trans(arma::log(theta0d));

            term1_d = (X_dk.row(d) + theta_p.row(d) - 1)*arma::trans(arma::log(theta0d));
            term2_d = 0.5*(1/sig2)*(-2*resp(d)*(theta0d*b2) + (theta0d*b2)*(theta0d*b2));
            term3_d = log(prod(g_theta0));
            term4_d = (theta0d-1)*arma::trans(arma::log(theta1d));

            crit = term1_n - term2_n - term3_n + term4_n - term1_d + term2_d + term3_d - term4_d;
        }


        if (crit(0) > log(U(d))){
            theta0.row(d) = theta1d;
        } else {
            theta0.row(d) = theta0d;
        }

  }

  return theta0;

}

