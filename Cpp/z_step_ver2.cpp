#include <Rcpp.h>
using namespace Rcpp;

// Code adapted from Glynn et al. (2019)
// https://github.com/G-Lynn/DLTM


// [[Rcpp::export]]
Rcpp::List Zstep(NumericMatrix phi0, NumericMatrix W, NumericVector Z0, NumericMatrix theta0, NumericMatrix Y_kv, NumericMatrix X_dk ){

  int K = theta0.ncol();
//  int D = theta0.nrow();
//  int V = phi0.ncol();

  int nr = W.nrow();

  Rcpp::NumericVector U(nr);
  U = Rcpp::runif(nr);

  Rcpp::List Retx;





  for(int i=0; i<nr; i++){
    int term_index = W(i,0) - 1;
    int doc_index =  W(i,1) - 1;

    int top_index = Z0(i) - 1;

    Y_kv(top_index,term_index)--;
    X_dk(doc_index,top_index)--;

    Rcpp::NumericVector Prob(K);
    Rcpp::NumericVector CumProb(K);



    for(int k=0; k<K; k++){
        Prob(k) = theta0(doc_index, k)*phi0(k, term_index);
        if(k==0){
            CumProb(k) = Prob(k);
        }
        else{
            CumProb(k) = CumProb( k - 1) + Prob(k);
        }
    }



    // this is my own R::sample function based on uniform draws and the inverse CDF

    for(int k = 0; k<K; k++){
      CumProb(k) = CumProb(k)/CumProb(K - 1);

      if( (k == 0) && ( U(i)<CumProb(k) ) ){
        Z0(i) = k+1;
        break;
      }
      else if( k>0 && ( U(i) > CumProb(k-1) ) && ( U(i)<CumProb(k) ) ){
        Z0(i) = k+1;
        break;
      }
    }

    int top_index_new = Z0(i) - 1;
    Y_kv(top_index_new,term_index)++;
    X_dk(doc_index,top_index_new)++;
  }

  Retx.push_back(Z0);
  Retx.push_back(Y_kv);
  Retx.push_back(X_dk);


  return Retx;
}


