#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List rw_MetropolisC(double x0, double sigma, int N) {
  NumericVector x(N);
  as<DoubleVector>(x)[0] = x0;
  NumericVector u(N);
  u = as<DoubleVector>(runif(N));
  List out(2);
  int k = 1;
  for(int i=1;i<(N-1);i++){
    double y = as<double>(rnorm(1,x[i-1],sigma));
    if(u[i] <= exp(abs(x[i-1])-abs(y))){
      x[i] = y;
      k+=1;
    }
    else{
      x[i] = x[i-1];
    }
  }  
  out[0] = x;
  out[1] = k;
  return(out);
}