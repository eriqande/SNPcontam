#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix P_likelihood(NumericMatrix gc, IntegerMatrix genos, double lambda) {
   int N = genos.ncol();
   int L = genos.nrow();
   int P = gc.nrow();
   double C = (2*N + 2*lambda)*(2*N + 1 + 2*lambda);
   double p;
   NumericVector like(P,1.0);
   NumericMatrix out(P,N);
   for(int i=0; i<N; i++){
     int j=0;
     for(j=0; j<P; j++) {
       int k = 0;
        for(k= 0; k<L; k++){
          int x1 = gc(j,k)*2*N;
          int x0 = 2*N - x1;
          int g = genos(k,i);
          if (g == 0) {p = (x0 + 1 + lambda)*(x0 + lambda);}
          if (g == 1) {p = 2*(x0 + lambda)*(x1 + lambda);}
          if (g == 2) {p = (x1 + 1 + lambda)*(x1 + lambda);}
          else {p = 1;}
          like[j] *= p/C;
        }
        out(j,i) = like[j];
      }
   }
  return out;
}