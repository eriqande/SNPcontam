#include <Rcpp.h>
using namespace Rcpp;

//' Compute likelihoods for populations of origin 
//' 
//' More later
//' @param snp_zeroes Matrix of counts of the "zero" allele at L SNPs in P populations.  
//' This is an L x P integer matrix where L is the number of SNP
//' loci and P is the number of populations in the baseline. he observed number of 0 alleles
//' at locus j in population p is snp_zeroes[j,p]
//' @param snp_ones Same as above, but for the 1 allele at each SNP.
//' @param genos  Fill in.
//' @param lambda Fill in.
//' @export
// [[Rcpp::export]]
NumericMatrix P_likelihood(IntegerMatrix snp_zeroes, IntegerMatrix snp_ones, IntegerMatrix genos, double lambda) {
   int N = genos.ncol();
   int L = genos.nrow();
   int P = snp_zeroes.ncol();
   double C = (2*N + 2*lambda)*(2*N + 1 + 2*lambda);
   double p;
   NumericVector like(P,1.0);
   NumericMatrix out(P,N);
   for(int i=0; i<N; i++){
     int j=0;
     for(j=0; j<P; j++) {
       int k = 0;
        for(k= 0; k<L; k++){
          int x1 = snp_ones(k,j);
          int x0 = snp_zeroes(k,j);
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