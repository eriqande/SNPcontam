#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// Vector x0 and x1 each of lenght (L x # of populations) these are
// the 0 allele and 1 alleles, resp. observed in each population at each
// locus.
// Let g = a vector of 0,1,2, or NAs of length L x # of individiuals
// Compute the likelihood of y(il) for an individual givin it came from
// a certain population P for all P.
// A matrix that nrow = # of population and ncol = number of individulas
// Computed using Beta-Binomial (P.M.F.)

// Use below code in R to get t
// library(fullsniplings)
// snp_genos <- get_snp_genos(sample_data)
// snp_indices <- genos_to_indicators(g = snp_genos$mat)
// geno_counts <- count_genos(snp_indices)


// [[Rcpp::export]]
NumericMatrix P_likelihood(NumericMatrix gc, NumericMatrix genos, double lambda) {
   int N = genos.ncol();
   int L = genos.nrow();
   int P = gc.nrow();
   double C = (2*N + 2*lambda)*(2*N + 1 + 2*lambda);
   double p
   NumericVector like(P,1)
   NumericMatric out
   for(int i=0; i<N; i++){
     for(int j=0; j<P; j++) {
        for(int k= 0; k<L; k++){
          int x1 = gc(j,k);
          int x0 = N - x1;
          if (genos(k,i) == 0) {p = (x0 + 1 + lambda)*(x0 + lambda);}
          if (genos(k,i) == 1) {p = 2*(x0 + lambda)*(x1 + lambda);}
          if (genos(k,i) == 2) {p = (x1 + 1 + lambda)*(x1 + lambda);}
          else {p = 1};
          like[j] *= p;
        }
        out(j,i) = like[j];
      }
   }
  return out
}
