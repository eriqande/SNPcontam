#include <Rcpp.h>
using namespace Rcpp;

//' Compute likelihoods for populations of origin 
//' 
//' @param snp_zeroes Matrix of counts of the "zero" allele at L SNPs in P populations.  
//' This is an L x P integer matrix where L is the number of SNP
//' loci and P is the number of populations in the baseline. The observed number of 0 alleles
//' at locus j in population p is snp_zeroes[j,p]
//' @param snp_ones Same as above, but for the 1 allele at each SNP.
//' @param genos  Matrix of the genotypes of N individuals at each of the L loci.  Genotypes are represented
//' as the number of "1" alleles.  This is an L x N integer matrix.  The genotype of the ith individual at 
//' the locus j is genos(j,i).
//' @param lambda The parameter of the beta distribution for the frequency of the "1" allele.
//' 
//' @return Returns a P x N numeric matrix which contains the the probabilities the individuals originate from
//' certain populations.  The probability that the ith originated from population p is ouput(p,i).   
//' @export
// [[Rcpp::export]]
NumericMatrix P_likelihood(IntegerMatrix snp_zeroes, IntegerMatrix snp_ones, IntegerMatrix genos, double lambda) {
   int N = genos.ncol();
   int L = genos.nrow();
   int P = snp_zeroes.ncol();
   double C;
   double p;
   double like;
   double new_p;
   int g;
   NumericMatrix out(P,N);
   for(int i=0; i<N; i++){
     int j=0;
     for(j=0; j<P; j++) {
       int k = 0;
       like = 1;
        for(k= 0; k<L; k++){
          int x1 = snp_ones(k,j);
          int x0 = snp_zeroes(k,j);
          int total = x1 + x0;
          C = (total + 2*lambda)*(total + 1 + 2*lambda);
          g = genos(k,i);
          if (g == 0) {p = (x0 + 1 + lambda)*(x0 + lambda);}
          else if (g == 1) {p = 2*(x0 + lambda)*(x1 + lambda);}
          else if (g == 2) {p = (x1 + 1 + lambda)*(x1 + lambda);}
          else {p = 1;}
          new_p = p/C;
          like *= new_p;
        }
        out(j,i) = like;
      }
   }
  return out;
}