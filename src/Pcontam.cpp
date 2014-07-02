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
//' @return Returns a P*P x N numeric matrix which contains the the probabilities that contaminated individuals 
//' originated from certain combintations of two populations.  The probability that the ith individual's contamination 
//' originated from population p1 and p2 is ouput(N*(p1-1) + p2,i) + output(N*(p2-1) + p1,i) if p1 and p2
//' are different populations and output(N*(p1-1) + p1,i) if p1 and p2 are the same population.    
//' @export
// [[Rcpp::export]]
NumericMatrix Pcontam(IntegerMatrix snp_zeroes, IntegerMatrix snp_ones, IntegerMatrix genos, double lambda) {
  int N = genos.ncol();
  int L = genos.nrow();
  int P = snp_zeroes.ncol();
  double C;
  double p;
  double like;
  double new_p;
  int g;
  NumericMatrix out(P*P,N);
  for(int i=0; i<N; i++){
     int j=0;
     for(j=0; j<P; j++) {
       int l = 0;
       for(l=0; l<P; l++){
         int k= 0;
         like = 1;
          for(k= 0; k<L; k++){
            int x1 = snp_ones(k,j);
            int x0 = snp_zeroes(k,j);
            int y1 = snp_ones(k,l);
            int y0 = snp_zeroes(k,l);
            int totalx = x1 + x0;
            int totaly = y1 + y0;
            g = genos(k,i);
            if(j!=l){
              C = (totalx + 2*lambda)*(totalx + 1 + 2*lambda)*(totaly + 2*lambda)*(totaly + 1 + 2*lambda);
              if (g == 0) {p = (x0 + 1 + lambda)*(x0 + lambda)*(y0 + 1 + lambda)*(y0 + lambda)/C;}
              else if (g == 1) {p = 1 - (x0 + 1 + lambda)*(x0 + lambda)*(y0 + 1 + lambda)*(y0 + lambda)/C - (x1 + 1 + lambda)*(x1 + lambda)*(y1 + lambda)*(y1 + 1 + lambda)/C;}
              else if (g == 2) {p = (x1 + 1 + lambda)*(x1 + lambda)*(y1 + lambda)*(y1 + 1 + lambda)/C;}
              else {p = 1;}
              }else{
              C = (totalx + 2*lambda)*(totalx + 1 + 2*lambda)*(totalx + 2 + 2*lambda)*(totalx + 3 + 2*lambda);
              if (g == 0) {p = (x0 + 1 + lambda)*(x0 + lambda)*(x0 + 2 + lambda)*(y0 + 3 + lambda)/C;}
              else if (g == 1) {p = 1 - (x0 + 1 + lambda)*(x0 + lambda)*(x0 + 2 + lambda)*(x0 + 3 + lambda)/C - (x1 + 1 + lambda)*(x1 + lambda)*(x1 + 2 + lambda)*(y1 + 3 + lambda)/C;}
              else if (g == 2) {p = (x1 + 1 + lambda)*(x1 + lambda)*(x1 + 2 + lambda)*(x1 + 3 + lambda)/C;}
              else {p = 1;}
              }  
            like *= p;
             }
          out(j*P + l,i) = like;
        }
     }
   }
  return out;
}