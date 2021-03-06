% Generated by roxygen2 (4.0.1): do not edit by hand
\name{mixed_MCMC}
\alias{mixed_MCMC}
\title{MCMC function for determing contamination proportions and mixing proportions of a mixture sample}
\usage{
mixed_MCMC(data, contam_data, clean_data, alpha = 0.5, beta = 0.5,
  contamination = TRUE, inters)
}
\arguments{
\item{data}{A L x N matrix containing the genotype of data of individuals in the form of 0s,1s, and 2s.
N is the number of individuals, and L is the number of loci.}

\item{contam_data}{A P*P x N matrix containing the likelihood that each individual orginates from each
combination of the 2 populations from P, assuming the individual is contaminated.  P is the number of
populations and is the number of individuals.}

\item{clean_data}{A P x N matrix containing the likelihood that each individual originates from each
of P populations, assuming the individual is not contaminated.}

\item{alpha}{The alpha parameter for the prior distribution of rho, the probability
of contamination.  The prior of rho is assumed to be a beta distribution.}

\item{beta}{The beta parameter for the prior distribution of rho, which is assumed to be
a beta distribution.}

\item{inters}{Total number of sweeps used in the MCMC model.}

\item{contamination}{Logical variable.  If contamination = FALSE, then the MCMC model will only update
the mixing proportions and the population identity, assuming no contamination is present.}
}
\value{
Returns a list of four named components:
\describe{
 \item{prob_contam}{A vector containing the rho value, which is the proportion of
 contaminated samples for each interation. The vector is a length of 1 plus the total number of iterations.}
 \item{pops}{A matrix the u, which denote the population or populations of origin, for each sweep.
 If contamintion=TRUE, u is a 2*inters x N matrix, and each pair of rows represents the population identification for each sweep.
 If an individual has u values of (any P, 0), than it is a non contaminated individual and the first value is
 its population of origin.  If an individual has two u values, it is contaminated and the two values
 represent the source of contamination.  If contamination=FALSE, u is a intersxN matrix, and each row represents the popultion
 identification for each sweep.}
 \item{z}{A matrix containing the z value, which denotes contaminated status, for each individual
 and every iteration.  The matrix has columns equal to the number of individuals and rows equal to
 the total number of iterations.}
 \item{mixing}{A matrix containing the mixing proportions of each population for each sweep.  The matrix
 has number of rows equal to the total interations plus one and columns equal to the number of populations.}
}
}
\description{
This function runs an MCMC tailored for a situation in which one is sampling from a mixture of individuals
from different populaitons.  The MCMC function uses the SNP genotype data and baseline data to estimate the
contamination proportion and the mixing proportions as well as identify contaminated individuals and
preform population assignment.
}
\examples{
# Creates data for the mixed_MCMC function
zeros <- t(matrix(c(40,160,50,80,100,20,30,15,140),nrow=3))
ones <- 200 - zeros
genos <- t(matrix(c(1,2,0,2,1,0,1,0,0,1,2,2,1,1,1),nrow=5))
data <- matrix(c(1,2,1,0,2,0,2,0,0,0,1,2,0,0,1),nrow=3)
clean_data <- P_likelihood(zeros,ones,genos,.5)
contam_data <- Pcontam(zeros,ones,genos,.5)
# Runs the MCMC that checks for contamination on the data for 5 interations
mixed_MCMC(data, contam_data, clean_data, contamination = TRUE, inters = 5)

# Runs the MCMC on same data, but without evaluating contamination
mixed_MCMC(data, contam_data, clean_data, contamination=FALSE, inters = 5)
}

