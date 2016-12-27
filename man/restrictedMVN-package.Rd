\name{restrictedMVN-package}
\alias{restrictedMVN-package}
\alias{restrictedMVN}
\docType{package}
\title{ 
  Sampler from multivariate normal with affine constraints
}
\description{
  The package implements a fast gibbs sampler for the multivariate nomral with affine constraints. For the d-dimensional Z~Normal(mu,sigma), the linear_part matrix A in d x r, and offset vector b in 1 x r define a multivariate normal with affine constraints in \{Z| A*Z<= b\}.
}
\details{
  Sampling is implemented in the main function, \code{sample_from_constraints}. It is parameterized by the parameters of the normal (\code{mean_param} and \code{covariance}), parameters of the restriction (\code{linear_part} and \code{offset}), and the number of samples \code{ndraw}. The user also 
  needs to specify an initial point that satisfies the constraints. 
  \code{thresh2constraints} is a helper function that translates coordinate-wise truncations into the affine form.
}
\author{
Jonathan Taylor and Yuval Benjamini.

Maintainer: Yuval Benjamini <yuval.benjamini@mail.huji.ac.il>
}
\keyword{ package }
\seealso{
  The package was originally part of the github selective-inference code base.
}
\examples{
constr = thresh2constraints(3, lower = c(0.2,0.2,0.2))
covariance = matrix(c(1,0.5,0,0.5,1,0.5,0,0.5,1),nc=3)
 
samp = sample_from_constraints(linear_part = constr$linear_part,
                                    offset= constr$offset,
                                    mean_param = c(0,0,0),
                                    covariance = covariance,
                                    initial_point = c(1,1,1), 
                                    ndraw=20000,
                                    burnin=2000)


# all points should be >= 0.2 
stopifnot(all(samp>=0.2))

mean_restricted = colMeans(samp)

# compare to rejection of multivariate normals
library("MASS")
full_samp = mvrnorm(n=100000,mu = c(0,0,0),Sigma = covariance)
# Add restrictions:
pass_restrictions = apply(full_samp, 1, 
          function(x){all(constr$linear_part\%*\% x - constr$offset <=0 )})
          
cond_samp = full_samp[pass_restrictions,]
mean_restricted_rej = colMeans(cond_samp)

stopifnot(all(abs(mean_restricted - mean_restricted_rej)<=0.05))
}
