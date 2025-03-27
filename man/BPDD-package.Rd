\name{BPDD-package}
\alias{BPDD-package}
\docType{package}
\title{A Bayesian Phase I/II Platform Design with Data Augmentation Accounting for Delayed Outcomes}
\description{
  This package implements Bayesian data augmentation methods for handling
  delayed outcomes in Phase I/II clinical trials. The code is specifically designed to facilitate
  platform designs, where data augmentation is employed to improve trial efficiency by accounting for
  delayed responses. The package includes functions for parallel simulations, model fitting, and
  statistical analysis, providing an easy-to-use solution for researchers to conduct Bayesian analysis in
  clinical trial scenarios. With integrated code, users can easily perform parallel computations for
  simulations after loading the package.
}
\author{
  Wentao Yang \email{jasper.yang@sjtu.edu.cn}
}
\seealso{
  \code{\link{parallel_simulation}}, \code{\link{main_BPDD}}
}
\examples{
  # Example usage of the BPDD package
  result <- parallel_simulation(num_simulations = 4, num_cores = 4)
  print(result)
}
