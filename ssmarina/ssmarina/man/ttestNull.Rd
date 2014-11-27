\name{ttestNull}
\alias{ttestNull}
\title{Null model by sample permutation testing}
\usage{
  ttestNull(x, y, per = 1000, repos = T, seed = 1,
    verbose = T)
}
\arguments{
  \item{x}{Matrix containing the test dataset}

  \item{y}{Matrix containing the reference dataset}

  \item{per}{Integer indicating the number of permutations}

  \item{repos}{Logical, whether the permutations should be
  performed with reposition}

  \item{seed}{Integer indicating the seed for the
  permutations, 0 for disable it}

  \item{verbose}{Logical, whether progression messages
  should be printed in the terminal}
}
\value{
  Matrix of z-scores with genes in rows and permutations in
  columns
}
\description{
  This function perfomrs sample permutation and t-test to
  generate a null model
}
\seealso{
  \code{\link{marina}}, \code{\link{ssmarina}}
}

