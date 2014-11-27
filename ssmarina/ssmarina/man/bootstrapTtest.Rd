\name{bootstrapTtest}
\alias{bootstrapTtest}
\title{Bootstrapped signature by t-test}
\usage{
  bootstrapTtest(x, y, per = 100, seed = 1, verbose = T)
}
\arguments{
  \item{x}{Matrix containing the test dataset}

  \item{y}{Matrix containing the reference dataset}

  \item{per}{Integer indicating the number of permutations}

  \item{seed}{Integer indicating the seed for the
  permutations, 0 for disable it #'}

  \item{verbose}{Logical, whether progression messages
  should be printed in the terminal}
}
\value{
  Matrix of z-scores with genes in rows and permutations in
  columns
}
\description{
  This function generates a bootstrapped signature matrix
  by t-test
}
\seealso{
  \code{\link{marina}}
}

