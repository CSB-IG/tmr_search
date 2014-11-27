\name{pwea3NULLgroups}
\alias{pwea3NULLgroups}
\title{Regulon-specific NULL model}
\usage{
  pwea3NULLgroups(pwnull, groups, tnorm = T, tw = 1,
    lw = 1, simsig = T, verbose = T)
}
\arguments{
  \item{pwnull}{Numerical matrix representing the null
  model, with genes as rows (geneID as rownames) and
  permutations as columns}

  \item{groups}{List containing the regulons}

  \item{tnorm}{Logical, whether the ranking vector should
  be normally transformed after the copula transformation}

  \item{tw}{Number indicating the exponent for the target
  pleotropy weight}

  \item{lw}{Number indicating the exponent for the edge
  likelihood weight}

  \item{simsig}{Logical, whether a symetrical distribution
  of the signature is used}

  \item{verbose}{Logical, whether progression messages
  should be printed in the terminal}
}
\value{
  A list containing two elements: \describe{
  \item{groups}{Regulon-specific NULL model containing the
  enrichment scores} \item{ss}{Direction of the
  regulon-specific NULL model} }
}
\description{
  This function generates the regulon-specific NULL models
}

