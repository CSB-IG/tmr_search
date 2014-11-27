\name{groupPwea3}
\alias{groupPwea3}
\title{Proportionally Weighted Enrichment Analysis for gene-set groups}
\usage{
  groupPwea3(rlist, groups, nullpw = NULL,
    alternative = c("two.sided", "less", "greater")[1],
    per = 0, tnorm = T, minsize = 5, tw = 1, lw = 1,
    simsig = T, verbose = T)
}
\arguments{
  \item{rlist}{Named vector containing the scores to rank
  the expression profile or matrix where columns contains
  bootstraped signatures}

  \item{groups}{List of gene-sets (regulons), each
  component is a list of two vectors: \emph{TFmode}
  containing the TFMoA index (-1; 1) and \emph{likelihood}
  containing the interaction relative likelihood}

  \item{nullpw}{Numerical matrix representing the null
  model, with genes as rows (geneID as rownames) and
  permutations as columns}

  \item{alternative}{Character string indicating the
  alternative hypothesis, either two.sided, greater or
  less}

  \item{per}{Integer indicating the number of permutations
  for the genes in case ???nullpw??? is ommited}

  \item{tnorm}{Logical, whether the ranking vector should
  be normally transformed after the copula transformation}

  \item{minsize}{Integer indicating the minimum size for
  the regulons}

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
  A list containing four matrices: \describe{
  \item{es}{Enrichment score} \item{nes}{Normalized
  Enrichment Score} \item{size}{Regulon size}
  \item{p.value}{Enrichment p.value} }
}
\description{
  This function performs a Proportionally Weighted
  Enrichment Analysis on groups of gene-sets
}

