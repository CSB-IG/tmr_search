\name{scaleGroups}
\alias{scaleGroups}
\title{Signatures with grouping variable}
\usage{
  scaleGroups(x, groups)
}
\arguments{
  \item{x}{Numerical matrix with genes in rows and samples
  in columns}

  \item{groups}{Vector of same length as columns has the
  dset containing the labels for grouping the samples}
}
\value{
  Numeric matrix of signatures (z-scores) with genes in
  rows and groups in columns
}
\description{
  \code{scaleGroups} compares each group vs. the remaining
  groups using a Student's t-test
}
\details{
  This function compute signatures using groups information
}

