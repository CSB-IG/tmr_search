\name{signatureDistance}
\alias{signatureDistance}
\title{Signature Distance}
\usage{
  signatureDistance(dset1, dset2 = NULL, nn = NULL,
    groups = NULL, scale. = T, two.tails = T, ws = 2)
}
\arguments{
  \item{dset1}{Dataset of any type in matrix format, with
  features in rows and samples in columns}

  \item{dset2}{Optional Dataset. If provided, distance
  between columns of dset and dset2 are computed and
  reported as rows and columns, respectively; if not,
  distance between all possible pairs of columns from dset
  are computed}

  \item{nn}{Optional size for the signature, default is
  either the full signature or 10% of it, depending or
  whether \code{ws}=0 or not)}

  \item{groups}{Optional vector indicating the group ID of
  the samples}

  \item{scale.}{Logical, whether the data should be scaled}

  \item{two.tails}{Logical, whether a two tails, instead of
  1 tail test should be performed}

  \item{ws}{Number indicating the exponent for the
  weighting the signatures, the default of 0 is uniform
  weighting, 1 is weighting by SD}
}
\value{
  Object of class \code{signatureDistance} as a matrix of
  normalized enrichment scores
}
\description{
  This function computes the similarity between columns of
  a data matrix
}

