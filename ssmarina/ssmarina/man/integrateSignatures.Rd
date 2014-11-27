\name{integrateSignatures}
\alias{integrateSignatures}
\title{Integrate signatures}
\usage{
  integrateSignatures(signature, score = 1)
}
\arguments{
  \item{signature}{Numeric matrix containing the signatures
  as z-scores or NES, genes in rows and signatures in
  columns}

  \item{score}{Number indicating the exponent score for the
  weight}
}
\value{
  Vector containing the integrated signatures
}
\description{
  This function integrates signatures represented as
  columns in the input matrix using self-weighting average
}

