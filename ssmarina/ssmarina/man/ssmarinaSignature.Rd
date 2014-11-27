\name{ssmarinaSignature}
\alias{ssmarinaSignature}
\title{Signature and sample-permutation null model for ssmarina}
\usage{
  ssmarinaSignature(eset, ref, per = 1000, seed = 1,
    verbose = T)
}
\arguments{
  \item{eset}{Numeric matrix containing the test dataset,
  with genes in rows and samples in columns}

  \item{ref}{Numeric matrix containing the reference
  samples (columns) and genes in rows}

  \item{per}{Integer indicating the number of sample
  permutations}

  \item{seed}{Integer indicating the seed for the random
  sample generation. The system default is used when set to
  zero}

  \item{verbose}{Logical, whether progression messages
  should be printed in the terminal}
}
\value{
  ssmarinaSignature S3 object containing the signature and
  null model
}
\description{
  This function generates a ssmarinaSignature object from a
  test dataset based on a set of samples to use as
  reference
}

