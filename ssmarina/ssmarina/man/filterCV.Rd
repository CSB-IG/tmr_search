\name{filterCV}
\alias{filterCV}
\title{Coefficient of variation filter}
\usage{
  filterCV(expset)
}
\arguments{
  \item{expset}{Matrix containing the gene expression data,
  with samples in columns and probes in rows. The
  \code{colnames} attribute should contain the sample names
  and the \code{rownames} attribute should contain the
  unique geneIDs}
}
\value{
  CV filtered dataset
}
\description{
  This function filter redundant probes based on the
  highest coefficient of variation
}

