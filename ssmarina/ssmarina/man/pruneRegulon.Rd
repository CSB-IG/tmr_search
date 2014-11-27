\name{pruneRegulon}
\alias{pruneRegulon}
\title{Prune Regulons}
\usage{
  pruneRegulon(regulon, cutoff = 50, eliminate = F)
}
\arguments{
  \item{regulon}{Object of class regulon}

  \item{cutoff}{Number indicating the maximum size for the
  regulons (maximum number of target genes)}

  \item{eliminate}{Logical whether regulons smalles than
  \code{cutoff} should be eliminated}
}
\value{
  Prunned regulon
}
\description{
  This function limits the maximum size of the regulons
}
\seealso{
  \code{\link{ssmarina}}, \code{\link{marina}}
}

