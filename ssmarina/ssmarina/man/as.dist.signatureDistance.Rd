\name{as.dist.signatureDistance}
\alias{as.dist.signatureDistance}
\title{Distance matrix from signatureDistance objects}
\usage{
  \method{as.dist}{signatureDistance} (m, diag = F,
    upper = F)
}
\arguments{
  \item{m}{signatureDistance object}

  \item{diag}{parameter included for compatibility}

  \item{upper}{parameter included for compatibility}
}
\value{
  Object of class dist
}
\description{
  This function transforms a signatureDistance object into
  a dist object
}

