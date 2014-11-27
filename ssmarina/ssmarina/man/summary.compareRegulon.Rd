\name{summary.compareRegulon}
\alias{summary.compareRegulon}
\title{Table of results for regulon comparison analysis}
\usage{
  \method{summary}{compareRegulon} (object, pval = 1,
    annot = c("Ntw1", "Ntw2"), ...)
}
\arguments{
  \item{object}{compareRegulon object}

  \item{pval}{Number indicating the p-value cutoff for
  regulon reliability}

  \item{annot}{Vector of 2 character strings indicating the
  names of the networks}

  \item{...}{Given for compatibility with the summary
  generic function}
}
\value{
  Data.frame of results
}
\description{
  This function generates a table of results for the
  regulon comparison analysis
}

