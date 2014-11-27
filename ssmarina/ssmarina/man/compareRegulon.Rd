\name{compareRegulon}
\alias{compareRegulon}
\title{Compare regulons}
\usage{
  compareRegulon(regulon1, regulon2, dset = NULL,
    subreg1 = NULL, subreg2 = NULL, method = "pearson",
    groups = NULL, nn = NULL, minsize = 25, ws = 2,
    rep = 1)
}
\arguments{
  \item{regulon1}{Regulon to be compare agains the
  reference regulon \code{regulon2}}

  \item{regulon2}{Regulon to use as reference to compare
  \code{regulon1}}

  \item{dset}{Numeric matrix or list of numeric matrixes
  containing the datasets, with genes in rows and samples
  in columns}

  \item{subreg1}{Single regulon or list of two regulons to
  use as positive conservation control for \code{regulon1},
  usually networks obtained from context-matching
  independent datasets. If a single regulon is given, then
  \code{regulon1} will be used as second regulon}

  \item{subreg2}{Single regulon or list of two regulons to
  use as positive conservation control for \code{regulon2},
  usually networks obtained from context-matching
  independent datasets. If a single regulon is given, then
  \code{regulon2} will be used as second regulon}

  \item{method}{Character string indicating the similarity
  metric, either signatureDistance, pearson, spearman,
  kendall, euclidean, jaccard, oddsratio or fet}

  \item{groups}{Vector of factors indicating the groups of
  datasets used during the integration of results}

  \item{nn}{Integer indicating the number of genes used to
  compute the similarity between regulator's signatures}

  \item{minsize}{Integer indicating the minimum allowed
  size for the regulons}

  \item{ws}{Number indicating the weight score for
  \code{signatureDistance}}

  \item{rep}{Integer indicating the number of iterations
  for the negative controls}
}
\value{
  List of six components containing the normalized
  enrichment scores for the comparison between regulons,
  for the internal comparion of regulons from network 1,
  and for the internal comparison of regulons from network
  2, the maximum difference in the normalized enrichment
  score and the associated p-value, a similarity matrix
  with regulon1 in rows and regulon2 in columns
}
\description{
  Thins funciton compare regulatory networks based on
  regulon activity
}

