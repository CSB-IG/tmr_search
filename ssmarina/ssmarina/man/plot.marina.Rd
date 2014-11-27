\name{plot.marina}
\alias{plot.marina}
\title{Plot MARINA results}
\usage{
  \method{plot}{marina} (x, mrs = 10,
    color = c("cornflowerblue", "salmon"), pval = NULL,
    bins = 500, cex = 0, density = 0, smooth = 0,
    sep = 0.2, hybrid = T,
    include = c("expression", "activity"), gama = 2, ...)
}
\arguments{
  \item{x}{MARINa object produced by \code{marina}
  function}

  \item{mrs}{Either an integer indicating the number of
  master regulators to include in the plot, or a character
  vector containing the names of the master regulators to
  include in the plot}

  \item{color}{Vector of two components indicating the
  colors for the negative and positive parts of the
  regulon}

  \item{pval}{Optional matrix of p-values to include in the
  plot}

  \item{bins}{Number of bins to split the vector of scores
  in order to compute the density color of the bars}

  \item{cex}{Number indicating the text size scaling, 0
  indicates automatic scaling}

  \item{density}{Integrer indicating the number of steps
  for the kernel density. Zero for not ploting it}

  \item{smooth}{Number indicating the proportion of point
  for smoothing the density distribution. Zero for not
  using the smoother}

  \item{sep}{Number indicating the separation from figure
  and text}

  \item{hybrid}{Logical, whether the 3-tail approach used
  for computingthe enrichment should be reflected in the
  plot}

  \item{include}{Vector indicating the information to
  include as heatmap to the right of the marina plot:
  expression and activity}

  \item{gama}{Positive number indicating the exponential
  transformation for the activity and expression color
  scale}

  \item{...}{Given for compatibility to the plot generic
  function}
}
\value{
  Nothing, a plot is generated in the default output device
}
\description{
  This function generate a plot for MARINA results showing
  the enrichment of the target genes for each significant
  master regulator on the gene expression signature
}
\seealso{
  \code{\link{marina}}
}

