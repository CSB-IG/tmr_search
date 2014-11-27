\name{plot.compareRegulon}
\alias{plot.compareRegulon}
\title{Plot the regulon comparison analysis}
\usage{
  \method{plot}{compareRegulon} (x, type = "raw",
    annot = c("Ntw1", "Ntw2"), pval = 0.05,
    color = c("blue", "red"), legend = "topleft", ...)
}
\arguments{
  \item{x}{compareRegulon object}

  \item{type}{Character string indicating the type of plot,
  either 'raw' or 'zscore'}

  \item{annot}{Vector of 2 character strings indicating the
  name of the networks to compare}

  \item{pval}{Numberr indicating the p-value threshold for
  Zscore type}

  \item{color}{Vecor of two character strings indicating
  the color for the positive controls, no positive ontrol
  is ploted when color is set to NULL}

  \item{legend}{Character string indicating the position of
  the legend, which is ommited when legend is set to none}

  \item{...}{Additional parameters to pass to the generic
  plot function}
}
\value{
  Nothing, a plot is generated in the default output device
}
\description{
  This function plots the regulon comparison analysis
  results
}

