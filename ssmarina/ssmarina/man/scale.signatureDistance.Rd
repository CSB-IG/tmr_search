\name{scale.signatureDistance}
\alias{scale.signatureDistance}
\title{Scaling of signatureDistance objects}
\usage{
  \method{scale}{signatureDistance} (x, center = TRUE,
    scale = TRUE)
}
\arguments{
  \item{x}{signatureDistance object}

  \item{center}{Not used, given for compatibility with the
  generic funtion scale}

  \item{scale}{Not used, given for compatibility with the
  generic function scale}
}
\value{
  Scaled signatureDistance object
}
\description{
  This function scales the signatureDistance so its range
  is (-1, 1)
}

