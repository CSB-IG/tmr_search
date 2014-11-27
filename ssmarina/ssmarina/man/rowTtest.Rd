\name{rowTtest}
\alias{rowTtest}
\title{Student's t-test for rows}
\usage{
  rowTtest(x, y = NULL, mu = 0, alternative = "two.sided")
}
\arguments{
  \item{x}{Numerical matrix containing the test samples}

  \item{y}{Optional numerical matrix containing the
  reference samples. If ommited \code{x} will be tested
  against mean = \code{mu}}

  \item{mu}{Number indicating the alternative hypothesis
  when \code{y} is ommited}

  \item{alternative}{Character string indicating the tail
  for the test, either two.sided, greater or lower}
}
\value{
  List of Student-t-statistic (\code{statistic}) and
  p-values (\code{p.value})
}
\description{
  This function performs a Student's t-test on each row of
  a matrix
}

