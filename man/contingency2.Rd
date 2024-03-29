% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epi_functions.R
\name{contingency2}
\alias{contingency2}
\title{Measures of association from two by two contingency tables (direct input).}
\usage{
contingency2(aa, bb, cc, dd, ...)
}
\arguments{
\item{aa}{Number of cases where both exposure and outcome are present.}

\item{bb}{Number of cases where exposure is present but outcome is absent.}

\item{cc}{Number of cases where exposure is absent but outcome is present.}

\item{dd}{Number of cases where both exposure and outcome are absent.}

\item{...}{Further arguments passed to \code{\link[epiR]{epi.2by2}}.}
}
\description{
\code{contingency2} is a wrap that calls \code{epi.2by2} from package \code{epiR}.
}
\examples{
## A case-control study on the effect of alcohol on oesophageal cancer.
Freq <- c(386, 29, 389, 171)
status <- gl(2, 1, 4, labels = c("Control", "Case"))
alcohol <- gl(2, 2, labels = c("0-39", "40+"))
cancer <- data.frame(Freq, status, alcohol)
cancer <- expand_df(cancer)

contingency2(171, 389, 29, 386, method = "case.control")
}
\seealso{
\code{\link[epiR]{epi.2by2}}.
}
