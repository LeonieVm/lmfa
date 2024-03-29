% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0PlotTransitions.r
\name{plot.lmfa_step3}
\alias{plot.lmfa_step3}
\title{Transition plot for a specific subject}
\usage{
\method{plot}{lmfa_step3}(x, identifier, id, ...)
}
\arguments{
\item{x}{Output from the main function step3() (must be of class lmfa_step3).}

\item{identifier}{The name of the column containing the subject identifiers (must be a single character).}

\item{id}{Identification number of the subject for which the transition pattern should be plotted.}

\item{...}{Further arguments for the default S3 plot method.}
}
\description{
Transition plot for a specific subject
}
\examples{
\dontrun{
plot(transitionmodel)
}
}
