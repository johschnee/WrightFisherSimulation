% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{graphic_variance}
\alias{graphic_variance}
\title{Graphic variance}
\usage{
graphic_variance(pop)
}
\arguments{
\item{pop}{Data from the function \code{simulation()}}
}
\value{
The described graphic
}
\description{
Creates a graphic that shows the variance of the number of alleles over generations.
}
\examples{
p <- simulation()
h <- graphic_variance(pop = p)
}
