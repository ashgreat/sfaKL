\name{plot_efficiencies}
\alias{plot_efficiencies}
\title{Plot efficiency scores}
\usage{
plot_efficiencies(object, type = "hist")
}
\arguments{
\item{object}{An object of class \code{sfaKL}.}

\item{type}{Plot type: \code{"hist"} or \code{"density"}.}
}
\value{
A \code{ggplot} object showing efficiency distributions.
}
\description{
Creates histograms or densities for input and output efficiencies, converting predicted \code{mu} and \code{delta} values to the efficiency scale.
}
