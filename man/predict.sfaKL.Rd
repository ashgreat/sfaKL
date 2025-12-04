\name{predict.sfaKL}
\alias{predict.sfaKL}
\title{Predict inefficiency components}
\usage{
\method{predict}{sfaKL}(object, newdata = NULL, n_cores = 1, ...)
}
\arguments{
\item{object}{An object of class \code{sfaKL} produced by \code{sfaKL_estimate}.}

\item{newdata}{Optional new data frame; original column names are mapped automatically.}

\item{n_cores}{Number of cores for the conditional expectation calculations.}

\item{...}{Additional arguments (unused).}
}
\value{
Data frame with predicted \code{mu} (inputs) and \code{delta} (outputs) components.
}
\description{
Computes conditional expectations of input and output inefficiency using the estimated closed skew-normal structure, reusing column mappings saved during estimation.
}
