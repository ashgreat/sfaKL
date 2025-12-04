\name{sfaKL_simulate}
\alias{sfaKL_simulate}
\title{Simulate data for the sfaKL model}
\usage{
sfaKL_simulate(n = 1000, J = 2, M = 2, seed = 123, Sigma_mu = NULL, Sigma_delta = NULL, Omega = NULL)
}
\arguments{
\item{n}{Number of observations.}

\item{J}{Number of inputs.}

\item{M}{Number of outputs.}

\item{seed}{Random seed.}

\item{Sigma_mu}{Optional covariance matrix for input inefficiency.}

\item{Sigma_delta}{Optional covariance matrix for output inefficiency.}

\item{Omega}{Optional covariance matrix for noise across equations.}
}
\value{
List with simulated data and the true parameters used to generate it.
}
\description{
Generates share equations with input- and output-specific inefficiency, allowing correlated inefficiency and noise components through user-supplied or random covariance matrices.
}
