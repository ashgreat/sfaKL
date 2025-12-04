\name{sfaKL_check}
\alias{sfaKL_check}
\title{Diagnostic check for sfaKL fits}
\usage{
sfaKL_check(object)
}
\arguments{
\item{object}{An object of class \code{sfaKL}.}
}
\value{
List with positive-definiteness flags and condition numbers for key covariance matrices.
}
\description{
Validates that the estimated covariance matrices are positive definite and reports condition numbers, helping diagnose convergence or identifiability issues.
}
