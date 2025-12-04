\name{sfaKL_estimate}
\alias{sfaKL_estimate}
\title{Estimate the sfaKL model}
\usage{
sfaKL_estimate(
  data,
  share_input = "S2",
  share_output = c("R1", "R2"),
  price_input_ratio = "ln_w2_w1",
  price_output_ratios = c("ln_p1_w1", "ln_p2_w1"),
  start_params = NULL,
  method = "L-BFGS-B",
  control = list(maxit = 10000, reltol = 1e-8),
  n_cores = 1,
  verbose = FALSE,
  restarts = 1,
  use_gradient = TRUE,
  cov_restart_scale = 0.15,
  lower = NULL,
  upper = NULL,
  parscale = NULL
)
}
\arguments{
\item{data}{Data frame containing share variables and price ratios.}

\item{share_input}{Name(s) of input share columns (default \code{"S2"}).}

\item{share_output}{Names of output share columns (default \code{c("R1", "R2")}).}

\item{price_input_ratio}{Name(s) of the log input price ratio columns.}

\item{price_output_ratios}{Names of the log output-to-input price ratio columns.}

\item{start_params}{Optional named vector of starting parameters; defaults to OLS-derived values.}

\item{method}{Optimization method passed to \code{optim}.}

\item{control}{Control list passed to \code{optim}.}

\item{n_cores}{Number of cores for likelihood and prediction evaluations.}

\item{verbose}{Logical; if TRUE prints optimizer diagnostics.}

\item{restarts}{Number of random restarts; best solution is returned.}

\item{use_gradient}{Logical; if TRUE uses finite-difference gradients for gradient-capable optimizers.}

\item{cov_restart_scale}{Jitter scale (sd) applied to covariance parameters during restarts.}

\item{lower}{Optional lower bounds for box-constrained optimization.}

\item{upper}{Optional upper bounds for box-constrained optimization.}

\item{parscale}{Optional parameter scaling vector passed to \code{optim}.}
}
\value{
An object of class \code{sfaKL} containing parameter estimates, likelihood, and mapped data.
}
\description{
Fits the Kumbhakar--Lai stochastic frontier system with homogeneity and adding-up restrictions, estimating full covariance structures for inefficiency and noise.
}
