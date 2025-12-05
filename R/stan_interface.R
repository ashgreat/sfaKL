#' Stan Estimation for sfaKL Model
#'
#' This file provides functions to run the Bayesian estimation of the
#' Kumbhakar-Lai (2021) model using Stan via cmdstanr.
#'
#' @section Prerequisites:
#' \preformatted{
#' install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
#' cmdstanr::install_cmdstan()
#' }
#'
#' @name stan_interface
NULL

#' Compile the Stan model
#'
#' @return A CmdStanModel object
compile_sfa_kl_model <- function() {
    stan_file <- system.file("stan", "sfa_kl.stan", package = "sfaKL")
    if (stan_file == "") {
        # Fallback to local path for development
        stan_file <- file.path(dirname(dirname(getwd())), "sfaKL", "stan", "sfa_kl.stan")
        if (!file.exists(stan_file)) {
            stan_file <- "sfaKL/stan/sfa_kl.stan"
        }
    }
    
    if (!file.exists(stan_file)) {
        stop("Stan model file not found. Expected at: ", stan_file)
    }
    
    message("Compiling Stan model from: ", stan_file)
    cmdstan_model(stan_file)
}

#' Prepare data for Stan
#'
#' @param data Data frame with share and price ratio columns
#' @param J Number of inputs
#' @param M Number of outputs
#' @return A list suitable for Stan
prepare_stan_data <- function(data, J = 2, M = 2) {
    N <- nrow(data)
    n_eq <- (J - 1) + M
    
    # Extract matrices
    S <- as.matrix(data[, paste0("S", 2:J), drop = FALSE])
    R <- as.matrix(data[, paste0("R", 1:M), drop = FALSE])
    ln_w <- as.matrix(data[, paste0("ln_w", 2:J, "_w1"), drop = FALSE])
    ln_p <- as.matrix(data[, paste0("ln_p", 1:M, "_w1"), drop = FALSE])
    
    list(
        N = N,
        J = J,
        M = M,
        S = S,
        R = R,
        ln_w = ln_w,
        ln_p = ln_p
    )
}

#' Run MCMC sampling
#'
#' @param data Data frame with share and price ratio columns
#' @param J Number of inputs (default 2)
#' @param M Number of outputs (default 2)
#' @param chains Number of MCMC chains (default 4)
#' @param iter_warmup Warmup iterations per chain (default 1000)
#' @param iter_sampling Sampling iterations per chain (default 1000)
#' @param ... Additional arguments passed to $sample()
#' @return A CmdStanMCMC object
sfaKL_stan <- function(data, J = 2, M = 2,
                       chains = 4,
                       iter_warmup = 1000,
                       iter_sampling = 1000,
                       ...) {
    
    # Compile model
    model <- compile_sfa_kl_model()
    
    # Prepare data
    stan_data <- prepare_stan_data(data, J, M)
    
    message("Running MCMC with ", chains, " chains...")
    
    # Sample
    fit <- model$sample(
        data = stan_data,
        chains = chains,
        parallel_chains = min(chains, parallel::detectCores()),
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        refresh = 500,
        ...
    )
    
    return(fit)
}

#' Extract efficiency estimates from Stan fit
#'
#' @param fit CmdStanMCMC object from sfaKL_stan()
#' @return A list with TE_input and TE_output matrices (posterior means)
extract_efficiencies <- function(fit) {
    draws <- fit$draws(format = "draws_matrix")
    
    # Get dimensions from variable names
    te_input_vars <- grep("^TE_input\\[", colnames(draws), value = TRUE)
    te_output_vars <- grep("^TE_output\\[", colnames(draws), value = TRUE)
    
    # Calculate posterior means
    te_input_means <- colMeans(draws[, te_input_vars, drop = FALSE])
    te_output_means <- colMeans(draws[, te_output_vars, drop = FALSE])
    
    # Parse dimensions
    N <- max(as.integer(gsub(".*\\[(\\d+),.*", "\\1", te_input_vars)))
    J <- max(as.integer(gsub(".*,(\\d+)\\]", "\\1", te_input_vars)))
    M <- max(as.integer(gsub(".*,(\\d+)\\]", "\\1", te_output_vars)))
    
    # Reshape to matrices
    TE_input <- matrix(te_input_means, nrow = N, ncol = J, byrow = TRUE)
    TE_output <- matrix(te_output_means, nrow = N, ncol = M, byrow = TRUE)
    
    list(
        TE_input = TE_input,
        TE_output = TE_output,
        TE_input_full = draws[, te_input_vars],   # Full posterior
        TE_output_full = draws[, te_output_vars]
    )
}

#' Summary of parameter estimates
#'
#' @param fit CmdStanMCMC object
#' @return Summary data frame
summarize_params <- function(fit) {
    fit$summary(variables = c("alpha", "beta", "L_Sigma_mu", "L_Sigma_delta", "L_Omega"))
}

# ============================================================================
# Example Usage
# ============================================================================

if (FALSE) {  # Set to TRUE to run example
    # Load package functions for simulation
    source("sfaKL/R/simulate_data.R")
    source("sfaKL/R/utils.R")
    
    # Simulate data
    set.seed(123)
    Sigma_mu <- matrix(c(0.04, 0.02, 0.02, 0.05), 2, 2)
    Sigma_delta <- matrix(c(0.03, 0.01, 0.01, 0.04), 2, 2)
    Omega <- diag(c(0.04, 0.05, 0.06))
    
    sim <- sfaKL_simulate(
        n = 200, J = 2, M = 2, seed = 123,
        Sigma_mu = Sigma_mu,
        Sigma_delta = Sigma_delta,
        Omega = Omega
    )
    
    # Run Stan
    fit <- sfaKL_stan(sim$data, J = 2, M = 2,
                      chains = 4,
                      iter_warmup = 500,
                      iter_sampling = 500)
    
    # Check diagnostics
    fit$cmdstan_diagnose()
    
    # Parameter summary
    print(summarize_params(fit))
    
    # Extract efficiencies
    eff <- extract_efficiencies(fit)
    cat("Mean input efficiency (input 1):", mean(eff$TE_input[, 1]), "\n")
    cat("Mean output efficiency (output 1):", mean(eff$TE_output[, 1]), "\n")
}
