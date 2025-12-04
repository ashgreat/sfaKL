#' Log-Likelihood for sfaKL model
#'
#' @param params Named vector of parameters
#' @param data Data frame containing the variables
#' @param J Number of inputs
#' @param M Number of outputs
#' @param n_cores Number of cores for parallel execution
#' @return Negative log-likelihood
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats dnorm
#' @importFrom parallel mclapply
sfaKL_loglik <- function(params, data, J = 2, M = 2, n_cores = 1) {
    # Dimensions
    n_input_shares <- J - 1
    n_output_shares <- M
    n_eq <- n_input_shares + n_output_shares
    n_ineff <- J + M

    # Extract Data
    # Assuming data has standard names from sfaKL_estimate
    S <- as.matrix(data[, paste0("S", 2:J), drop = FALSE])
    R <- as.matrix(data[, paste0("R", 1:M), drop = FALSE])
    ln_w <- as.matrix(data[, paste0("ln_w", 2:J, "_w1"), drop = FALSE])
    ln_p <- as.matrix(data[, paste0("ln_p", 1:M, "_w1"), drop = FALSE])

    # Parameters with homogeneity/add-up enforcement
    profit_params <- assemble_profit_parameters(params, J, M)
    alpha_0 <- profit_params$alpha
    beta_0 <- profit_params$beta
    A <- profit_params$A
    B <- profit_params$B
    Gamma <- profit_params$Gamma

    # Covariance matrices (positive definite via Cholesky)
    # Construct Residuals
    eps_S <- matrix(0, nrow(data), n_input_shares)
    for (j in 1:n_input_shares) {
        det_part <- alpha_0[j] + as.vector(A[j, , drop = FALSE] %*% t(ln_w)) + as.vector(Gamma[j, , drop = FALSE] %*% t(ln_p))
        eps_S[, j] <- -S[, j] - det_part
    }

    eps_R <- matrix(0, nrow(data), n_output_shares)
    for (m in 1:n_output_shares) {
        det_part <- beta_0[m] + as.vector(B[m, , drop = FALSE] %*% t(ln_p)) + as.vector(t(Gamma[, m, drop = FALSE]) %*% t(ln_w))
        eps_R[, m] <- R[, m] - det_part
    }

    epsilon <- cbind(eps_S, eps_R)

    # Robust Matrix Inversion
    tryCatch(
        {
            mats <- get_matrices(params, J, M)
            Theta <- mats$Theta
            Delta <- mats$Delta
            W <- mats$Psi

            # Term 1: PDF of N(0, Theta)
            term1 <- mvtnorm::dmvnorm(epsilon, mean = rep(0, n_eq), sigma = Theta, log = TRUE)

            # Term 2: CDF term
            args <- epsilon %*% t(W) # N x (J+M)

            if (any(is.na(args))) {
                return(1e10)
            }

            calc_term2 <- function(i) {
                prob <- suppressWarnings(mvtnorm::pmvnorm(lower = rep(-Inf, n_ineff), upper = args[i, ], mean = rep(0, n_ineff), sigma = Delta))
                val <- as.numeric(prob)
                if (val <= 0 || is.na(val)) val <- 1e-100
                return(log(val))
            }

            term2_list <- safe_parallel_lapply(1:nrow(data), calc_term2, n_cores = n_cores)
            term2 <- unlist(term2_list)

            # Term 3: Denominator
            term3 <- n_ineff * log(0.5)

            ll_i <- term1 + term2 - term3

            if (any(is.na(ll_i)) || any(is.infinite(ll_i))) {
                return(1e10)
            }

            return(-sum(ll_i))
        },
        error = function(e) {
            return(1e10)
        }
    )
}
