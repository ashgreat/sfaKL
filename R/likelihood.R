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

    # Parse Parameters
    # Alpha
    alpha_0 <- params[paste0("alpha", 2:J)]

    # A Matrix (Symmetric)
    A <- matrix(0, n_input_shares, n_input_shares)
    for (j in 1:n_input_shares) {
        for (k in j:n_input_shares) {
            val <- params[paste0("A_", j, "_", k)]
            A[j, k] <- val
            A[k, j] <- val
        }
    }

    # Beta
    beta_0 <- params[paste0("beta", 1:M)]

    # B Matrix (Symmetric)
    B <- matrix(0, n_output_shares, n_output_shares)
    for (m in 1:n_output_shares) {
        for (k in m:n_output_shares) {
            val <- params[paste0("B_", m, "_", k)]
            B[m, k] <- val
            B[k, m] <- val
        }
    }

    # Gamma Matrix
    Gamma <- matrix(0, n_input_shares, n_output_shares)
    for (j in 1:n_input_shares) {
        for (m in 1:n_output_shares) {
            Gamma[j, m] <- params[paste0("Gamma_", j, "_", m)]
        }
    }

    # Variances
    # Limit extreme values
    sigma_mu <- exp(pmin(params[paste0("log_sigma_mu", 1:J)], 10))
    sigma_delta <- exp(pmin(params[paste0("log_sigma_delta", 1:M)], 10))

    sigma_v_S <- exp(pmin(params[paste0("log_sigma_v_S", 2:J)], 10))
    sigma_v_R <- exp(pmin(params[paste0("log_sigma_v_R", 1:M)], 10))

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

    # Matrix Construction for Likelihood
    l_J_1 <- matrix(1, n_input_shares, 1)
    l_M <- matrix(1, n_output_shares, 1)
    I_J_1 <- diag(1, n_input_shares)
    O_M_J_1 <- matrix(0, n_output_shares, n_input_shares)

    # H1 = -A(-l, I) + Gamma(l, 0)
    # -l is J-1 x 1, I is J-1 x J-1. Combined is J-1 x J.
    neg_l_I <- cbind(-l_J_1, I_J_1)
    l_O <- cbind(l_M, O_M_J_1)

    part1 <- -A %*% neg_l_I
    part2 <- Gamma %*% l_O
    H1 <- part1 + part2

    # H2 = -Gamma'(-l, I) + B(l, 0)
    part3 <- -t(Gamma) %*% neg_l_I
    part4 <- B %*% l_O
    H2 <- part3 + part4

    # H = (-H1, -Gamma; -H2, -B)
    top_row <- cbind(-H1, -Gamma)
    bottom_rows <- cbind(-H2, -B)
    H <- rbind(top_row, bottom_rows)

    # Sigma = diag(sigma_mu^2, sigma_delta^2)
    Sigma <- diag(c(sigma_mu^2, sigma_delta^2))

    # Omega = diag(sigma_v_S^2, sigma_v_R^2)
    Omega <- diag(c(sigma_v_S^2, sigma_v_R^2))

    # Robust Matrix Inversion
    tryCatch(
        {
            Theta <- Omega + H %*% Sigma %*% t(H)
            Theta <- (Theta + t(Theta)) / 2
            Theta_inv <- solve(Theta)

            # Delta
            Omega_inv <- solve(Omega)
            Sigma_inv <- solve(Sigma)
            Delta_inv <- t(H) %*% Omega_inv %*% H + Sigma_inv
            Delta <- solve(Delta_inv)
            Delta <- (Delta + t(Delta)) / 2

            # W
            W <- Sigma %*% t(H) %*% Theta_inv

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

            if (n_cores > 1) {
                term2_list <- parallel::mclapply(1:nrow(data), calc_term2, mc.cores = n_cores)
                term2 <- unlist(term2_list)
            } else {
                term2 <- numeric(nrow(data))
                for (i in 1:nrow(data)) term2[i] <- calc_term2(i)
            }

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
