#' Simulate data for sfaKL model
#'
#' @param n Number of observations
#' @param J Number of inputs (default 2)
#' @param M Number of outputs (default 2)
#' @param seed Random seed
#' @param Sigma_mu Optional JxJ covariance matrix for input inefficiency (default: random SPD)
#' @param Sigma_delta Optional MxM covariance matrix for output inefficiency (default: random SPD)
#' @param Omega Optional (J+M-1)x(J+M-1) covariance matrix for noise (default: random SPD)
#' @return A list containing the simulated data and true parameters
#' @export
sfaKL_simulate <- function(n = 1000, J = 2, M = 2, seed = 123, Sigma_mu = NULL, Sigma_delta = NULL, Omega = NULL) {
    set.seed(seed)

    # Dimensions
    n_input_shares <- J - 1
    n_output_shares <- M
    n_eq <- n_input_shares + n_output_shares

    # 1. Exogenous prices
    # ln(w_j/w_1) for j=2..J
    ln_w <- matrix(rnorm(n * n_input_shares, 0, 0.5), nrow = n)
    colnames(ln_w) <- paste0("ln_w", 2:J, "_w1")

    # ln(p_m/w_1) for m=1..M
    ln_p <- matrix(rnorm(n * n_output_shares, 0, 0.5), nrow = n)
    colnames(ln_p) <- paste0("ln_p", 1:M, "_w1")

    # 2. Parameters (free) -> enforce restrictions via assemble_profit_parameters
    params <- c()
    params[paste0("alpha", 2:J)] <- runif(n_input_shares, -0.4, -0.1)
    params[paste0("beta", 1:M)] <- runif(n_output_shares, 0.3, 0.6)

    for (j in 1:(J - 1)) {
        for (k in j:max(j, J - 2)) {
            if (k <= (J - 2)) params[paste0("A_", j, "_", k)] <- runif(1, 0.05, 0.15)
        }
    }

    for (m in 1:M) {
        for (k in m:max(m, M - 1)) {
            if (k <= (M - 1)) params[paste0("B_", m, "_", k)] <- runif(1, 0.05, 0.15)
        }
    }

    for (j in 1:(J - 1)) {
        for (m in 1:M) params[paste0("Gamma_", j, "_", m)] <- runif(1, 0.03, 0.08)
    }

    prof <- assemble_profit_parameters(params, J, M)
    alpha_0 <- prof$alpha
    A <- prof$A
    beta_0 <- prof$beta
    B <- prof$B
    Gamma <- prof$Gamma

    # Covariance matrices (random SPD if not supplied)
    if (is.null(Sigma_mu)) {
        L_mu <- matrix(rnorm(J^2, sd = 0.05), J, J)
        Sigma_mu <- crossprod(L_mu) + diag(0.05, J)
    }

    if (is.null(Sigma_delta)) {
        L_delta <- matrix(rnorm(M^2, sd = 0.05), M, M)
        Sigma_delta <- crossprod(L_delta) + diag(0.05, M)
    }

    if (is.null(Omega)) {
        L_omega <- matrix(rnorm(n_eq^2, sd = 0.02), n_eq, n_eq)
        Omega <- crossprod(L_omega) + diag(0.05, n_eq)
    }

    # 3. Inefficiencies
    mu_raw <- mvtnorm::rmvnorm(n, sigma = Sigma_mu)
    mu <- -abs(mu_raw)

    delta_raw <- mvtnorm::rmvnorm(n, sigma = Sigma_delta)
    delta <- abs(delta_raw)

    # 4. Composite Errors
    # u (length J-1)
    mu_diff <- mu[, 2:J, drop = FALSE] - mu[, 1]
    delta_diff <- delta - mu[, 1]

    u <- matrix(0, n, n_input_shares)
    for (j in 1:n_input_shares) {
        term1 <- as.vector(mu_diff %*% A[j, ])
        term2 <- as.vector(delta_diff %*% Gamma[j, ])
        u[, j] <- -(term1 + term2)
    }

    omega <- matrix(0, n, M)
    for (m_idx in 1:M) {
        term1 <- as.vector(delta_diff %*% B[m_idx, ])
        term2 <- as.vector(mu_diff %*% Gamma[, m_idx])
        omega[, m_idx] <- -(term1 + term2)
    }

    # 5. Noise
    v <- mvtnorm::rmvnorm(n, sigma = Omega)

    # 6. Shares
    # S (length J-1)
    S <- matrix(0, n, n_input_shares)
    colnames(S) <- paste0("S", 2:J)

    for (j in 1:n_input_shares) {
        # S_j = - (alpha_j + sum A_jk ln_w + sum G_jm ln_p + u_j + v_j)

        det_part <- alpha_0[j]
        det_part <- det_part + as.vector(A[j, , drop = FALSE] %*% t(ln_w))
        det_part <- det_part + as.vector(Gamma[j, , drop = FALSE] %*% t(ln_p))

        S[, j] <- -(det_part + u[, j] + v[, j])
    }

    # R (length M)
    R <- matrix(0, n, M)
    colnames(R) <- paste0("R", 1:M)

    for (m in 1:M) {
        # R_m = beta_m + sum B_mk ln_p + sum G_jm ln_w + omega_m + v_{M+m}

        det_part <- beta_0[m]
        det_part <- det_part + as.vector(B[m, , drop = FALSE] %*% t(ln_p))
        det_part <- det_part + as.vector(t(Gamma[, m, drop = FALSE]) %*% t(ln_w))

        R[, m] <- det_part + omega[, m] + v[, n_input_shares + m]
    }

    data <- cbind(as.data.frame(S), as.data.frame(R), as.data.frame(ln_w), as.data.frame(ln_p))

    true_params <- list(
        J = J,
        M = M,
        alpha_0 = alpha_0,
        A = A,
        beta_0 = beta_0,
        B = B,
        Gamma = Gamma,
        Sigma_mu = Sigma_mu,
        Sigma_delta = Sigma_delta,
        Omega = Omega
    )

    return(list(data = data, true_params = true_params))
}
