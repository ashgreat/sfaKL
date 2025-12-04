#' Simulate data for sfaKL model
#'
#' @param n Number of observations
#' @param J Number of inputs (default 2)
#' @param M Number of outputs (default 2)
#' @param seed Random seed
#' @return A list containing the simulated data and true parameters
#' @export
sfaKL_simulate <- function(n = 1000, J = 2, M = 2, seed = 123) {
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

    # 2. Parameters
    # Alpha (Input share params) - Symmetric (J-1)x(J-1)
    alpha_0 <- runif(n_input_shares, -0.4, -0.2)
    names(alpha_0) <- paste0("alpha", 2:J)

    A <- matrix(runif(n_input_shares^2, 0.05, 0.15), n_input_shares, n_input_shares)
    A <- (A + t(A)) / 2 # Symmetry

    # Beta (Revenue share params) - Symmetric MxM
    beta_0 <- runif(n_output_shares, 0.3, 0.5)
    names(beta_0) <- paste0("beta", 1:M)

    B <- matrix(runif(n_output_shares^2, 0.05, 0.15), n_output_shares, n_output_shares)
    B <- (B + t(B)) / 2 # Symmetry

    # Gamma (Cross params) - (J-1)xM
    Gamma <- matrix(runif(n_input_shares * n_output_shares, 0.03, 0.07), n_input_shares, n_output_shares)

    # Variances
    sigma_mu <- rep(0.1, J)
    names(sigma_mu) <- paste0("mu", 1:J)

    sigma_delta <- rep(0.1, M)
    names(sigma_delta) <- paste0("delta", 1:M)

    sigma_v <- rep(0.05, n_eq)
    names(sigma_v) <- c(paste0("v_S", 2:J), paste0("v_R", 1:M))

    # 3. Inefficiencies
    mu <- matrix(0, n, J)
    for (j in 1:J) mu[, j] <- -abs(rnorm(n, 0, sigma_mu[j]))

    delta <- matrix(0, n, M)
    for (m in 1:M) delta[, m] <- abs(rnorm(n, 0, sigma_delta[m]))

    # 4. Composite Errors
    # u (length J-1)
    u <- matrix(0, n, n_input_shares)
    for (j in 1:n_input_shares) {
        # u_j corresponds to input j+1
        # Formula: - sum_k=2^J A_jk (mu_k - mu_1) - sum_m=1^M G_jm (delta_m - mu_1)

        term1 <- 0
        for (k in 1:n_input_shares) {
            term1 <- term1 + A[j, k] * (mu[, k + 1] - mu[, 1])
        }

        term2 <- 0
        for (m in 1:M) {
            term2 <- term2 + Gamma[j, m] * (delta[, m] - mu[, 1])
        }
        u[, j] <- -term1 - term2
    }

    # omega (length M)
    omega <- matrix(0, n, M)
    for (m in 1:M) {
        # Formula: - sum_k=1^M B_mk (delta_k - mu_1) - sum_j=2^J G_jm (mu_j - mu_1)
        # Note G_jm is Gamma[j-1, m]

        term1 <- 0
        for (k in 1:M) {
            term1 <- term1 + B[m, k] * (delta[, k] - mu[, 1])
        }

        term2 <- 0
        for (j in 1:n_input_shares) {
            term2 <- term2 + Gamma[j, m] * (mu[, j + 1] - mu[, 1])
        }
        omega[, m] <- -term1 - term2
    }

    # 5. Noise
    v <- matrix(0, n, n_eq)
    for (i in 1:n_eq) v[, i] <- rnorm(n, 0, sigma_v[i])

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
        sigma_mu = sigma_mu,
        sigma_delta = sigma_delta,
        sigma_v = sigma_v
    )

    return(list(data = data, true_params = true_params))
}
