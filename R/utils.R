get_matrices <- function(params, J = 2, M = 2) {
    # Dimensions
    n_input_shares <- J - 1
    n_output_shares <- M

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
    sigma_mu <- numeric(J)
    sigma_delta <- numeric(M)
    sigma_v_S <- numeric(n_input_shares)
    sigma_v_R <- numeric(n_output_shares)

    if (any(grepl("log_sigma_", names(params)))) {
        sigma_mu <- exp(params[paste0("log_sigma_mu", 1:J)])
        sigma_delta <- exp(params[paste0("log_sigma_delta", 1:M)])
        sigma_v_S <- exp(params[paste0("log_sigma_v_S", 2:J)])
        sigma_v_R <- exp(params[paste0("log_sigma_v_R", 1:M)])
    } else {
        sigma_mu <- params[paste0("sigma_mu", 1:J)]
        sigma_delta <- params[paste0("sigma_delta", 1:M)]
        sigma_v_S <- params[paste0("sigma_v_S", 2:J)]
        sigma_v_R <- params[paste0("sigma_v_R", 1:M)]
    }

    # Matrix Construction
    l_J_1 <- matrix(1, n_input_shares, 1)
    l_M <- matrix(1, n_output_shares, 1)
    I_J_1 <- diag(1, n_input_shares)
    O_M_J_1 <- matrix(0, n_output_shares, n_input_shares)

    neg_l_I <- cbind(-l_J_1, I_J_1)
    l_O <- cbind(l_M, O_M_J_1)

    part1 <- -A %*% neg_l_I
    part2 <- Gamma %*% l_O
    H1 <- part1 + part2

    part3 <- -t(Gamma) %*% neg_l_I
    part4 <- B %*% l_O
    H2 <- part3 + part4

    top_row <- cbind(-H1, -Gamma)
    bottom_rows <- cbind(-H2, -B)
    H <- rbind(top_row, bottom_rows)

    Sigma <- diag(c(sigma_mu^2, sigma_delta^2))
    Omega <- diag(c(sigma_v_S^2, sigma_v_R^2))

    Theta <- Omega + H %*% Sigma %*% t(H)
    Theta <- (Theta + t(Theta)) / 2

    Omega_inv <- solve(Omega)
    Sigma_inv <- solve(Sigma)
    Delta_inv <- t(H) %*% Omega_inv %*% H + Sigma_inv
    Delta <- solve(Delta_inv)
    Delta <- (Delta + t(Delta)) / 2

    Psi <- Sigma %*% t(H) %*% solve(Theta)

    return(list(H = H, Sigma = Sigma, Omega = Omega, Theta = Theta, Delta = Delta, Psi = Psi))
}
