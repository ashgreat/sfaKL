safe_parallel_lapply <- function(X, FUN, n_cores = 1) {
    if (.Platform$OS.type == "windows" || n_cores <= 1) {
        if (.Platform$OS.type == "windows" && n_cores > 1) {
            warning("parallel::mclapply is not available on Windows. Falling back to serial execution.")
        }
        return(lapply(X, FUN))
    }
    parallel::mclapply(X, FUN, mc.cores = n_cores)
}

assemble_profit_parameters <- function(params, J, M) {
    n_input_shares <- J - 1
    n_output_shares <- M

    # Betas (outputs)
    beta <- params[paste0("beta", 1:M)]
    beta[is.na(beta)] <- 0

    # Alpha intercepts for input shares (last one enforces adding-up)
    alpha_names <- paste0("alpha", 2:J)
    alpha <- params[alpha_names]
    alpha[is.na(alpha)] <- 0
    if (n_input_shares >= 1) {
        alpha[n_input_shares] <- 1 - sum(beta, na.rm = TRUE) - sum(alpha[-n_input_shares], na.rm = TRUE)
    }

    # Gamma matrix (input x output)
    Gamma <- matrix(0, n_input_shares, n_output_shares)
    for (j in 1:n_input_shares) {
        for (m in 1:n_output_shares) {
            nm <- paste0("Gamma_", j, "_", m)
            if (nm %in% names(params)) Gamma[j, m] <- params[[nm]]
        }
    }

    # A matrix (input Hessian) with homogeneity restriction
    A <- matrix(0, n_input_shares, n_input_shares)
    if (n_input_shares > 0) {
        for (j in 1:n_input_shares) {
            for (k in j:min(n_input_shares - 1, n_input_shares)) {
                nm <- paste0("A_", j, "_", k)
                if (nm %in% names(params)) {
                    A[j, k] <- params[[nm]]
                    A[k, j] <- params[[nm]]
                }
            }
        }

        if (n_input_shares == 1) {
            A[1, 1] <- -sum(Gamma[1, ])
        } else {
            for (k in 1:(n_input_shares - 1)) {
                A[n_input_shares, k] <- -sum(A[1:(n_input_shares - 1), k]) - sum(Gamma[k, ])
                A[k, n_input_shares] <- A[n_input_shares, k]
            }
            A[n_input_shares, n_input_shares] <- -sum(A[1:(n_input_shares - 1), n_input_shares]) - sum(Gamma[n_input_shares, ])
        }
    }

    # B matrix (output Hessian) with homogeneity restriction
    B <- matrix(0, n_output_shares, n_output_shares)
    if (n_output_shares > 0) {
        for (m in 1:n_output_shares) {
            for (k in m:min(n_output_shares - 1, n_output_shares)) {
                nm <- paste0("B_", m, "_", k)
                if (nm %in% names(params)) {
                    B[m, k] <- params[[nm]]
                    B[k, m] <- params[[nm]]
                }
            }
        }

        if (n_output_shares == 1) {
            B[1, 1] <- -sum(Gamma[, 1])
        } else {
            gamma_col_sums <- apply(Gamma, 2, sum)
            for (n in 1:(n_output_shares - 1)) {
                B[n_output_shares, n] <- -sum(B[1:(n_output_shares - 1), n]) - gamma_col_sums[n]
                B[n, n_output_shares] <- B[n_output_shares, n]
            }
            B[n_output_shares, n_output_shares] <- -sum(B[1:(n_output_shares - 1), n_output_shares]) - gamma_col_sums[n_output_shares]
        }
    }

    list(alpha = alpha, beta = beta, A = A, B = B, Gamma = Gamma)
}

build_cov_matrix <- function(params, prefix, dim, default_sd = 0.1) {
    L <- matrix(0, dim, dim)
    for (i in 1:dim) {
        for (j in 1:i) {
            nm <- paste0("chol_", prefix, "_", i, "_", j)
            val <- params[[nm]]
            if (is.null(val) || is.na(val)) val <- if (i == j) log(default_sd) else 0
            L[i, j] <- if (i == j) exp(val) else val
        }
    }
    cov <- L %*% t(L)
    cov <- (cov + t(cov)) / 2
    list(L = L, cov = cov)
}

assemble_covariance_matrices <- function(params, J, M) {
    n_eq <- (J - 1) + M

    Sigma_mu <- build_cov_matrix(params, "Sigma_mu", J)$cov
    Sigma_delta <- build_cov_matrix(params, "Sigma_delta", M)$cov

    Sigma <- matrix(0, J + M, J + M)
    Sigma[1:J, 1:J] <- Sigma_mu
    Sigma[(J + 1):(J + M), (J + 1):(J + M)] <- Sigma_delta

    Omega <- build_cov_matrix(params, "Omega", n_eq)$cov

    list(Sigma_mu = Sigma_mu, Sigma_delta = Sigma_delta, Sigma = Sigma, Omega = Omega)
}

get_matrices <- function(params, J = 2, M = 2) {
    n_input_shares <- J - 1
    n_output_shares <- M

    profit_params <- assemble_profit_parameters(params, J, M)
    cov_params <- assemble_covariance_matrices(params, J, M)

    A <- profit_params$A
    B <- profit_params$B
    Gamma <- profit_params$Gamma

    Sigma <- cov_params$Sigma
    Omega <- cov_params$Omega

    l_J_1 <- matrix(1, n_input_shares, 1)
    l_M <- matrix(1, n_output_shares, 1)
    I_J_1 <- diag(1, n_input_shares)
    O_M_J_1 <- matrix(0, n_output_shares, n_input_shares)

    neg_l_I <- cbind(-l_J_1, I_J_1)
    l_O <- cbind(l_M, O_M_J_1)

    H1 <- -A %*% neg_l_I + Gamma %*% l_O
    H2 <- -t(Gamma) %*% neg_l_I + B %*% l_O

    H <- rbind(cbind(-H1, -Gamma), cbind(-H2, -B))

    Theta <- Omega + H %*% Sigma %*% t(H)
    Theta <- (Theta + t(Theta)) / 2

    # Add small regularization for numerical stability
    ridge <- 1e-10
    Omega_inv <- solve(Omega + ridge * diag(nrow(Omega)))
    Sigma_inv <- solve(Sigma + ridge * diag(nrow(Sigma)))
    Delta_inv <- t(H) %*% Omega_inv %*% H + Sigma_inv
    Delta <- solve(Delta_inv + ridge * diag(nrow(Delta_inv)))
    Delta <- (Delta + t(Delta)) / 2

    Psi <- Sigma %*% t(H) %*% solve(Theta + ridge * diag(nrow(Theta)))

    list(
        H = H,
        Sigma = Sigma,
        Omega = Omega,
        Theta = Theta,
        Delta = Delta,
        Psi = Psi,
        alpha = profit_params$alpha,
        beta = profit_params$beta,
        A = A,
        B = B,
        Gamma = Gamma
    )
}
