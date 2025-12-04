get_matrices <- function(params) {
    # Parse parameters (handling both raw log and exponentiated forms if needed)
    # Assuming params passed here are already exponentiated for sigmas if coming from result$par
    # But likelihood uses log_sigmas.
    # Let's standardize: The input params should be the raw optimization parameters (with log sigmas)

    alpha2 <- params["alpha2"]
    alpha22 <- params["alpha22"]

    beta1 <- params["beta1"]
    beta2 <- params["beta2"]
    beta11 <- params["beta11"]
    beta22 <- params["beta22"]
    beta12 <- params["beta12"]

    gamma21 <- params["gamma21"]
    gamma22 <- params["gamma22"]

    # Check if log_sigma or sigma is present
    if ("log_sigma_mu1" %in% names(params)) {
        sigma_mu1 <- exp(params["log_sigma_mu1"])
        sigma_mu2 <- exp(params["log_sigma_mu2"])
        sigma_delta1 <- exp(params["log_sigma_delta1"])
        sigma_delta2 <- exp(params["log_sigma_delta2"])

        sigma_v12 <- exp(params["log_sigma_v12"])
        sigma_v21 <- exp(params["log_sigma_v21"])
        sigma_v22 <- exp(params["log_sigma_v22"])
    } else {
        # Assume already exponentiated (e.g. from print output), but predict should use raw_par
        sigma_mu1 <- params["sigma_mu1"]
        sigma_mu2 <- params["sigma_mu2"]
        sigma_delta1 <- params["sigma_delta1"]
        sigma_delta2 <- params["sigma_delta2"]

        sigma_v12 <- params["sigma_v12"]
        sigma_v21 <- params["sigma_v21"]
        sigma_v22 <- params["sigma_v22"]
    }


    # Matrix construction
    A <- matrix(alpha22, 1, 1)
    B <- matrix(c(beta11, beta12, beta12, beta22), 2, 2)
    Gamma <- matrix(c(gamma21, gamma22), 1, 2)

    l_J_1 <- matrix(1, 1, 1)
    l_M <- matrix(1, 2, 1)
    I_J_1 <- diag(1)
    O_M_J_1 <- matrix(0, 2, 1)

    part1 <- -A %*% cbind(-l_J_1, I_J_1)
    part2 <- Gamma %*% cbind(l_M, O_M_J_1)
    H1 <- part1 + part2

    part3 <- -t(Gamma) %*% cbind(-l_J_1, I_J_1)
    part4 <- B %*% cbind(l_M, O_M_J_1)
    H2 <- part3 + part4

    top_row <- cbind(-H1, -Gamma)
    bottom_rows <- cbind(-H2, -B)
    H <- rbind(top_row, bottom_rows)

    Sigma <- diag(c(sigma_mu1^2, sigma_mu2^2, sigma_delta1^2, sigma_delta2^2))
    Omega <- diag(c(sigma_v12^2, sigma_v21^2, sigma_v22^2))

    Theta <- Omega + H %*% Sigma %*% t(H)

    Omega_inv <- solve(Omega)
    Sigma_inv <- solve(Sigma)
    Delta_inv <- t(H) %*% Omega_inv %*% H + Sigma_inv
    Delta <- solve(Delta_inv)

    Psi <- Sigma %*% t(H) %*% solve(Theta)

    return(list(H = H, Sigma = Sigma, Omega = Omega, Theta = Theta, Delta = Delta, Psi = Psi))
}
