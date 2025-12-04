#' Log-Likelihood for sfaKL model
#'
#' @param params Named vector of parameters
#' @param data Data frame containing the variables
#' @return Negative log-likelihood
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats dnorm
sfaKL_loglik <- function(params, data) {
    # Extract data
    # Dependent variables (residuals from the share equations)

    # Parse parameters
    alpha2 <- params["alpha2"]
    alpha22 <- params["alpha22"]

    beta1 <- params["beta1"]
    beta2 <- params["beta2"]
    beta11 <- params["beta11"]
    beta22 <- params["beta22"]
    beta12 <- params["beta12"]

    gamma21 <- params["gamma21"]
    gamma22 <- params["gamma22"]

    # Standard deviations (exponentiated to ensure positivity)
    # Limit extreme values to prevent numerical issues
    sigma_mu1 <- exp(min(params["log_sigma_mu1"], 10))
    sigma_mu2 <- exp(min(params["log_sigma_mu2"], 10))
    sigma_delta1 <- exp(min(params["log_sigma_delta1"], 10))
    sigma_delta2 <- exp(min(params["log_sigma_delta2"], 10))

    sigma_v12 <- exp(min(params["log_sigma_v12"], 10))
    sigma_v21 <- exp(min(params["log_sigma_v21"], 10))
    sigma_v22 <- exp(min(params["log_sigma_v22"], 10))

    # Construct Residuals (epsilon)
    eps1 <- -data$S2 - (alpha2 + alpha22 * data$ln_w2_w1 + gamma21 * data$ln_p1_w1 + gamma22 * data$ln_p2_w1)
    eps2_1 <- data$R1 - (beta1 + beta11 * data$ln_p1_w1 + beta12 * data$ln_p2_w1 + gamma21 * data$ln_w2_w1)
    eps2_2 <- data$R2 - (beta2 + beta12 * data$ln_p1_w1 + beta22 * data$ln_p2_w1 + gamma22 * data$ln_w2_w1)

    epsilon <- cbind(eps1, eps2_1, eps2_2) # N x 3 matrix

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

    # Robust matrix inversion
    tryCatch(
        {
            Theta <- Omega + H %*% Sigma %*% t(H)
            Theta <- (Theta + t(Theta)) / 2 # Force symmetry
            Theta_inv <- solve(Theta)

            Omega_inv <- solve(Omega)
            Sigma_inv <- solve(Sigma)
            Delta_inv <- t(H) %*% Omega_inv %*% H + Sigma_inv
            Delta <- solve(Delta_inv)
            Delta <- (Delta + t(Delta)) / 2 # Force symmetry

            # W = Sigma H' Theta^-1
            W <- Sigma %*% t(H) %*% Theta_inv

            # Term 1: PDF of N(0, Theta)
            term1 <- mvtnorm::dmvnorm(epsilon, mean = rep(0, 3), sigma = Theta, log = TRUE)

            # Term 2: CDF term
            args <- epsilon %*% t(W) # N x 4 matrix

            if (any(is.na(args))) {
                return(1e10)
            }

            term2 <- numeric(nrow(data))
            lower <- rep(-Inf, 4)
            mean_delta <- rep(0, 4)

            for (i in 1:nrow(data)) {
                prob <- suppressWarnings(mvtnorm::pmvnorm(lower = lower, upper = args[i, ], mean = mean_delta, sigma = Delta))
                val <- as.numeric(prob)
                if (val <= 0 || is.na(val)) val <- 1e-100
                term2[i] <- log(val)
            }

            # Term 3: Denominator
            term3 <- log(0.0625)

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
