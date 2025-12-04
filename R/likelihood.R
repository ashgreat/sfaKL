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
    # We need to calculate residuals inside the function

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
    sigma_mu1 <- exp(params["log_sigma_mu1"])
    sigma_mu2 <- exp(params["log_sigma_mu2"])
    sigma_delta1 <- exp(params["log_sigma_delta1"])
    sigma_delta2 <- exp(params["log_sigma_delta2"])

    sigma_v12 <- exp(params["log_sigma_v12"])
    sigma_v21 <- exp(params["log_sigma_v21"])
    sigma_v22 <- exp(params["log_sigma_v22"])

    # Construct Residuals (epsilon)
    # -S2 - (...)
    # R1 - (...)
    # R2 - (...)

    eps1 <- -data$S2 - (alpha2 + alpha22 * data$ln_w2_w1 + gamma21 * data$ln_p1_w1 + gamma22 * data$ln_p2_w1)
    eps2_1 <- data$R1 - (beta1 + beta11 * data$ln_p1_w1 + beta12 * data$ln_p2_w1 + gamma21 * data$ln_w2_w1)
    eps2_2 <- data$R2 - (beta2 + beta12 * data$ln_p1_w1 + beta22 * data$ln_p2_w1 + gamma22 * data$ln_w2_w1) # Note: beta21 = beta12

    # Matrix construction
    # J=2, M=2
    # Dimensions:
    # epsilon has length J+M-1 = 3
    # mu has length J = 2
    # delta has length M = 2
    # eta has length J+M = 4

    # Matrices A, B, Gamma
    # A is (J-1)x(J-1) => 1x1
    A <- matrix(alpha22, 1, 1)

    # B is MxM => 2x2
    B <- matrix(c(beta11, beta12, beta12, beta22), 2, 2)

    # Gamma is (J-1)xM => 1x2
    Gamma <- matrix(c(gamma21, gamma22), 1, 2)

    # Helper matrices
    l_J_1 <- matrix(1, 1, 1) # J-1 = 1
    l_M <- matrix(1, 2, 1) # M = 2
    I_J_1 <- diag(1)
    O_M_J_1 <- matrix(0, 2, 1)

    # H1 = -A(-l_J-1, I_J-1) + Gamma(l_M, O_M_(J-1))
    # -A(-1, 1) = (alpha22, -alpha22)
    part1 <- -A %*% cbind(-l_J_1, I_J_1)
    # Gamma(1, 0; 1, 0)
    part2 <- Gamma %*% cbind(l_M, O_M_J_1)
    H1 <- part1 + part2

    # H2 = -Gamma'(-l_J-1, I_J-1) + B(l_M, O_M_(J-1))
    part3 <- -t(Gamma) %*% cbind(-l_J_1, I_J_1)
    part4 <- B %*% cbind(l_M, O_M_J_1)
    H2 <- part3 + part4

    # H = (-H1, -Gamma; -H2, -B)
    # Dimensions: (1+2) x (2+2) = 3 x 4
    top_row <- cbind(-H1, -Gamma)
    bottom_rows <- cbind(-H2, -B)
    H <- rbind(top_row, bottom_rows)

    # Sigma = diag(sigma_mu, sigma_delta)
    Sigma <- diag(c(sigma_mu1^2, sigma_mu2^2, sigma_delta1^2, sigma_delta2^2))

    # Omega = diag(sigma_v)
    # Assuming independent noise for now as per paper application (diagonal Omega)
    Omega <- diag(c(sigma_v12^2, sigma_v21^2, sigma_v22^2))

    # Theta = Omega + H Sigma H'
    Theta <- Omega + H %*% Sigma %*% t(H)

    # Delta = (H' Omega^-1 H + Sigma^-1)^-1
    # Use solve() for inverse
    Omega_inv <- solve(Omega)
    Sigma_inv <- solve(Sigma)
    Delta_inv <- t(H) %*% Omega_inv %*% H + Sigma_inv
    Delta <- solve(Delta_inv)

    # Calculate Likelihood terms
    # 1. Phi_3(epsilon; 0, Theta)
    # We can use dmvnorm for the PDF part
    # But we need to loop over observations or use vectorized approach

    epsilon <- cbind(eps1, eps2_1, eps2_2) # N x 3 matrix

    # Term 1: PDF of N(0, Theta) evaluated at epsilon
    term1 <- mvtnorm::dmvnorm(epsilon, mean = rep(0, 3), sigma = Theta, log = TRUE)

    # Term 2: CDF term
    # Phi_4(Sigma H' Theta^-1 epsilon; 0, Delta)
    # Let W = Sigma H' Theta^-1
    W <- Sigma %*% t(H) %*% solve(Theta) # 4 x 3 matrix

    # Argument for CDF: W %*% epsilon[i,]
    # This results in a 4x1 vector for each observation
    # We need to calculate P(X <= arg) where X ~ N(0, Delta)

    # This part is computationally intensive.
    # We can try to optimize by checking if Delta is diagonal (unlikely)

    # Loop over observations (can be slow in R, but unavoidable with pmvnorm)
    # To speed up, maybe we can use a C++ implementation later, but for now R loop.

    term2 <- numeric(nrow(data))

    # Pre-calculate arguments
    args <- epsilon %*% t(W) # N x 4 matrix

    # Define lower and upper bounds
    lower <- rep(-Inf, 4)
    mean_delta <- rep(0, 4)

    for (i in 1:nrow(data)) {
        # pmvnorm returns an object with attributes, we need just the value
        # suppress warnings about precision
        prob <- suppressWarnings(mvtnorm::pmvnorm(lower = lower, upper = args[i, ], mean = mean_delta, sigma = Delta))
        term2[i] <- log(max(as.numeric(prob), 1e-100)) # Avoid log(0)
    }

    # Term 3: Denominator Phi_4(0; 0, Sigma)
    # Since Sigma is diagonal, this is 0.5^4 = 0.0625
    # log(0.0625)
    term3 <- log(0.0625)

    # Log Likelihood per observation
    ll_i <- term1 + term2 - term3

    # Return negative sum
    return(-sum(ll_i))
}
