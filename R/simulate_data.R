#' Simulate data for sfaKL model
#'
#' @param n Number of observations
#' @param seed Random seed
#' @return A list containing the simulated data and true parameters
#' @export
sfaKL_simulate <- function(n = 1000, seed = 123) {
    set.seed(seed)

    # 1. Exogenous prices (normalized by w1)
    # ln(w2/w1), ln(p1/w1), ln(p2/w1)
    ln_w2_w1 <- rnorm(n, 0, 0.5)
    ln_p1_w1 <- rnorm(n, 0, 0.5)
    ln_p2_w1 <- rnorm(n, 0, 0.5)

    # 2. True Parameters
    # Alpha (Input share params)
    alpha2 <- -0.3
    alpha22 <- 0.1

    # Beta (Revenue share params)
    beta1 <- 0.4
    beta2 <- 0.3
    beta11 <- 0.1
    beta22 <- 0.1
    beta12 <- -0.05 # Symmetry: beta12 = beta21
    beta21 <- beta12

    # Gamma (Cross params)
    gamma21 <- 0.05
    gamma22 <- 0.04

    # Variance parameters for inefficiencies (Half-Normal)
    sigma_mu1 <- 0.1
    sigma_mu2 <- 0.1
    sigma_delta1 <- 0.1
    sigma_delta2 <- 0.1

    # Variance parameters for noise (Normal)
    sigma_v12 <- 0.05 # Error in S2
    sigma_v21 <- 0.05 # Error in R1
    sigma_v22 <- 0.05 # Error in R2

    # 3. Simulate Inefficiencies (Half-Normal)
    # mu <= 0, so we simulate |N| and take negative
    mu1 <- -abs(rnorm(n, 0, sigma_mu1))
    mu2 <- -abs(rnorm(n, 0, sigma_mu2))
    delta1 <- abs(rnorm(n, 0, sigma_delta1))
    delta2 <- abs(rnorm(n, 0, sigma_delta2))

    # 4. Calculate Composite Inefficiency Terms (u and omega)
    # u2 = -alpha22*(mu2 - mu1) - gamma21*(delta1 - mu1) - gamma22*(delta2 - mu1)
    u2 <- -alpha22 * (mu2 - mu1) - gamma21 * (delta1 - mu1) - gamma22 * (delta2 - mu1)

    # omega1 = -beta11*(delta1 - mu1) - beta12*(delta2 - mu1) - gamma21*(mu2 - mu1)
    omega1 <- -beta11 * (delta1 - mu1) - beta12 * (delta2 - mu1) - gamma21 * (mu2 - mu1)

    # omega2 = -beta21*(delta1 - mu1) - beta22*(delta2 - mu1) - gamma22*(mu2 - mu1)
    omega2 <- -beta21 * (delta1 - mu1) - beta22 * (delta2 - mu1) - gamma22 * (mu2 - mu1)

    # 5. Simulate Noise
    v12 <- rnorm(n, 0, sigma_v12)
    v21 <- rnorm(n, 0, sigma_v21)
    v22 <- rnorm(n, 0, sigma_v22)

    # 6. Calculate Shares
    # -S2 = alpha2 + alpha22*ln_w2_w1 + gamma21*ln_p1_w1 + gamma22*ln_p2_w1 + u2 + v12
    # S2 = -(...)
    neg_S2 <- alpha2 + alpha22 * ln_w2_w1 + gamma21 * ln_p1_w1 + gamma22 * ln_p2_w1 + u2 + v12
    S2 <- -neg_S2

    # R1 = beta1 + beta11*ln_p1_w1 + beta12*ln_p2_w1 + gamma21*ln_w2_w1 + omega1 + v21
    R1 <- beta1 + beta11 * ln_p1_w1 + beta12 * ln_p2_w1 + gamma21 * ln_w2_w1 + omega1 + v21

    # R2 = beta2 + beta21*ln_p1_w1 + beta22*ln_p2_w1 + gamma22*ln_w2_w1 + omega2 + v22
    R2 <- beta2 + beta21 * ln_p1_w1 + beta22 * ln_p2_w1 + gamma22 * ln_w2_w1 + omega2 + v22

    data <- data.frame(
        S2 = S2,
        R1 = R1,
        R2 = R2,
        ln_w2_w1 = ln_w2_w1,
        ln_p1_w1 = ln_p1_w1,
        ln_p2_w1 = ln_p2_w1
    )

    true_params <- list(
        alpha = c(alpha2 = alpha2, alpha22 = alpha22),
        beta = c(beta1 = beta1, beta2 = beta2, beta11 = beta11, beta22 = beta22, beta12 = beta12),
        gamma = c(gamma21 = gamma21, gamma22 = gamma22),
        sigma_ineff = c(mu1 = sigma_mu1, mu2 = sigma_mu2, delta1 = sigma_delta1, delta2 = sigma_delta2),
        sigma_noise = c(v12 = sigma_v12, v21 = sigma_v21, v22 = sigma_v22)
    )

    return(list(data = data, true_params = true_params))
}
