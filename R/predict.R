#' Predict inefficiencies from sfaKL model
#'
#' @param object An object of class sfaKL
#' @param newdata Optional new data frame
#' @param ... Additional arguments
#' @param n_cores Number of cores for parallel execution (default 1)
#' @return A data frame with predicted inefficiencies
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats dnorm
#' @importFrom parallel mclapply
#' @export
predict.sfaKL <- function(object, newdata = NULL, n_cores = 1, ...) {
    if (is.null(newdata)) {
        data <- object$data
    } else {
        data <- newdata
    }

    params <- object$raw_par
    mats <- get_matrices(params)

    alpha2 <- params["alpha2"]
    alpha22 <- params["alpha22"]
    beta1 <- params["beta1"]
    beta2 <- params["beta2"]
    beta11 <- params["beta11"]
    beta22 <- params["beta22"]
    beta12 <- params["beta12"]
    gamma21 <- params["gamma21"]
    gamma22 <- params["gamma22"]

    # Calculate residuals
    eps1 <- -data$S2 - (alpha2 + alpha22 * data$ln_w2_w1 + gamma21 * data$ln_p1_w1 + gamma22 * data$ln_p2_w1)
    eps2_1 <- data$R1 - (beta1 + beta11 * data$ln_p1_w1 + beta12 * data$ln_p2_w1 + gamma21 * data$ln_w2_w1)
    eps2_2 <- data$R2 - (beta2 + beta12 * data$ln_p1_w1 + beta22 * data$ln_p2_w1 + gamma22 * data$ln_w2_w1)

    epsilon <- cbind(eps1, eps2_1, eps2_2)

    Psi <- mats$Psi
    Delta <- mats$Delta

    # E(eta | epsilon) = Psi * epsilon + Delta * (Phi* / Phi)

    n <- nrow(data)
    eta_pred <- matrix(0, n, 4)

    # Loop over observations
    # Use parallel::mclapply if n_cores > 1

    calc_ineff <- function(i) {
        eps_i <- matrix(as.numeric(epsilon[i, ]), ncol = 1)
        mu_cond <- Psi %*% eps_i # 4x1 vector

        # Denominator: Phi(mu_cond; 0, Delta)
        denom <- mvtnorm::pmvnorm(lower = rep(-Inf, 4), upper = as.numeric(mu_cond), mean = rep(0, 4), sigma = Delta)
        denom <- max(as.numeric(denom), 1e-100)

        # Numerator: Gradient vector Phi*
        num <- numeric(4)
        for (j in 1:4) {
            # Partial derivative w.r.t s_j
            # phi(s_j) * Phi(s_-j | s_j)

            s_j <- mu_cond[j]
            sigma_jj <- Delta[j, j]

            # Marginal PDF of s_j
            pdf_val <- dnorm(s_j, mean = 0, sd = sqrt(sigma_jj))

            # Conditional distribution of s_-j | s_j
            # Indices other than j
            idx <- (1:4)[-j]

            Sigma_11 <- Delta[idx, idx]
            Sigma_12 <- matrix(Delta[idx, j], ncol = 1)
            Sigma_21 <- matrix(Delta[j, idx], nrow = 1)
            Sigma_22 <- Delta[j, j]

            mu_c <- Sigma_12 %*% (1 / Sigma_22) * s_j
            Sigma_c <- Sigma_11 - Sigma_12 %*% (1 / Sigma_22) %*% Sigma_21

            upper_c <- mu_cond[idx]

            cdf_cond <- mvtnorm::pmvnorm(lower = rep(-Inf, 3), upper = as.numeric(upper_c), mean = as.numeric(mu_c), sigma = Sigma_c)

            num[j] <- pdf_val * as.numeric(cdf_cond)
        }

        # E[eta] = mu_cond + Delta %*% (num / denom)
        return(mu_cond + Delta %*% (num / denom))
    }

    if (n_cores > 1) {
        eta_pred_list <- parallel::mclapply(1:n, calc_ineff, mc.cores = n_cores)
        eta_pred <- do.call(rbind, lapply(eta_pred_list, t))
    } else {
        eta_pred <- matrix(0, n, 4)
        for (i in 1:n) {
            eta_pred[i, ] <- t(calc_ineff(i))
        }
    }

    # eta = (-mu1, -mu2, delta1, delta2)
    result <- data.frame(
        mu1 = -eta_pred[, 1],
        mu2 = -eta_pred[, 2],
        delta1 = eta_pred[, 3],
        delta2 = eta_pred[, 4]
    )

    return(result)
}
