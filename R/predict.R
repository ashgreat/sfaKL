#' Predict inefficiencies from sfaKL model
#'
#' @param object An object of class sfaKL
#' @param newdata Optional new data frame
#' @param n_cores Number of cores for parallel execution (default 1)
#' @param ... Additional arguments
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

    # Extract dimensions
    J <- object$J
    M <- object$M
    if (is.null(J)) J <- 2 # Fallback
    if (is.null(M)) M <- 2

    n_input_shares <- J - 1
    n_output_shares <- M
    n_ineff <- J + M

    params <- object$raw_par
    mats <- get_matrices(params, J = J, M = M)

    # Parse parameters needed for residuals
    alpha_0 <- params[paste0("alpha", 2:J)]

    A <- matrix(0, n_input_shares, n_input_shares)
    for (j in 1:n_input_shares) {
        for (k in j:n_input_shares) {
            val <- params[paste0("A_", j, "_", k)]
            A[j, k] <- val
            A[k, j] <- val
        }
    }

    beta_0 <- params[paste0("beta", 1:M)]

    B <- matrix(0, n_output_shares, n_output_shares)
    for (m in 1:n_output_shares) {
        for (k in m:n_output_shares) {
            val <- params[paste0("B_", m, "_", k)]
            B[m, k] <- val
            B[k, m] <- val
        }
    }

    Gamma <- matrix(0, n_input_shares, n_output_shares)
    for (j in 1:n_input_shares) {
        for (m in 1:n_output_shares) {
            Gamma[j, m] <- params[paste0("Gamma_", j, "_", m)]
        }
    }

    # Extract Data
    S <- as.matrix(data[, paste0("S", 2:J), drop = FALSE])
    R <- as.matrix(data[, paste0("R", 1:M), drop = FALSE])
    ln_w <- as.matrix(data[, paste0("ln_w", 2:J, "_w1"), drop = FALSE])
    ln_p <- as.matrix(data[, paste0("ln_p", 1:M, "_w1"), drop = FALSE])

    # Calculate residuals
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

    Psi <- mats$Psi
    Delta <- mats$Delta

    n <- nrow(data)

    # Loop over observations
    calc_ineff <- function(i) {
        eps_i <- matrix(as.numeric(epsilon[i, ]), ncol = 1)
        mu_cond <- Psi %*% eps_i

        # Denominator: Phi(mu_cond; 0, Delta)
        denom <- mvtnorm::pmvnorm(lower = rep(-Inf, n_ineff), upper = as.numeric(mu_cond), mean = rep(0, n_ineff), sigma = Delta)
        denom <- max(as.numeric(denom), 1e-100)

        # Numerator: Gradient vector Phi*
        num <- numeric(n_ineff)
        for (j in 1:n_ineff) {
            s_j <- mu_cond[j]
            sigma_jj <- Delta[j, j]

            pdf_val <- dnorm(s_j, mean = 0, sd = sqrt(sigma_jj))

            idx <- (1:n_ineff)[-j]

            Sigma_11 <- Delta[idx, idx]
            Sigma_12 <- matrix(Delta[idx, j], ncol = 1)
            Sigma_21 <- matrix(Delta[j, idx], nrow = 1)
            Sigma_22 <- Delta[j, j]

            mu_c <- Sigma_12 %*% (1 / Sigma_22) * s_j
            Sigma_c <- Sigma_11 - Sigma_12 %*% (1 / Sigma_22) %*% Sigma_21

            upper_c <- mu_cond[idx]

            cdf_cond <- mvtnorm::pmvnorm(lower = rep(-Inf, n_ineff - 1), upper = as.numeric(upper_c), mean = as.numeric(mu_c), sigma = Sigma_c)

            num[j] <- pdf_val * as.numeric(cdf_cond)
        }

        return(mu_cond + Delta %*% (num / denom))
    }

    if (n_cores > 1) {
        eta_pred_list <- parallel::mclapply(1:n, calc_ineff, mc.cores = n_cores)
        eta_pred <- do.call(rbind, lapply(eta_pred_list, t))
    } else {
        eta_pred <- matrix(0, n, n_ineff)
        for (i in 1:n) {
            eta_pred[i, ] <- t(calc_ineff(i))
        }
    }

    # eta = (-mu1, ..., -muJ, delta1, ..., deltaM)
    # Return data frame with named columns
    result <- as.data.frame(eta_pred)
    names(result) <- c(paste0("mu", 1:J), paste0("delta", 1:M))

    # Flip signs for mu (since eta_j = -mu_j for inputs)
    for (j in 1:J) result[, j] <- -result[, j]

    return(result)
}
