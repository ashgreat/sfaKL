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
    # Extract dimensions
    J <- object$J
    M <- object$M
    if (is.null(J)) J <- 2
    if (is.null(M)) M <- 2

    required_std <- c(paste0("S", 2:J), paste0("R", 1:M), paste0("ln_w", 2:J, "_w1"), paste0("ln_p", 1:M, "_w1"))
    data_internal <- NULL

    if (is.null(newdata)) {
        data_internal <- object$data
    } else {
        # Map newdata columns to internal names using stored mapping when available
        if (!is.null(object$cols_map)) {
            if (all(names(object$cols_map) %in% names(newdata))) {
                data_internal <- newdata[, names(object$cols_map), drop = FALSE]
            } else if (all(object$cols_map %in% names(newdata))) {
                data_internal <- data.frame(matrix(NA, nrow = nrow(newdata), ncol = length(object$cols_map)))
                colnames(data_internal) <- names(object$cols_map)
                for (std_name in names(object$cols_map)) data_internal[[std_name]] <- newdata[[object$cols_map[[std_name]]]]
            } else {
                stop("newdata is missing required columns.")
            }
        } else {
            if (!all(required_std %in% names(newdata))) stop("newdata must include: ", paste(required_std, collapse = ", "))
            data_internal <- newdata[, required_std, drop = FALSE]
        }
    }

    # Ensure columns ordered
    data_internal <- data_internal[, required_std, drop = FALSE]

    n_input_shares <- J - 1
    n_output_shares <- M
    n_ineff <- J + M

    params <- object$raw_par
    profit_params <- assemble_profit_parameters(params, J = J, M = M)
    mats <- get_matrices(params, J = J, M = M)

    alpha_0 <- profit_params$alpha
    beta_0 <- profit_params$beta
    A <- profit_params$A
    B <- profit_params$B
    Gamma <- profit_params$Gamma

    # Extract Data
    S <- as.matrix(data_internal[, paste0("S", 2:J), drop = FALSE])
    R <- as.matrix(data_internal[, paste0("R", 1:M), drop = FALSE])
    ln_w <- as.matrix(data_internal[, paste0("ln_w", 2:J, "_w1"), drop = FALSE])
    ln_p <- as.matrix(data_internal[, paste0("ln_p", 1:M, "_w1"), drop = FALSE])

    # Calculate residuals
    eps_S <- matrix(0, nrow(data_internal), n_input_shares)
    for (j in 1:n_input_shares) {
        det_part <- alpha_0[j] + as.vector(A[j, , drop = FALSE] %*% t(ln_w)) + as.vector(Gamma[j, , drop = FALSE] %*% t(ln_p))
        eps_S[, j] <- -S[, j] - det_part
    }

    eps_R <- matrix(0, nrow(data_internal), n_output_shares)
    for (m in 1:n_output_shares) {
        det_part <- beta_0[m] + as.vector(B[m, , drop = FALSE] %*% t(ln_p)) + as.vector(t(Gamma[, m, drop = FALSE]) %*% t(ln_w))
        eps_R[, m] <- R[, m] - det_part
    }

    epsilon <- cbind(eps_S, eps_R)

    Psi <- mats$Psi
    Delta <- mats$Delta

    n <- nrow(data_internal)

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

    eta_pred_list <- safe_parallel_lapply(1:n, calc_ineff, n_cores = n_cores)
    eta_pred <- do.call(rbind, lapply(eta_pred_list, t))

    # eta = (-mu1, ..., -muJ, delta1, ..., deltaM)
    # Return data frame with named columns
    result <- as.data.frame(eta_pred)
    names(result) <- c(paste0("mu", 1:J), paste0("delta", 1:M))

    # Flip signs for mu (since eta_j = -mu_j for inputs)
    for (j in 1:J) result[, j] <- -result[, j]

    return(result)
}
