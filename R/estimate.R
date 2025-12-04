#' Estimate sfaKL model
#'
#' @param data Data frame containing the variables
#' @param share_input Name of the input share column (default "S2") or vector of names
#' @param share_output Names of the output share columns (default c("R1", "R2"))
#' @param price_input_ratio Name of the log input price ratio column (default "ln_w2_w1") or vector
#' @param price_output_ratios Names of the log output-to-input price ratio columns (default c("ln_p1_w1", "ln_p2_w1"))
#' @param start_params Optional named vector of starting parameters
#' @param method Optimization method (default "Nelder-Mead")
#' @param control List of control parameters for optim
#' @param n_cores Number of cores for parallel execution (default 1)
#' @return An object of class sfaKL
#' @export
sfaKL_estimate <- function(data,
                           share_input = "S2",
                           share_output = c("R1", "R2"),
                           price_input_ratio = "ln_w2_w1",
                           price_output_ratios = c("ln_p1_w1", "ln_p2_w1"),
                           start_params = NULL,
                           method = "Nelder-Mead",
                           control = list(maxit = 5000),
                           n_cores = 1) {
    # Infer Dimensions
    n_input_shares <- length(share_input)
    n_output_shares <- length(share_output)
    J <- n_input_shares + 1
    M <- n_output_shares

    # Validate price ratio lengths
    if (length(price_input_ratio) != n_input_shares) stop("Length of price_input_ratio must match share_input")
    if (length(price_output_ratios) != n_output_shares) stop("Length of price_output_ratios must match share_output")

    # Map user columns to internal names
    cols_map <- c()
    for (i in 1:n_input_shares) cols_map[paste0("S", i + 1)] <- share_input[i]
    for (i in 1:n_output_shares) cols_map[paste0("R", i)] <- share_output[i]
    for (i in 1:n_input_shares) cols_map[paste0("ln_w", i + 1, "_w1")] <- price_input_ratio[i]
    for (i in 1:n_output_shares) cols_map[paste0("ln_p", i, "_w1")] <- price_output_ratios[i]

    # Check if columns exist
    if (!all(cols_map %in% names(data))) {
        missing <- cols_map[!cols_map %in% names(data)]
        stop("Missing columns in data: ", paste(missing, collapse = ", "))
    }

    # Create internal data frame with standard names
    data_internal <- data.frame(matrix(NA, nrow = nrow(data), ncol = length(cols_map)))
    colnames(data_internal) <- names(cols_map)
    for (std_name in names(cols_map)) {
        data_internal[[std_name]] <- data[[cols_map[[std_name]]]]
    }

    # Default starting parameters if not provided
    if (is.null(start_params)) {
        start_params <- c()

        # Alpha vector
        for (j in 2:J) start_params[paste0("alpha", j)] <- -0.1

        # A matrix (upper triangular + diagonal)
        for (j in 1:(J - 1)) {
            for (k in j:(J - 1)) {
                start_params[paste0("A_", j, "_", k)] <- if (j == k) 0.1 else 0.0
            }
        }

        # Beta vector
        for (m in 1:M) start_params[paste0("beta", m)] <- 0.4

        # B matrix (upper triangular + diagonal)
        for (m in 1:M) {
            for (k in m:M) {
                start_params[paste0("B_", m, "_", k)] <- if (m == k) 0.1 else 0.0
            }
        }

        # Gamma matrix (full)
        for (j in 1:(J - 1)) {
            for (m in 1:M) {
                start_params[paste0("Gamma_", j, "_", m)] <- 0.05
            }
        }

        # Sigmas
        for (j in 1:J) start_params[paste0("log_sigma_mu", j)] <- log(0.1)
        for (m in 1:M) start_params[paste0("log_sigma_delta", m)] <- log(0.1)

        # Noise
        for (j in 2:J) start_params[paste0("log_sigma_v_S", j)] <- log(0.05)
        for (m in 1:M) start_params[paste0("log_sigma_v_R", m)] <- log(0.05)
    }

    message("Starting optimization...")
    opt <- optim(
        par = start_params,
        fn = sfaKL_loglik,
        data = data_internal,
        J = J,
        M = M,
        method = method,
        control = control,
        hessian = TRUE,
        n_cores = n_cores
    )

    if (opt$convergence != 0) {
        warning("Optimization did not converge. Code: ", opt$convergence)
    } else {
        message("Optimization converged.")
    }

    # Transform log-sigmas back to sigmas for reporting
    est_params <- opt$par
    log_sigma_indices <- grep("log_sigma_", names(est_params))
    if (length(log_sigma_indices) > 0) {
        est_params[log_sigma_indices] <- exp(est_params[log_sigma_indices])
        names(est_params)[log_sigma_indices] <- gsub("log_", "", names(est_params)[log_sigma_indices])
    }

    result <- list(
        par = est_params,
        raw_par = opt$par,
        value = opt$value,
        counts = opt$counts,
        convergence = opt$convergence,
        hessian = opt$hessian,
        data = data_internal,
        J = J,
        M = M
    )

    class(result) <- "sfaKL"
    return(result)
}

#' Print method for sfaKL
#' @export
print.sfaKL <- function(x, ...) {
    cat("sfaKL Model Estimation\n")
    cat("Log-Likelihood:", -x$value, "\n")
    cat("Convergence:", x$convergence, "\n\n")
    cat("Coefficients:\n")
    print(x$par)
}
