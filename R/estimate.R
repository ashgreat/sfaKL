#' Estimate sfaKL model
#'
#' @param data Data frame containing the variables
#' @param start_params Optional named vector of starting parameters
#' @param method Optimization method (default "Nelder-Mead")
#' @param control List of control parameters for optim
#' @param n_cores Number of cores for parallel execution (default 1)
#' @return An object of class sfaKL
#' @export
#' @param share_input Name of the input share column (default "S2")
#' @param share_output Names of the output share columns (default c("R1", "R2"))
#' @param price_input_ratio Name of the log input price ratio column (default "ln_w2_w1")
#' @param price_output_ratios Names of the log output-to-input price ratio columns (default c("ln_p1_w1", "ln_p2_w1"))
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
    # Map user columns to internal names
    cols_map <- c(
        S2 = share_input,
        R1 = share_output[1],
        R2 = share_output[2],
        ln_w2_w1 = price_input_ratio,
        ln_p1_w1 = price_output_ratios[1],
        ln_p2_w1 = price_output_ratios[2]
    )

    # Check if columns exist
    if (!all(cols_map %in% names(data))) {
        missing <- cols_map[!cols_map %in% names(data)]
        stop("Missing columns in data: ", paste(missing, collapse = ", "))
    }

    # Create internal data frame with standard names
    data_internal <- data
    for (std_name in names(cols_map)) {
        data_internal[[std_name]] <- data[[cols_map[[std_name]]]]
    }

    # Default starting parameters if not provided
    if (is.null(start_params)) {
        # Simple OLS to get rough starting values for alphas and betas
        # S2 ~ ln_w2_w1 + ln_p1_w1 + ln_p2_w1
        m_S2 <- lm(-S2 ~ ln_w2_w1 + ln_p1_w1 + ln_p2_w1, data = data_internal)
        m_R1 <- lm(R1 ~ ln_p1_w1 + ln_p2_w1 + ln_w2_w1, data = data_internal)
        m_R2 <- lm(R2 ~ ln_p1_w1 + ln_p2_w1 + ln_w2_w1, data = data_internal)

        coef_S2 <- coef(m_S2)
        coef_R1 <- coef(m_R1)
        coef_R2 <- coef(m_R2)

        start_params <- c(
            alpha2 = as.numeric(coef_S2["(Intercept)"]),
            alpha22 = as.numeric(coef_S2["ln_w2_w1"]),
            beta1 = as.numeric(coef_R1["(Intercept)"]),
            beta2 = as.numeric(coef_R2["(Intercept)"]),
            beta11 = as.numeric(coef_R1["ln_p1_w1"]),
            beta22 = as.numeric(coef_R2["ln_p2_w1"]),
            beta12 = as.numeric((coef_R1["ln_p2_w1"] + coef_R2["ln_p1_w1"]) / 2), # Average for symmetry

            gamma21 = as.numeric(coef_S2["ln_p1_w1"]),
            gamma22 = as.numeric(coef_S2["ln_p2_w1"]),

            # Log variances (start small/moderate)
            log_sigma_mu1 = log(0.1),
            log_sigma_mu2 = log(0.1),
            log_sigma_delta1 = log(0.1),
            log_sigma_delta2 = log(0.1),
            log_sigma_v12 = log(sd(residuals(m_S2))),
            log_sigma_v21 = log(sd(residuals(m_R1))),
            log_sigma_v22 = log(sd(residuals(m_R2)))
        )
    }

    message("Starting optimization...")
    opt <- optim(
        par = start_params,
        fn = sfaKL_loglik,
        data = data_internal,
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
    est_params[grep("log_", names(est_params))] <- exp(est_params[grep("log_", names(est_params))])
    names(est_params) <- gsub("log_", "", names(est_params))

    result <- list(
        par = est_params,
        raw_par = opt$par,
        value = opt$value,
        counts = opt$counts,
        convergence = opt$convergence,
        hessian = opt$hessian,
        data = data
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
