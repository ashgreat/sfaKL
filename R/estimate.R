#' Estimate sfaKL model
#'
#' @param data Data frame containing the variables
#' @param start_params Optional named vector of starting parameters
#' @param method Optimization method (default "Nelder-Mead")
#' @param control List of control parameters for optim
#' @return An object of class sfaKL
#' @export
sfaKL_estimate <- function(data, start_params = NULL, method = "Nelder-Mead", control = list(maxit = 5000)) {
    # Check if data has required columns
    required_cols <- c("S2", "R1", "R2", "ln_w2_w1", "ln_p1_w1", "ln_p2_w1")
    if (!all(required_cols %in% names(data))) {
        stop("Data must contain columns: ", paste(required_cols, collapse = ", "))
    }

    # Default starting parameters if not provided
    if (is.null(start_params)) {
        # Simple OLS to get rough starting values for alphas and betas
        # S2 ~ ln_w2_w1 + ln_p1_w1 + ln_p2_w1
        m_S2 <- lm(-S2 ~ ln_w2_w1 + ln_p1_w1 + ln_p2_w1, data = data)
        m_R1 <- lm(R1 ~ ln_p1_w1 + ln_p2_w1 + ln_w2_w1, data = data)
        m_R2 <- lm(R2 ~ ln_p1_w1 + ln_p2_w1 + ln_w2_w1, data = data)

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
        data = data,
        method = method,
        control = control,
        hessian = TRUE
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
