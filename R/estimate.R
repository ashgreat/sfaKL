#' Estimate sfaKL model
#'
#' @param data Data frame containing the variables
#' @param share_input Name of the input share column (default "S2") or vector of names
#' @param share_output Names of the output share columns (default c("R1", "R2"))
#' @param price_input_ratio Name of the log input price ratio column (default "ln_w2_w1") or vector
#' @param price_output_ratios Names of the log output-to-input price ratio columns (default c("ln_p1_w1", "ln_p2_w1"))
#' @param start_params Optional named vector of starting parameters
#' @param optimizer Which optimizer to use ("optim" or "nlminb")
#' @param method Optimization method (default "L-BFGS-B")
#' @param control List of control parameters for optim
#' @param n_cores Number of cores for parallel execution (default 1)
#' @param verbose Logical; if TRUE prints optimizer progress/diagnostics
#' @param restarts Number of random restarts (>=1); best solution is returned
#' @param use_gradient Logical; if TRUE uses finite-difference gradients with a gradient-capable optimizer
#' @param cov_restart_scale Standard deviation for restart jitter on covariance parameters (Cholesky entries)
#' @param lower Optional lower bounds (use with box methods such as L-BFGS-B)
#' @param upper Optional upper bounds (use with box methods such as L-BFGS-B)
#' @param parscale Optional parscale vector passed to optim for scaling
#' @return An object of class sfaKL
#' @export
sfaKL_estimate <- function(data,
                           share_input = "S2",
                           share_output = c("R1", "R2"),
                           price_input_ratio = "ln_w2_w1",
                           price_output_ratios = c("ln_p1_w1", "ln_p2_w1"),
                           start_params = NULL,
                           optimizer = c("optim", "nlminb"),
                           method = "L-BFGS-B",
                           control = list(maxit = 10000, reltol = 1e-8),
                           n_cores = 1,
                           verbose = FALSE,
                           restarts = 1,
                           use_gradient = TRUE,
                           cov_restart_scale = 0.15,
                           lower = NULL,
                           upper = NULL,
                           parscale = NULL) {
    optimizer <- match.arg(optimizer)
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
        n_input_shares <- J - 1

        # Initialize containers
        alpha_init <- rep(0, n_input_shares)
        beta_init <- rep(0, M)
        A_init <- matrix(0, n_input_shares, n_input_shares)
        B_init <- matrix(0, M, M)
        Gamma_from_inputs <- matrix(0, n_input_shares, M)
        Gamma_from_outputs <- matrix(0, n_input_shares, M)

        ln_w_names <- paste0("ln_w", 2:J, "_w1")
        ln_p_names <- paste0("ln_p", 1:M, "_w1")

        # OLS for input share equations (-S_j as dependent variable)
        for (j in 1:n_input_shares) {
            dep <- -data_internal[[paste0("S", j + 1)]]
            df_j <- data.frame(dep, data_internal[, c(ln_w_names, ln_p_names), drop = FALSE])
            fit <- try(lm(dep ~ ., data = df_j), silent = TRUE)
            if (!inherits(fit, "try-error")) {
                coefs <- coef(fit)
                alpha_init[j] <- coefs["(Intercept)"]
                A_init[j, ] <- coefs[ln_w_names]
                Gamma_from_inputs[j, ] <- coefs[ln_p_names]
            }
        }

        # OLS for output share equations (R_m)
        for (m_idx in 1:M) {
            dep <- data_internal[[paste0("R", m_idx)]]
            df_m <- data.frame(dep, data_internal[, c(ln_p_names, ln_w_names), drop = FALSE])
            fit <- try(lm(dep ~ ., data = df_m), silent = TRUE)
            if (!inherits(fit, "try-error")) {
                coefs <- coef(fit)
                beta_init[m_idx] <- coefs["(Intercept)"]
                B_init[m_idx, ] <- coefs[ln_p_names]
                Gamma_from_outputs[, m_idx] <- coefs[ln_w_names]
            }
        }

        # Symmetrize A and B, average Gamma from both systems
        A_init[is.na(A_init)] <- 0
        B_init[is.na(B_init)] <- 0
        Gamma_from_inputs[is.na(Gamma_from_inputs)] <- 0
        Gamma_from_outputs[is.na(Gamma_from_outputs)] <- 0

        A_init <- (A_init + t(A_init)) / 2
        B_init <- (B_init + t(B_init)) / 2

        Gamma_init <- Gamma_from_inputs
        # Average when both input and output regressions deliver information
        shared_gamma <- (Gamma_from_inputs + Gamma_from_outputs) / 2
        shared_gamma[is.na(shared_gamma)] <- 0
        Gamma_init[!is.na(shared_gamma)] <- shared_gamma[!is.na(shared_gamma)]

        # Populate parameter vector (note: constraints handled downstream)
        for (j in 2:J) start_params[paste0("alpha", j)] <- alpha_init[j - 1]

        if (n_input_shares > 1) {
            for (j in 1:(n_input_shares - 1)) {
                for (k in j:(n_input_shares - 1)) start_params[paste0("A_", j, "_", k)] <- A_init[j, k]
            }
        }

        for (m in 1:M) start_params[paste0("beta", m)] <- beta_init[m]

        if (M > 1) {
            for (m in 1:(M - 1)) {
                for (k in m:(M - 1)) start_params[paste0("B_", m, "_", k)] <- B_init[m, k]
            }
        }

        for (j in 1:(J - 1)) {
            for (m in 1:M) {
                start_params[paste0("Gamma_", j, "_", m)] <- Gamma_init[j, m]
            }
        }

        # Inefficiency covariance (diagonal start, log-scale on diag)
        # Use larger starting values to avoid collapsing toward zero
        for (j in 1:J) start_params[paste0("chol_Sigma_mu_", j, "_", j)] <- log(0.15)
        if (J > 1) {
            for (j in 2:J) {
                for (k in 1:(j - 1)) start_params[paste0("chol_Sigma_mu_", j, "_", k)] <- 0
            }
        }
        for (m in 1:M) start_params[paste0("chol_Sigma_delta_", m, "_", m)] <- log(0.15)
        if (M > 1) {
            for (m in 2:M) {
                for (k in 1:(m - 1)) start_params[paste0("chol_Sigma_delta_", m, "_", k)] <- 0
            }
        }

        # Noise covariance from residuals (captures correlation)
        residual_mat <- c()
        if (nrow(data_internal) > 1) {
            res_cols <- list()
            X_w <- as.matrix(data_internal[, ln_w_names, drop = FALSE])
            X_p <- as.matrix(data_internal[, ln_p_names, drop = FALSE])

            for (j in 1:n_input_shares) {
                res_cols[[paste0("S", j + 1)]] <- -data_internal[[paste0("S", j + 1)]] -
                    (alpha_init[j] + X_w %*% A_init[j, ] + X_p %*% Gamma_init[j, ])
            }
            for (m_idx in 1:M) {
                res_cols[[paste0("R", m_idx)]] <- data_internal[[paste0("R", m_idx)]] -
                    (beta_init[m_idx] + X_p %*% B_init[m_idx, ] + X_w %*% Gamma_init[, m_idx])
            }
            residual_mat <- as.matrix(as.data.frame(res_cols))
        }

        if (!is.null(residual_mat) && ncol(residual_mat) > 0) {
            cov_start <- try(cov(residual_mat), silent = TRUE)
        } else {
            cov_start <- diag(0.05, n_input_shares + M)
        }

        if (inherits(cov_start, "try-error")) cov_start <- diag(0.05, n_input_shares + M)
        chol_start <- try(chol(cov_start), silent = TRUE)
        if (inherits(chol_start, "try-error")) chol_start <- diag(sqrt(diag(cov_start)), n_input_shares + M)
        L <- t(chol_start)
        for (i in 1:nrow(L)) {
            for (k in 1:i) {
                nm <- paste0("chol_Omega_", i, "_", k)
                if (i == k) {
                    start_params[nm] <- log(pmax(abs(L[i, k]), 1e-3))
                } else {
                    start_params[nm] <- L[i, k]
                }
            }
        }
    }

    # Build default bounds (applied only if user has not supplied)
    n_par <- length(start_params)
    lower_bounds <- if (is.null(lower)) rep(-Inf, n_par) else lower
    upper_bounds <- if (is.null(upper)) rep(Inf, n_par) else upper
    if (length(lower_bounds) != n_par || length(upper_bounds) != n_par) stop("lower/upper must have length equal to number of parameters")

    # Apply sensible bounds to Cholesky diagonals (log-sd scale) to avoid collapse
    diag_idx <- which(grepl("^chol_.*_\\d+_\\d+$", names(start_params)))
    if (length(diag_idx) > 0) {
        for (idx in diag_idx) {
            parts <- strsplit(names(start_params)[idx], "_")[[1]]
            last_two <- tail(parts, 2)
            if (length(last_two) == 2 && last_two[1] == last_two[2]) {
                # Raise minimum to prevent variance collapse (SD >= 0.05)
                lower_bounds[idx] <- ifelse(is.infinite(lower_bounds[idx]), log(0.05), lower_bounds[idx])
                upper_bounds[idx] <- ifelse(is.infinite(upper_bounds[idx]), log(5), upper_bounds[idx])
            }
        }
    }

    # Incorporate optional scaling
    if (!is.null(parscale)) {
        if (optimizer == "optim") {
            control$parscale <- parscale
        } else if (optimizer == "nlminb") {
            control$scale <- parscale
        }
    }

    # Gradient function (finite differences)
    gr_fun <- NULL
    if (use_gradient && optimizer == "optim" && method %in% c("BFGS", "L-BFGS-B")) {
        gr_fun <- function(par, ...) {
            numDeriv::grad(
                func = function(p) sfaKL_loglik(p, data_internal, J, M, n_cores),
                x = par
            )
        }
    } else if (use_gradient && optimizer == "nlminb") {
        gr_fun <- function(par) {
            numDeriv::grad(
                func = function(p) sfaKL_loglik(p, data_internal, J, M, n_cores),
                x = par
            )
        }
    }

    message("Starting optimization...")
    do_optim <- function(par_start) {
        if (optimizer == "nlminb") {
            nlminb(
                start = par_start,
                objective = function(p) sfaKL_loglik(p, data_internal, J, M, n_cores),
                gradient = gr_fun,
                lower = lower_bounds,
                upper = upper_bounds,
                control = control
            )
        } else if (is.null(lower) && is.null(upper)) {
            optim(
                par = par_start,
                fn = sfaKL_loglik,
                data = data_internal,
                J = J,
                M = M,
                method = method,
                control = control,
                hessian = TRUE,
                n_cores = n_cores,
                gr = gr_fun
            )
        } else {
            lower_opt <- lower_bounds
            upper_opt <- upper_bounds

            optim(
                par = par_start,
                fn = sfaKL_loglik,
                data = data_internal,
                J = J,
                M = M,
                method = method,
                control = control,
                hessian = TRUE,
                n_cores = n_cores,
                lower = lower_opt,
                upper = upper_opt,
                gr = gr_fun
            )
        }
    }

    # First run
    opt <- do_optim(start_params)
    best_opt <- opt

    # Random restarts (jittered)
    if (restarts > 1) {
        chol_idx <- grep("^chol_", names(start_params))
        base_sd <- 0.05
        for (r in 2:restarts) {
            jitter <- rnorm(length(start_params), 0, base_sd)
            if (length(chol_idx) > 0) {
                jitter[chol_idx] <- rnorm(length(chol_idx), 0, cov_restart_scale)
            }
            names(jitter) <- names(start_params)
            par_restart <- start_params + jitter
            opt_r <- do_optim(par_restart)
            if (opt_r$value < best_opt$value) best_opt <- opt_r
        }
    }

    opt <- best_opt

    # Harmonize optimizer outputs
    opt_par <- if (!is.null(opt$par)) opt$par else opt$parameters
    opt_value <- if (!is.null(opt$value)) opt$value else opt$objective
    opt_conv <- if (!is.null(opt$convergence)) opt$convergence else if (!is.null(opt$message)) 0 else 1
    opt_counts <- if (!is.null(opt$counts)) opt$counts else if (!is.null(opt$evaluations)) opt$evaluations else NA
    opt_message <- if (!is.null(opt$message)) opt$message else ""
    opt_hessian <- if (!is.null(opt$hessian)) opt$hessian else NULL

    if (opt_conv != 0) {
        warning("Optimization did not converge. Code: ", opt_conv, " ", opt_message)
    } else if (verbose) {
        message("Optimization converged.")
    }

    # Transform Cholesky diagonals for readability (raw_par retains optimization scale)
    est_params <- opt_par
    chol_names <- names(est_params)[grepl("^chol_", names(est_params))]
    for (nm in chol_names) {
        idx <- as.integer(tail(strsplit(nm, "_")[[1]], 2))
        if (length(idx) == 2 && idx[1] == idx[2]) {
            est_params[nm] <- exp(est_params[nm])
        }
    }

    cond_numbers <- list(
        Theta = NA_real_,
        Delta = NA_real_
    )
    if (verbose) {
        mats_check <- try(get_matrices(opt_par, J, M), silent = TRUE)
        if (!inherits(mats_check, "try-error")) {
            cond_numbers$Theta <- kappa(mats_check$Theta)
            cond_numbers$Delta <- kappa(mats_check$Delta)
        }
    }

    result <- list(
        par = est_params,
        raw_par = opt_par,
        value = opt_value,
        counts = opt_counts,
        convergence = opt_conv,
        hessian = opt_hessian,
        data = data_internal,
        J = J,
        M = M,
        cols_map = cols_map,
        method = method,
        control = control,
        restarts = restarts,
        condition_numbers = cond_numbers
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
