test_that("parameters are roughly recovered on correlated simulation", {
    skip_on_cran()

    J <- 2
    M <- 2
    Sigma_mu <- matrix(c(0.04, 0.02, 0.02, 0.05), 2, 2)
    Sigma_delta <- matrix(c(0.03, 0.01, 0.01, 0.04), 2, 2)
    Omega <- matrix(
        c(0.04, 0.01, 0.00,
          0.01, 0.05, 0.01,
          0.00, 0.01, 0.06),
        3, 3,
        byrow = TRUE
    )

    sim <- sfaKL_simulate(n = 80, J = J, M = M, seed = 7, Sigma_mu = Sigma_mu, Sigma_delta = Sigma_delta, Omega = Omega)
    fit <- suppressWarnings(sfaKL_estimate(sim$data, control = list(maxit = 180, reltol = 1e-6), n_cores = 1, verbose = FALSE, use_gradient = FALSE))

    # Check that estimated inefficiency covariance diagonals are in the right ballpark
    est <- fit$par
    sigma_mu_hat <- matrix(0, J, J)
    sigma_delta_hat <- matrix(0, M, M)

    for (i in 1:J) sigma_mu_hat[i, i] <- est[paste0("chol_Sigma_mu_", i, "_", i)]^2
    for (i in 1:M) sigma_delta_hat[i, i] <- est[paste0("chol_Sigma_delta_", i, "_", i)]^2

    expect_lt(abs(sigma_mu_hat[1, 1] - Sigma_mu[1, 1]), 0.05)
    expect_lt(abs(sigma_delta_hat[1, 1] - Sigma_delta[1, 1]), 0.05)

    preds <- predict(fit, n_cores = 1)
    expect_equal(nrow(preds), nrow(sim$data))
})
