test_that("estimation, prediction, and plotting work with correlated structure", {
    sim <- sfaKL_simulate(n = 25, J = 2, M = 2, seed = 99)
    result <- suppressWarnings(sfaKL_estimate(sim$data, control = list(maxit = 200), n_cores = 1, verbose = FALSE, use_gradient = FALSE))

    expect_s3_class(result, "sfaKL")
    preds <- predict(result, n_cores = 1)
    expect_equal(nrow(preds), nrow(sim$data))
    expect_true(all(grepl("^mu|^delta", names(preds))))

    p <- plot_efficiencies(result, type = "density")
    expect_s3_class(p, "ggplot")

    diag_res <- sfaKL_check(result)
    expect_true(diag_res$ok)
    expect_true(diag_res$Theta_pd)
})

test_that("predict honors original column names via mapping", {
    sim <- sfaKL_simulate(n = 15, J = 2, M = 2, seed = 123)
    renamed <- sim$data
    names(renamed) <- sub("S2", "share_in", names(renamed))
    names(renamed) <- sub("R1", "rev1", names(renamed))
    names(renamed) <- sub("R2", "rev2", names(renamed))
    names(renamed) <- sub("ln_w2_w1", "lw", names(renamed))
    names(renamed) <- sub("ln_p1_w1", "lp1", names(renamed))
    names(renamed) <- sub("ln_p2_w1", "lp2", names(renamed))

    fit <- suppressWarnings(sfaKL_estimate(
        renamed,
        share_input = "share_in",
        share_output = c("rev1", "rev2"),
        price_input_ratio = "lw",
        price_output_ratios = c("lp1", "lp2"),
        control = list(maxit = 30),
        n_cores = 1,
        verbose = FALSE,
        use_gradient = FALSE
    ))

    preds <- predict(fit, newdata = renamed, n_cores = 1)
    expect_equal(nrow(preds), nrow(renamed))
})
