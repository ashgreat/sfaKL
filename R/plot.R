#' Plot efficiency scores
#'
#' @param object An object of class sfaKL
#' @param type Type of plot ("hist" or "density")
#' @return A ggplot object
#' @import ggplot2
#' @importFrom stats reshape
#' @export
plot_efficiencies <- function(object, type = "hist") {
    preds <- predict(object)

    # Convert to efficiency scores (inputs: exp(mu) <= 1, outputs: exp(-delta) <= 1)
    eff <- preds
    input_cols <- grep("^mu", names(preds), value = TRUE)
    output_cols <- grep("^delta", names(preds), value = TRUE)

    if (length(input_cols) > 0) eff[input_cols] <- lapply(eff[input_cols], exp)
    if (length(output_cols) > 0) eff[output_cols] <- lapply(eff[output_cols], function(x) exp(-x))

    # Convert to long format
    preds_long <- stack(eff)
    names(preds_long) <- c("Efficiency", "Variable")

    p <- ggplot(preds_long, aes(x = Efficiency, fill = Variable)) +
        theme_minimal() +
        labs(title = "Distribution of Efficiency Scores", x = "Efficiency", y = "Count/Density")

    if (type == "hist") {
        p <- p + geom_histogram(bins = 30, alpha = 0.7, position = "identity")
    } else if (type == "density") {
        p <- p + geom_density(alpha = 0.5)
    }

    p <- p + facet_wrap(~Variable, scales = "free")

    return(p)
}
