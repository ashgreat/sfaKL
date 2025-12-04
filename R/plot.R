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

    # Convert to long format
    preds_long <- stack(preds)
    names(preds_long) <- c("Efficiency", "Variable")

    p <- ggplot(preds_long, aes(x = Efficiency, fill = Variable)) +
        theme_minimal() +
        labs(title = "Distribution of Efficiency Scores", x = "Log Efficiency", y = "Count/Density")

    if (type == "hist") {
        p <- p + geom_histogram(bins = 30, alpha = 0.7, position = "identity")
    } else if (type == "density") {
        p <- p + geom_density(alpha = 0.5)
    }

    p <- p + facet_wrap(~Variable, scales = "free")

    return(p)
}
