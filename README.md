# sfaKL: Multi-Output Multi-Input Stochastic Frontier Analysis

**Author:** Ashwin Malshe

`sfaKL` implements the multi-output multi-input stochastic frontier system with input- and output-specific inefficiency proposed by Kumbhakar and Lai (2021).

## Installation

You can install the development version of sfaKL from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("ashgreat/sfaKL")
```

## Example

This example demonstrates how to simulate data, estimate the model, and recover parameters.

```r
library(sfaKL)

# 1. Simulate Data
set.seed(123)
sim <- sfaKL_simulate(n = 1000)
data <- sim$data
true_params <- sim$true_params

head(data)

# 2. Estimate Model
# The estimation uses the share equations derived from the translog profit function.
# You can specify custom column names if your data doesn't match the defaults.
result <- sfaKL_estimate(
  data, 
  share_input = "S2", 
  share_output = c("R1", "R2"), 
  price_input_ratio = "ln_w2_w1", 
  price_output_ratios = c("ln_p1_w1", "ln_p2_w1"),
  n_cores = 1
)

print(result)

### Arguments

*   `data`: A data frame containing the variables.
*   `share_input`: Name of the input share column (default `"S2"`).
*   `share_output`: Vector of names for the output share columns (default `c("R1", "R2")`).
*   `price_input_ratio`: Name of the log input price ratio column $\ln(w_2/w_1)$ (default `"ln_w2_w1"`).
*   `price_output_ratios`: Vector of names for the log output-to-input price ratio columns $\ln(p_m/w_1)$ (default `c("ln_p1_w1", "ln_p2_w1")`).
*   `n_cores`: Number of cores to use for parallel execution (default `1`).
*   `start_params`: Optional named vector of starting parameters. If `NULL`, OLS estimates are used.
*   `method`: Optimization method passed to `optim` (default `"Nelder-Mead"`).
*   `control`: List of control parameters passed to `optim`.

### Model Dimensions

The `sfaKL` package now supports **arbitrary numbers of inputs ($J$) and outputs ($M$)**. 
The dimensions are automatically inferred from the lengths of the `share_input` and `share_output` arguments.

*   $J = \text{length(share\_input)} + 1$
*   $M = \text{length(share\_output)}$

Ensure that `price_input_ratio` has length $J-1$ and `price_output_ratios` has length $M$.

# Compare estimated parameters with true parameters
print("True Parameters:")
print(unlist(true_params))

print("Estimated Parameters:")
print(result$par)

# 3. Predict Inefficiencies
ineff <- predict(result)
head(ineff)

# 4. Visualize Results
library(ggplot2)
plot_efficiencies(result, type = "density")
```

## References

Kumbhakar, S. C., & Lai, H. P. (2021). A multi-output multi-input stochastic frontier system with input- and output-specific inefficiency. *Economics Letters*, 201, 109807.
