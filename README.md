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
library(ggplot2)

# 1. Simulate Data (with correlated inefficiency/noise)
set.seed(123)
Sigma_mu <- matrix(c(0.04, 0.02, 0.02, 0.05), 2, 2)
Sigma_delta <- matrix(c(0.03, 0.01, 0.01, 0.04), 2, 2)
Omega <- matrix(c(
  0.04, 0.01, 0.00,
  0.01, 0.05, 0.01,
  0.00, 0.01, 0.06
), 3, 3, byrow = TRUE)

sim <- sfaKL_simulate(n = 500, J = 2, M = 2, Sigma_mu = Sigma_mu, Sigma_delta = Sigma_delta, Omega = Omega)
data <- sim$data
true_params <- sim$true_params

# 2. Estimate Model (custom columns supported)
result <- sfaKL_estimate(
  data,
  share_input = "S2",
  share_output = c("R1", "R2"),
  price_input_ratio = "ln_w2_w1",
  price_output_ratios = c("ln_p1_w1", "ln_p2_w1"),
  method = "BFGS",
  control = list(maxit = 8000, reltol = 1e-8),
  n_cores = 1,
  restarts = 2,
  verbose = TRUE,
  # Optional: box constraints / scaling
  lower = NULL,
  upper = NULL,
  parscale = NULL
)

print(result)

# Optional diagnostics
sfaKL_check(result)
```

### Arguments

-   `data`: A data frame containing the variables.
-   `share_input`: Name of the input share column (default `"S2"`).
-   `share_output`: Vector of names for the output share columns (default `c("R1", "R2")`).
-   `price_input_ratio`: Name of the log input price ratio column $\ln(w_2/w_1)$ (default `"ln_w2_w1"`).
-   `price_output_ratios`: Vector of names for the log output-to-input price ratio columns $\ln(p_m/w_1)$ (default `c("ln_p1_w1", "ln_p2_w1")`).
-   `n_cores`: Number of cores to use for parallel execution (default `1`).
-   `start_params`: Optional named vector of starting parameters. If `NULL`, OLS estimates are used.
-   `method`: Optimization method passed to `optim` (default `"Nelder-Mead"`).
-   `control`: List of control parameters passed to `optim`.

### Model Dimensions

The `sfaKL` package now supports **arbitrary numbers of inputs ($J$) and outputs ($M$)**. 
The dimensions are automatically inferred from the lengths of the `share_input` and `share_output` arguments.

-   $J = \text{length(share\_input)} + 1$
-   $M = \text{length(share\_output)}$

Ensure that `price_input_ratio` has length $J-1$ and `price_output_ratios` has length $M$.

### Distributional Assumptions

The model relies on the following distributional assumptions (Kumbhakar & Lai, 2021):

1. **[A1]** The vector of input-specific inefficiency $\mu \sim N^-(0_J, \Sigma_\mu)$, i.e., $\mu \le 0$ is multivariate half-normal.
2. **[A2]** The vector of output-specific inefficiency $\delta \sim N^+(0_M, \Sigma_\delta)$, i.e., $\delta \ge 0$ is multivariate half-normal.
3. **[A3]** The symmetric random error vector $v \sim N(0_{J+M-1}, \Omega)$.
4. **[A4]** The inefficiency components $\mu$ and $\delta$ are independent of each other and of the noise $v$.

Under these assumptions, the composite error follows a **closed skew normal (CSN)** distribution, enabling closed-form likelihood evaluation.

```r
# Compare estimated parameters with true parameters
print("True Parameters:")
print(unlist(true_params))

print("Estimated Parameters:")
print(result$par)

# 3. Predict Inefficiencies (works on new data with original column names)
ineff <- predict(result, n_cores = 1)
head(ineff)

# 4. Visualize Results (efficiencies on the correct scale)
plot_efficiencies(result, type = "density")

# 5. Using custom column names
colnames(data) <- c("share_in", "rev1", "rev2", "lw", "lp1", "lp2")
fit_custom <- sfaKL_estimate(
  data,
  share_input = "share_in",
  share_output = c("rev1", "rev2"),
  price_input_ratio = "lw",
  price_output_ratios = c("lp1", "lp2")
)

predict(fit_custom, newdata = data)
```

## References

Kumbhakar, S. C., & Lai, H. P. (2021). A multi-output multi-input stochastic frontier system with input- and output-specific inefficiency. *Economics Letters*, 201, 109807.
