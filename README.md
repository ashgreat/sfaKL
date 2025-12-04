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
