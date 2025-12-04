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
# Note: Optimization can be sensitive to starting values.
result <- sfaKL_estimate(data, method = "Nelder-Mead", control = list(maxit = 10000))

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
