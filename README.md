# sfaKL: Multi-Output Multi-Input Stochastic Frontier Analysis

**Author:** Ashwin Malshe

`sfaKL` implements the multi-output multi-input stochastic frontier system with input- and output-specific inefficiency proposed by Kumbhakar and Lai (2021).

## Key Features

- **Multi-Output, Multi-Input**: Supports arbitrary numbers of inputs (J) and outputs (M)
- **Correlated Inefficiencies**: Models correlation between input-specific and output-specific inefficiencies
- **Two Estimation Methods**:
  - Maximum Likelihood (R-based, fast for small samples)
  - **Bayesian MCMC via Stan** (recommended for robust inference)

## Installation

### Basic Installation

```r
# install.packages("remotes")
remotes::install_github("ashgreat/sfaKL")
```

### For Bayesian Estimation (Recommended)

The package includes a Stan model for Bayesian estimation, which provides better parameter recovery and uncertainty quantification.

```r
# Install cmdstanr
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# Install CmdStan (one-time setup)
cmdstanr::install_cmdstan()
```

## Quick Start

### Bayesian Estimation with Stan (Recommended)

```r
library(cmdstanr)
library(mvtnorm)

# Source package functions
source("sfaKL/R/utils.R")
source("sfaKL/R/simulate_data.R")
source("sfaKL/R/stan_interface.R")

# Simulate data
set.seed(123)
Sigma_mu <- matrix(c(0.04, 0.02, 0.02, 0.05), 2, 2)
Sigma_delta <- matrix(c(0.03, 0.01, 0.01, 0.04), 2, 2)
Omega <- diag(c(0.04, 0.05, 0.06))

sim <- sfaKL_simulate(n = 200, J = 2, M = 2, seed = 123,
                      Sigma_mu = Sigma_mu, Sigma_delta = Sigma_delta, Omega = Omega)

# Compile and run Stan model
model <- cmdstan_model("sfaKL/stan/sfa_kl.stan")
stan_data <- prepare_stan_data(sim$data, J = 2, M = 2)

fit <- model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 5000,
    adapt_delta = 0.9
)

# View results
fit$summary(variables = c("alpha", "beta", "L_Sigma_mu", "L_Sigma_delta"))
```

### Maximum Likelihood Estimation

```r
library(sfaKL)

sim <- sfaKL_simulate(n = 500, J = 2, M = 2)

result <- sfaKL_estimate(
  sim$data,
  share_input = "S2",
  share_output = c("R1", "R2"),
  method = "BFGS",
  restarts = 2
)

print(result)
```

## Parameter Recovery Example

The Bayesian Stan model successfully recovers the true parameters from simulated data. Here are results from a simulation with N=200 observations:

### True vs Estimated Parameters

| Parameter | True Value | Posterior Mean | 95% CI | R-hat |
|-----------|------------|----------------|--------|-------|
| L_Σ_μ[1,1] | 0.20 | 0.17 | [0.04, 0.32] | 1.02 ✓ |
| L_Σ_μ[2,2] | 0.22 | 0.24 | [0.12, 0.38] | 1.01 ✓ |
| L_Σ_δ[1,1] | 0.17 | 0.19 | [0.06, 0.33] | 1.01 ✓ |
| L_Σ_δ[2,2] | 0.20 | 0.25 | [0.12, 0.39] | 1.02 ✓ |

The estimates are close to the true values, and all R-hat values are below 1.05, indicating good convergence.

### Technical Efficiency Estimates

| Component | Mean TE |
|-----------|---------|
| Input 1 (labor) | 88.1% |
| Input 2 (materials) | 47.9% |
| Output 1 | 86.9% |
| Output 2 | 63.6% |

## Model Specification

### Distributional Assumptions (Kumbhakar & Lai, 2021)

1. **[A1]** Input-specific inefficiency: μ ~ N⁻(0, Σ_μ), where μ ≤ 0 (multivariate half-normal)
2. **[A2]** Output-specific inefficiency: δ ~ N⁺(0, Σ_δ), where δ ≥ 0 (multivariate half-normal)
3. **[A3]** Symmetric noise: v ~ N(0, Ω)
4. **[A4]** Independence between μ, δ, and v

Under these assumptions, the composite error follows a **Closed Skew Normal (CSN)** distribution.

### Model Dimensions

- J = number of inputs
- M = number of outputs
- The package automatically infers dimensions from the input data

## Package Structure

```
sfaKL/
├── R/
│   ├── simulate_data.R   # Data simulation
│   ├── estimate.R        # ML estimation
│   ├── likelihood.R      # Likelihood function
│   ├── utils.R           # Utilities
│   └── stan_interface.R  # Stan wrapper functions
├── stan/
│   └── sfa_kl.stan       # Bayesian model
└── inst/extdata/
    └── example_results.rds  # Pre-computed example
```

## Why Bayesian Estimation?

The ML approach can suffer from:
1. **Variance collapse**: Inefficiency variances may collapse to zero
2. **Identification issues**: Flat likelihood surfaces make optimization difficult
3. **No uncertainty quantification**: Standard errors from the Hessian may be unreliable

The Stan model addresses these issues through:
- **Informative priors** that prevent parameter collapse
- **NUTS sampler** that explores complex posterior geometries
- **Full posterior distributions** for proper uncertainty quantification

## Arguments Reference

### sfaKL_estimate()

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Data frame with share and price variables | Required |
| `share_input` | Column name(s) for input shares | `"S2"` |
| `share_output` | Column names for output shares | `c("R1", "R2")` |
| `price_input_ratio` | Column for log input price ratio | `"ln_w2_w1"` |
| `price_output_ratios` | Columns for log output price ratios | `c("ln_p1_w1", "ln_p2_w1")` |
| `method` | Optimization method | `"Nelder-Mead"` |
| `n_cores` | Parallel cores for likelihood | `1` |
| `restarts` | Number of random restarts | `1` |

### prepare_stan_data()

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Data frame from sfaKL_simulate() | Required |
| `J` | Number of inputs | `2` |
| `M` | Number of outputs | `2` |

## References

Kumbhakar, S. C., & Lai, H. P. (2021). A multi-output multi-input stochastic frontier system with input- and output-specific inefficiency. *Economics Letters*, 201, 109807. https://doi.org/10.1016/j.econlet.2021.109807

## License

MIT
