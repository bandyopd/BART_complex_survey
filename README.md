# complexsurvey

R package "complexsurvey" implementing Soft Bayesian Additive Regression Trees (SBART) 
for clustered complex survey data with non-Gaussian errors.

## Reference
Mandal A, Linero AR, Bandyopadhyay D, Sinha D (2026). 
Soft Bayesian Additive Regression Trees (SBART) for correlated survey 
response with non-Gaussian error, 
*Journal of Nonparametric Statistics*, 38(1), 235–253. 
DOI: https://doi.org/10.1080/10485252.2025.2574708

## Required Packages
```r
install.packages(c("rlang", "vctrs", "fs", "devtools", "remotes", "SoftBart", "invgamma", "zeallot", "truncnorm",
                   "sn", "glmnet", "extraDistr", "ald", "loo", "diversitree"),
                  repos = "https://cloud.r-project.org",
                  dependencies = TRUE)
```
## Installation
You need to have Rtools 45 pre-installed in your computer
```r
devtools::install_github("amandal-stat/complex_survey", force = TRUE)
```

## Models
- `SNSBARTW()` — SBART with Skew-Normal errors
- `TSBARTW()`  — SBART with Skew-t errors (heavy-tailed)
- `QRSBARTW()` — SBART Quantile Regression via Asymmetric Laplace Distribution



## Quick Example
```r
library(complexsurvey)
library(SoftBart)

# Simulate data from skew-normal model
dat <- sim_fried_SN(N = 100, P = 5, alpha = 5, sigma = 1, sigma1 = 0.2)

# Fit model
hypers <- Hypers(dat$X, dat$Y, num_tree = 50)
opts   <- Opts(num_burn = 500, num_save = 500)

fit <- SNSBARTW(X = dat$X, Y = dat$Y, test_X = dat$X,
                w = dat$w, ni = dat$ni,
                hypers = hypers, opts = opts)



## Contact
Abhishek Mandal — amandal2@fsu.edu  
Department of Statistics, Florida State University

