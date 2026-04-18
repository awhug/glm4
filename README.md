# glm4

[![CRAN status](https://www.r-pkg.org/badges/version/glm4)](https://CRAN.R-project.org/package=glm4)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/glm4)](https://CRAN.R-project.org/package=glm4)
[![License: GPL (>= 2)](https://img.shields.io/badge/License-GPL%20(%3E%3D%202)-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.html)

## Overview

`glm4` is an R package to help fit and evaluate generalised linear models (GLMs) on large, high-dimensional datasets where standard approaches run out of memory.

This package adds functionality to `MatrixModels::glm4()`, which does the core work of fiting GLMs using iteratively reweighted least squares (IRLS) with support for sparse matrix operations via the `Matrix` package. When a model includes many factor levels (for example, thousands of fixed effects) the usual dense format model matrix can require tens or hundreds of gigabytes of RAM. Sparse storage keeps only the non-zero elements, making fitting and working with these models more feasible.

The raw output of `MatrixModels::glm4()` however is a bare-bones S4 object with no print, summary, or inference methods familiar to R users of the standard `stats` package. `glm4` fills these gaps by wrapping the fit in a familiar S3 object with the same general structure of a `stats::glm` fit, and adding methods and functions including (but not limited to):

- `print()` and `summary()` - to examine coefficients, standard errors, and test statistics
- `logLik()`, `deviance()`, `AIC()` - to extract information criteria components
- `confint()` - to calculate Wald confidence intervals on model parameters
- `vcov()` - to return a variance-covariance matrix (returned as a `Matrix`)
- `predict()`, `residuals()` - to inspect various forms of model predictions and prediction error
- `anova()` - for sequential deviance analysis and multi-model comparison
- `hatvalues()`, `influence()`, `rstandard()`, `cooks.distance()` - for model diagnostics based on leverage and influence measures

The interface is intentionally close to `stats::glm()` and associated methods and functions. Switching from `glm()` to `glm4()` should give virtually identical results, and fitting very sparse models requires only setting `sparse = TRUE`.

## Installation

Install the released version from CRAN:

```r
install.packages("glm4")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("awhug/glm4")
```

## Contributing

Contributions are more than welcome. If you find a bug or have a feature request, please [open an issue](https://github.com/awhug/glm4/issues). If you would like to contribute code, feel free to submit a pull request.
