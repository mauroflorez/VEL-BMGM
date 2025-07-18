## VEL.BMGM <img src="man/figures/logo.png" align="right" width="100"/>

**Varying-Effects Logistic Regression & Bayesian Mixed Graphical Structure**

This R package implements Bayesian varying-effects logistic regression for binary outcomes, with mixed-type predictors (continuous, discrete, zero-inflated, and categorical). Predictor effects can vary nonlinearly with external covariates, modeled via Gaussian Process priors. The package also estimates a sparse conditional independence graph among predictors, enabling structure discovery and variable selection.

## Key Features

-   **Flexible logistic regression**: Models complex predictor effects using Gaussian Processes.
-   **Mixed data support**: Handles continuous, discrete, zero-inflated, and categorical predictors.
-   **Spike-and-slab priors**: For variable and covariate selection, and for graphical model edges.
-   **Efficient MCMC**: Uses Polya-Gamma augmentation for fast posterior sampling.
-   **Graph estimation**: Jointly estimates a sparse conditional independence graph among predictors.

## Installation

Clone/download this repository, then install with:

``` r
devtools::install("path/to/VEL.BMGM")
```

## Basic Usage

library(VEL.BMGM) fit \<- bmgm_gp(X, Y, Z, type_y = 'b_new', type = type_vector, ...)

## Authors

Mauro Florez

This package is under active development. Please see the paper for full model and algorithm details.
