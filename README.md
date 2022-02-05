
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`SuperNOVA` <img src="man/figures/SuperNOVA_sticker.png" height="300" align="right"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/blind-contours/SuperNOVA/workflows/R-CMD-check/badge.svg)](https://github.com/blind-contours/SuperNOVA/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/blind-contours/SuperNOVA/master.svg)](https://codecov.io/github/blind-contours/SuperNOVA?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/SupernOVA)](https://www.r-pkg.org/pkg/SuperNOVA)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/SuperNOVA)](https://CRAN.R-project.org/package=SuperNOVA)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/SuperNOVA)](https://CRAN.R-project.org/package=SuperNOVA)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042) -->
<!-- [![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447) -->
<!-- badges: end -->

> Efficient Estimation of the Causal Effects of Non-Parametric
> Interactions and Effect Modifications using Stochastic Interventions
> **Authors:** [David McCoy](https://davidmccoy.org)

------------------------------------------------------------------------

## What’s `SuperNOVA`?

The `SuperNOVA` R package is designed to provide facilities to
data-adaptively identify variable sets that are most predictive of an
outcome of interest and construct efficient estimators for the
counterfactual mean of an outcome under stochastic interventions on
these variables. Stochastic interventions are shifts to the exposures
that depend on the naturally observed values (Dı́az and van der Laan
2012; Haneuse and Rotnitzky 2013). `SuperNOVA` builds off of the
`txshift` package which implements the targeted maximum likelihood (TML)
estimator of a stochastic shift causal parameter originally proposed by
Dı́az and van der Laan (2018). `SuperNOVA` extends the original
stochastic intervention methodology from one treatment/exposure variable
to include joint stochastic interventions on two variables which allows
for the construction of a non-parametric interaction parameter.
Likewise, `SuperNOVA` also estimates both individual stochastic
intervention outcomes under some delta shift compared to outcome under
no intervention and a target parameter for effect modification which is
the difference in mean outcomes under intervention compared to no
intervention across strata of an effect modifier.

Of course, it is not known a priori in most cases what variables are
interacting or modifying effects and therefore it is necessary to
identify these variable sets first. As such `SuperNOVA` uses V-fold
cross-validation framework to estimate a data-adaptive parameter in
training folds and a non-parametric interaction target parameter in
estimation folds. Our data-adaptive parameters are variable sets used in
basis functions in the best fitting multivariate adaptive regression
spline model. The best fitting model is determined using a Super Learner
which selects the model from an ensemble with the lowest cross-validated
MSE. Variable sets are considered important based on ANOVA-like variance
decompositions for the basis functions in the best fitting model.
Individual variables and variable sets used in all the training folds
are considered consistent predictors. The interaction target parameter
is applied to variable sets composed of two variables in the mixed
exposure. This target parameter is the expected outcome under a dual
shift of both variables by some delta compared to the sum of individual
shifts. Other parameters exist for effect modification and individual
variable shifts. Cross-validated targeted minimum loss-based estimation
(TMLE) is used to update the initial expected outcomes given stochastic
shift interventions. This method, called SuperNOVA, guarantees
consistency, efficiency, and multiple robustness. SuperNOVA provides
researchers with V-fold specific and pooled results for each target
parameter. Additional information is provided in the vignette.

`SuperNOVA` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) (Coyle, Hejazi, Malenica, et
al. 2022) to allow for ensemble machine learning to be leveraged in the
estimation procedure for each nuisance parameter and estimation of the
data-adaptive parameters in the iterative backfitting procedure used to
identify basis functions in the best fittin model. There are several
stacks of machine learning algorithms used that are constructed from
`sl3`.

------------------------------------------------------------------------

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=SuperNOVA) via

``` r
install.packages("SuperNOVA")
```

*Note:* If `SuperNOVA` is installed from
[CRAN](https://CRAN.R-project.org/package=SuperNOVA), the `sl3`, an
enhancing dependency that allows ensemble machine learning to be used
for nuisance parameter estimation, won’t be included. We highly
recommend additionally  
installing `sl3` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("tlverse/sl3@devel")
```

For the latest features, install the most recent *stable version* of
`SuperNOVA` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("blind-contours/SuperNOVA@main")
```

To contribute, install the *development version* of `SuperNOVA` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("blind-contours/SuperNOVA@devel")
```

------------------------------------------------------------------------

## Example

To illustrate how `SuperNOVA` may be used to ascertain the effect of a
mixed exposure, consider the following example:

``` r
library(SuperNOVA)
library(devtools)
#> Loading required package: usethis
load_all('~/sl3')
#> ℹ Loading sl3
set.seed(429153)
# simulate simple data
n_obs <- 400
```

The `simulate_data` function creates simulated data with a multivariate
exposure, covariates (confounders), and a continuous outcome.

``` r
data_info <- simulate_data(n_obs = n_obs)
data <- data_info$data
head(data)
#>          M1        M2       M3 W1       W2       W3        Y
#> 1 1.8118767 0.5594007 3.452442  0 5.805874 6.010026 4.365341
#> 2 1.3757667 3.1725122 3.442254  1 4.826966 7.257671 4.319712
#> 3 0.5968083 2.0818659 3.237153  1 6.010977 7.565084 5.182973
#> 4 2.6877639 1.6314774 3.201915  1 6.333759 7.441065 5.678622
#> 5 3.0809005 2.8610332 3.951894  0 3.643420 7.530927 3.714323
#> 6 0.8583468 0.8110409 3.091384  0 5.035758 7.665307 4.592644
```

``` r
## this is the stack of learners to learn the density of our exposures, or the g mechanism
sl_density_lrnr <- make_density_superlearner()

# this stack will be used to find basis functions in the mixture during an iterative backfitting procedure, details in vignette
Lrnr_earth_1 <- Lrnr_earth$new(linpreds = FALSE, degree = 1)
Lrnr_earth_2 <- Lrnr_earth$new(linpreds = FALSE, degree = 2)
Lrnr_earth_3 <- Lrnr_earth$new(linpreds = FALSE, degree = 2, pmethod = "none")

learners <- c(
  Lrnr_earth_1,
  Lrnr_earth_2,
  Lrnr_earth_3
)

names(learners) <- c(
  "full earth 1",
  "full earth 2",
  "full earth 3"
)

Exposures_stack <- make_learner(Stack, learners)

## this is the stack of learners used to estimate Y given covariates in the iterative backfitting procedure

Lrnr_glm_basic <- Lrnr_glm$new()
Lrnr_mean_base <- Lrnr_mean$new()

Lrnr_ridge <- Lrnr_glmnet$new(alpha = 0)
Lrnr_lasso <- Lrnr_glmnet$new(alpha = 1)

learners <- c(
  Lrnr_earth_1,
  Lrnr_earth_2,
  Lrnr_earth_3,
  Lrnr_glm_basic,
  Lrnr_mean_base,
  Lrnr_ridge,
  Lrnr_lasso
)

names(learners) <- c(
  "full earth 1",
  "full earth 2",
  "full earth 3",
  "Lrnr_glm",
  "Lrnr_mean",
  "Lrnr_ridge",
  "Lrnr_lasso"
)

Covariate_stack <- make_learner(Stack, learners)

## and now we make the estimators for our Q or outcome mechanism: 

mean_lrnr <- Lrnr_mean$new()
fglm_lrnr <- Lrnr_glm_fast$new()
rf_lrnr <- Lrnr_ranger$new()
lasso_learner <- Lrnr_glmnet$new(alpha = 1)
ridge_learner <- Lrnr_glmnet$new(alpha = 0)
lrn_polspline <- Lrnr_polspline$new()
lrn_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)
hal_lrnr <- Lrnr_hal9001$new(max_degree = 3, n_folds = 3)

Outcome_stack <- make_learner(
  Stack, mean_lrnr, fglm_lrnr, rf_lrnr, lasso_learner, ridge_learner, lrn_polspline, lrn_ranger100, Lrnr_earth_1, Lrnr_earth_2, Lrnr_earth_3
)
```

``` r
W <- data[, c("W2", "W3")]
A <- data[, c("M1", "M2", "M3")]
V <- data[, c("W1")]
Y <- data[, c("Y")]


sim_results <- SuperNOVA(W = W,
                         V = V,
                         A = A,
                         Y = Y,
                         delta = 1,
                         LOD_val = 0,
                         Density_stack = sl_density_lrnr,
                         Exposures_stack = Exposures_stack,
                         Covariate_stack = Covariate_stack,
                         Outcome_stack = Outcome_stack,
                         n_folds = 3,
                         family = "continuous",
                         quantile_thresh = 0) 

indiv_shift_results <- sim_results$`Indiv Shift Results`
em_results <- sim_results$`Effect Mod Results`
joint_shift_results <- sim_results$`Joint Shift Results`
```

Let’s first look at the results for individual stochastic shifts by
delta compared to no shift:

``` r
indiv_shift_results
#> # A tibble: 6 × 11
#>   Condition    Psi Variance      SE `Lower CI` `Upper CI` `P-value`  Fold Type  
#>   <chr>      <dbl>    <dbl>   <dbl>      <dbl>      <dbl>     <dbl> <int> <chr> 
#> 1 M1        0.147  1.34e- 3 3.66e-2    7.57e-2    2.19e-1  5.69e- 5     1 Indiv…
#> 2 V         0.0359 6.33e+ 8 2.52e+4   -4.93e+4    4.93e+4  1.00e+ 0     1 Indiv…
#> 3 M1        0.187  4.66e- 4 2.16e-2    1.45e-1    2.30e-1  4.08e-18     2 Indiv…
#> 4 V         0.0403 1.07e+18 1.03e+9   -2.03e+9    2.03e+9  1.00e+ 0     2 Indiv…
#> 5 M1        0.172  4.34e- 3 6.59e-2    4.24e-2    3.01e-1  9.22e- 3     3 Indiv…
#> 6 V         0.0647 1.23e+18 1.11e+9   -2.17e+9    2.17e+9  1.00e+ 0     3 Indiv…
#> # … with 2 more variables: Variables <chr>, N <int>
```

Next we can look at effect modifications:

``` r
em_results
#> [1] NA
```

And finally results for the joint shift

``` r
joint_shift_results
#> [1] NA
```

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/SuperNOVA/issues). Further
details on filing issues are provided in our [contribution
guidelines](https://github.com/nhejazi/SuperNOVA/blob/master/CONTRIBUTING.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/SuperNOVA/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `SuperNOVA` R package, please cite the following:

        @article{hejazi2020efficient,
          author = {Hejazi, Nima S and {van der Laan}, Mark J and Janes, Holly
            E and Gilbert, Peter B and Benkeser, David C},
          title = {Efficient nonparametric inference on the effects of
            stochastic interventions under two-phase sampling, with
            applications to vaccine efficacy trials},
          year = {2020},
          doi = {10.1111/biom.13375},
          url = {https://doi.org/10.1111/biom.13375},
          journal = {Biometrics},
          publisher = {Wiley Online Library}
        }

        @article{hejazi2020SuperNOVA-joss,
          author = {Hejazi, Nima S and Benkeser, David C},
          title = {{SuperNOVA}: Efficient estimation of the causal effects of
            stochastic interventions in {R}},
          year  = {2020},
          doi = {10.21105/joss.02447},
          url = {https://doi.org/10.21105/joss.02447},
          journal = {Journal of Open Source Software},
          publisher = {The Open Journal}
        }

        @software{hejazi2022SuperNOVA-rpkg,
          author = {Hejazi, Nima S and Benkeser, David C},
          title = {{SuperNOVA}: Efficient Estimation of the Causal Effects of
            Stochastic Interventions},
          year  = {2022},
          doi = {10.5281/zenodo.4070042},
          url = {https://CRAN.R-project.org/package=SuperNOVA},
          note = {R package version 0.3.7}
        }

------------------------------------------------------------------------

## Related

-   [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R
    package providing an independent implementation of the same core
    routines for the TML estimation procedure and statistical
    methodology as is made available here, through reliance on a unified
    interface for Targeted Learning provided by the
    [`tmle3`](https://github.com/tlverse/tmle3) engine of the [`tlverse`
    ecosystem](https://github.com/tlverse).

-   [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
    providing facilities to estimate the causal effect of stochastic
    treatment regimes in the mediation setting, including classical
    (IPW) and augmented double robust (one-step) estimators. This is an
    implementation of the methodology explored by Dı́az and
    Hejazi (2020).

-   [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
    package for estimating the conditional density treatment mechanism
    component of this parameter based on using the [highly adaptive
    lasso](https://github.com/tlverse/hal9001) (Coyle, Hejazi, Phillips,
    et al. 2022; Hejazi, Coyle, and van der Laan 2020) in combination
    with a pooled hazard regression. This package implements a variant
    of the approach advocated by Dı́az and van der Laan (2011).

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through grants
from the National Library of Medicine (award no. [T32
LM012417](https://reporter.nih.gov/project-details/9248418)) and the
National Institute of Allergy and Infectious Diseases (award no. [R01
AI074345](https://reporter.nih.gov/project-details/9926564)) of the
National Institutes of Health, as well as by the National Science
Foundation (award no. [DMS
2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

------------------------------------------------------------------------

## License

© 2017-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2017-2022 Nima S. Hejazi
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-coyle-sl3-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2022. *<span class="nocase">sl3</span>: Modern Machine
Learning Pipelines for Super Learning*.
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-coyle-hal9001-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars W van der Laan,
and Mark J van der Laan. 2022. *<span class="nocase">hal9001</span>: The
Scalable Highly Adaptive Lasso*.
<https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal" class="csl-entry">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super" class="csl-entry">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-diaz2012population" class="csl-entry">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2): 541–49.

</div>

<div id="ref-diaz2018stochastic" class="csl-entry">

———. 2018. “Stochastic Treatment Regimes.” In *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*, 167–80.
Springer Science & Business Media.

</div>

<div id="ref-haneuse2013estimation" class="csl-entry">

Haneuse, Sebastian, and Andrea Rotnitzky. 2013. “Estimation of the
Effect of Interventions That Modify the Received Treatment.” *Statistics
in Medicine* 32 (30): 5260–77.

</div>

<div id="ref-hejazi2020hal9001-joss" class="csl-entry">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020. “<span
class="nocase">hal9001</span>: Scalable Highly Adaptive Lasso Regression
in R.” *Journal of Open Source Software* 5 (53): 2526.
<https://doi.org/10.21105/joss.02526>.

</div>

</div>
