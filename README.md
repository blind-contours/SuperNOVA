
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
for the construction of a non-parametric interaction parameter and
mediation. Likewise, `SuperNOVA` also estimates both individual
stochastic intervention outcomes under some delta shift compared to
outcome under no intervention and a target parameter for effect
modification which is the difference in mean outcomes under intervention
compared to no intervention across strata of an effect modifier which is
data-adaptively determined.

Of course, it is not known a priori in most cases what variables are
interacting or modifying effects (especially in the case of a mixed
exposure) and therefore it is necessary to identify these variable sets
first. As such `SuperNOVA` uses k-fold cross-validation framework to
estimate a data-adaptive parameter in training folds and a
non-parametric interaction target parameter in estimation folds. Our
data-adaptive parameters are variable sets used in basis functions in
the best fitting multivariate adaptive regression spline model or highly
adaptive lasso model. The best fitting model is determined using a Super
Learner which selects the model from an ensemble with the lowest
cross-validated MSE. Variable sets are considered important based on
ANOVA-like variance decompositions for the basis functions in the best
fitting model. Individual variables and variable sets used in all the
training folds are considered consistent predictors. The interaction
target parameter is applied to variable sets composed of two variables
in the mixed exposure. This target parameter is the expected outcome
under a dual shift of both variables by some delta compared to the sum
of individual shifts. Other parameters exist for effect modification and
individual variable shifts. Cross-validated targeted minimum loss-based
estimation (TMLE) is used to update the initial expected outcomes given
stochastic shift interventions. This method, called SuperNOVA,
guarantees consistency, efficiency, and multiple robustness. SuperNOVA
provides researchers with k-fold specific and pooled results for each
target parameter. Additional information is provided in the vignette.

`SuperNOVA` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) (Coyle et al. 2022) to allow
for ensemble machine learning to be leveraged in the estimation
procedure for each nuisance parameter and estimation of the
data-adaptive parameters. There are several stacks of machine learning
algorithms used that are constructed from `sl3` automatically. If the
stack parameters are NULL, SuperNOVA automatically builds ensembles of
machine learning algorithms that are flexible yet not overly
computationally taxing.

------------------------------------------------------------------------

## Installation

*Note:* Because `SuperNOVA` package (currently) depends on `sl3` that
allows ensemble machine learning to be used for nuisance parameter
estimation and `sl3` is not on CRAN the `SuperNOVA` package is not
available on CRAN and must be downloaded here.

``` r
remotes::install_github("blind-contours/SuperNOVA@main")
```

`SuperNOVA` uses `sl3` so please download sl3 from:

``` r
remotes::install_github("tlverse/sl3@devel")
```

To contribute, install the *mediation* version of `SuperNOVA` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("blind-contours/SuperNOVA@mediation")
```

This is under development and issues should be written in the issues
page.

------------------------------------------------------------------------

## Example

To illustrate how `SuperNOVA` may be used to ascertain the effect of a
mixed exposure, consider the following example:

``` r
library(SuperNOVA)
library(devtools)
#> Loading required package: usethis
library(kableExtra)

load_all("~/sl3")
#> ℹ Loading sl3
set.seed(429153)
# simulate simple data
n_obs <- 100000
```

The `simulate_data` function creates simulated data with a multivariate
exposure, covariates (confounders), and a continuous outcome.

``` r
data <- simulate_data(n_obs = n_obs, shift_var_index = c(3))
effect <- data$effect
data <- data$data
head(data) %>%
  kbl(caption = "Simulated Data") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Simulated Data
</caption>
<thead>
<tr>
<th style="text-align:right;">
M1
</th>
<th style="text-align:right;">
M2
</th>
<th style="text-align:right;">
M3
</th>
<th style="text-align:right;">
W1
</th>
<th style="text-align:right;">
W2
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
Y_shifted
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
5.122589
</td>
<td style="text-align:right;">
0.6810488
</td>
<td style="text-align:right;">
3.098372
</td>
<td style="text-align:right;">
7.873068
</td>
<td style="text-align:right;">
7.362137
</td>
<td style="text-align:right;">
8.727373
</td>
<td style="text-align:right;">
8.796372
</td>
</tr>
<tr>
<td style="text-align:right;">
5.752568
</td>
<td style="text-align:right;">
3.9186685
</td>
<td style="text-align:right;">
4.723786
</td>
<td style="text-align:right;">
6.677126
</td>
<td style="text-align:right;">
7.302253
</td>
<td style="text-align:right;">
8.549657
</td>
<td style="text-align:right;">
8.552309
</td>
</tr>
<tr>
<td style="text-align:right;">
4.515052
</td>
<td style="text-align:right;">
0.9926542
</td>
<td style="text-align:right;">
2.979040
</td>
<td style="text-align:right;">
7.598614
</td>
<td style="text-align:right;">
7.647076
</td>
<td style="text-align:right;">
8.290229
</td>
<td style="text-align:right;">
8.369829
</td>
</tr>
<tr>
<td style="text-align:right;">
3.816552
</td>
<td style="text-align:right;">
-0.1666900
</td>
<td style="text-align:right;">
1.502374
</td>
<td style="text-align:right;">
7.321127
</td>
<td style="text-align:right;">
6.907626
</td>
<td style="text-align:right;">
6.539904
</td>
<td style="text-align:right;">
7.549933
</td>
</tr>
<tr>
<td style="text-align:right;">
4.797729
</td>
<td style="text-align:right;">
1.8740593
</td>
<td style="text-align:right;">
4.024303
</td>
<td style="text-align:right;">
6.203951
</td>
<td style="text-align:right;">
6.547683
</td>
<td style="text-align:right;">
7.192016
</td>
<td style="text-align:right;">
7.200583
</td>
</tr>
<tr>
<td style="text-align:right;">
4.380270
</td>
<td style="text-align:right;">
0.2098730
</td>
<td style="text-align:right;">
2.961153
</td>
<td style="text-align:right;">
7.677073
</td>
<td style="text-align:right;">
8.250861
</td>
<td style="text-align:right;">
8.473848
</td>
<td style="text-align:right;">
8.556300
</td>
</tr>
</tbody>
</table>

The `shift_var_index` parameter above shifts a variable set and gets the
expected outcome under this shift. Here, we shift our first exposure
variable and true effect for this DGP is:

``` r
effect
#> [1] 0.2057182
```

And therefore, in `SuperNOVA` we would expect most of the fold CIs to
cover this number and the pooled estimate to also cover this true
effect. Let’s run `SuperNOVA` to see if it correctly identifies the
exposures that drive the outcome and any interaction/effect modification
that exists in the DGP.

Of note, there are three exposures M1, M2, M3 - M1 and M3 have
individual effects and interactions that drive the outcome. There is
also effect modification between M3 and W1.

``` r
data_sample <- data[sample(nrow(data), 1000), ]

w <- data_sample[, c("W1", "W2")]
a <- data_sample[, c("M1", "M2", "M3")]
y <- data_sample$Y

deltas <- list("M1" = 1, "M2" = 1, "M3" = 1)

ptm <- proc.time()
sim_results <- SuperNOVA(
  w = w,
  a = a,
  y = y,
  delta = deltas,
  n_folds = 2,
  num_cores = 6,
  family = "continuous",
  quantile_thresh = 0,
  seed = 294580
)
#> 
#> Iter: 1 fn: 736.1309  Pars:  0.999992466 0.000007534
#> Iter: 2 fn: 736.1309  Pars:  0.99999789 0.00000211
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 731.4051  Pars:  0.9997671 0.0002329
#> Iter: 2 fn: 731.4051  Pars:  0.99996729 0.00003271
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 736.7229  Pars:  0.99996647 0.00003353
#> Iter: 2 fn: 736.7229  Pars:  0.999992346 0.000007654
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 731.7438  Pars:  0.97476 0.02524
#> Iter: 2 fn: 731.5662  Pars:  0.74017 0.25983
#> Iter: 3 fn: 731.5662  Pars:  0.74017 0.25983
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 733.0992  Pars:  0.99996393 0.00003607
#> Iter: 2 fn: 733.0992  Pars:  0.99998573 0.00001427
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 738.3274  Pars:  0.99996537 0.00003463
#> Iter: 2 fn: 738.3274  Pars:  0.99998015 0.00001985
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 734.4877  Pars:  0.999991798 0.000008202
#> Iter: 2 fn: 734.4877  Pars:  0.999994932 0.000005068
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 737.4888  Pars:  0.99998951 0.00001049
#> Iter: 2 fn: 737.4888  Pars:  0.999997809 0.000002191
#> solnp--> Completed in 2 iterations
proc.time() - ptm
#>     user   system  elapsed 
#>   62.457    6.183 1035.798

indiv_shift_results <- sim_results$`Indiv Shift Results`
em_results <- sim_results$`Effect Mod Results`
joint_shift_results <- sim_results$`Joint Shift Results`
```

Let’s first look at the results for individual stochastic shifts by
delta compared to no shift:

``` r
indiv_shift_results$M3 %>%
  kbl(caption = "Individual Stochastic Intervention Results for M1") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Individual Stochastic Intervention Results for M1
</caption>
<thead>
<tr>
<th style="text-align:left;">
Condition
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:left;">
Fold
</th>
<th style="text-align:left;">
Type
</th>
<th style="text-align:left;">
Variables
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
0.4781274
</td>
<td style="text-align:right;">
0.0187579
</td>
<td style="text-align:right;">
0.1369596
</td>
<td style="text-align:right;">
0.2097
</td>
<td style="text-align:right;">
0.7466
</td>
<td style="text-align:right;">
0.0004812
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
0.2560727
</td>
<td style="text-align:right;">
0.1685680
</td>
<td style="text-align:right;">
0.4105704
</td>
<td style="text-align:right;">
-0.5486
</td>
<td style="text-align:right;">
1.0608
</td>
<td style="text-align:right;">
0.5328247
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
0.3098651
</td>
<td style="text-align:right;">
0.0347369
</td>
<td style="text-align:right;">
0.1863783
</td>
<td style="text-align:right;">
-0.0554
</td>
<td style="text-align:right;">
0.6752
</td>
<td style="text-align:right;">
0.0964005
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
1000
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

Next we can look at effect modifications:

``` r
em_results$M3W1 %>%
  kbl(caption = "Effect Modification Stochastic Intervention Results for M3 and W1") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Effect Modification Stochastic Intervention Results for M3 and W1
</caption>
<thead>
<tr>
<th style="text-align:left;">
Condition
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:left;">
Fold
</th>
<th style="text-align:left;">
Type
</th>
<th style="text-align:left;">
Variables
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 6.49113162474391
</td>
<td style="text-align:right;">
7.127136
</td>
<td style="text-align:right;">
0.0328169
</td>
<td style="text-align:right;">
0.1811544
</td>
<td style="text-align:right;">
6.7721
</td>
<td style="text-align:right;">
7.4822
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 6.49113162474391
</td>
<td style="text-align:right;">
7.997249
</td>
<td style="text-align:right;">
0.0874509
</td>
<td style="text-align:right;">
0.2957210
</td>
<td style="text-align:right;">
7.4176
</td>
<td style="text-align:right;">
8.5769
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 7.05053831347608
</td>
<td style="text-align:right;">
7.025923
</td>
<td style="text-align:right;">
0.3024329
</td>
<td style="text-align:right;">
0.5499390
</td>
<td style="text-align:right;">
5.9481
</td>
<td style="text-align:right;">
8.1038
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 7.05053831347608
</td>
<td style="text-align:right;">
9.351465
</td>
<td style="text-align:right;">
0.0739362
</td>
<td style="text-align:right;">
0.2719121
</td>
<td style="text-align:right;">
8.8185
</td>
<td style="text-align:right;">
9.8844
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 6.02939925010252
</td>
<td style="text-align:right;">
6.686384
</td>
<td style="text-align:right;">
0.0872849
</td>
<td style="text-align:right;">
0.2954401
</td>
<td style="text-align:right;">
6.1073
</td>
<td style="text-align:right;">
7.2654
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
1000
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 6.02939925010252
</td>
<td style="text-align:right;">
8.391628
</td>
<td style="text-align:right;">
0.1062345
</td>
<td style="text-align:right;">
0.3259364
</td>
<td style="text-align:right;">
7.7528
</td>
<td style="text-align:right;">
9.0305
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
1000
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

And finally results for the joint shift which is a joint shift compared
to additive individual shifts.

``` r
joint_shift_results$M1M3 %>%
  kbl(caption = "Interactions Stochastic Intervention Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Interactions Stochastic Intervention Results
</caption>
<tbody>
<tr>
</tr>
</tbody>
</table>

    ---

    ## Issues

    If you encounter any bugs or have any specific feature requests, please [file an
    issue](https://github.com/blind-contours/SuperNOVA/issues). Further details on filing
    issues are provided in our [contribution
    guidelines](https://github.com/blind-contours/SuperNOVA/blob/master/CONTRIBUTING.md).

    ---

    ## Contributions

    Contributions are very welcome. Interested contributors should consult our
    [contribution
    guidelines](https://github.com/blind-contours/SuperNOVA/blob/master/CONTRIBUTING.md)
    prior to submitting a pull request.

    ---

    ## Citation

    After using the `SuperNOVA` R package, please cite the following:


    ---

    ## Related

    * [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R package
      providing an independent implementation of the same core routines for the TML
      estimation procedure and statistical methodology as is made available here,
      through reliance on a unified interface for Targeted Learning provided by the
      [`tmle3`](https://github.com/tlverse/tmle3) engine of the [`tlverse`
      ecosystem](https://github.com/tlverse).

    * [R/`medshift`](https://github.com/nhejazi/medshift) - An R package providing
      facilities to estimate the causal effect of stochastic treatment regimes in
      the mediation setting, including classical (IPW) and augmented double robust
      (one-step) estimators. This is an implementation of the methodology explored
      by @diaz2020causal.

    * [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal package
      for estimating the conditional density treatment mechanism component of this
      parameter based on using the [highly adaptive
      lasso](https://github.com/tlverse/hal9001) [@coyle-hal9001-rpkg;
      @hejazi2020hal9001-joss] in combination with a pooled hazard regression. This
      package implements a variant of the approach advocated by @diaz2011super.

    ---

    ## Funding

    The development of this software was supported in part through grants from the


    ---

    ## License

    &copy; 2020-2022 [David B. McCoy](https://davidmccoy.org)

    The contents of this repository are distributed under the MIT license. See below
    for details:

MIT License Copyright (c) 2020-2022 David B. McCoy Permission is hereby
granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the “Software”), to deal in the
Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions: The above
copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
“AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE. \`\`\`

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-coyle-sl3-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2022. *<span class="nocase">sl3</span>: Modern Machine
Learning Pipelines for Super Learning*.
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-diaz2012population" class="csl-entry">

Dı́az, Iván, and Mark J van der Laan. 2012. “Population Intervention
Causal Effects Based on Stochastic Interventions.” *Biometrics* 68 (2):
541–49.

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

</div>
