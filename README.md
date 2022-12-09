
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
intervention across strata of an effect modifier which is
data-adaptively determined. Future work will also include stochastic
shift mediation analysis which is under development.

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
data <- simulate_data(n_obs = n_obs, shift_var_index = c(1))
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
6.122589
</td>
<td style="text-align:right;">
0.6810488
</td>
<td style="text-align:right;">
2.0983722
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
9.881214
</td>
</tr>
<tr>
<td style="text-align:right;">
6.752568
</td>
<td style="text-align:right;">
3.9186685
</td>
<td style="text-align:right;">
3.7237856
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
10.729448
</td>
</tr>
<tr>
<td style="text-align:right;">
5.515052
</td>
<td style="text-align:right;">
0.9926542
</td>
<td style="text-align:right;">
1.9790401
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
8.926546
</td>
</tr>
<tr>
<td style="text-align:right;">
4.816552
</td>
<td style="text-align:right;">
-0.1666900
</td>
<td style="text-align:right;">
0.5023735
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
6.791823
</td>
</tr>
<tr>
<td style="text-align:right;">
5.797729
</td>
<td style="text-align:right;">
1.8740593
</td>
<td style="text-align:right;">
3.0243031
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
8.042248
</td>
</tr>
<tr>
<td style="text-align:right;">
5.380270
</td>
<td style="text-align:right;">
0.2098730
</td>
<td style="text-align:right;">
1.9611529
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
9.032642
</td>
</tr>
</tbody>
</table>

The `shift_var_index` parameter above shifts a variable set and gets the
expected outcome under this shift. Here, we shift our first exposure
variable and true effect for this DGP is:

``` r
effect
#> [1] 1.149002
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

ptm <- proc.time()
sim_results <- SuperNOVA(
  w = w,
  a = a,
  y = y,
  delta = 1,
  n_folds = 5,
  num_cores = 6,
  family = "continuous",
  quantile_thresh = 0,
  seed = 294580
)
#> 
#> Iter: 1 fn: 1166.4198     Pars:  0.76564 0.23436
#> Iter: 2 fn: 1166.4198     Pars:  0.76569 0.23431
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1177.9483     Pars:  0.99997794 0.00002207
#> Iter: 2 fn: 1177.9483     Pars:  0.99998569 0.00001431
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1187.4375     Pars:  0.99998354 0.00001646
#> Iter: 2 fn: 1187.4375     Pars:  0.999995465 0.000004535
#> Iter: 3 fn: 1187.4375     Pars:  0.999998193 0.000001807
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 1170.4194     Pars:  0.68990 0.31010
#> Iter: 2 fn: 1170.4194     Pars:  0.68991 0.31009
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1166.6846     Pars:  0.992424 0.007576
#> Iter: 2 fn: 1166.6846     Pars:  0.992995 0.007005
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1180.3693     Pars:  0.00002575 0.99997425
#> Iter: 2 fn: 1180.3693     Pars:  0.00001267 0.99998733
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1163.1972     Pars:  0.999998791 0.000001209
#> Iter: 2 fn: 1163.1972     Pars:  0.9999992831 0.0000007169
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1175.6696     Pars:  0.85537 0.14463
#> Iter: 2 fn: 1175.6696     Pars:  0.85549 0.14451
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1176.4228     Pars:  0.83788 0.16212
#> Iter: 2 fn: 1176.3986     Pars:  0.59815 0.40185
#> Iter: 3 fn: 1176.3986     Pars:  0.59814 0.40186
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 1157.5911     Pars:  0.999998904 0.000001096
#> Iter: 2 fn: 1157.5911     Pars:  0.9999998731 0.0000001269
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1164.1017     Pars:  0.9999990908 0.0000009093
#> Iter: 2 fn: 1164.1017     Pars:  0.9999994669 0.0000005331
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1172.1758     Pars:  0.22288 0.77712
#> Iter: 2 fn: 1172.1690     Pars:  0.13125 0.86875
#> Iter: 3 fn: 1172.1690     Pars:  0.13125 0.86875
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 1172.7308     Pars:  0.72105 0.27895
#> Iter: 2 fn: 1172.6851     Pars:  0.99992727 0.00007273
#> Iter: 3 fn: 1172.6851     Pars:  0.99995324 0.00004676
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 1170.7932     Pars:  0.16966 0.83034
#> Iter: 2 fn: 1170.6595     Pars:  0.77188 0.22812
#> Iter: 3 fn: 1170.6595     Pars:  0.77188 0.22812
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 1157.6767     Pars:  0.99993859 0.00006141
#> Iter: 2 fn: 1157.6767     Pars:  0.999990081 0.000009919
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1167.1661     Pars:  0.999997510 0.000002488
#> Iter: 2 fn: 1167.1661     Pars:  0.9999996965 0.0000003035
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1040.0668     Pars:  0.57310 0.42690
#> Iter: 2 fn: 1040.0668     Pars:  0.57310 0.42690
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1174.4439     Pars:  0.99997832 0.00002168
#> Iter: 2 fn: 1174.4439     Pars:  0.999992464 0.000007536
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1162.2626     Pars:  0.999991615 0.000008386
#> Iter: 2 fn: 1162.2626     Pars:  0.999997643 0.000002357
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1157.9708     Pars:  0.99998993 0.00001006
#> Iter: 2 fn: 1157.9708     Pars:  0.999996087 0.000003913
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1162.9490     Pars:  0.999995269 0.000004731
#> Iter: 2 fn: 1162.9490     Pars:  0.9999993592 0.0000006408
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1165.8277     Pars:  0.64072 0.35928
#> Iter: 2 fn: 1165.8276     Pars:  0.63379 0.36621
#> Iter: 3 fn: 1165.8276     Pars:  0.63379 0.36621
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 1159.5643     Pars:  0.999998931 0.000001069
#> Iter: 2 fn: 1159.5643     Pars:  0.99999977 0.00000023
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 1025.5563     Pars:  0.999998917 0.000001082
#> Iter: 2 fn: 1025.5563     Pars:  0.9999993613 0.0000006387
#> solnp--> Completed in 2 iterations
proc.time() - ptm
#>     user   system  elapsed 
#>  186.260   13.332 2592.653

indiv_shift_results <- sim_results$`Indiv Shift Results`
em_results <- sim_results$`Effect Mod Results`
joint_shift_results <- sim_results$`Joint Shift Results`
```

Let’s first look at the results for individual stochastic shifts by
delta compared to no shift:

``` r
indiv_shift_results$M1 %>%
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
0.8118000
</td>
<td style="text-align:right;">
0.0421613
</td>
<td style="text-align:right;">
0.2053321
</td>
<td style="text-align:right;">
0.4094
</td>
<td style="text-align:right;">
1.2142
</td>
<td style="text-align:right;">
7.7e-05
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
1.1249803
</td>
<td style="text-align:right;">
0.0175037
</td>
<td style="text-align:right;">
0.1323017
</td>
<td style="text-align:right;">
0.8657
</td>
<td style="text-align:right;">
1.3843
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
0.9113407
</td>
<td style="text-align:right;">
0.0418939
</td>
<td style="text-align:right;">
0.2046801
</td>
<td style="text-align:right;">
0.5102
</td>
<td style="text-align:right;">
1.3125
</td>
<td style="text-align:right;">
8.5e-06
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
0.9725846
</td>
<td style="text-align:right;">
0.0182654
</td>
<td style="text-align:right;">
0.1351497
</td>
<td style="text-align:right;">
0.7077
</td>
<td style="text-align:right;">
1.2375
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
1.2015947
</td>
<td style="text-align:right;">
0.0421926
</td>
<td style="text-align:right;">
0.2054084
</td>
<td style="text-align:right;">
0.7990
</td>
<td style="text-align:right;">
1.6042
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
1.0515509
</td>
<td style="text-align:right;">
0.0058070
</td>
<td style="text-align:right;">
0.0762038
</td>
<td style="text-align:right;">
0.9022
</td>
<td style="text-align:right;">
1.2009
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
1000
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 6.5051649726032
</td>
<td style="text-align:right;">
7.396953
</td>
<td style="text-align:right;">
0.0216254
</td>
<td style="text-align:right;">
0.1470559
</td>
<td style="text-align:right;">
7.1087
</td>
<td style="text-align:right;">
7.6852
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
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 6.5051649726032
</td>
<td style="text-align:right;">
8.286543
</td>
<td style="text-align:right;">
0.0409230
</td>
<td style="text-align:right;">
0.2022943
</td>
<td style="text-align:right;">
7.8901
</td>
<td style="text-align:right;">
8.6830
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
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 6.49354965732094
</td>
<td style="text-align:right;">
7.025393
</td>
<td style="text-align:right;">
0.0042238
</td>
<td style="text-align:right;">
0.0649907
</td>
<td style="text-align:right;">
6.8980
</td>
<td style="text-align:right;">
7.1528
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 6.49354965732094
</td>
<td style="text-align:right;">
8.274424
</td>
<td style="text-align:right;">
0.2327762
</td>
<td style="text-align:right;">
0.4824689
</td>
<td style="text-align:right;">
7.3288
</td>
<td style="text-align:right;">
9.2200
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 6.78732393155181
</td>
<td style="text-align:right;">
7.370033
</td>
<td style="text-align:right;">
0.0560709
</td>
<td style="text-align:right;">
0.2367929
</td>
<td style="text-align:right;">
6.9059
</td>
<td style="text-align:right;">
7.8341
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 6.78732393155181
</td>
<td style="text-align:right;">
8.457109
</td>
<td style="text-align:right;">
0.0097595
</td>
<td style="text-align:right;">
0.0987901
</td>
<td style="text-align:right;">
8.2635
</td>
<td style="text-align:right;">
8.6507
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 6.5051649726032
</td>
<td style="text-align:right;">
7.102586
</td>
<td style="text-align:right;">
0.1482954
</td>
<td style="text-align:right;">
0.3850914
</td>
<td style="text-align:right;">
6.3478
</td>
<td style="text-align:right;">
7.8574
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 6.5051649726032
</td>
<td style="text-align:right;">
7.628454
</td>
<td style="text-align:right;">
0.0289861
</td>
<td style="text-align:right;">
0.1702529
</td>
<td style="text-align:right;">
7.2948
</td>
<td style="text-align:right;">
7.9621
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W1
</td>
<td style="text-align:right;">
200
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W1 \<= 6.5051649726032
</td>
<td style="text-align:right;">
6.982645
</td>
<td style="text-align:right;">
0.0050493
</td>
<td style="text-align:right;">
0.0710581
</td>
<td style="text-align:right;">
6.8434
</td>
<td style="text-align:right;">
7.1219
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
800
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W1 \<= 6.5051649726032
</td>
<td style="text-align:right;">
8.316520
</td>
<td style="text-align:right;">
0.0229209
</td>
<td style="text-align:right;">
0.1513967
</td>
<td style="text-align:right;">
8.0198
</td>
<td style="text-align:right;">
8.6133
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
800
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
