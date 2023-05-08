
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

The `SuperNOVA` R package provides users with the tools necessary to
identify the most predictive variable sets for a given outcome and
develop efficient estimators for the counterfactual mean of the outcome
under stochastic interventions on those variables as described fully in our 
paper: https://arxiv.org/abs/2305.01849. 

These interventions
are shifts to the exposures that are dependent on naturally observed
values (Dı́az and van der Laan 2012; Haneuse and Rotnitzky 2013).
Building on the `txshift` package, which implements the TML estimator
for a stochastic shift causal parameter, `SuperNOVA` extends this
methodology to include joint stochastic interventions on two variables,
allowing for the construction of a non-parametric interaction parameter
and mediation (in development). Additionally, `SuperNOVA` estimates
individual stochastic intervention outcomes under some delta shift
compared to the outcome under no intervention, and a target parameter
for effect modification, which is the mean outcome under intervention in
regions of the covariate space that are also data-adaptively determined.

The `SuperNOVA` package provides a comprehensive solution for
identifying variable sets that interact or modify effects in the context
of mixed exposures. To achieve this, we use a k-fold cross-validation
framework to estimate a data-adaptive parameter, namely, stochastic
shift target parameters for variable sets that are discovered to be
predictive of the outcome. We ensure unbiased estimation of the
data-adaptive parameter by employing cross-validated targeted maximum
likelihood estimation. In this approach, we split the data into
parameter-generating and estimation samples. In the parameter-generating
sample, we fit an ensemble of basis function estimators to the data and
select the estimator with the lowest cross-validated mean squared error.
We then extract important variable sets using ANOVA-like variance
decomposition of the linear combination of basis functions. In the
estimation fold, we use targeted learning to estimate causal target
parameters for interaction, effect modification, and individual variable
shifts. That is, we estimate the counterfactual mean different of these
discovered variable sets under a shift of $\delta$ (an amount to shift
an exposure by) compared to the observed outcome under observed
exposure. `SuperNOVA` is a versatile tool that provides researchers with
k-fold specific and pooled results for each target parameter. More
details are available in the accompanying vignette.

Users simply input a vector for exposures, covariates, and an outcome.
The user also specifies the respective $\delta$ for each exposure (the
amount to shift by) and if this delta should be adaptive based on
positivity violations (see vignette). `SuperNOVA` comes with flexible
default machine learning algorithms used in the data-adaptive procedure
and for estimation of each of the nuisance parameters. Given these
inputs `SuperNOVA` outputs tables for fold specific results (say
exposure 1 under shifts for each fold it was found predictive) and
pooled results which uses a pooled TMLE fluctuation across the folds to
esimate an average effect.

`SuperNOVA` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) (Coyle, Hejazi, Malenica, et
al. 2022) to allow for ensemble machine learning to be leveraged in the
estimation procedure for each nuisance parameter and estimation of the
data-adaptive parameters. There are several stacks of machine learning
algorithms used that are constructed from `sl3` automatically. If the
stack parameters are NULL, SuperNOVA automatically builds ensembles of
machine learning algorithms that are flexible yet not overly
computationally taxing.

------------------------------------------------------------------------

## Installation

*Note:* Because the `SuperNOVA` package (currently) depends on `sl3`
that allows ensemble machine learning to be used for nuisance parameter
estimation and `sl3` is not on CRAN the `SuperNOVA` package is not
available on CRAN and must be downloaded here.

There are many depedencies for `SuperNOVA` so it’s easier to break up
installation of the various packages to ensure proper installation.

First install the basis estimators used in the data-adaptive variable
discovery of the exposure and covariate space:

``` r
install.packages("earth")
install.packages("hal9001")
```

`SuperNOVA` uses the `sl3` package to build ensemble machine learners
for each nuisance parameter. We have to install off the development
branch, first download these two packages for `sl3`

``` r
install.packages(c("ranger", "arm", "xgboost", "nnls"))
```

Now install `sl3` on devel:

``` r
remotes::install_github("tlverse/sl3@devel")
```

Make sure `sl3` installs correctly then install `SuperNOVA`

``` r
remotes::install_github("blind-contours/SuperNOVA@main")
```

`SuperNOVA` has some other miscellaneous dependencies that are used in
the examples as well as in the plotting functions.

``` r
install.packages(c("kableExtra", "hrbrthemes", "viridis"))
```

------------------------------------------------------------------------

## Example

To illustrate how `SuperNOVA` may be used to ascertain the effect of a
mixed exposure, consider the following example:

``` r
library(SuperNOVA)
library(devtools)
#> Loading required package: usethis
library(kableExtra)
library(sl3)

set.seed(429153)
# simulate simple data
n_obs <- 100000
```

The `simulate_data` is a function for generating synthetic data with a
complex structure to study the causal effects of shifting values in the
mixtures of exposures The primary purpose of this simulation is to
provide a controlled environment for testing and validating estimates
given by the SuperNOVA package.

The simulate_data function generates synthetic data for n_obs
observations with a pre-specified covariance structure (sigma_mod) and a
shift parameter (delta). It simulates four mixture components (M1, M2,
M3, M4) and three covariates (W1, W2, W3) with specific relationships
between them. The outcome variable Y is generated as a function of these
mixtures and covariates.

After generating the data, the function applies a shift (delta) to each
mixture component separately and calculates the average treatment effect
for each component. Additionally, it calculates the interaction effect
of shifting two mixture components simultaneously (m14_intxn). These
ground truth effects can be used for validating and comparing the
performance of various causal inference methods on this synthetic
dataset.

The function returns a list containing the generated data and the ground
truth effects of the shifts applied to each mixture component, the
interaction effect, and the modified effect results based on a specific
level in the W3 covariate.

``` r
sim_out <- simulate_data(n_obs = n_obs)
data <- sim_out$data
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
M4
</th>
<th style="text-align:right;">
W1
</th>
<th style="text-align:right;">
W2
</th>
<th style="text-align:right;">
W3
</th>
<th style="text-align:right;">
Y
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
23.01099
</td>
<td style="text-align:right;">
3.864569
</td>
<td style="text-align:right;">
4.461534
</td>
<td style="text-align:right;">
4.747169
</td>
<td style="text-align:right;">
7.873068
</td>
<td style="text-align:right;">
7.362137
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
72.76739
</td>
</tr>
<tr>
<td style="text-align:right;">
23.70848
</td>
<td style="text-align:right;">
4.169044
</td>
<td style="text-align:right;">
4.990670
</td>
<td style="text-align:right;">
5.568766
</td>
<td style="text-align:right;">
6.677126
</td>
<td style="text-align:right;">
7.302253
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
87.81536
</td>
</tr>
<tr>
<td style="text-align:right;">
22.44274
</td>
<td style="text-align:right;">
2.922274
</td>
<td style="text-align:right;">
4.885933
</td>
<td style="text-align:right;">
3.347086
</td>
<td style="text-align:right;">
7.598614
</td>
<td style="text-align:right;">
7.647076
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
42.33824
</td>
</tr>
<tr>
<td style="text-align:right;">
22.99681
</td>
<td style="text-align:right;">
2.507461
</td>
<td style="text-align:right;">
4.227508
</td>
<td style="text-align:right;">
1.719390
</td>
<td style="text-align:right;">
7.321127
</td>
<td style="text-align:right;">
6.907626
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
38.72763
</td>
</tr>
<tr>
<td style="text-align:right;">
22.88541
</td>
<td style="text-align:right;">
3.349974
</td>
<td style="text-align:right;">
4.369185
</td>
<td style="text-align:right;">
6.618387
</td>
<td style="text-align:right;">
6.203951
</td>
<td style="text-align:right;">
6.547683
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
90.86418
</td>
</tr>
<tr>
<td style="text-align:right;">
22.88781
</td>
<td style="text-align:right;">
3.890133
</td>
<td style="text-align:right;">
5.370329
</td>
<td style="text-align:right;">
5.736624
</td>
<td style="text-align:right;">
7.677073
</td>
<td style="text-align:right;">
8.250861
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
68.59041
</td>
</tr>
</tbody>
</table>

And therefore, in `SuperNOVA` we would expect most of the fold CIs to
cover this number and the pooled estimate to also cover this true
effect. Let’s run `SuperNOVA` to see if it correctly identifies the
exposures that drive the outcome and any interaction/effect modification
that exists in the DGP.

Of note, there are three exposures M1, M2, M3 - M1 and M3 have
individual effects and interactions that drive the outcome. There is
also effect modification between M3 and W1.

``` r
data_sample <- data[sample(nrow(data), 4000), ]

w <- data_sample[, c("W1", "W2", "W3")]
a <- data_sample[, c("M1", "M2", "M3", "M4")]
y <- data_sample$Y

deltas <- list("M1" = 1, "M2" = 1, "M3" = 1, "M4" = 1)

ptm <- proc.time()
sim_results <- SuperNOVA(
  w = w,
  a = a,
  y = y,
  delta = deltas,
  n_folds = 3,
  num_cores = 6,
  outcome_type = "continuous",
  quantile_thresh = 0,
  seed = 294580
)
#> 
#> Iter: 1 fn: 5617.4930     Pars:  0.83767 0.16233
#> Iter: 2 fn: 5617.1485     Pars:  0.999991224 0.000008777
#> Iter: 3 fn: 5617.1484     Pars:  0.999994372 0.000005628
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 3814.5345     Pars:  0.99998864 0.00001136
#> Iter: 2 fn: 3814.5345     Pars:  0.999996147 0.000003853
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3732.3492     Pars:  0.999994933 0.000005067
#> Iter: 2 fn: 3732.3492     Pars:  0.999996634 0.000003366
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3731.9982     Pars:  0.88391 0.11609
#> Iter: 2 fn: 3731.8248     Pars:  0.999991263 0.000008737
#> Iter: 3 fn: 3731.8248     Pars:  0.999994386 0.000005614
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5615.8454     Pars:  0.43869 0.56131
#> Iter: 2 fn: 5615.8454     Pars:  0.43869 0.56131
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 5620.9486     Pars:  0.44141 0.55859
#> Iter: 2 fn: 5620.4981     Pars:  0.99997607 0.00002393
#> Iter: 3 fn: 5620.4981     Pars:  0.99998463 0.00001537
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5614.1167     Pars:  0.82113 0.17887
#> Iter: 2 fn: 5614.1077     Pars:  0.9994545 0.0005455
#> Iter: 3 fn: 5614.1077     Pars:  0.9996802 0.0003198
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 3845.2259     Pars:  0.99998184 0.00001816
#> Iter: 2 fn: 3845.2259     Pars:  0.999993039 0.000006961
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3815.6067     Pars:  0.69644 0.30356
#> Iter: 2 fn: 3815.6067     Pars:  0.69650 0.30350
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3840.8999     Pars:  0.89353 0.10647
#> Iter: 2 fn: 3840.8999     Pars:  0.89353 0.10647
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 5615.8566     Pars:  0.73031 0.26969
#> Iter: 2 fn: 5615.7797     Pars:  0.9992627 0.0007373
#> Iter: 3 fn: 5615.7796     Pars:  0.9995925 0.0004075
#> Iter: 4 fn: 5615.7795     Pars:  0.9998947 0.0001053
#> Iter: 5 fn: 5615.7795     Pars:  0.99995525 0.00004475
#> solnp--> Completed in 5 iterations
#> 
#> Iter: 1 fn: 5610.6494     Pars:  0.44630 0.55370
#> Iter: 2 fn: 5610.6494     Pars:  0.44627 0.55373
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 5599.9737     Pars:  0.99997813 0.00002187
#> Iter: 2 fn: 5599.9736     Pars:  0.999992193 0.000007807
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3829.3779     Pars:  0.89231 0.10769
#> Iter: 2 fn: 3829.1156     Pars:  0.61038 0.38962
#> Iter: 3 fn: 3829.1156     Pars:  0.61031 0.38969
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 3802.6635     Pars:  0.999992896 0.000007104
#> Iter: 2 fn: 3802.6635     Pars:  0.99999763 0.00000237
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3831.5452     Pars:  0.47081 0.52919
#> Iter: 2 fn: 3831.5452     Pars:  0.47082 0.52918
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 5598.6825     Pars:  0.99998629 0.00001371
#> Iter: 2 fn: 5598.6825     Pars:  0.999994535 0.000005465
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 5600.0257     Pars:  0.82318 0.17682
#> Iter: 2 fn: 5599.8034     Pars:  0.99998028 0.00001972
#> Iter: 3 fn: 5599.8034     Pars:  0.99998732 0.00001268
#> solnp--> Completed in 3 iterations
proc.time() - ptm
#>    user  system elapsed 
#>  27.248   1.843 460.624

basis_in_folds <- sim_results$`Basis Fold Proportions`
indiv_shift_results <- sim_results$`Indiv Shift Results`
em_results <- sim_results$`Effect Mod Results`
joint_shift_results <- sim_results$`Joint Shift Results`
```

Let’s first look at the variable relationships used in the folds:

``` r
basis_in_folds
#> 
#>   M1 M1M4 M3M4 M3W3   M4   W3 
#> 1.00 0.33 0.67 1.00 1.00 1.00
```

The above list shows that marginal effects for exposures M1 and M4 were
found, an interaction for M1 and M4, and effect modification for M3 and
W3 - all of which are correct as the outcome is generated from these
relationships. There is no effect for M2 which we correctly reject.

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
<th style="text-align:right;">
Delta
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
2.594282
</td>
<td style="text-align:right;">
0.7612363
</td>
<td style="text-align:right;">
0.8724886
</td>
<td style="text-align:right;">
0.8842
</td>
<td style="text-align:right;">
4.3043
</td>
<td style="text-align:right;">
0.0029449
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
1334
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
1.907373
</td>
<td style="text-align:right;">
0.5215788
</td>
<td style="text-align:right;">
0.7222041
</td>
<td style="text-align:right;">
0.4919
</td>
<td style="text-align:right;">
3.3229
</td>
<td style="text-align:right;">
0.0082651
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
1333
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
2.256934
</td>
<td style="text-align:right;">
0.5683514
</td>
<td style="text-align:right;">
0.7538909
</td>
<td style="text-align:right;">
0.7793
</td>
<td style="text-align:right;">
3.7345
</td>
<td style="text-align:right;">
0.0027560
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
1333
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
1.522048
</td>
<td style="text-align:right;">
0.1679864
</td>
<td style="text-align:right;">
0.4098615
</td>
<td style="text-align:right;">
0.7187
</td>
<td style="text-align:right;">
2.3254
</td>
<td style="text-align:right;">
0.0002044
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
4000
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

The true effect for a shifted M1 vs observed M1 is:

``` r
sim_out$m1_effect
#> [1] 1.603984
```

And so we see that we have proper coverage.

Next we can look at effect modifications:

``` r
em_results$M3W3 %>%
  kbl(caption = "Effect Modification Stochastic Intervention Results for M3 and W3") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Effect Modification Stochastic Intervention Results for M3 and W3
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
Lower_CI
</th>
<th style="text-align:right;">
Upper_CI
</th>
<th style="text-align:right;">
P_value
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
Level 1 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
-3.248960
</td>
<td style="text-align:right;">
3.3836568
</td>
<td style="text-align:right;">
1.8394719
</td>
<td style="text-align:right;">
-6.8543
</td>
<td style="text-align:right;">
0.3563
</td>
<td style="text-align:right;">
0.0773546
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
4.308103
</td>
<td style="text-align:right;">
8.9139246
</td>
<td style="text-align:right;">
2.9856196
</td>
<td style="text-align:right;">
-1.5436
</td>
<td style="text-align:right;">
10.1598
</td>
<td style="text-align:right;">
0.1490342
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
3.458559
</td>
<td style="text-align:right;">
4.2797587
</td>
<td style="text-align:right;">
2.0687578
</td>
<td style="text-align:right;">
-0.5961
</td>
<td style="text-align:right;">
7.5132
</td>
<td style="text-align:right;">
0.0945629
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
1333
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
11.950912
</td>
<td style="text-align:right;">
3.4042433
</td>
<td style="text-align:right;">
1.8450592
</td>
<td style="text-align:right;">
8.3347
</td>
<td style="text-align:right;">
15.5672
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
1333
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
4.819357
</td>
<td style="text-align:right;">
38.8875342
</td>
<td style="text-align:right;">
6.2359870
</td>
<td style="text-align:right;">
-7.4030
</td>
<td style="text-align:right;">
17.0417
</td>
<td style="text-align:right;">
0.4396231
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
1333
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
11.272408
</td>
<td style="text-align:right;">
6.9938039
</td>
<td style="text-align:right;">
2.6445801
</td>
<td style="text-align:right;">
6.0891
</td>
<td style="text-align:right;">
16.4557
</td>
<td style="text-align:right;">
0.0000202
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
1333
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 1 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
1.704760
</td>
<td style="text-align:right;">
0.9422027
</td>
<td style="text-align:right;">
0.9706713
</td>
<td style="text-align:right;">
-0.1977
</td>
<td style="text-align:right;">
3.6072
</td>
<td style="text-align:right;">
0.0790424
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
4000
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Level 0 Shift Diff in W3 \<= 0
</td>
<td style="text-align:right;">
10.180014
</td>
<td style="text-align:right;">
1.3675158
</td>
<td style="text-align:right;">
1.1694083
</td>
<td style="text-align:right;">
7.8880
</td>
<td style="text-align:right;">
12.4720
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Effect Mod
</td>
<td style="text-align:left;">
M3W3
</td>
<td style="text-align:right;">
4000
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

Let’s first look at the truth:

``` r
sim_out$effect_mod
#> $`Level 0 Shift Diff in W3 <= 0`
#> [1] 10.99749
#> 
#> $`Level 1 Shift Diff in W3 <= 0`
#> [1] 1
```

When W3 is greater than 0 the truth effect is 11, our estimates are 10
with CI coverage. When W3 is 0 the truth is 1, our estimate is 1.9 with
CIs that cover the truth as well.

And finally results for the joint shift which is a joint shift compared
to additive individual shifts.

``` r
joint_shift_results$M1M4 %>%
  kbl(caption = "Interactions Stochastic Intervention Results for M1 and M4") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Interactions Stochastic Intervention Results for M1 and M4
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
Delta M1
</th>
<th style="text-align:right;">
Delta M4
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
2.599696
</td>
<td style="text-align:right;">
0.7610882
</td>
<td style="text-align:right;">
0.8724037
</td>
<td style="text-align:right;">
0.8898
</td>
<td style="text-align:right;">
4.3096
</td>
<td style="text-align:right;">
0.0053805
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M4
</td>
<td style="text-align:right;">
10.263304
</td>
<td style="text-align:right;">
0.3342128
</td>
<td style="text-align:right;">
0.5781114
</td>
<td style="text-align:right;">
9.1302
</td>
<td style="text-align:right;">
11.3964
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
11.587589
</td>
<td style="text-align:right;">
0.3370926
</td>
<td style="text-align:right;">
0.5805968
</td>
<td style="text-align:right;">
10.4496
</td>
<td style="text-align:right;">
12.7255
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Psi
</td>
<td style="text-align:right;">
-1.275411
</td>
<td style="text-align:right;">
0.7560428
</td>
<td style="text-align:right;">
0.8695072
</td>
<td style="text-align:right;">
-2.9796
</td>
<td style="text-align:right;">
0.4288
</td>
<td style="text-align:right;">
0.1713836
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
2.599696
</td>
<td style="text-align:right;">
0.7610882
</td>
<td style="text-align:right;">
0.8724037
</td>
<td style="text-align:right;">
0.8898
</td>
<td style="text-align:right;">
4.3096
</td>
<td style="text-align:right;">
0.0053805
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M4
</td>
<td style="text-align:right;">
10.263304
</td>
<td style="text-align:right;">
0.3342128
</td>
<td style="text-align:right;">
0.5781114
</td>
<td style="text-align:right;">
9.1302
</td>
<td style="text-align:right;">
11.3964
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
11.587589
</td>
<td style="text-align:right;">
0.3370926
</td>
<td style="text-align:right;">
0.5805968
</td>
<td style="text-align:right;">
10.4496
</td>
<td style="text-align:right;">
12.7255
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Psi
</td>
<td style="text-align:right;">
-1.275411
</td>
<td style="text-align:right;">
0.7560428
</td>
<td style="text-align:right;">
0.8695072
</td>
<td style="text-align:right;">
-2.9796
</td>
<td style="text-align:right;">
0.4288
</td>
<td style="text-align:right;">
0.1713836
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Interaction
</td>
<td style="text-align:left;">
M1&M4
</td>
<td style="text-align:right;">
1334
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

Let’s look at the truth again:

``` r
sim_out$m1_effect
#> [1] 1.603984
sim_out$m4_effect
#> [1] 10.41265
sim_out$m14_effect
#> [1] 12.41664
sim_out$m14_intxn
#> [1] 0.4
```

So comparing the results to the above table in the pooled section we can
see all our estimates for the marginal shifts, dual shift, and
difference between dual and sum of marginals have CIs that cover the
truth.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/SuperNOVA/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/%20SuperNOVA/main/contributing.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/SuperNOVA/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `SuperNOVA` R package, please cite the following:

------------------------------------------------------------------------

## Related

- [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R package
  providing an independent implementation of the same core routines for
  the TML estimation procedure and statistical methodology as is made
  available here, through reliance on a unified interface for Targeted
  Learning provided by the [`tmle3`](https://github.com/tlverse/tmle3)
  engine of the [`tlverse` ecosystem](https://github.com/tlverse).

- [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
  providing facilities to estimate the causal effect of stochastic
  treatment regimes in the mediation setting, including classical (IPW)
  and augmented double robust (one-step) estimators. This is an
  implementation of the methodology explored by Dı́az and Hejazi (2020).

- [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
  package for estimating the conditional density treatment mechanism
  component of this parameter based on using the [highly adaptive
  lasso](https://github.com/tlverse/hal9001) (Coyle, Hejazi, Phillips,
  et al. 2022; Hejazi, Coyle, and van der Laan 2020) in combination with
  a pooled hazard regression. This package implements a variant of the
  approach advocated by Dı́az and van der Laan (2011).

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through grants
from the

------------------------------------------------------------------------

## License

© 2020-2022 [David B. McCoy](https://davidmccoy.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2020-2022 David B. McCoy
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
