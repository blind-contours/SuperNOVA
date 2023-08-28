
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
> Interactions, Effect Modifications and Mediation using Stochastic
> Interventions **Authors:** [David McCoy](https://davidmccoy.org)

------------------------------------------------------------------------

## What’s `SuperNOVA`?

The `SuperNOVA` R package offers a comprehensive toolset for identifying
predictive variable sets, be it subsets of exposures,
exposure-covariates, or exposure-mediators, for a specified outcome. It
further assists in creating efficient estimators for the counterfactual
mean of the outcome under stochastic interventions on these variable
sets. This means making exposure changes that depend on naturally
observed values, as described in past literature (Dı́az and van der Laan
2012; Haneuse and Rotnitzky 2013).

`SuperNOVA` introduces several estimators, constructed based on the
patterns found in the data. At present, semi-parametric estimators are
available for various contexts: interaction, effect modification,
marginal impacts, and mediation. The target parameters that this package
calculates are grounded in previously published works: one regarding
interaction and effect modification (McCoy et al. 2023) and another
dedicated to mediation \[McCoy2023mediation\].

The `SuperNOVA` package builds upon the capabilities of the `txshift`
package, which implements the TML estimator for a stochastic shift
causal parameter (Dı́az and van der Laan 2012). A notable extension in
SuperNOVA is its support for joint stochastic interventions on two
exposures. This allows for the creation of a non-parametric interaction
parameter, shedding light on the combined effect of concurrently
shifting two variables relative to the cumulative effect of individual
shifts. At this stage, the focus is on two-way shifts.

Various parameters are calculated based on identified patterns:

- Interaction Parameter: When two exposures show signs of interaction,
  this parameter is calculated.

- Individual Stochastic Interventions: When marginal effects are
  identified, this estimates the difference in outcomes upon making a
  specific shift compared to no intervention.

- Effect Modification: In scenarios where effect modification is
  evident, the target parameter for it is the mean outcome under
  intervention within specific regions of the covariate space that the
  data suggests. This contrasts with the marginal effect by diving
  deeper into strata of covariates, for example, gauging the effects of
  exposure shifts among distinct genders.

- Mediation: In scenarios where mediation is found, the target parameter
  are the natural direct and indirect effects as defined by stochastic
  intervention. That is, for the natural indirect effect, `SuperNOVA`
  estimates the impact of shifting an exposure through a mediator and
  the direct effect which is not through the mediator.

Mediation in `SuperNOVA` is a new feature. Here, `SuperNOVA` brings to
the table estimates initially conceived in another work (Díaz and Hejazi
2020), while also supporting mediation estimates for continuous
exposures.

The package ensures robustness by employing a k-fold cross-validation
framework. This framework helps in estimating a data-adaptive parameter,
which is the stochastic shift target parameters for the variable sets
identified as influential for the outcome. The process begins by
partitioning the data into parameter-generating and estimation samples.
The former sample assists in fitting a collection of basis function
estimators to the data. The one with the lowest cross-validated mean
squared error is selected. Key variable sets are then extracted using an
ANOVA-like variance decomposition methodology. For the estimation
sample, targeted learning is harnessed to gauge causal target parameters
across different contexts: interaction, mediation, effect modification,
and individual variable shifts.

By using SuperNOVA, users get access to a tool that offers both k-fold
specific and aggregated results for each target parameter, ensuring that
researchers can glean the most information from their data. For a more
in-depth exploration, there’s an accompanying vignette.

To utilize the package, users need to provide vectors for exposures,
covariates, mediators, and outcomes. They also specify the respective
$\delta$ for each exposure (indicating the degree of shift) and if this
delta should be adaptive in response to positivity violations. A
detailed guide is provided in the vignette. With these inputs,
`SuperNOVA` processes the data and delivers tables showcasing
fold-specific results and aggregated outcomes, allowing users to glean
insights effectively.

`SuperNOVA` also incorporates features from the `sl3` package (Coyle,
Hejazi, Malenica, et al. 2022), facilitating ensemble machine learning
in the estimation process. If the user does not specify any stack
parameters, `SuperNOVA` will automatically create an ensemble of machine
learning algorithms that strike a balance between flexibility and
computational efficiency.

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
mixtures of exposures. The primary purpose of this simulation is to
provide a controlled environment for testing and validating estimates
given by the `SuperNOVA` package.

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
#> Iter: 1 fn: 5616.2821     Pars:  0.83768 0.16232
#> Iter: 2 fn: 5615.9460     Pars:  0.99996649 0.00003351
#> Iter: 3 fn: 5615.9460     Pars:  0.99997849 0.00002151
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 3815.7486     Pars:  0.99998993 0.00001007
#> Iter: 2 fn: 3815.7486     Pars:  0.999998834 0.000001166
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3737.1413     Pars:  0.999996084 0.000003915
#> Iter: 2 fn: 3737.1412     Pars:  0.9999991151 0.0000008849
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3731.8788     Pars:  0.88038 0.11962
#> Iter: 2 fn: 3731.8383     Pars:  0.99995226 0.00004774
#> Iter: 3 fn: 3731.8383     Pars:  0.99997338 0.00002662
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5615.8427     Pars:  0.37837 0.62163
#> Iter: 2 fn: 5615.7080     Pars:  0.61125 0.38875
#> Iter: 3 fn: 5615.7080     Pars:  0.61254 0.38746
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5621.0789     Pars:  0.84929 0.15071
#> Iter: 2 fn: 5620.8684     Pars:  0.00007712 0.99992288
#> Iter: 3 fn: 5620.8684     Pars:  0.000002771 0.999997229
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5612.7158     Pars:  0.82116 0.17884
#> Iter: 2 fn: 5612.5096     Pars:  0.9999838 0.0000162
#> Iter: 3 fn: 5612.5096     Pars:  0.99998959 0.00001041
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 3842.7726     Pars:  0.999856 0.000144
#> Iter: 2 fn: 3842.7726     Pars:  0.99995859 0.00004141
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3811.0844     Pars:  0.67596 0.32404
#> Iter: 2 fn: 3811.0844     Pars:  0.67601 0.32399
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3841.7125     Pars:  0.997439 0.002561
#> Iter: 2 fn: 3841.7125     Pars:  0.9993128 0.0006872
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 5612.7076     Pars:  0.43156 0.56844
#> Iter: 2 fn: 5612.3977     Pars:  0.999891 0.000109
#> Iter: 3 fn: 5612.3977     Pars:  0.99992987 0.00007013
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5621.1314     Pars:  0.51060 0.48940
#> Iter: 2 fn: 5621.1314     Pars:  0.51078 0.48922
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 5600.7646     Pars:  0.99997834 0.00002166
#> Iter: 2 fn: 5600.7645     Pars:  0.999995592 0.000004408
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3831.8447     Pars:  0.45887 0.54113
#> Iter: 2 fn: 3831.8447     Pars:  0.45886 0.54114
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 3799.4818     Pars:  0.02808 0.97192
#> Iter: 2 fn: 3797.3950     Pars:  0.9998439 0.0001561
#> Iter: 3 fn: 3797.3948     Pars:  0.99998699 0.00001301
#> Iter: 4 fn: 3797.3948     Pars:  0.999993236 0.000006764
#> solnp--> Completed in 4 iterations
#> 
#> Iter: 1 fn: 3830.8913     Pars:  0.84699 0.15301
#> Iter: 2 fn: 3830.8719     Pars:  0.75597 0.24403
#> Iter: 3 fn: 3830.8719     Pars:  0.75597 0.24403
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5600.6732     Pars:  0.13089 0.86911
#> Iter: 2 fn: 5600.6602     Pars:  0.007388 0.992612
#> Iter: 3 fn: 5600.6602     Pars:  0.007384 0.992616
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 5600.7959     Pars:  0.9998956 0.0001044
#> Iter: 2 fn: 5600.7958     Pars:  0.99998238 0.00001762
#> solnp--> Completed in 2 iterations
proc.time() - ptm
#>    user  system elapsed 
#>  51.213   4.756 824.470

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
2.664497
</td>
<td style="text-align:right;">
0.7389329
</td>
<td style="text-align:right;">
0.8596120
</td>
<td style="text-align:right;">
0.9797
</td>
<td style="text-align:right;">
4.3493
</td>
<td style="text-align:right;">
0.0019375
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
1.752571
</td>
<td style="text-align:right;">
0.5247758
</td>
<td style="text-align:right;">
0.7244141
</td>
<td style="text-align:right;">
0.3327
</td>
<td style="text-align:right;">
3.1724
</td>
<td style="text-align:right;">
0.0155506
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
2.214215
</td>
<td style="text-align:right;">
0.5609105
</td>
<td style="text-align:right;">
0.7489396
</td>
<td style="text-align:right;">
0.7463
</td>
<td style="text-align:right;">
3.6821
</td>
<td style="text-align:right;">
0.0031119
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
1.496780
</td>
<td style="text-align:right;">
0.1667218
</td>
<td style="text-align:right;">
0.4083158
</td>
<td style="text-align:right;">
0.6965
</td>
<td style="text-align:right;">
2.2971
</td>
<td style="text-align:right;">
0.0002466
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
-3.335721
</td>
<td style="text-align:right;">
3.3306066
</td>
<td style="text-align:right;">
1.8249950
</td>
<td style="text-align:right;">
-6.9126
</td>
<td style="text-align:right;">
0.2412
</td>
<td style="text-align:right;">
0.0675799
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
4.242965
</td>
<td style="text-align:right;">
8.7593547
</td>
<td style="text-align:right;">
2.9596207
</td>
<td style="text-align:right;">
-1.5578
</td>
<td style="text-align:right;">
10.0437
</td>
<td style="text-align:right;">
0.1516814
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
3.485335
</td>
<td style="text-align:right;">
4.2329425
</td>
<td style="text-align:right;">
2.0574116
</td>
<td style="text-align:right;">
-0.5471
</td>
<td style="text-align:right;">
7.5178
</td>
<td style="text-align:right;">
0.0902579
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
11.881071
</td>
<td style="text-align:right;">
3.4080300
</td>
<td style="text-align:right;">
1.8460850
</td>
<td style="text-align:right;">
8.2628
</td>
<td style="text-align:right;">
15.4993
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
4.945589
</td>
<td style="text-align:right;">
41.1813216
</td>
<td style="text-align:right;">
6.4172675
</td>
<td style="text-align:right;">
-7.6320
</td>
<td style="text-align:right;">
17.5232
</td>
<td style="text-align:right;">
0.4409031
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
11.089110
</td>
<td style="text-align:right;">
6.3564564
</td>
<td style="text-align:right;">
2.5212014
</td>
<td style="text-align:right;">
6.1476
</td>
<td style="text-align:right;">
16.0306
</td>
<td style="text-align:right;">
0.0000109
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
1.689869
</td>
<td style="text-align:right;">
0.9935718
</td>
<td style="text-align:right;">
0.9967807
</td>
<td style="text-align:right;">
-0.2638
</td>
<td style="text-align:right;">
3.6435
</td>
<td style="text-align:right;">
0.0900134
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
10.114451
</td>
<td style="text-align:right;">
1.4270177
</td>
<td style="text-align:right;">
1.1945785
</td>
<td style="text-align:right;">
7.7731
</td>
<td style="text-align:right;">
12.4558
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

When W3 is 1 the truth effect is 11, our estimates are 11 with CI
coverage. When W3 is 0 the truth is 1, our estimate is 1.9 with CIs that
cover the truth as well.

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
2.701885
</td>
<td style="text-align:right;">
0.7471848
</td>
<td style="text-align:right;">
0.8643985
</td>
<td style="text-align:right;">
1.0077
</td>
<td style="text-align:right;">
4.3961
</td>
<td style="text-align:right;">
0.0036597
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
10.236543
</td>
<td style="text-align:right;">
0.3365067
</td>
<td style="text-align:right;">
0.5800920
</td>
<td style="text-align:right;">
9.0996
</td>
<td style="text-align:right;">
11.3735
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
11.832085
</td>
<td style="text-align:right;">
0.3484322
</td>
<td style="text-align:right;">
0.5902815
</td>
<td style="text-align:right;">
10.6752
</td>
<td style="text-align:right;">
12.9890
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
-1.106343
</td>
<td style="text-align:right;">
0.7334625
</td>
<td style="text-align:right;">
0.8564242
</td>
<td style="text-align:right;">
-2.7849
</td>
<td style="text-align:right;">
0.5722
</td>
<td style="text-align:right;">
0.2318963
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
2.701885
</td>
<td style="text-align:right;">
0.7471848
</td>
<td style="text-align:right;">
0.8643985
</td>
<td style="text-align:right;">
1.0077
</td>
<td style="text-align:right;">
4.3961
</td>
<td style="text-align:right;">
0.0036597
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
10.236543
</td>
<td style="text-align:right;">
0.3365067
</td>
<td style="text-align:right;">
0.5800920
</td>
<td style="text-align:right;">
9.0996
</td>
<td style="text-align:right;">
11.3735
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
11.832085
</td>
<td style="text-align:right;">
0.3484322
</td>
<td style="text-align:right;">
0.5902815
</td>
<td style="text-align:right;">
10.6752
</td>
<td style="text-align:right;">
12.9890
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
-1.106343
</td>
<td style="text-align:right;">
0.7334625
</td>
<td style="text-align:right;">
0.8564242
</td>
<td style="text-align:right;">
-2.7849
</td>
<td style="text-align:right;">
0.5722
</td>
<td style="text-align:right;">
0.2318963
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
Oleg Sofrygin. 2022. “<span class="nocase">sl3</span>: Modern Machine
Learning Pipelines for Super Learning.”
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-coyle-hal9001-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars W van der Laan,
and Mark J van der Laan. 2022. “<span class="nocase">hal9001</span>: The
Scalable Highly Adaptive Lasso.”
<https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-Diaz2020a" class="csl-entry">

Díaz, Iván, and Nima S. Hejazi. 2020. “<span class="nocase">Causal
mediation analysis for stochastic interventions</span>.” *Journal of the
Royal Statistical Society. Series B: Statistical Methodology* 82 (3):
661–83. <https://doi.org/10.1111/rssb.12362>.

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

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020.
“<span class="nocase">hal9001</span>: Scalable Highly Adaptive Lasso
Regression in R.” *Journal of Open Source Software* 5 (53): 2526.
<https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-mccoy2023semiparametric" class="csl-entry">

McCoy, David B., Alan E. Hubbard, Alejandro Schuler, and Mark J. van der
Laan. 2023. “Semi-Parametric Identification and Estimation of
Interaction and Effect Modification in Mixed Exposures Using Stochastic
Interventions.” <https://arxiv.org/abs/2305.01849>.

</div>

</div>
