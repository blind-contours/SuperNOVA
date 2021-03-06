% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SuperNOVA.R
\name{SuperNOVA}
\alias{SuperNOVA}
\title{Data-Adaptive Efficient Estimation of Interactions and Effect Modification using Stochastic Shift Interventions}
\usage{
SuperNOVA(
  W,
  V,
  A,
  Y,
  delta,
  LOD_vals = NULL,
  estimator = c("tmle", "onestep"),
  fluctuation = c("standard", "weighted"),
  max_iter = 10,
  sl_density_lrnr,
  Q1_stack,
  Q2_stack,
  n_folds,
  family,
  quantile_thresh,
  verbose,
  parallel
)
}
\arguments{
\item{W}{A \code{matrix}, \code{data.frame}, or similar containing a set of
baseline covariates.}

\item{V}{\code{matrix}, \code{data.frame}, or similar containing a set of
baseline covariates considered to be possible effect modifiers of exposure}

\item{A}{\code{matrix}, \code{data.frame}, or similar containing an individual or set of
exposures}

\item{Y}{\code{numeric} vector of the observed outcomes.}

\item{delta}{A \code{numeric} value indicating the shift in the exposures to
be used in defining the target parameter. This is defined with respect to
the scale of the exposures (A).}

\item{estimator}{The type of estimator to be fit, either \code{"tmle"} for
targeted maximum likelihood or \code{"onestep"} for a one-step estimator.}

\item{fluctuation}{The method to be used in the submodel fluctuation step
(targeting step) to compute the TML estimator. The choices are "standard"
and "weighted" for where to place the auxiliary covariate in the logistic
tilting regression.}

\item{max_iter}{A \code{numeric} integer giving the maximum number of steps
to be taken in iterating to a solution of the efficient influence function.}

\item{sl_density_lrnr}{Learners used to fit Super Learner ensembles to densities via \pkg{sl3}}

\item{Q1_stack}{Learners used to fit Super Learner ensembles to the outcome model via \pkg{sl3}}

\item{n_folds}{Number of folds to use in cross-validation}

\item{family}{Outcome type family}

\item{quantile_thresh}{Threshold based on quantiles of the f-statistic used to identify "important" basis functions in the data-adaptive procedure}

\item{verbose}{Whether to run verbosely}

\item{parallel}{TRUE/FALSE parallelize across cores}
}
\value{

}
\description{
Data-Adaptive Efficient Estimation of Interactions and Effect Modification using Stochastic Shift Interventions
}
\details{
Treat variables identified in the same basis functions in ensemble b-spline models as
a data-adaptive parameter. Given variable sets identified, depending on if basis functions contain
variables for \code{A} or \code{AV} construct targeted minimum loss estimations of
the counterfactual mean differences under various target parameters for individual variables, effect modifying variables,
or interacting variables. Ensemble machine learning may be used to construct the initial
estimates of nuisance functions using \pkg{sl3}.
}
