---
title: "`SuperNOVA`: Semi-Parametric Identification and Estimation of Interaction and Effect Modification in Mixed Exposures using Stochastic Interventions in `R`"
tags:
  - causal inference
  - machine learning
  - stochastic interventions
  - efficient estimation
  - targeted learning
  - mixed exposures
  - R
authors:
  - name: David McCoy
    orcid: 0000-0002-5515-6307
    affiliation: 1
  - name: Alan Hubbard
    orcid: 0000-0002-3769-0127
    affiliation: 2
  - name: Mark Van der Laan
    orcid: 0000-0003-1432-5511
    affiliation: 2
affiliations:
  - name: Division of Environmental Health Sciences, University of California, Berkeley
    index: 1
  - name: Department of Biostatistics, University of California, Berkeley
    index: 2
date: 05 January 2022
bibliography: paper.bib
---

# Summary

In environmental epidemiology studies, analysts are interested in understanding the joint impact of mixed exposures, i.e. a vector of exposures, on health outcomes while robustly adjusting for covariates. However, traditional statistical methods often make overly simplistic assumptions that lead to statistical quantities not be directly applicable to public policy decisions. For example, ultimately we are interested in causal questions such as how a decrease in exposures (like toxic chemicals) lead to a decrease in deleterious health outcomes (like cancer). Within the context of mixtures, we are also interested estimating statistical quantities of interaction and effect modification but doing so outside of a parametric context as nonlinear and non-additive relationships likely exist in mixed exposures. New methods in statistics are necessary to bridge this gap by delivering estimates on modified exposure policies. To address these limitations, the package SuperNOVA has been developed to use data-adaptive machine learning methods to identify the variables and variable sets that have the most explanatory power on an outcome of interest and applies non-parametric definitions of interaction and effect modification to these variable sets in a mixed exposure. 

# Statement of Need

In public health and medical research, it is crucial to have reliable and accurate methods for estimating the effects of treatments and interventions. However, traditional parametric models have limitations that can result in biased estimates and inconsistent findings. This is particularly true when dealing with complex exposure scenarios, such as mixed exposures or treatments. To address these limitations, there is a growing need for semi-parametric statistical methods that are both theoretically proven to converge to the truth with minimum bias and accessible to researchers. This is where SuperNOVA, an open-source R package, comes in. SuperNOVA provides a powerful framework for estimating modified treatment policies of a mixed exposure using stochastic interventions [@diaz2018stochastic] and target maximum likelihood estimation [@vdl2011targeted, @vdl2018targeted]. This software estimates non-parametric definitions of interaction and effect modification target parameters, reducing the risk of model bias. Furthermore, by being open source, SuperNOVA enables researchers who are path-dependent on parametric models to adopt these new methods more easily. This is especially important for researchers who may not have the technical expertise or resources to develop and implement these methods from scratch. SuperNOVA provides a much-needed solution to the limitations of traditional parametric models in public health and medical research by offering a powerful, interpretable and accessible framework for semi-parametric statistical analysis of joint impacts, interaction and effect modification. SuperNOVA can help drive more consistent findings and faster public health decisions by removing human bias due to model selection. 

# Background

The package SuperNOVA was developed to address the limitations of traditional statistical methods in environmental epidemiology studies. These traditional methods often make overly simplistic assumptions, such as linear and additive relationships, and the resulting statistical quantities may not be directly applicable to public policy decisions. SuperNOVA addresses these limitations by using data-adaptive machine learning methods to identify the variables and variable sets that have the most explanatory power on an outcome of interest. In the variable set discovery, the package builds a discrete Super Learner [@coyle-sl3-rpkg] which is a library of machine learning estimators which uses cross-validation to select the best fitting estimator. This Super Learner is composed of flexible basis function estimators, the best of which is analyzed using ANOVA style analysis to determine the variables that contribute most to the model fit through basis functions. The variable sets used in the basis functions drive the target parameters estimated. In the event of basis functions for an individual exposure $A$, the effects of an individual shift are estimated, for basis function with $A$ and $W$ (a baseline covariate), the effect modification parameter is estimated, which is an individual shift in a covariate region and if two exposures are included in a basis function $A_1, A_2$ the interaction target parameter is estimated, which is the expected outcome under dual shift of both exposures compared to the sum of expected outcomes given individual shifts independently. For each target parameter we use ensemble machine learning to ascertain the expected outcome under a shift and we use cross-validated targeted maximum likelihood estimation [@Hubbard2016] to debias our initial estimates thereby creating an asymptotically unbiased estimator with minimum variance. In this way, SuperNOVA allows analysts to explore modified treatment policies and ask causal questions (under assumptions) about the impact of mixed exposures on health outcomes. SuperNOVA uses V-fold cross-validation procedures to avoid over-fitting and incorrect model assumptions by creating parameter generating samples wherein the variable sets are determined and estimators for nuisance parameters are trained, an estimation sample is then used to estimate the target parameters of interest. Additionally, to avoid positivity violations (user inputs a shift amount that there isn't enough experimentation in the data to estimate) the shift amount can also be input as a data-adaptive parameter which finds the maximum shift possible for each exposure.

# `SuperNOVA`'s Scope

The SuperNOVA package is a novel tool for estimating the impact of a mixed exposure on an outcome of interest. It uses cross-validated targeted minimum loss-based estimation to guarantee consistency, efficiency, and robustness, despite using highly flexible basis function estimators in an ensemble. The package is based on prior work related to data-adaptive parameters and CV-TMLE, and builds on these concepts to create a unique and innovative approach for causal inference for interaction and effect modification in complex exposure situations.

The SuperNOVA software package is built for the R language and implements our proposed methodology. As input, SuperNOVA takes in variable sets $A$, $W$, $Y$ and a vector of deltas for each exposure in $A$. Default ensemble machine learning estimators are built using the sl3 package to estimate the nuisance parameters for each target parameter. The output of SuperNOVA is a dose-response analyses for variable sets data-adaptively identified in the mixed exposure. Estimates are for the expected outcome under a change in exposure compared to the observed outcome under the observed exposure. By making use of cross-validated targeted minimum loss-based estimation and ensemble basis function estimators, as well as its implementation of data-adaptive modified treatment policies, SuperNOVA is a valuable tool for researchers in many fields who need an interpretable and robust statistical approach to answer modified treatment policy questions.

`SuperNOVA` is designed to provide analysts with both V-fold specific and pooled results for stochastic intervention causal effects. It integrates with the [`sl3` package](https://github.com/tlverse/sl3) [@coyle2020sl3] to allow for ensemble machine learning to be leveraged in the estimation of nuisance parameters. 

SuperNOVA comes with synthetic mixtures data used by the NIEHS to assess new statistical methods for validiity. In this synthetic data, SuperNOVA identifies the correct marginal and interaction impacts and delivers expected outcome changes given a modification in these exposure sets. SuperNOVA also comes with real-world National Health and Nutrition Examination Survey Data (NHANES) of mixed metal exposure on telomere length and delivers interpretable outputs for what toxic metals are associated with telomere length. 

SuperNOVA is currently limited to the assessment of two-way interactions. This is due to the fact that in most public health setting there is not enough data to estimate the conditional joint density of three exposures without leading to positivity violations. Developments are also being made to include mediation target parameters as well.

# Availability

The `SuperNOVA` package has been made publicly available  [via GitHub](https://github.com/blind-contours/SuperNOVA). Use of the `SuperNOVA`package has been extensively documented in the package's `README` and a vignette. 


# Acknowledgments

David McCoy's contributions to this work were supported in part by Core E of the NIEHS Superfund Center at Berkeley funded by NIH grant P42ES004705.

# References

