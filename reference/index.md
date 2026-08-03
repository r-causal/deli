# Package index

## Core

Primary estimator classes and estimation interface.

- [`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
  : One-step M-estimation
- [`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
  : One-step GMM estimation
- [`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md)
  : M-Estimator
- [`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md)
  : GMM Estimator
- [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
  : Estimate parameters and sandwich variance

## Statistical Inference

Tools for statistical inference from fitted estimators.

- [`deli-generics`](https://r-causal.github.io/deli/reference/deli-generics.md)
  : Standard S3 generics for deli estimators
- [`deli-display`](https://r-causal.github.io/deli/reference/deli-display.md)
  : Display methods for deli estimators
- [`deli-tidiers`](https://r-causal.github.io/deli/reference/deli-tidiers.md)
  : Broom tidiers for deli estimators
- [`reexports`](https://r-causal.github.io/deli/reference/reexports.md)
  [`tidy`](https://r-causal.github.io/deli/reference/reexports.md)
  [`glance`](https://r-causal.github.io/deli/reference/reexports.md)
  [`augment`](https://r-causal.github.io/deli/reference/reexports.md) :
  Objects exported from other packages
- [`confidence_intervals()`](https://r-causal.github.io/deli/reference/confidence_intervals.md)
  : Confidence intervals for M-Estimator parameters
- [`z_scores()`](https://r-causal.github.io/deli/reference/z_scores.md)
  : Z-scores for M-Estimator parameters
- [`p_values()`](https://r-causal.github.io/deli/reference/p_values.md)
  : P-values for M-Estimator parameters
- [`s_values()`](https://r-causal.github.io/deli/reference/s_values.md)
  : S-values (surprisal) for M-Estimator parameters
- [`confidence_bands()`](https://r-causal.github.io/deli/reference/confidence_bands.md)
  : Confidence bands for parameter vectors
- [`compute_confidence_bands()`](https://r-causal.github.io/deli/reference/compute_confidence_bands.md)
  : Compute confidence bands from theta and covariance
- [`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md)
  : Delta method for variance of transformed parameters
- [`influence_functions()`](https://r-causal.github.io/deli/reference/influence_functions.md)
  : Influence functions for M-Estimator
- [`compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.md)
  : Compute the empirical sandwich variance estimator

## Basic Estimating Equations

Estimating equations for means, variances, and percentiles.

- [`ee_mean()`](https://r-causal.github.io/deli/reference/ee_mean.md) :
  Estimating equation for the mean
- [`ee_mean_variance()`](https://r-causal.github.io/deli/reference/ee_mean_variance.md)
  : Estimating equations for the mean and variance
- [`ee_mean_geometric()`](https://r-causal.github.io/deli/reference/ee_mean_geometric.md)
  : Estimating equation for the geometric mean
- [`ee_mean_robust()`](https://r-causal.github.io/deli/reference/ee_mean_robust.md)
  : Estimating equation for the robust mean
- [`ee_percentile()`](https://r-causal.github.io/deli/reference/ee_percentile.md)
  : Estimating equation for the percentile
- [`ee_positive_mean_deviation()`](https://r-causal.github.io/deli/reference/ee_positive_mean_deviation.md)
  : Estimating equations for the positive mean deviation

## Regression Estimating Equations

Estimating equations for generalized linear models and penalized
regression.

- [`ee_regression()`](https://r-causal.github.io/deli/reference/ee_regression.md)
  : Estimating equation for regression
- [`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md) :
  Estimating equation for generalized linear models
- [`ee_robust_regression()`](https://r-causal.github.io/deli/reference/ee_robust_regression.md)
  : Estimating equation for robust regression
- [`ee_ridge_regression()`](https://r-causal.github.io/deli/reference/ee_ridge_regression.md)
  : Estimating equation for ridge regression
- [`ee_lasso_regression()`](https://r-causal.github.io/deli/reference/ee_lasso_regression.md)
  : Estimating equation for approximate LASSO regression
- [`ee_dlasso_regression()`](https://r-causal.github.io/deli/reference/ee_dlasso_regression.md)
  : Estimating equation for differentiable LASSO regression
- [`ee_elasticnet_regression()`](https://r-causal.github.io/deli/reference/ee_elasticnet_regression.md)
  : Estimating equation for elastic net regression
- [`ee_bridge_regression()`](https://r-causal.github.io/deli/reference/ee_bridge_regression.md)
  : Estimating equation for bridge penalized regression
- [`ee_beta_regression()`](https://r-causal.github.io/deli/reference/ee_beta_regression.md)
  : Estimating equation for beta regression
- [`ee_tobit()`](https://r-causal.github.io/deli/reference/ee_tobit.md)
  : Estimating equation for Tobit regression (Type I)
- [`ee_mlogit()`](https://r-causal.github.io/deli/reference/ee_mlogit.md)
  : Estimating equation for multinomial logistic regression
- [`ee_additive_regression()`](https://r-causal.github.io/deli/reference/ee_additive_regression.md)
  : Estimating equation for additive regression (GAM)

## Causal Estimating Equations

Estimating equations for causal inference.

- [`ee_gformula()`](https://r-causal.github.io/deli/reference/ee_gformula.md)
  : Estimating equations for the g-formula (g-computation)
- [`ee_ipw()`](https://r-causal.github.io/deli/reference/ee_ipw.md) :
  Estimating equations for inverse probability weighting (IPW)
- [`ee_aipw()`](https://r-causal.github.io/deli/reference/ee_aipw.md) :
  Estimating equations for augmented inverse probability weighting
  (AIPW)
- [`ee_ipw_msm()`](https://r-causal.github.io/deli/reference/ee_ipw_msm.md)
  : Estimating equations for IPW marginal structural model
- [`ee_gestimation_snmm()`](https://r-causal.github.io/deli/reference/ee_gestimation_snmm.md)
  : Estimating equations for g-estimation of structural nested mean
  models
- [`ee_iv_causal()`](https://r-causal.github.io/deli/reference/ee_iv_causal.md)
  : Estimating equations for instrumental variable (IV) estimation
- [`ee_2sls()`](https://r-causal.github.io/deli/reference/ee_2sls.md) :
  Estimating equations for Two-Stage Least Squares (2SLS)
- [`ee_mean_sensitivity_analysis()`](https://r-causal.github.io/deli/reference/ee_mean_sensitivity_analysis.md)
  : Estimating equations for weighted sensitivity analysis of the mean

## Survival Estimating Equations

Estimating equations for survival analysis.

- [`ee_aft()`](https://r-causal.github.io/deli/reference/ee_aft.md) :
  Estimating equation for accelerated failure time models
- [`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md)
  : Estimating equation for pooled logistic regression
- [`ee_survival_model()`](https://r-causal.github.io/deli/reference/ee_survival_model.md)
  : Estimating equation for parametric survival models

## Measurement Estimating Equations

Estimating equations for measurement error correction.

- [`ee_rogan_gladen()`](https://r-causal.github.io/deli/reference/ee_rogan_gladen.md)
  : Estimating equation for Rogan-Gladen correction
- [`ee_rogan_gladen_extended()`](https://r-causal.github.io/deli/reference/ee_rogan_gladen_extended.md)
  : Estimating equation for extended Rogan-Gladen correction
- [`ee_regression_calibration()`](https://r-causal.github.io/deli/reference/ee_regression_calibration.md)
  : Estimating equation for regression calibration

## Pharmacokinetics Estimating Equations

Estimating equations for dose-response models.

- [`ee_emax()`](https://r-causal.github.io/deli/reference/ee_emax.md) :
  Estimating equation for E-max dose-response model
- [`ee_emax_ed()`](https://r-causal.github.io/deli/reference/ee_emax_ed.md)
  : Estimating equation for delta-effective dose (E-max)
- [`ee_loglogistic()`](https://r-causal.github.io/deli/reference/ee_loglogistic.md)
  : Estimating equation for 4-parameter log-logistic dose-response model
- [`ee_loglogistic_ed()`](https://r-causal.github.io/deli/reference/ee_loglogistic_ed.md)
  : Estimating equation for delta-effective dose (log-logistic)

## Predictions

Prediction helpers for fitted models.

- [`deli-predict`](https://r-causal.github.io/deli/reference/deli-predict.md)
  : Predictions from a fitted deli estimator
- [`deli-augment`](https://r-causal.github.io/deli/reference/deli-augment.md)
  : Augment data with predictions from a fitted deli estimator
- [`regression_predictions()`](https://r-causal.github.io/deli/reference/regression_predictions.md)
  : Generate predicted values from a regression model
- [`survival_predictions()`](https://r-causal.github.io/deli/reference/survival_predictions.md)
  : Generate predicted survival measures from a parametric survival
  model
- [`aft_predictions_individual()`](https://r-causal.github.io/deli/reference/aft_predictions_individual.md)
  : Predicted survival measures from an AFT model
- [`aft_predictions_function()`](https://r-causal.github.io/deli/reference/aft_predictions_function.md)
  : Function-level predicted survival measures from an AFT model
- [`plogit_predict()`](https://r-causal.github.io/deli/reference/plogit_predict.md)
  : Predicted survival measures from a pooled logistic regression model
- [`convert_survival_measures()`](https://r-causal.github.io/deli/reference/convert_survival_measures.md)
  : Convert between survival analysis measures

## Transformations and Math

Transformations and special functions that mirror base R counterparts
such as [`plogis()`](https://rdrr.io/r/stats/Logistic.html),
[`qlogis()`](https://rdrr.io/r/stats/Logistic.html),
[`pnorm()`](https://rdrr.io/r/stats/Normal.html),
[`dnorm()`](https://rdrr.io/r/stats/Normal.html), and
[`psigamma()`](https://rdrr.io/r/base/Special.html). Several are
required in place of those counterparts under exact differentiation,
which the base R versions do not support.

- [`logit()`](https://r-causal.github.io/deli/reference/logit.md) :
  Logistic transformation
- [`inverse_logit()`](https://r-causal.github.io/deli/reference/inverse_logit.md)
  : Inverse logistic transformation
- [`identity_transform()`](https://r-causal.github.io/deli/reference/identity_transform.md)
  : Identity transformation
- [`standard_normal_cdf()`](https://r-causal.github.io/deli/reference/standard_normal_cdf.md)
  : Standard normal CDF
- [`standard_normal_pdf()`](https://r-causal.github.io/deli/reference/standard_normal_pdf.md)
  : Standard normal PDF
- [`deli_polygamma()`](https://r-causal.github.io/deli/reference/deli_polygamma.md)
  : Polygamma function
- [`deli_digamma()`](https://r-causal.github.io/deli/reference/deli_digamma.md)
  : Digamma function

## Model Building Utilities

Helpers for building estimating equations and design matrices. None has
a base R counterpart that returns the same values.

- [`robust_loss_functions()`](https://r-causal.github.io/deli/reference/robust_loss_functions.md)
  : Robust loss function derivatives
- [`aggregate_efuncs()`](https://r-causal.github.io/deli/reference/aggregate_efuncs.md)
  : Aggregate estimating function contributions by group
- [`deli_spline()`](https://r-causal.github.io/deli/reference/deli_spline.md)
  : Generate polynomial spline basis terms
- [`additive_design_matrix()`](https://r-causal.github.io/deli/reference/additive_design_matrix.md)
  : Build an additive design matrix for GAMs

## Conditions

The condition classes deli’s errors and warnings carry, so that a caller
can answer one by class rather than by matching its message.

- [`deli-conditions`](https://r-causal.github.io/deli/reference/deli-conditions.md)
  : The condition classes deli raises

## Datasets

Bundled example datasets.

- [`shaq_free_throws`](https://r-causal.github.io/deli/reference/shaq_free_throws.md)
  : Shaquille O'Neal free throw data
- [`inderjit`](https://r-causal.github.io/deli/reference/inderjit.md) :
  Inderjit dose-response data
- [`robust_regress`](https://r-causal.github.io/deli/reference/robust_regress.md)
  : Robust regression example data
- [`breast_cancer`](https://r-causal.github.io/deli/reference/breast_cancer.md)
  : Breast cancer survival data
- [`mroz`](https://r-causal.github.io/deli/reference/mroz.md) : Mroz
  labor force participation data
- [`crime`](https://r-causal.github.io/deli/reference/crime.md) : US
  state crime data
- [`collett_bladder`](https://r-causal.github.io/deli/reference/collett_bladder.md)
  : Collett bladder cancer recurrence data
- [`cutler1995`](https://r-causal.github.io/deli/reference/cutler1995.md)
  : Cutler (1995) pharmacodynamic data
- [`bonate_adverse`](https://r-causal.github.io/deli/reference/bonate_adverse.md)
  : Bonate adverse events data
- [`lau_wihs`](https://r-causal.github.io/deli/reference/lau_wihs.md) :
  Lau WIHS HIV/CD4 data
- [`sdss_quasar`](https://r-causal.github.io/deli/reference/sdss_quasar.md)
  : SDSS quasar survey data
- [`nsduh`](https://r-causal.github.io/deli/reference/nsduh.md) : NSDUH
  substance use survey data
- [`actg175`](https://r-causal.github.io/deli/reference/actg175.md) :
  ACTG 175 clinical trial data
- [`pparg`](https://r-causal.github.io/deli/reference/pparg.md) : PPARg
  cheminformatics data
- [`get_tested`](https://r-causal.github.io/deli/reference/get_tested.md)
  : GetTested randomized trial data
