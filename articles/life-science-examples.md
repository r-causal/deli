# Zivich et al. (2022): Life-Science Examples

> **Note**
>
> This article is translated from the [Zivich et al. (2022):
> Life-Science
> Examples](https://deli.readthedocs.io/en/latest/Examples/LifeScienceExamples.html)
> in the documentation of [delicatessen](https://deli.readthedocs.io/),
> deli’s Python counterpart.

``` r

library(deli)
```

The following replicates the life-science case studies provided in
Zivich et al. (2022). The original paper accompanies the Delicatessen
Python library and demonstrates M-estimation for three applied examples:
linear vs robust regression, dose-response modeling, and standardization
to an external population.

## Case Study 1: Linear vs Robust Regression

Height and weight data with a simulated outlier are used to compare
simple linear regression to robust linear regression. Three models are
fit: (a) linear regression without the outlier (benchmark), (b) linear
regression with the outlier, and (c) robust linear regression with the
outlier using Huber loss (k = 1.345).

### Data

``` r

# Height is the same for both scenarios
x <- robust_regress$height

# Weight without outlier (benchmark)
y_no <- robust_regress$weight_no_outlier

# Weight with outlier (observation 9 has +3 added)
y_out <- robust_regress$weight

# Design matrix (intercept + height)
X <- cbind(1, x)
```

### (a) Linear Regression without Outlier (Benchmark)

``` r

estr_bench <- m_estimate(
  stacked_equations = function(theta) {
    ee_regression(theta, X = cbind(1, x), y = y_no, model = "linear")
  },
  init = c(0, 0)
)

data.frame(
  Param = c("Intercept", "Height"),
  Coef = round(estr_bench@theta, 3),
  LCL = round(confint(estr_bench)[, 1], 3),
  UCL = round(confint(estr_bench)[, 2], 3)
)
#>             Param    Coef     LCL     UCL
#> theta_1 Intercept -53.740 -80.816 -26.665
#> theta_2    Height   0.721   0.559   0.883
```

### (b) Linear Regression with Outlier

``` r

estr_outlier <- m_estimate(
  stacked_equations = function(theta) {
    ee_regression(theta, X = X, y = y_out, model = "linear")
  },
  init = c(0, 0)
)

data.frame(
  Param = c("Intercept", "Height"),
  Coef = round(estr_outlier@theta, 3),
  LCL = round(confint(estr_outlier)[, 1], 3),
  UCL = round(confint(estr_outlier)[, 2], 3)
)
#>             Param    Coef     LCL    UCL
#> theta_1 Intercept -19.407 -80.527 41.714
#> theta_2    Height   0.515   0.147  0.882
```

### (c) Robust Regression with Outlier (Huber, k=1.345)

Robust linear regression reduces the influence of the outlier. We use
Huber loss with tuning constant k = 1.345.

``` r

estr_robust <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y_out,
                         model = "linear", k = 1.345)
  },
  # Use OLS estimates as initial values for stability
  init = estr_outlier@theta
)

data.frame(
  Param = c("Intercept", "Height"),
  Coef = round(estr_robust@theta, 3),
  LCL = round(confint(estr_robust)[, 1], 3),
  UCL = round(confint(estr_robust)[, 2], 3)
)
#>             Param    Coef      LCL    UCL
#> theta_1 Intercept -36.791 -146.312 72.729
#> theta_2    Height   0.619   -0.037  1.276
```

The robust regression estimates are closer to the benchmark (without
outlier) than the simple linear regression with the outlier,
demonstrating how Huber loss down-weights the influence of extreme
observations.

## Case Study 2: Dose-Response (Log-Logistic)

Dose-response data from Inderjit et al. (2002) on the effect of ferulic
acid on ryegrass root length. A 3-parameter log-logistic model is fit
(the lower limit is fixed at zero) along with an estimating equation for
the 20% effective dose (ED20).

### Data

``` r

# Load Inderjit data from the deli package
response <- inderjit$response
dose <- inderjit$dose
```

### 3-Parameter Log-Logistic + ED20

The lower limit is fixed at zero. We estimate three log-logistic
parameters (upper limit, ED50, steepness) and the ED20 simultaneously.

``` r

estr_dr <- m_estimate(
  stacked_equations = function(theta) {
    lower_limit <- 0

    # 3-parameter log-logistic (lower limit fixed at 0)
    pl3 <- ee_loglogistic(c(lower_limit, theta[1:3]),
                          dose = dose, response = response)

    # Effective dose at 20%
    ed20 <- ee_loglogistic_ed(theta[4], dose = dose, delta = 0.20,
                               lower = lower_limit, upper = theta[1],
                               steepness = theta[3], ed50 = theta[2])

    # Drop the first row (fixed lower limit) and stack
    rbind(pl3[2:4, , drop = FALSE], ed20)
  },
  init = c(8, 3, 2, 2),
  solver = "nleqslv"
)

data.frame(
  Param = c("Upper", "ED50", "Steepness", "ED20"),
  Coef = round(estr_dr@theta, 3),
  LCL = round(confint(estr_dr)[, 1], 3),
  UCL = round(confint(estr_dr)[, 2], 3)
)
#>             Param  Coef   LCL   UCL
#> theta_1     Upper 7.855 7.554 8.157
#> theta_2      ED50 3.263 2.743 3.784
#> theta_3 Steepness 2.470 1.897 3.043
#> theta_4      ED20 1.862 1.581 2.143
```

The first parameter is the upper limit of the dose-response curve, the
second is the ED50, the third is the steepness, and the fourth is the
ED20. As the dose of ferulic acid increases, root length decreases.

## Case Study 3: Standardization to External Population

> **Note**
>
> This example requires the Kamat et al. (2012) biomarker CSV data,
> which is not included in the repository. The code below is shown but
> not executed.

Inverse odds weights are used to standardize biomarker results from
Kamat et al. (2012) to the Women’s Interagency HIV Study (WIHS)
population. This approach generalizes study results beyond a convenience
sample by reweighting based on drug use prevalence.

### Data

``` r

# Load Kamat et al. biomarker data
d1 <- read.csv(file.path(getwd(), "..", "docs", "Examples", "data",
                          "kamat.et.al.2012_biomarkers.csv"))
d1$drug_use <- ifelse(d1$Cocaine + d1$Opiate > 0, 1, 0)
d1$S <- 1
biomarkers <- c("IFN_alpha", "CXCL9", "CXCL10", "sIL.2R", "IL12")
d1 <- d1[, c("drug_use", "S", biomarkers)]

# Log-transform biomarkers
for (bm in biomarkers) {
  d1[[bm]] <- log(d1[[bm]])
}

# WIHS population: drug use prevalence
d0 <- data.frame(
  drug_use = c(rep(1, 300), rep(0, 4016 - 300)),
  S = 0
)

# Stack datasets
d <- rbind(
  cbind(d0, setNames(as.data.frame(matrix(NA, nrow(d0), 5)), biomarkers)),
  d1
)
d$constant <- 1

# Fill missing biomarkers with placeholder
for (bm in biomarkers) {
  d[[bm]] <- ifelse(is.na(d[[bm]]), 9999, d[[bm]])
}
```

### Naive Means (Study Population)

First, compute simple means of the log-transformed biomarkers in the
study sample (Kamat et al.). These may not generalize to the target
population.

``` r

# Subset to study data only
d1_r <- d[d$S == 1, ]

estr_naive <- m_estimate(

  stacked_equations = function(theta) {
    rbind(
      matrix(d1_r$IFN_alpha - theta[1], nrow = 1),
      matrix(d1_r$CXCL9    - theta[2], nrow = 1),
      matrix(d1_r$CXCL10   - theta[3], nrow = 1),
      matrix(d1_r$sIL.2R   - theta[4], nrow = 1),
      matrix(d1_r$IL12     - theta[5], nrow = 1)
    )
  },
  init = rep(1, 5)
)

data.frame(
  Biomarker = c("IFN-a", "CXCL9", "CXCL10", "sIL-2R", "IL-12"),
  Mean = round(estr_naive@theta, 3),
  LCL = round(confint(estr_naive)[, 1], 3),
  UCL = round(confint(estr_naive)[, 2], 3)
)
```

### Standardized Means (IOSW)

Now generalize the biomarker means to the WIHS population using inverse
odds of sampling weights. A logistic regression model for study
membership (S) given drug use is fit, and inverse odds weights reweight
the study sample to match the target population’s drug use distribution.

``` r

# Design matrix and study indicator
x_mat <- as.matrix(d[, c("constant", "drug_use")])
s <- d$S

estr_iosw <- m_estimate(
  stacked_equations = function(theta) {
    # Logistic regression for study membership
    nuisance <- ee_regression(theta[1:2], X = x_mat, y = s,
                               model = "logistic")

    # Inverse odds weights
    pi_s <- inverse_logit(as.numeric(x_mat %*% theta[1:2]))
    wt <- ifelse(s == 1, (1 - pi_s) / pi_s, 0)

    # Weighted biomarker means
    ee_ifn   <- matrix(s * wt * (d$IFN_alpha - theta[3]), nrow = 1)
    ee_cxcl9 <- matrix(s * wt * (d$CXCL9    - theta[4]), nrow = 1)
    ee_cxcl10 <- matrix(s * wt * (d$CXCL10  - theta[5]), nrow = 1)
    ee_sil2r <- matrix(s * wt * (d$sIL.2R   - theta[6]), nrow = 1)
    ee_il12  <- matrix(s * wt * (d$IL12     - theta[7]), nrow = 1)

    rbind(nuisance, ee_ifn, ee_cxcl9, ee_cxcl10, ee_sil2r, ee_il12)
  },
  init = c(0, 0, rep(1, 5))
)

data.frame(
  Biomarker = c("IFN-a", "CXCL9", "CXCL10", "sIL-2R", "IL-12"),
  Mean = round(estr_iosw@theta[3:7], 3),
  LCL = round(confint(estr_iosw)[3:7, 1], 3),
  UCL = round(confint(estr_iosw)[3:7, 2], 3)
)
```

After standardization, some biomarker means shift noticeably
(particularly IL-12 and sIL-2R), illustrating the importance of
accounting for differences in drug use prevalence when generalizing
results.

## References

Inderjit, Streibig JC, & Olofsdotter M. (2002). Joint action of phenolic
acid mixtures and its significance in allelopathy research. *Physiologia
Plantarum*, 114(3), 422-428.

Kamat A, et al. (2012). A Plasma Biomarker Signature of Immune
Activation in HIV Patients on Antiretroviral Therapy. *PLoS ONE*, 7(2),
e30881.

Ritz C, Baty F, Streibig JC, & Gerhard D. (2015). Dose-Response Analysis
Using R. *PLoS ONE*, 10(12), e0146021.

Zivich PN, Klose M, Cole SR, Edwards JK, & Shook-Sa BE. (2022).
Delicatessen: M-Estimation in Python. *arXiv preprint arXiv:2203.11300*.
