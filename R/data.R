#' Shaquille O'Neal free throw data
#'
#' Example data from Boos and Stefanski (2013) on Shaquille O'Neal free throws
#' in the 2000 NBA playoffs (Table 7.1 on pg 324).
#'
#' @format A data frame with 23 rows and 3 columns:
#' \describe{
#'   \item{game}{Game number}
#'   \item{ft_success}{Free throws made during game}
#'   \item{ft_attempt}{Free throws attempted during game}
#' }
#'
#' @references
#' Boos DD, & Stefanski LA. (2013). M-estimation (estimating equations). In
#' Essential Statistical Inference (pp. 297-337). Springer, New York, NY.
"shaq_free_throws"

#' Inderjit dose-response data
#'
#' Example data from Inderjit et al. (2002) on the dose-response of herbicide
#' on perennial ryegrass growth.
#'
#' @format A data frame with 24 rows and 2 columns:
#' \describe{
#'   \item{response}{Ryegrass root length}
#'   \item{dose}{Herbicide dose}
#' }
#'
#' @references
#' Inderjit, Streibig JC, & Olofsdotter M. (2002). Joint action of phenolic
#' acid mixtures and its significance in allelopathy research. *Physiologia
#' Plantarum*, 114(3), 422-428.
"inderjit"

#' Robust regression example data
#'
#' Illustrative example of robust linear regression from Zivich et al. (2022).
#'
#' @format A data frame with 15 rows and 3 columns:
#' \describe{
#'   \item{height}{Height measurement}
#'   \item{weight}{Weight measurement (with outlier induced at observation 9)}
#'   \item{weight_no_outlier}{Weight measurement (without outlier)}
#' }
#'
#' @references
#' Zivich PN, Klose M, Cole SR, Edwards JK, & Shook-Sa BE. (2022).
#' Delicatessen: M-estimation in Python. *arXiv:2203.11300*.
"robust_regress"

#' Breast cancer survival data
#'
#' Survival data from Collett (2015) on the survival times of 45 women with
#' breast cancer in Middlesex Hospital July 1987 (Table 1.2).
#'
#' @format A data frame with 45 rows and 3 columns:
#' \describe{
#'   \item{delta}{Event indicator (1 = event, 0 = censored)}
#'   \item{times}{Observation time}
#'   \item{stain}{Tumor HPA stain indicator (1 = positive, 0 = negative)}
#' }
#'
#' @references
#' Collett, D. (2015). Survival analysis. In Modelling survival data in medical
#' research. 3rd Ed. Chapman and Hall/CRC. pg 6-7
"breast_cancer"

# ── Applied example datasets ────────────────────────────────────────────────

#' Mroz labor force participation data
#'
#' Data from Mroz (1987) on married women's labor force participation, used to
#' demonstrate OLS, 2SLS, and instrumental variable estimation.
#'
#' @format A data frame with 753 rows.
#'
#' @references
#' Mroz TA. (1987). The sensitivity of an empirical model of married women's
#' hours of work to economic and statistical assumptions. *Econometrica*,
#' 55(4), 765-799.
"mroz"

#' US state crime data
#'
#' Data on violent crime rates by US state from the 2005 Statistical Abstract
#' of the United States, used to demonstrate robust regression.
#'
#' @format A data frame with 51 rows (50 states + DC).
#'
#' @references
#' Agresti A & Finlay B. (2009). *Statistical Methods for the Social Sciences*,
#' 4th Edition. Pearson.
"crime"

#' Collett bladder cancer recurrence data
#'
#' Bladder cancer recurrence data from Collett (2015), used to demonstrate
#' pooled logistic regression for survival analysis.
#'
#' @format A data frame with 85 rows and 6 columns:
#' \describe{
#'   \item{patient}{Patient ID}
#'   \item{time}{Follow-up time}
#'   \item{delta}{Event indicator}
#'   \item{treat}{Treatment group}
#'   \item{init}{Number of initial tumors}
#'   \item{size}{Size of largest initial tumor}
#' }
#'
#' @references
#' Collett, D. (2015). Modelling survival data in medical research. 3rd Ed.
#' Chapman and Hall/CRC.
"collett_bladder"

#' Cutler (1995) pharmacodynamic data
#'
#' Pharmacodynamic data from Cutler et al. (1995) on acetylcholinesterase
#' inhibition, used to demonstrate E-max dose-response models.
#'
#' @format A data frame.
#'
#' @references
#' Cutler NR et al. (1995). Dose-dependent CSF acetylcholinesterase inhibition
#' by SDZ ENA 713 in Alzheimer's disease. *Acta Neurologica Scandinavica*,
#' 91(5), 347-351.
"cutler1995"

#' Bonate adverse events data
#'
#' Adverse events case study data from Bonate (2011), used to demonstrate
#' generalized linear models (logistic, probit, complementary log-log).
#'
#' @format A data frame.
#'
#' @references
#' Bonate PL. (2011). *Pharmacokinetic-Pharmacodynamic Modeling and
#' Simulation*. 2nd Ed. Springer.
"bonate_adverse"

#' Lau WIHS HIV/CD4 data
#'
#' CD4 T cell count data from the Women's Interagency HIV Study (WIHS),
#' used to demonstrate sensitivity analysis for missing data.
#'
#' @format A data frame with columns: id, black, age, cd4, cd41-cd44.
#'
#' @references
#' Cole SR, Zivich PN, Edwards JK, Shook-Sa BE, & Hudgens MG. (2023).
#' Sensitivity analyses for means or proportions with missing outcome data.
#' *Epidemiology*, 34(5), 645-651.
"lau_wihs"

#' SDSS quasar survey data
#'
#' Quasar photometric data from the Sloan Digital Sky Survey 3rd Data Release
#' (n=46,420 quasars), used to demonstrate generalized additive models.
#'
#' @format A data frame with columns for redshift (z) and magnitude bands.
#'
#' @references
#' Weinstein MA, et al. (2004). An empirical algorithm for broadband
#' photometric redshifts of quasars from the Sloan Digital Sky Survey.
#' *The Astrophysical Journal Supplement Series*, 155(2), 243.
"sdss_quasar"

#' NSDUH substance use survey data
#'
#' Data from the National Survey on Drug Use and Health, used to demonstrate
#' causal mediation analysis.
#'
#' @format A data frame with approximately 40,000 rows.
#'
#' @references
#' Coffman DL, Zhong W. (2021). Assessing mediation using marginal structural
#' models in the presence of confounding and moderation. *Psychological
#' Methods*, 17(4), 532.
"nsduh"

#' ACTG 175 clinical trial data
#'
#' Data from AIDS Clinical Trials Group Study 175, used to demonstrate
#' confidence bands for treatment effects.
#'
#' @format A data frame.
#'
#' @references
#' Hammer SM, et al. (1996). A trial comparing nucleoside monotherapy with
#' combination therapy in HIV-infected adults with CD4 cell counts from 200
#' to 500 per cubic millimeter. *New England Journal of Medicine*, 335(15),
#' 1081-1090.
"actg175"

#' PPARg cheminformatics data
#'
#' Virtual screening data for the PPARg receptor, used to demonstrate hit
#' enrichment curves and confidence bands.
#'
#' @format A data frame.
#'
#' @references
#' Ash JR & Hughes-Oliver JM. (2022). Confidence bands for hit enrichment
#' curves. *Journal of Cheminformatics*, 14(1), 1-15.
"pparg"

#' GetTested randomized trial data
#'
#' Data from the GetTested randomized trial, used to demonstrate
#' standardization and inverse probability weighting for treatment effects.
#'
#' @format A data frame.
#'
#' @references
#' Morris TP, Walker AS, Williamson EJ, & White IR. (2022). Planning a method
#' for covariate adjustment in randomised trials. *BMJ*, 376.
"get_tested"
