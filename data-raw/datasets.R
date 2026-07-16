# Script to generate bundled datasets for the deli package
# Source: delicatessen/data.py

# Shaq free throws from Boos & Stefanski (2013), Table 7.1 pg 324
shaq_free_throws <- data.frame(
  game = 1:23,
  ft_success = c(4, 5, 5, 5, 2, 7, 6, 9, 4, 1, 13, 5, 6, 9, 7, 3, 8, 1, 18, 3, 10, 1, 3),
  ft_attempt = c(5, 11, 14, 12, 7, 10, 14, 15, 12, 4, 27, 17, 12, 9, 12, 10, 12, 6, 39, 13, 17, 6, 12)
)

# Inderjit et al. (2002) dose-response data
inderjit <- data.frame(
  response = c(7.58, 8.0, 8.3285714, 7.25, 7.375, 7.9625,
               8.3555556, 6.9142857, 7.75,
               6.8714286, 6.45, 5.9222222,
               1.925, 2.8857143, 4.2333333,
               1.1875, 0.8571429, 1.0571429,
               0.6875, 0.525, 0.825,
               0.25, 0.22, 0.44),
  dose = c(0, 0, 0, 0, 0, 0,
           0.94, 0.94, 0.94,
           1.88, 1.88, 1.88,
           3.75, 3.75, 3.75,
           7.5, 7.5, 7.5,
           15, 15, 15,
           30, 30, 30)
)

# Robust regression data from Zivich et al. (2022)
# Note: outlier=TRUE version (weight[9] += 3)
height <- c(168.519, 166.944, 164.327, 164.058, 166.212, 167.358,
            165.244, 169.352, 159.386, 166.953, 163.876,
            164.245, 165.061, 162.876, 164.185)
weight_outlier <- c(67.634, 67.418, 63.394, 66.18, 65.573, 67.66, 64.592,
                    68.4, 62.043 + 3, 67.093, 65.202, 64.328, 64.754,
                    62.179, 64.716)
weight_no_outlier <- c(67.634, 67.418, 63.394, 66.18, 65.573, 67.66, 64.592,
                       68.4, 62.043, 67.093, 65.202, 64.328, 64.754,
                       62.179, 64.716)
robust_regress <- data.frame(
  height = height,
  weight = weight_outlier,
  weight_no_outlier = weight_no_outlier
)

# Breast cancer survival data from Collett (2015), Table 1.2
breast_cancer <- data.frame(
  delta = c(1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0,
            0, 0, 0, 0, 0),
  times = c(23, 47, 69, 70, 71, 100, 101, 148, 181, 198, 208, 212, 224,
            5, 8, 10, 13, 18, 24, 26, 26, 31, 35, 40, 41, 48, 50, 59, 61,
            68, 71, 76, 105, 107, 109, 113, 116, 118, 143, 154, 162, 188,
            212, 217, 225),
  stain = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1)
)

usethis::use_data(shaq_free_throws, inderjit, robust_regress, breast_cancer,
                  overwrite = TRUE)

# ── Applied example datasets ────────────────────────────────────────────────
# Source: data-raw/data/ (originally from docs/Examples/data/). This directory
# is build-ignored, so the raw sources are versioned but are excluded from the
# built package. Run this script from the package root.

data_dir <- file.path("data-raw", "data")

# Mroz (1987) labor force participation
mroz <- read.csv(file.path(data_dir, "mroz.csv"))
usethis::use_data(mroz, overwrite = TRUE)

# Agresti crime data (converted from Stata .dta)
crime <- haven::read_dta(file.path(data_dir, "crime.dta"))
crime <- as.data.frame(crime)
usethis::use_data(crime, overwrite = TRUE)

# Collett bladder cancer recurrence (whitespace-delimited, no header)
collett_bladder <- read.table(
  file.path(data_dir, "collett.dat"),
  header = FALSE,
  col.names = c("patient", "time", "delta", "treat", "init", "size")
)
usethis::use_data(collett_bladder, overwrite = TRUE)

# Cutler (1995) pharmacodynamic data
cutler1995 <- read.csv(file.path(data_dir, "cutler1995_pd.csv"))
usethis::use_data(cutler1995, overwrite = TRUE)

# Bonate adverse events case study
bonate_adverse <- read.csv(file.path(data_dir, "bonate.csv"))
usethis::use_data(bonate_adverse, overwrite = TRUE)

# Lau WIHS missing data (whitespace-delimited)
lau_wihs <- read.table(
  file.path(data_dir, "lau_wihs.dat"),
  header = FALSE,
  col.names = c("id", "black", "age", "cd4", "cd41", "cd42", "cd43", "cd44")
)
lau_wihs$cd41 <- as.numeric(lau_wihs$cd41)
lau_wihs$cd42 <- as.numeric(lau_wihs$cd42)
lau_wihs$cd43 <- as.numeric(lau_wihs$cd43)
lau_wihs$cd44 <- as.numeric(lau_wihs$cd44)
usethis::use_data(lau_wihs, overwrite = TRUE)

# The two survey datasets below are the largest bundled objects. They are stored
# with xz compression, which is the smallest lossless option for their mostly
# numeric columns (xz beats the bzip2 default by about 12 percent and gzip by
# about a third). The values are unchanged by the compressor choice. After xz,
# nsduh.rda is about 503 KB and sdss_quasar.rda about 1.50 MB. sdss_quasar is
# kept at full size on purpose: it is the Sloan Digital Sky Survey 3rd Data
# Release table (n = 46,420 quasars) that the generalized additive model article
# uses to recreate Weinstein et al. (2004) Figure 1, so subsetting it or dropping
# columns would change the documented sample size and the article's results. The
# article that consumes it lives under vignettes/articles, which is excluded from
# the built package, so the only footprint is the data directory itself.

# SDSS quasar survey (whitespace-delimited)
sdss_quasar <- read.table(
  file.path(data_dir, "SDSS_quasar.dat"),
  header = TRUE
)
usethis::use_data(sdss_quasar, overwrite = TRUE, compress = "xz")

# NSDUH substance use survey
nsduh <- read.csv(file.path(data_dir, "nsduh.csv"))
usethis::use_data(nsduh, overwrite = TRUE, compress = "xz")

# ACTG 175 clinical trial
actg175 <- read.csv(file.path(data_dir, "actg175.csv"))
usethis::use_data(actg175, overwrite = TRUE)

# PPARg cheminformatics data
pparg <- read.csv(file.path(data_dir, "pparg.csv"))
usethis::use_data(pparg, overwrite = TRUE)

# GetTested randomized trial (Excel): read from the "data" sheet and recode
get_tested <- readxl::read_xls(file.path(data_dir, "GetTestedData.xls"), sheet = "data")
get_tested <- as.data.frame(get_tested)

# Subset desired columns
get_tested <- get_tested[, c("group", "anytest", "anydiag", "gender", "msm",
                             "age", "partners", "ethnicgrp")]

# Recode categorical variables to numeric (matching Python preprocessing)
get_tested$group <- ifelse(get_tested$group == "SH:24", 1,
                    ifelse(get_tested$group == "Control", 0, NA))

get_tested$gender <- ifelse(get_tested$gender == "Female", 0,
                    ifelse(get_tested$gender == "Male", 1,
                    ifelse(get_tested$gender == "Transgender", 2, NA)))

get_tested$msm <- ifelse(get_tested$msm == "other", 0,
                  ifelse(get_tested$msm == "msm", 1,
                  ifelse(get_tested$msm == "99", 2, NA)))

partner_cats <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10+")
get_tested$partners <- match(get_tested$partners, partner_cats) - 1

ethnic_cats <- c("White/ White British", "Black/ Black British",
                 "Mixed/ Multiple ethnicity", "Asian/ Asian British", "Other")
get_tested$ethnicgrp <- match(get_tested$ethnicgrp, ethnic_cats) - 1

usethis::use_data(get_tested, overwrite = TRUE)
