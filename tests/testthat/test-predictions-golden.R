# Golden-value pins for the five prediction helpers in R/predictions.R.
#
# Each block below records the output of one of those functions across every
# branch it has, so that any later reworking of the file's internals is held to
# the values it produces today. The tolerance-based tests in
# test-predictions.R check parity with Python Delicatessen at 1e-6 and the
# qualitative properties of the results; these check that the numbers do not
# move at all.
#
# The values are recorded at 17 significant digits, the full precision of a
# double, and compared at a relative tolerance of 1e-14, about fifty units in
# the last place. Exact equality is deliberately not asserted, for two reasons.
# R's parser does not correctly round a decimal literal to the nearest double,
# so a literal in this file can sit a unit or two in the last place away from
# the value it was printed from. And the matrix products, exponentials, and
# logarithms underneath these functions come from BLAS and the system math
# library, either of which may differ in the last place from one platform to
# another. A tolerance of 1e-14 is eight orders of magnitude tighter than the
# parity tests elsewhere in the suite, and tight enough that any rearrangement
# of the arithmetic beyond the last place or two fails here.
#
# The fixtures are written out rather than fitted, so nothing here depends on a
# solver, on a random seed, or on an optional package. The parameter vectors are
# plausible fits and the covariance matrices are positive definite, but no value
# below is claimed to be the answer to a statistical question; they are pins.

# ---- Helpers ----------------------------------------------------------------

golden_tolerance <- 1e-14

expect_golden <- function(actual, expected, label) {
  expect_equal(
    as.numeric(actual),
    expected,
    tolerance = golden_tolerance,
    label = label
  )
}

# A design whose rows are named, so the row names of the returned data frame
# are pinned alongside its values.
golden_regression_design <- function() {
  X <- cbind(1, c(-1.5, 0.25, 2), c(0, 1, 1))
  rownames(X) <- c("r1", "r2", "r3")
  X
}

golden_regression_theta <- function() c(0.8, -0.4, 1.25)

# Built row by row so that the symmetry stays visible; the formatter puts a
# matrix() of nine numbers on nine lines, where the rows cannot be read off.
golden_regression_covariance <- function() {
  rbind(
    c(0.040, 0.008, -0.006),
    c(0.008, 0.030, 0.004),
    c(-0.006, 0.004, 0.050)
  )
}

golden_regression_offset <- function() c(0.25, -1.5, 0.75)

golden_survival_theta <- function() c(0.35, 1.4)

golden_survival_covariance <- function() {
  matrix(c(0.0125, 0.0032, 0.0032, 0.0410), nrow = 2)
}

golden_survival_times <- function() c(0.5, 2)

# The AFT parameters are two regression coefficients followed by log(1 / sigma),
# the layout ee_aft() solves for. The exponential distribution ignores the last
# entry, which is why the same vector serves every distribution below.
golden_aft_theta <- function() c(2.1, 0.45, 0.3)

golden_aft_covariance <- function() {
  rbind(
    c(0.0220, 0.0041, -0.0013),
    c(0.0041, 0.0310, 0.0025),
    c(-0.0013, 0.0025, 0.0180)
  )
}

golden_aft_times <- function() c(0.5, 3)

golden_aft_individual_design <- function() {
  rbind(c(1, -0.5), c(1, 0.6), c(1, 1.2))
}

golden_aft_pattern <- function() matrix(c(1, 0.6), nrow = 1)

# Five people, three unique event times, one covariate.
golden_plogit_data <- function() {
  list(
    X = cbind(c(0, 1, 0, 1, 0)),
    time = c(1, 2, 2, 3, 3),
    event = c(1, 1, 0, 1, 0)
  )
}

golden_measures <- function() {
  c("survival", "risk", "cumulative_hazard", "hazard", "density")
}

# ---- regression_predictions -------------------------------------------------

test_that("regression_predictions reproduces its recorded values", {
  golden <- list(
    "no offset" = c(
      1.4000000000000001,
      1.9500000000000002,
      1.25,
      0.083499999999999991,
      0.085875000000000007,
      0.246,
      0.83364162270876663,
      1.3756435982535691,
      0.2778894765044696,
      1.9663583772912336,
      2.5243564017464313,
      2.2221105234955303
    ),
    "offset" = c(
      1.6500000000000001,
      0.45000000000000018,
      2,
      0.083499999999999991,
      0.085875000000000007,
      0.246,
      1.0836416227087666,
      -0.12435640174643092,
      1.0278894765044697,
      2.2163583772912334,
      1.0243564017464313,
      2.9721105234955303
    ),
    "alpha 0.01" = c(
      1.4000000000000001,
      1.9500000000000002,
      1.25,
      0.083499999999999991,
      0.085875000000000007,
      0.246,
      0.65567892265147365,
      1.1951677367701508,
      -0.027569788250809246,
      2.1443210773485264,
      2.7048322632298496,
      2.5275697882508092
    ),
    "alpha 0.5" = c(
      1.4000000000000001,
      1.9500000000000002,
      1.25,
      0.083499999999999991,
      0.085875000000000007,
      0.246,
      1.2050969694168809,
      1.7523445792916541,
      0.91546396294656851,
      1.5949030305831193,
      2.1476554207083463,
      1.5845360370534314
    ),
    "offset and alpha 0.2" = c(
      1.6500000000000001,
      0.45000000000000018,
      2,
      0.083499999999999991,
      0.085875000000000007,
      0.246,
      1.2796781620468358,
      0.07444853495028092,
      1.3643710940125087,
      2.0203218379531647,
      0.82555146504971944,
      2.6356289059874913
    )
  )

  X <- golden_regression_design()
  theta <- golden_regression_theta()
  covariance <- golden_regression_covariance()
  offset <- golden_regression_offset()

  scenarios <- list(
    "no offset" = list(),
    "offset" = list(offset = offset),
    "alpha 0.01" = list(alpha = 0.01),
    "alpha 0.5" = list(alpha = 0.5),
    "offset and alpha 0.2" = list(offset = offset, alpha = 0.2)
  )

  for (nm in names(scenarios)) {
    result <- do.call(
      regression_predictions,
      c(list(X = X, theta = theta, covariance = covariance), scenarios[[nm]])
    )
    expect_golden(unlist(result, use.names = FALSE), golden[[nm]], nm)
  }
})

test_that("regression_predictions reproduces its recorded shape and names", {
  result <- regression_predictions(
    golden_regression_design(),
    golden_regression_theta(),
    golden_regression_covariance()
  )

  expect_s3_class(result, "data.frame")
  expect_identical(names(result), c("predicted", "variance", "lower", "upper"))
  # The variance is formed with rowSums(), which carries the row names of the
  # design into the data frame. Pinning them keeps that incidental but visible
  # part of the return value from changing.
  expect_identical(rownames(result), c("r1", "r2", "r3"))
})

test_that("regression_predictions reproduces its recorded values for one row", {
  golden <- c(0, 0.192, -0.8588131889842352, 0.8588131889842352)

  result <- regression_predictions(
    matrix(c(1, 2, 0), nrow = 1),
    golden_regression_theta(),
    golden_regression_covariance()
  )

  expect_identical(nrow(result), 1L)
  expect_golden(unlist(result, use.names = FALSE), golden, "one row")
})

# ---- survival_predictions ---------------------------------------------------

test_that("survival_predictions reproduces its recorded values", {
  golden <- list(
    "exponential/survival/capprox" = c(
      0.83945702076920736,
      0.49658530379140953,
      0.0022021499006691482,
      0.012329848730618957,
      0.7474816536298684,
      0.27895119303654548,
      0.93143238790854632,
      0.71421941454627358
    ),
    "exponential/survival/exact" = c(
      0.83945702076920736,
      0.49658530379140953,
      0.0022021502803709798,
      0.012329848197080325,
      0.74748164570052311,
      0.27895119774528965,
      0.93143239583789161,
      0.71421940983752941
    ),
    "exponential/risk/capprox" = c(
      0.16054297923079264,
      0.50341469620859047,
      0.0022021499006691482,
      0.012329848041468404,
      0.06856761209145372,
      0.28578059153582336,
      0.25251834637013154,
      0.72104880088135759
    ),
    "exponential/risk/exact" = c(
      0.16054297923079264,
      0.50341469620859047,
      0.0022021502803709798,
      0.012329848197080325,
      0.068567604162108442,
      0.28578059016247059,
      0.25251835429947683,
      0.72104880225471035
    ),
    "exponential/cumulative_hazard/capprox" = c(
      0.17499999999999999,
      0.69999999999999996,
      0.0031249994762932547,
      0.050000002722922016,
      0.065434691608743059,
      0.26173871777819685,
      0.28456530839125693,
      1.138261282221803
    ),
    "exponential/cumulative_hazard/exact" = c(
      0.17499999999999999,
      0.69999999999999996,
      0.0031250000000000002,
      0.050000000000000003,
      0.06543468242792734,
      0.26173872971170936,
      0.28456531757207265,
      1.1382612702882906
    ),
    "exponential/hazard/capprox" = c(
      0.34999999999999998,
      0.34999999999999998,
      0.012500000680730504,
      0.012500000680730504,
      0.13086935888909842,
      0.13086935888909842,
      0.5691306411109015,
      0.5691306411109015
    ),
    "exponential/hazard/exact" = c(
      0.34999999999999998,
      0.34999999999999998,
      0.012500000000000001,
      0.012500000000000001,
      0.13086936485585468,
      0.13086936485585468,
      0.5691306351441453,
      0.5691306351441453
    ),
    "exponential/density/capprox" = c(
      0.29380995726922254,
      0.17380485632699333,
      0.0059953546423461672,
      0.00027742160160755579,
      0.1420505820266032,
      0.14115973940965887,
      0.44556933251184189,
      0.20644997324432779
    ),
    "exponential/density/exact" = c(
      0.29380995726922254,
      0.17380485632699333,
      0.0059953541383099924,
      0.00027742158443430734,
      0.1420505884058936,
      0.14115974042007534,
      0.44556932613255151,
      0.20644997223391132
    ),
    "weibull/survival/capprox" = c(
      0.87579327830036235,
      0.39706489922339999,
      0.0014714296283216321,
      0.018079560244967003,
      0.8006105919640335,
      0.13352764155894892,
      0.95097596463669121,
      0.66060215688785107
    ),
    "weibull/survival/exact" = c(
      0.87579327830036235,
      0.39706489922339999,
      0.0014714297064164323,
      0.018079559214071277,
      0.80061058996890688,
      0.1335276490723904,
      0.95097596663181783,
      0.66060214937440964
    ),
    "weibull/risk/capprox" = c(
      0.12420672169963765,
      0.60293510077659995,
      0.0014714296283216321,
      0.018079561009681389,
      0.049024035363308807,
      0.33939783753870773,
      0.1993894080359665,
      0.86647236401449224
    ),
    "weibull/risk/exact" = c(
      0.12420672169963765,
      0.60293510077659995,
      0.0014714297064164323,
      0.018079559214071277,
      0.049024033368182202,
      0.33939785062559036,
      0.19938941003109309,
      0.86647235092760955
    ),
    "weibull/cumulative_hazard/capprox" = c(
      0.13262519956965987,
      0.92365553754102592,
      0.001918387262323973,
      0.11467397211943439,
      0.046779957008632203,
      0.25994223972167896,
      0.21847044213068756,
      1.587368835360373
    ),
    "weibull/cumulative_hazard/exact" = c(
      0.13262519956965987,
      0.92365553754102592,
      0.0019183873604314188,
      0.1146739690205333,
      0.0467799548135443,
      0.25994224868963323,
      0.21847044432577545,
      1.5873688263924186
    ),
    "weibull/hazard/capprox" = c(
      0.37135055879504753,
      0.64655887627871811,
      0.014127386199546149,
      0.087366809767342457,
      0.13839182354436205,
      0.067235130850772595,
      0.60430929404573297,
      1.2258826217066636
    ),
    "weibull/hazard/exact" = c(
      0.37135055879504753,
      0.64655887627871811,
      0.014127384651522168,
      0.087366799957133354,
      0.13839183630771804,
      0.067235163376209051,
      0.60430928128237704,
      1.2258825891812273
    ),
    "weibull/density/capprox" = c(
      0.32522632328578616,
      0.25672583505160396,
      0.0083651779194993132,
      0.0017003273131658079,
      0.14596537799781681,
      0.17590667052221814,
      0.50448726857375548,
      0.33754499958098977
    ),
    "weibull/density/exact" = c(
      0.32522632328578616,
      0.25672583505160396,
      0.0083651770147434024,
      0.0017003268247517541,
      0.14596538769201611,
      0.17590668212975347,
      0.50448725887955614,
      0.33754498797345445
    ),
    "gompertz/survival/capprox" = c(
      0.77612772608910308,
      0.021043539190673869,
      0.0042688293493298064,
      0.001386641001147202,
      0.64807097242877898,
      -0.051940872192806145,
      0.90418447974942717,
      0.094027950574153876
    ),
    "gompertz/survival/exact" = c(
      0.77612772608910308,
      0.021043539190673869,
      0.0042688287514033429,
      0.0013866407976756685,
      0.64807098139710706,
      -0.051940866838049354,
      0.90418447078109909,
      0.094027945219397085
    ),
    "gompertz/risk/capprox" = c(
      0.22387227391089692,
      0.97895646080932608,
      0.0042688293493298064,
      0.0013866409822931066,
      0.0958155202505728,
      0.90597204992202895,
      0.35192902757122102,
      1.0519408716966232
    ),
    "gompertz/risk/exact" = c(
      0.22387227391089692,
      0.97895646080932608,
      0.0042688287514033429,
      0.0013866407976756685,
      0.095815529218900908,
      0.90597205478060283,
      0.35192901860289294,
      1.0519408668380492
    ),
    "gompertz/cumulative_hazard/capprox" = c(
      0.2534381768676191,
      3.861161692774262,
      0.0070866733541853225,
      3.1313129434871549,
      0.088443741004426396,
      0.39290424455730344,
      0.41843261273081178,
      7.3294191409912202
    ),
    "gompertz/cumulative_hazard/exact" = c(
      0.2534381768676191,
      3.861161692774262,
      0.0070866720667690647,
      3.1313124710146925,
      0.088443755991467449,
      0.39290450621371242,
      0.41843259774377073,
      7.3294188793348116
    ),
    "gompertz/hazard/capprox" = c(
      0.7048134476146668,
      5.7556263698839665,
      0.060323639540816085,
      10.024708741375925,
      0.22342921728463128,
      -0.44997640659536753,
      1.1861976779447023,
      11.961229146363301
    ),
    "gompertz/hazard/exact" = c(
      0.7048134476146668,
      5.7556263698839665,
      0.060323634006793816,
      10.024706923287477,
      0.22342923936545339,
      -0.44997584386903178,
      1.1861976558638803,
      11.961228583636965
    ),
    "gompertz/density/capprox" = c(
      0.54702525841419247,
      0.12111874908152923,
      0.021033374540144239,
      0.022131541873389889,
      0.2627738869905874,
      -0.17045869645412073,
      0.8312766298377976,
      0.41269619461717921
    ),
    "gompertz/density/exact" = c(
      0.54702525841419247,
      0.12111874908152923,
      0.021033374252092115,
      0.022131538835675865,
      0.2627738889369991,
      -0.17045867644356452,
      0.8312766278913859,
      0.412696174606623
    )
  )

  times <- golden_survival_times()
  theta <- golden_survival_theta()
  covariance <- golden_survival_covariance()

  for (distribution in c("exponential", "weibull", "gompertz")) {
    for (measure in golden_measures()) {
      for (deriv_method in c("capprox", "exact")) {
        key <- paste(distribution, measure, deriv_method, sep = "/")
        result <- survival_predictions(
          times = times,
          theta = theta,
          covariance = covariance,
          distribution = distribution,
          measure = measure,
          deriv_method = deriv_method
        )
        expect_identical(result$time, times, label = paste0(key, ": time"))
        expect_golden(
          unlist(result[-1], use.names = FALSE),
          golden[[key]],
          key
        )
      }
    }
  }
})

test_that("survival_predictions reproduces its recorded values off the defaults", {
  golden <- list(
    "weibull/survival/fapprox" = c(
      0.87579327830036235,
      0.39706489922339999,
      0.0014714291963715297,
      0.01807955948025268,
      0.800610602999278,
      0.13352764713238968,
      0.95097595360144671,
      0.66060215131441025
    ),
    "weibull/survival/bapprox" = c(
      0.87579327830036235,
      0.39706489922339999,
      0.0014714300602718118,
      0.018079561009681389,
      0.80061058092878867,
      0.13352763598550776,
      0.95097597567193604,
      0.66060216246129222
    ),
    "weibull/survival/exact/alpha=0.01" = c(
      0.87579327830036235,
      0.39706489922339999,
      0.0014714297064164323,
      0.018079559214071277,
      0.77698647712886393,
      0.050718243432697652,
      0.97460007947186078,
      0.74341155501410228
    ),
    "weibull/survival/exact/alpha=0.5" = c(
      0.87579327830036235,
      0.39706489922339999,
      0.0014714297064164323,
      0.018079559214071277,
      0.84992037804722154,
      0.30637283782645025,
      0.90166617855350317,
      0.48775696062034973
    )
  )

  times <- golden_survival_times()
  theta <- golden_survival_theta()
  covariance <- golden_survival_covariance()

  for (deriv_method in c("fapprox", "bapprox")) {
    key <- paste("weibull", "survival", deriv_method, sep = "/")
    result <- survival_predictions(
      times,
      theta,
      covariance,
      "weibull",
      "survival",
      deriv_method = deriv_method
    )
    expect_golden(unlist(result[-1], use.names = FALSE), golden[[key]], key)
  }

  for (alpha in c(0.01, 0.5)) {
    key <- paste0("weibull/survival/exact/alpha=", alpha)
    result <- survival_predictions(
      times,
      theta,
      covariance,
      "weibull",
      "survival",
      alpha = alpha,
      deriv_method = "exact"
    )
    expect_golden(unlist(result[-1], use.names = FALSE), golden[[key]], key)
  }
})

test_that("survival_predictions reproduces its recorded values for one time", {
  golden <- c(
    0.53932253771734195,
    0.012884813003423104,
    0.31684449708188056,
    0.76180057835280335
  )

  result <- survival_predictions(
    1.5,
    golden_survival_theta(),
    golden_survival_covariance(),
    "weibull",
    "survival",
    deriv_method = "exact"
  )

  expect_identical(nrow(result), 1L)
  expect_golden(unlist(result[-1], use.names = FALSE), golden, "one time")
})

test_that("survival_predictions reproduces its recorded shape and names", {
  result <- survival_predictions(
    golden_survival_times(),
    golden_survival_theta(),
    golden_survival_covariance(),
    "weibull"
  )

  expect_s3_class(result, "data.frame")
  expect_identical(
    names(result),
    c("time", "predicted", "variance", "lower", "upper")
  )
  expect_identical(rownames(result), c("1", "2"))
})

# ---- aft_predictions_individual ---------------------------------------------

test_that("aft_predictions_individual reproduces its recorded values", {
  golden <- list(
    "exponential/survival" = c(
      0.92618851659254886,
      0.95433514600829228,
      0.96494841524348351,
      0.63124267619021091,
      0.7554496237573709,
      0.80728072506427084
    ),
    "exponential/risk" = c(
      0.073811483407451139,
      0.045664853991707721,
      0.035051584756516485,
      0.36875732380978909,
      0.2445503762426291,
      0.19271927493572916
    ),
    "exponential/cumulative_hazard" = c(
      0.076677483422464221,
      0.046740363139029246,
      0.035680634778193006,
      0.46006490053478549,
      0.28044217883417549,
      0.21408380866915819
    ),
    "exponential/hazard" = c(
      0.15335496684492844,
      0.093480726278058451,
      0.071361269556386039,
      0.15335496684492847,
      0.093480726278058479,
      0.071361269556386053
    ),
    "exponential/density" = c(
      0.14203560925420378,
      0.089211942561532115,
      0.068859943968197751,
      0.096804199678253711,
      0.070619979495325058,
      0.057608577428986209
    ),
    "weibull/survival" = c(
      0.96926048797109099,
      0.98412180238310054,
      0.98894454644916818,
      0.70424181851157974,
      0.8354798987035299,
      0.88263077642764043
    ),
    "weibull/risk" = c(
      0.030739512028909011,
      0.015878197616899459,
      0.011055453550831817,
      0.29575818148842026,
      0.1645201012964701,
      0.11736922357235957
    ),
    "weibull/cumulative_hazard" = c(
      0.03122188178178437,
      0.016005606680700817,
      0.011117019255824558,
      0.3506334896000875,
      0.17974899023846908,
      0.12484831256696781
    ),
    "weibull/hazard" = c(
      0.084290264224476844,
      0.043210618297082615,
      0.030012812712933622,
      0.15776856805592901,
      0.080878585875430148,
      0.056175864776507739
    ),
    "weibull/density" = c(
      0.081699222633428617,
      0.042524511560613126,
      0.02968100745605597,
      0.11110722327167537,
      0.067572432734489124,
      0.049582547144183163
    ),
    "log-normal/survival" = c(
      0.9997364925645843,
      0.99998223805846487,
      0.99999659076910852,
      0.85268387332285622,
      0.95693672034349442,
      0.98126728774730632
    ),
    "log-normal/risk" = c(
      0.00026350743541569965,
      1.7761941535132664e-05,
      3.4092308914823732e-06,
      0.14731612667714378,
      0.043063279656505582,
      0.018732712252693684
    ),
    "log-normal/cumulative_hazard" = c(
      0.00026354215960014714,
      1.7762099280284122e-05,
      3.4092367029232175e-06,
      0.15936640589936532,
      0.044018012658457999,
      0.018910391954288706
    ),
    "log-normal/hazard" = c(
      0.0026470467255962139,
      0.00020880158795490533,
      4.3292135925026773e-05,
      0.12155956757521445,
      0.043015551736348336,
      0.021000591746627006
    ),
    "log-normal/density" = c(
      0.0026463492091021263,
      0.00020879787923330762,
      4.3291988332139622e-05,
      0.10365188291948534,
      0.041163161002347085,
      0.020607193704301148
    ),
    "log-logistic/survival" = c(
      0.96972341032190079,
      0.98424653705111798,
      0.9890052100359199,
      0.74039331002823905,
      0.84763793677658894,
      0.88900875685001757
    ),
    "log-logistic/risk" = c(
      0.030276589678099208,
      0.015753462948882024,
      0.010994789964080098,
      0.25960668997176095,
      0.15236206322341106,
      0.11099124314998243
    ),
    "log-logistic/cumulative_hazard" = c(
      0.030744392144028507,
      0.01587886752758047,
      0.011055679389559663,
      0.30057373394187775,
      0.16530169571254533,
      0.11764819329175398
    ),
    "log-logistic/hazard" = c(
      0.081738242480693807,
      0.042529901422741233,
      0.029682828140923646,
      0.11681079232134478,
      0.068555757660857766,
      0.049940835709937834
    ),
    "log-logistic/density" = c(
      0.079263487252096865,
      0.041859908196458476,
      0.029356471679974305,
      0.086485929173821668,
      0.058110460977805312,
      0.044397840270542797
    )
  )

  X <- golden_aft_individual_design()
  times <- golden_aft_times()
  theta <- golden_aft_theta()

  distributions <- c("exponential", "weibull", "log-normal", "log-logistic")
  for (distribution in distributions) {
    for (measure in golden_measures()) {
      key <- paste(distribution, measure, sep = "/")
      result <- aft_predictions_individual(
        X = X,
        times = times,
        theta = theta,
        distribution = distribution,
        measure = measure
      )
      expect_golden(unlist(result, use.names = FALSE), golden[[key]], key)
    }
  }
})

test_that("aft_predictions_individual reproduces its recorded shape and names", {
  result <- aft_predictions_individual(
    golden_aft_individual_design(),
    golden_aft_times(),
    golden_aft_theta(),
    "weibull"
  )

  expect_s3_class(result, "data.frame")
  expect_identical(names(result), c("t_0.5", "t_3"))
  expect_identical(rownames(result), c("1", "2", "3"))
})

test_that("aft_predictions_individual reads the two spellings alike", {
  # tolower() plus the paired branch conditions mean each distribution has two
  # accepted spellings and a case-insensitive form. Comparing the results with
  # each other rather than with a recorded value keeps the check exact.
  args <- list(
    X = golden_aft_individual_design(),
    times = golden_aft_times(),
    theta = golden_aft_theta(),
    measure = "density"
  )
  spellings <- list(
    c("log-normal", "lognormal"),
    c("log-logistic", "loglogistic"),
    c("weibull", "WEIBULL"),
    c("exponential", "Exponential")
  )

  for (pair in spellings) {
    expect_identical(
      do.call(aft_predictions_individual, c(args, distribution = pair[1])),
      do.call(aft_predictions_individual, c(args, distribution = pair[2])),
      label = paste(pair, collapse = " and ")
    )
  }
})

# ---- aft_predictions_function -----------------------------------------------

test_that("aft_predictions_function reproduces its recorded values", {
  golden <- list(
    "exponential/survival/capprox" = c(
      0.95433514600829228,
      0.7554496237573709,
      7.576754557177082e-05,
      0.0017092063089318547,
      0.93727472672097634,
      0.67441971784706933,
      0.97139556529560822,
      0.83647952966767247
    ),
    "exponential/survival/exact" = c(
      0.95433514600829228,
      0.7554496237573709,
      7.5767497022132037e-05,
      0.0017092068450235073,
      0.93727473218688728,
      0.67441970513957683,
      0.97139555982969727,
      0.83647954237516497
    ),
    "exponential/risk/capprox" = c(
      0.045664853991707721,
      0.2445503762426291,
      7.576754557177082e-05,
      0.0017092063089318547,
      0.028604434704391748,
      0.16352047033232758,
      0.062725273279023691,
      0.32558028215293061
    ),
    "exponential/risk/exact" = c(
      0.045664853991707721,
      0.2445503762426291,
      7.5767497022132037e-05,
      0.0017092068450235073,
      0.028604440170302777,
      0.16352045762483508,
      0.062725267813112662,
      0.32558029486042311
    ),
    "exponential/cumulative_hazard/capprox" = c(
      0.046740363139029246,
      0.28044217883417549,
      8.3191966419345629e-05,
      0.0029949077455154551,
      0.028863604114740628,
      0.17318167922609287,
      0.064617122163317858,
      0.38770267844225814
    ),
    "exponential/cumulative_hazard/exact" = c(
      0.046740363139029246,
      0.28044217883417549,
      8.3191911685705686e-05,
      0.0029949088206854069,
      0.028863609995478247,
      0.17318165997286944,
      0.064617116282580239,
      0.38770269769548155
    ),
    "exponential/hazard/capprox" = c(
      0.093480726278058451,
      0.093480726278058479,
      0.00033276771158338069,
      0.00033276771158338069,
      0.057727216507629263,
      0.057727216507629291,
      0.12923423604848763,
      0.12923423604848766
    ),
    "exponential/hazard/exact" = c(
      0.093480726278058451,
      0.093480726278058479,
      0.00033276764674282274,
      0.00033276764674282301,
      0.057727219990956452,
      0.057727219990956466,
      0.12923423256516045,
      0.12923423256516048
    ),
    "exponential/density/capprox" = c(
      0.089211942561532115,
      0.070619979495325058,
      0.00027540093891392916,
      9.8329480666101797e-05,
      0.056685931954352571,
      0.051184737000334962,
      0.12173795316871167,
      0.090055221990315154
    ),
    "exponential/density/exact" = c(
      0.089211942561532115,
      0.070619979495325058,
      0.00027540089083776529,
      9.8329427390910491e-05,
      0.056685934793351392,
      0.051184742265370411,
      0.12173795032971284,
      0.090055216725279705
    ),
    "weibull/survival/capprox" = c(
      0.98412180238310054,
      0.8354798987035299,
      9.412245235717436e-05,
      0.0027814368954283135,
      0.96510687338421053,
      0.73211271195576566,
      1.0031367313819906,
      0.93884708545129414
    ),
    "weibull/survival/exact" = c(
      0.98412180238310054,
      0.8354798987035299,
      9.4122401404978217e-05,
      0.002781436498193849,
      0.96510687853097676,
      0.73211271933702438,
      1.0031367262352244,
      0.93884707807003542
    ),
    "weibull/risk/capprox" = c(
      0.015878197616899459,
      0.1645201012964701,
      9.412245235717436e-05,
      0.0027814368954283135,
      -0.0031367313819905922,
      0.061152914548705847,
      0.034893126615789513,
      0.26788728804423434
    ),
    "weibull/risk/exact" = c(
      0.015878197616899459,
      0.1645201012964701,
      9.4122401404978217e-05,
      0.002781436498193849,
      -0.0031367262352243359,
      0.061152921929964535,
      0.03489312146902325,
      0.26788728066297568
    ),
    "weibull/cumulative_hazard/capprox" = c(
      0.016005606680700817,
      0.17974899023846908,
      9.718416640339811e-05,
      0.003984714408601482,
      -0.0033161161529943993,
      0.056027058824232229,
      0.03532732951439603,
      0.30347092165270595
    ),
    "weibull/cumulative_hazard/exact" = c(
      0.016005606680700817,
      0.17974899023846908,
      9.7184116846509868e-05,
      0.0039847138412458299,
      -0.0033161112266539422,
      0.056027067632183361,
      0.035327324588055577,
      0.3034709128447548
    ),
    "weibull/hazard/capprox" = c(
      0.043210618297082615,
      0.080878585875430148,
      0.00046299202690414326,
      0.00051680432157952087,
      0.001037584191602707,
      0.036322076589979019,
      0.085383652402562515,
      0.12543509516088128
    ),
    "weibull/hazard/exact" = c(
      0.043210618297082615,
      0.080878585875430148,
      0.0004629919857246814,
      0.00051680421251873194,
      0.0010375860670808826,
      0.036322081291340944,
      0.085383650527084354,
      0.12543509045951934
    ),
    "weibull/density/capprox" = c(
      0.042524511560613126,
      0.067572432734489124,
      0.00043094551993066939,
      0.00022764678971017195,
      0.0018371736189845972,
      0.03800055811797487,
      0.083211849502241655,
      0.097144307351003378
    ),
    "weibull/density/exact" = c(
      0.042524511560613126,
      0.067572432734489124,
      0.00043094548848837184,
      0.00022764675872675107,
      0.0018371751032828995,
      0.03800056013038601,
      0.083211848017943346,
      0.097144305338592238
    ),
    "log-normal/survival/capprox" = c(
      0.99998223805846487,
      0.95693672034349442,
      2.2692182209560327e-09,
      0.0010321611164179321,
      0.99988887260258863,
      0.89396843924069502,
      1.000075603514341,
      1.0199050014462938
    ),
    "log-normal/survival/exact" = c(
      0.99998223805846487,
      0.95693672034349442,
      2.2691672009162199e-09,
      0.0010321610714903744,
      0.99988887365218715,
      0.89396844061112601,
      1.0000756024647426,
      1.0199050000758627
    ),
    "log-normal/risk/capprox" = c(
      1.7761941535132664e-05,
      0.043063279656505582,
      2.2692182209560327e-09,
      0.0010321611164179321,
      -7.560351434107124e-05,
      -0.019905001446293771,
      0.00011112739741133657,
      0.10603156075930494
    ),
    "log-normal/risk/exact" = c(
      1.7761941535132664e-05,
      0.043063279656505582,
      2.2691672009162199e-09,
      0.0010321610714903744,
      -7.560246474257793e-05,
      -0.019905000075862808,
      0.00011112634781284326,
      0.10603155938887397
    ),
    "log-normal/cumulative_hazard/capprox" = c(
      1.7762099280284122e-05,
      0.044018012658457999,
      2.2692988283496717e-09,
      0.0011271482745844063,
      -7.5605014849663605e-05,
      -0.021783914985285292,
      0.00011112921341023184,
      0.1098199403022013
    ),
    "log-normal/cumulative_hazard/exact" = c(
      1.7762099280284122e-05,
      0.044018012658457999,
      2.2692478126942581e-09,
      0.0011271482352283601,
      -7.5603965360007984e-05,
      -0.021783913836499573,
      0.00011112816392057622,
      0.10981993915341556
    ),
    "log-normal/hazard/capprox" = c(
      0.00020880158795490533,
      0.043015551736348336,
      2.5663222140747605e-07,
      0.00057346095973545303,
      -0.0007840942326627847,
      -0.0039197938157989609,
      0.0012016974085725953,
      0.089950897288495632
    ),
    "log-normal/hazard/exact" = c(
      0.00020880158795490533,
      0.043015551736348336,
      2.5663218183884279e-07,
      0.0005734608190054607,
      -0.00078409415611835874,
      -0.0039197880567223287,
      0.0012016973320281693,
      0.089950891529419
    ),
    "log-normal/density/capprox" = c(
      0.00020879787923330762,
      0.041163161002347085,
      2.5661303089154537e-07,
      0.00046474936557737542,
      -0.00078406081716985137,
      -0.0010898335638278231,
      0.0012016565756364666,
      0.083416155568521999
    ),
    "log-normal/density/exact" = c(
      0.00020879787923330762,
      0.041163161002347085,
      2.5661299094931436e-07,
      0.00046474926313521292,
      -0.00078406073989982202,
      -0.0010898289070292375,
      0.0012016564983664373,
      0.083416150911723413
    ),
    "log-logistic/survival/capprox" = c(
      0.98424653705111798,
      0.84763793677658894,
      9.120333132017537e-05,
      0.0020570206994135307,
      0.96552879569227401,
      0.75874496568149896,
      1.0029642784099619,
      0.93653090787167892
    ),
    "log-logistic/survival/exact" = c(
      0.98424653705111798,
      0.84763793677658894,
      9.1203367530470089e-05,
      0.0020570207905549922,
      0.96552879197653885,
      0.75874496371218592,
      1.0029642821256972,
      0.93653090984099197
    ),
    "log-logistic/risk/capprox" = c(
      0.015753462948882024,
      0.15236206322341106,
      9.120333132017537e-05,
      0.0020570206994135307,
      -0.002964278409961895,
      0.063469092128321081,
      0.034471204307725947,
      0.24125503431850104
    ),
    "log-logistic/risk/exact" = c(
      0.015753462948882024,
      0.15236206322341106,
      9.1203367530470089e-05,
      0.0020570207905549922,
      -0.0029642821256971319,
      0.063469090159008076,
      0.034471208023461181,
      0.24125503628781403
    ),
    "log-logistic/cumulative_hazard/capprox" = c(
      0.01587886752758047,
      0.16530169571254533,
      9.4146222680196903e-05,
      0.0028629771838190187,
      -0.0031384623993539508,
      0.060430303253840995,
      0.034896197454514888,
      0.27017308817124963
    ),
    "log-logistic/cumulative_hazard/exact" = c(
      0.01587886752758047,
      0.16530169571254533,
      9.4146262415862223e-05,
      0.0028629773341496984,
      -0.0031384664126120239,
      0.060430300500520506,
      0.034896201467772961,
      0.27017309092457015
    ),
    "log-logistic/hazard/capprox" = c(
      0.042529901422741233,
      0.068555757660857766,
      0.00043132718268858465,
      0.00025277520958549393,
      0.0018245502936303504,
      0.037394474655901101,
      0.083235252551852115,
      0.099717040665814438
    ),
    "log-logistic/hazard/exact" = c(
      0.042529901422741233,
      0.068555757660857766,
      0.00043132714148562575,
      0.00025277515760872425,
      0.0018245522378400048,
      0.037394477859662489,
      0.083235250607642461,
      0.099717037462053043
    ),
    "log-logistic/density/capprox" = c(
      0.041859908196458476,
      0.058110460977805312,
      0.00040151948875847122,
      0.00011577966627146267,
      0.0025862454977512925,
      0.037021060316441506,
      0.081133570895165652,
      0.079199861639169117
    ),
    "log-logistic/density/exact" = c(
      0.041859908196458476,
      0.058110460977805312,
      0.00040151944638605776,
      0.00011577963662563727,
      0.0025862475700291834,
      0.037021063016452443,
      0.081133568822887775,
      0.07919985893915818
    )
  )

  X <- golden_aft_pattern()
  times <- golden_aft_times()
  theta <- golden_aft_theta()
  covariance <- golden_aft_covariance()

  distributions <- c("exponential", "weibull", "log-normal", "log-logistic")
  for (distribution in distributions) {
    for (measure in golden_measures()) {
      for (deriv_method in c("capprox", "exact")) {
        key <- paste(distribution, measure, deriv_method, sep = "/")
        result <- aft_predictions_function(
          X = X,
          times = times,
          theta = theta,
          covariance = covariance,
          distribution = distribution,
          measure = measure,
          deriv_method = deriv_method
        )
        expect_identical(result$time, times, label = paste0(key, ": time"))
        expect_golden(
          unlist(result[-1], use.names = FALSE),
          golden[[key]],
          key
        )
      }
    }
  }
})

test_that("aft_predictions_function reproduces its values off the defaults", {
  golden <- list(
    "weibull/survival/fapprox" = c(
      0.98412180238310054,
      0.8354798987035299,
      9.412245235717436e-05,
      0.0027814374778530556,
      0.96510687338421053,
      0.73211270113337312,
      1.0031367313819906,
      0.93884709627368668
    ),
    "weibull/survival/bapprox" = c(
      0.98412180238310054,
      0.8354798987035299,
      9.412245235717436e-05,
      0.0027814363130037622,
      0.96510687338421053,
      0.73211272277815576,
      1.0031367313819906,
      0.93884707462890404
    ),
    "weibull/survival/exact/alpha=0.01" = c(
      0.98412180238310054,
      0.8354798987035299,
      9.4122401404978217e-05,
      0.002781436498193849,
      0.95913195642195381,
      0.6996323975567561,
      1.0091116483442473,
      0.9713273998503037
    ),
    "weibull/survival/exact/alpha=0.5" = c(
      0.98412180238310054,
      0.8354798987035299,
      9.4122401404978217e-05,
      0.002781436498193849,
      0.97757812538689992,
      0.79990776394017138,
      0.99066547937930116,
      0.87105203346688842
    )
  )

  X <- golden_aft_pattern()
  times <- golden_aft_times()
  theta <- golden_aft_theta()
  covariance <- golden_aft_covariance()

  for (deriv_method in c("fapprox", "bapprox")) {
    key <- paste("weibull", "survival", deriv_method, sep = "/")
    result <- aft_predictions_function(
      X,
      times,
      theta,
      covariance,
      "weibull",
      "survival",
      deriv_method = deriv_method
    )
    expect_golden(unlist(result[-1], use.names = FALSE), golden[[key]], key)
  }

  for (alpha in c(0.01, 0.5)) {
    key <- paste0("weibull/survival/exact/alpha=", alpha)
    result <- aft_predictions_function(
      X,
      times,
      theta,
      covariance,
      "weibull",
      "survival",
      alpha = alpha,
      deriv_method = "exact"
    )
    expect_golden(unlist(result[-1], use.names = FALSE), golden[[key]], key)
  }
})

test_that("aft_predictions_function reproduces its recorded values for one time", {
  golden <- c(
    0.8354798987035299,
    0.002781436498193849,
    0.73211271933702438,
    0.93884707807003542
  )

  result <- aft_predictions_function(
    golden_aft_pattern(),
    3,
    golden_aft_theta(),
    golden_aft_covariance(),
    "weibull",
    "survival",
    deriv_method = "exact"
  )

  expect_identical(nrow(result), 1L)
  expect_golden(unlist(result[-1], use.names = FALSE), golden, "one time")
})

test_that("aft_predictions_function agrees with the individual predictions", {
  # The two functions answer the same question for one covariate pattern, one
  # through a matrix product and the other through a scalar accumulation that
  # keeps derivatives attached. Their point estimates have to agree, and this
  # pins that they agree exactly rather than to a tolerance.
  pattern <- golden_aft_pattern()
  times <- golden_aft_times()
  theta <- golden_aft_theta()

  distributions <- c("exponential", "weibull", "log-normal", "log-logistic")
  for (distribution in distributions) {
    for (measure in golden_measures()) {
      individual <- aft_predictions_individual(
        X = pattern,
        times = times,
        theta = theta,
        distribution = distribution,
        measure = measure
      )
      curve <- aft_predictions_function(
        X = pattern,
        times = times,
        theta = theta,
        covariance = golden_aft_covariance(),
        distribution = distribution,
        measure = measure
      )
      expect_identical(
        unlist(individual[1, ], use.names = FALSE),
        curve$predicted,
        label = paste(distribution, measure, sep = "/")
      )
    }
  }
})

# ---- plogit_predict ---------------------------------------------------------

test_that("plogit_predict reproduces its recorded values", {
  golden <- list(
    "disjoint/survival" = c(
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987,
      0.75026010559511769,
      0.5176603270687391,
      0.40679559434386792,
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987,
      0.75026010559511769,
      0.5176603270687391,
      0.40679559434386792,
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987
    ),
    "disjoint/risk" = c(
      0.18242552380635635,
      0.3716737526889573,
      0.46872833303880013,
      0.24973989440488231,
      0.4823396729312609,
      0.59320440565613208,
      0.18242552380635635,
      0.3716737526889573,
      0.46872833303880013,
      0.24973989440488231,
      0.4823396729312609,
      0.59320440565613208,
      0.18242552380635635,
      0.3716737526889573,
      0.46872833303880013
    ),
    "disjoint/cumulative_hazard" = c(
      0.20141327798275241,
      0.46469574532078367,
      0.63248177470704969,
      0.2873353251154307,
      0.65843599106320827,
      0.89944444489620057,
      0.20141327798275241,
      0.46469574532078367,
      0.63248177470704969,
      0.2873353251154307,
      0.65843599106320827,
      0.89944444489620057,
      0.20141327798275241,
      0.46469574532078367,
      0.63248177470704969
    ),
    "disjoint/hazard" = c(
      0.18242552380635635,
      0.23147521650098238,
      0.1544652650835347,
      0.24973989440488234,
      0.31002551887238755,
      0.21416501695744142,
      0.18242552380635635,
      0.23147521650098238,
      0.1544652650835347,
      0.24973989440488234,
      0.31002551887238755,
      0.21416501695744142,
      0.18242552380635635,
      0.23147521650098238,
      0.1544652650835347
    ),
    "disjoint/density" = c(
      0.14914645207033286,
      0.14544195412957339,
      0.082063018868533108,
      0.18736987954752057,
      0.16048791149913569,
      0.087121385360866929,
      0.14914645207033286,
      0.14544195412957339,
      0.082063018868533108,
      0.18736987954752057,
      0.16048791149913569,
      0.087121385360866929,
      0.14914645207033286,
      0.14544195412957339,
      0.082063018868533108
    ),
    "unique_times/survival" = c(
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987,
      0.42617757168393511,
      0.75026010559511769,
      0.5176603270687391,
      0.40679559434386792,
      0.29739140899397615,
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987,
      0.42617757168393511,
      0.75026010559511769,
      0.5176603270687391,
      0.40679559434386792,
      0.29739140899397615,
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987,
      0.42617757168393511
    ),
    "S/survival" = c(
      0.77729986117469119,
      0.56825173167966991,
      0.38594447192636266,
      0.700567142473973,
      0.45232559347290596,
      0.26534214452542054,
      0.77729986117469119,
      0.56825173167966991,
      0.38594447192636266,
      0.700567142473973,
      0.45232559347290596,
      0.26534214452542054,
      0.77729986117469119,
      0.56825173167966991,
      0.38594447192636266
    ),
    "times_to_predict/survival" = c(
      1,
      1,
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987,
      1,
      1,
      0.75026010559511769,
      0.5176603270687391,
      0.40679559434386792,
      1,
      1,
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987,
      1,
      1,
      0.75026010559511769,
      0.5176603270687391,
      0.40679559434386792,
      1,
      1,
      0.81757447619364365,
      0.6283262473110427,
      0.53127166696119987
    ),
    "times_to_predict/hazard" = c(
      0,
      0,
      0.18242552380635635,
      0.23147521650098238,
      0.1544652650835347,
      0,
      0,
      0.24973989440488234,
      0.31002551887238755,
      0.21416501695744142,
      0,
      0,
      0.18242552380635635,
      0.23147521650098238,
      0.1544652650835347,
      0,
      0,
      0.24973989440488234,
      0.31002551887238755,
      0.21416501695744142,
      0,
      0,
      0.18242552380635635,
      0.23147521650098238,
      0.1544652650835347
    )
  )

  d <- golden_plogit_data()
  theta <- c(0.4, -1.5, 0.3, -0.2)

  for (measure in golden_measures()) {
    key <- paste0("disjoint/", measure)
    result <- plogit_predict(
      theta,
      time = d$time,
      event = d$event,
      X = d$X,
      measure = measure
    )
    expect_identical(dim(result), c(3L, 5L), label = paste0(key, ": dim"))
    expect_golden(result, golden[[key]], key)
  }

  # unique_times supplied: the time grid is the one given rather than the one
  # read off the event times, so theta carries one more time parameter.
  result <- plogit_predict(
    c(0.4, -1.5, 0.3, -0.2, 0.1),
    time = d$time,
    event = d$event,
    X = d$X,
    unique_times = c(1, 2, 3, 4),
    measure = "survival"
  )
  expect_identical(dim(result), c(4L, 5L))
  expect_golden(result, golden[["unique_times/survival"]], "unique_times")

  # S supplied: an intercept and a linear time term in place of the disjoint
  # indicators.
  result <- plogit_predict(
    c(0.4, -1.5, 0.25),
    time = d$time,
    event = d$event,
    X = d$X,
    S = cbind(1, seq_len(3)),
    measure = "survival"
  )
  expect_identical(dim(result), c(3L, 5L))
  expect_golden(result, golden[["S/survival"]], "S")

  # times_to_predict covers each arm of the lookup: zero, a time before the
  # first event time, an interior event time, a time between two event times,
  # and the last event time.
  for (measure in c("survival", "hazard")) {
    key <- paste0("times_to_predict/", measure)
    result <- plogit_predict(
      theta,
      time = d$time,
      event = d$event,
      X = d$X,
      times_to_predict = c(0, 0.5, 1, 2.5, 3),
      measure = measure
    )
    expect_identical(dim(result), c(5L, 5L), label = paste0(key, ": dim"))
    expect_golden(result, golden[[key]], key)
  }
})

# ---- the AFT per-time core under exact differentiation ----------------------
#
# aft_measure_at_time() is shared by aft_predictions_individual(), which reaches
# it with plain doubles, and aft_predictions_function(), which reaches it with
# tangent-carrying pairs when deriv_method = "exact". Whether a tangent survives
# the call is invisible to every value-based test in this file: a core that
# quietly stopped propagating tangents would return the same numbers on the
# numeric path, and the exact-mode variance it fed into would drop to zero only
# where the recorded values happen to notice. These tests read the tangent
# directly instead.

test_that("the AFT per-time core carries a tangent through every branch", {
  theta <- golden_aft_theta()
  x_vec <- c(1, 0.6)
  t_val <- 3

  # The tangent-carrying argument auto_differentiation() builds for direction
  # `j`: each parameter is a primal-tangent pair whose tangent is one in that
  # direction and zero in the others.
  direction <- function(j) {
    primal_tangent_vector(lapply(
      seq_along(theta),
      function(i) primal_tangent(theta[i], if (i == j) 1 else 0)
    ))
  }

  # The linear predictor and the scale exactly as aft_predictions_function()
  # forms them, then the shared core.
  measure_at <- function(th, distribution, measure) {
    xbeta <- 0
    for (k in seq_along(x_vec)) {
      xbeta <- xbeta + x_vec[k] * th[k]
    }
    sigma <- if (distribution == "exponential") 1 else exp(-th[length(theta)])
    aft_measure_at_time(t_val, xbeta, sigma, distribution, measure)
  }

  distributions <- c("exponential", "weibull", "log-normal", "log-logistic")
  for (distribution in distributions) {
    for (measure in golden_measures()) {
      for (j in seq_along(theta)) {
        label <- paste(distribution, measure, j, sep = "/")
        exact <- measure_at(direction(j), distribution, measure)

        expect_s3_class(exact, "PrimalTangent")

        # The tangent must be the derivative, not merely present. A central
        # difference of the same core on numeric arguments is accurate to about
        # 1e-9 here, which is far below any of the derivatives involved.
        step <- 1e-6
        up <- theta
        down <- theta
        up[j] <- up[j] + step
        down[j] <- down[j] - step
        finite <- (measure_at(up, distribution, measure) -
          measure_at(down, distribution, measure)) /
          (2 * step)

        expect_equal(exact$tangent, finite, tolerance = 1e-7, label = label)

        # sigma is fixed at one for the exponential distribution, so the third
        # parameter genuinely has no effect there. Everywhere else a zero
        # tangent would mean the derivative was lost rather than absent.
        if (distribution == "exponential" && j == length(theta)) {
          expect_identical(exact$tangent, 0, label = label)
        } else {
          expect_false(exact$tangent == 0, label = label)
        }
      }
    }
  }
})

test_that("the AFT per-time core returns plain doubles for plain arguments", {
  # The same call with numeric arguments must produce an ordinary numeric
  # vector, one entry per person, which is what aft_predictions_individual()
  # writes into a column of its result matrix.
  result <- aft_measure_at_time(
    3,
    c(1.875, 2.37, 2.64),
    exp(-0.3),
    "log-normal",
    "hazard"
  )

  expect_type(result, "double")
  expect_length(result, 3)
  expect_false(is_pt(result))
})
