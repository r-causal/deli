# print() unfitted MEstimator

    Code
      print(m)
    Message
      <MEstimator>
        Parameters: 1
      i Call `estimate()` to fit.

# print() fitted single-parameter MEstimator

    Code
      print(m)
    Message
      <MEstimator>
        Parameters: 1
        Observations: 14
      Coefficients:
      theta_1: 4.3571

# print() fitted multi-parameter MEstimator

    Code
      print(m)
    Message
      <MEstimator>
        Parameters: 2
        Observations: 14
      Coefficients:
      theta_1: 4.3571
      theta_2: 6.6582

# print() fitted named-parameter MEstimator

    Code
      print(m)
    Message
      <MEstimator>
        Parameters: 3
        Observations: 32
      Coefficients:
      intercept: 37.2273
      wt: -3.8778
      hp: -0.0318

# print() unfitted GMMEstimator

    Code
      print(g)
    Message
      <GMMEstimator>
        Parameters: 1
      i Call `estimate()` to fit.

# print() fitted GMMEstimator

    Code
      print(g)
    Message
      <GMMEstimator>
        Parameters: 1
        Observations: 5
      Coefficients:
      mean: 3.0000

# print(summary()) single-parameter output

    Code
      print(summary(m))
    Message
      -- MEstimator Results ----------------------------------------------------------
      Observations: 14
      Parameters: 1
      
                Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
      theta_1     4.3571     0.6896     6.3181     3.0055     5.7088   2.65e-10    <S-value>

# print(summary()) multi-parameter output

    Code
      print(summary(m))
    Message
      -- MEstimator Results ----------------------------------------------------------
      Observations: 14
      Parameters: 2
      
                Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
      theta_1     4.3571     0.6896     6.3181     3.0055     5.7088   2.65e-10    <S-value>
      theta_2     6.6582     2.9533     2.2545     0.8698    12.4466     0.0242     <S-value>

# print(summary()) named-parameter output

    Code
      print(s)
    Message
      -- MEstimator Results ----------------------------------------------------------
      Observations: 32
      Parameters: 3
      
                  Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
      intercept    37.2273     1.9389    19.2000    33.4271    41.0275     <2e-16   <S-value>
      wt           -3.8778     0.6199    -6.2553    -5.0929    -2.6628   3.97e-10    <S-value>
      hp           -0.0318     0.0066    -4.7807    -0.0448    -0.0187   1.75e-06    <S-value>

# print(summary()) GMMEstimator output

    Code
      print(summary(g))
    Message
      -- GMMEstimator Results --------------------------------------------------------
      Observations: 5
      Parameters: 1
      
              Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
      mean      3.0000     0.6325     4.7434     1.7604     4.2396    2.1e-06    <S-value>

# print(summary()) with alpha = 0.10

    Code
      print(summary(m, alpha = 0.1))
    Message
      -- MEstimator Results ----------------------------------------------------------
      Observations: 14
      Parameters: 1
      
                Estimate    Std.Err    Z-score    90% LCL    90% UCL    P-value    S-value
      theta_1     4.3571     0.6896     6.3181     3.2228     5.4915   2.65e-10    <S-value>

# summary() subset row labels track the original parameter index

    Code
      print(s)
    Message
      -- MEstimator Results ----------------------------------------------------------
      Observations: 14
      Parameters: 2
      
                Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
      theta_2     6.6582     2.9533     2.2545     0.8698    12.4466     0.0242     <S-value>

# print(summary()) honors subset

    Code
      print(summary(m, subset = c(1, 3)))
    Message
      -- MEstimator Results ----------------------------------------------------------
      Observations: 32
      Parameters: 3
      
                  Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
      intercept    37.2273     1.9389    19.2000    33.4271    41.0275     <2e-16   <S-value>
      hp           -0.0318     0.0066    -4.7807    -0.0448    -0.0187   1.75e-06    <S-value>

# print() honors subset on the coefficient list

    Code
      print(m, subset = c(1, 3))
    Message
      <MEstimator>
        Parameters: 3
        Observations: 32
      Coefficients:
      intercept: 37.2273
      hp: -0.0318

