# the grid mismatch abort names both counts it compared

    Code
      plogit_predict(f$theta, time = f$time, event = f$event, X = f$X, S = cbind(1, c(
        1, 2, 3)))
    Condition
      Error in `plogit_predict()`:
      ! Dimension mismatch between time intervals and `S`.
      x Found 8 unit-time intervals but `S` has 3 rows.
      i These values must match.

