# Forward-mode automatic differentiation via primal-tangent pairs
#
# PrimalTangent is a simple S3 class (not S7) because R's Ops/Math group
# generics dispatch on S3 classes.  This lets us override +, -, *, /, ^,
# exp, log, sin, cos, etc. without needing special infrastructure.
#
# Design: uses one forward pass per parameter direction with scalar tangents.
# This avoids shape mismatches when PrimalTangent interacts with data vectors.

# ---- PrimalTangent constructor (internal) -----------------------------------

#' @noRd
primal_tangent <- function(primal, tangent) {
  structure(
    list(primal = primal, tangent = tangent),
    class = "PrimalTangent"
  )
}

#' @noRd
is_pt <- function(x) inherits(x, "PrimalTangent")

# ---- PrimalTangentVector: indexable collection of PrimalTangent objects ------
# This allows user functions to use x[1], x[2] syntax.

#' @noRd
primal_tangent_vector <- function(elements) {
  structure(
    list(elements = elements),
    class = "PrimalTangentVector"
  )
}

#' @export
`[.PrimalTangentVector` <- function(x, i) {
  # Select from the underlying list, which resolves positive, negative, and
  # logical indices exactly as base R does. A negative index such as `theta[-p]`
  # (drop the last element, used by the GLM equations to split coefficients from
  # a nuisance parameter) selects several elements, so the single-index shortcut
  # cannot rely on `[[`.
  selected <- x$elements[i]
  # An out-of-range positive subscript selects a NULL slot rather than erroring
  # (list `[` semantics), which would silently produce a zero-row Jacobian under
  # exact differentiation. Reject it with the same base-style message the scalar
  # surface uses. Negative and in-range logical indices never yield NULL slots,
  # so the GLM `theta[-p]` idiom still passes. An over-length logical index or an
  # NA index does select a NULL slot, so those abort here, which is the intended
  # strict behavior.
  if (any(vapply(selected, is.null, logical(1)))) {
    cli::cli_abort("subscript out of bounds", call = NULL)
  }
  if (length(selected) == 1) {
    return(selected[[1]])
  }
  primal_tangent_vector(selected)
}

#' @export
`[[.PrimalTangentVector` <- function(x, i) {
  x$elements[[i]]
}

#' @export
length.PrimalTangentVector <- function(x) length(x$elements)

# ---- indexing, length, dim, and transpose of a PrimalTangent ----------------
# A genuinely scalar PrimalTangent is a single tangent-carrying value, so it
# indexes like a length-1 vector: the sole subscript selects the value itself
# and preserves the class. Without this `theta[1]` on a scalar pair falls
# through to base list indexing and returns the bare primal element, breaking
# downstream arithmetic. The stacked delta-effective-dose equations
# (`ee_emax_ed`, `ee_loglogistic_ed`) index their scalar theta argument as
# `theta[1]`, which is a scalar PrimalTangent under exact mode.
#
# A scalar pair can also carry a vector or matrix payload: `theta[k] * X` (the
# interaction and scaling idiom) multiplies a scalar tangent through a data
# vector or matrix, so both slots become non-scalar. Selecting an element or
# transposing such a pair must operate on the payloads, mirroring the array
# surface, rather than self-selecting the whole pair or falling through to
# `t.default` (which strips the tangent).

# Validate a length-1 self-select and return the pair unchanged, otherwise raise
# a base-style out-of-bounds error.
#' @noRd
pt_scalar_subscript <- function(x, i) {
  if (length(i) != 1L || is.na(i) || i != 1) {
    cli::cli_abort("subscript out of bounds", call = NULL)
  }
  x
}

# Whether a PrimalTangent carries a genuinely scalar payload (a single value
# with no dimensions) rather than a vector or matrix in its slots.
#' @noRd
pt_is_scalar <- function(x) length(x$primal) == 1L && is.null(dim(x$primal))

#' @export
`[.PrimalTangent` <- function(x, i, j, ...) {
  if (pt_is_scalar(x)) {
    return(pt_scalar_subscript(x, i))
  }
  # Non-scalar payload: recycle a broadcast tangent to the primal's shape, then
  # subset both slots identically, exactly as the array surface does. A single
  # index selects column-major; the two-dimensional form propagates a missing i
  # or j to select a whole row or column.
  primal <- x$primal
  tangent <- pt_recycle_tangent(primal, x$tangent)
  if (nargs() <= 2L) {
    p <- primal[i]
    t <- tangent[i]
  } else {
    p <- primal[i, j, drop = FALSE]
    t <- tangent[i, j, drop = FALSE]
  }
  if (length(p) == 1L) {
    return(primal_tangent(as.vector(p), as.vector(t)))
  }
  primal_tangent_array(p, t)
}

#' @export
`[[.PrimalTangent` <- function(x, i, ...) {
  if (pt_is_scalar(x)) {
    return(pt_scalar_subscript(x, i))
  }
  primal <- x$primal
  tangent <- pt_recycle_tangent(primal, x$tangent)
  primal_tangent(primal[[i]], tangent[[i]])
}

#' @export
length.PrimalTangent <- function(x) length(x$primal)

#' @export
dim.PrimalTangent <- function(x) dim(x$primal)

#' @export
t.PrimalTangent <- function(x) {
  tangent <- pt_recycle_tangent(x$primal, x$tangent)
  primal_tangent_array(t(x$primal), t(tangent))
}

# A PrimalTangentVector is a list of independent scalar pairs, so an elementwise
# Math member applies per element. Collapsing to parallel primal/tangent vectors
# and re-dispatching the same generic on the array representation reuses the
# vectorized array rules and preserves each element's own tangent. The extra
# arguments are forwarded so that `log(x, base)` keeps its base rather than
# dropping it on the way to the array method.
#' @export
Math.PrimalTangentVector <- function(x, ...) {
  parts <- pt_arrays(x)
  arr <- primal_tangent_array(parts$primal, parts$tangent)
  get(.Generic)(arr, ...)
}

# ---- c() method for combining PrimalTangent objects -------------------------

#' @export
c.PrimalTangent <- function(...) {
  list(...)
}

# ---- Ops group generic (arithmetic + comparison) ----------------------------

# Whether an Ops result carries a vector or matrix of tangents. It does when
# either operand is a PrimalTangentArray or a PrimalTangentVector; otherwise both
# operands are scalar pairs (or numeric constants) and the result stays a scalar
# PrimalTangent.
#' @noRd
pt_ops_is_array <- function(e1, e2) {
  is_pt_array(e1) ||
    inherits(e1, "PrimalTangentVector") ||
    (!missing(e2) && (is_pt_array(e2) || inherits(e2, "PrimalTangentVector")))
}

# Raise the abort for an Ops member that carries no tangent rule. The message
# names the operator through a local (cli rejects an interpolated expression that
# begins with a dot, so `{.Generic}` cannot be interpolated directly) and the
# wording differs by result shape: a tangent array reports "tangent arrays"
# while a scalar pair reports the {.cls PrimalTangent} class.
#' @noRd
pt_ops_unsupported <- function(generic, array_result, unary) {
  fname <- generic
  if (array_result) {
    if (unary) {
      cli::cli_abort(
        "Unary {.val {fname}} is not supported for tangent arrays."
      )
    } else {
      cli::cli_abort(
        "Operator {.val {fname}} is not supported for tangent arrays."
      )
    }
  } else if (unary) {
    cli::cli_abort(
      "Unary {.val {fname}} is not supported for a {.cls PrimalTangent}."
    )
  } else {
    cli::cli_abort(
      "Operator {.val {fname}} is not supported for a {.cls PrimalTangent}."
    )
  }
}

# Apply a unary Ops member. Reused by the scalar and array Ops methods, which
# are the same function object so that mixing a tangent array with a scalar pair
# resolves to a single compatible method rather than R's incompatible-methods
# fallback.
#' @noRd
pt_ops_unary <- function(e1, generic) {
  array_result <- pt_ops_is_array(e1)
  wrap <- if (array_result) primal_tangent_array else primal_tangent
  a <- pt_arrays(e1)
  switch(
    generic,
    "-" = wrap(-a$primal, -a$tangent),
    "+" = e1,
    pt_ops_unsupported(generic, array_result, unary = TRUE)
  )
}

# Apply a binary Ops member. Both operands are normalized to shape-matched
# primal/tangent slots (a numeric constant contributes a zero tangent), so every
# arithmetic rule is written once and broadcasts across the scalar-scalar,
# array-array, and mixed array-scalar cases identically. Comparisons return the
# primal comparison for control flow.
#' @noRd
pt_ops_binary <- function(e1, e2, generic) {
  array_result <- pt_ops_is_array(e1, e2)
  wrap <- if (array_result) primal_tangent_array else primal_tangent

  a <- pt_arrays(e1)
  b <- pt_arrays(e2)
  p1 <- a$primal
  t1 <- a$tangent
  p2 <- b$primal
  t2 <- b$tangent

  # Whether the exponent carries a tangent decides the power rule to apply
  pow_is_variable <- is_pt(e2) ||
    is_pt_array(e2) ||
    inherits(e2, "PrimalTangentVector")

  switch(
    generic,
    "+" = wrap(p1 + p2, t1 + t2),
    "-" = wrap(p1 - p2, t1 - t2),
    "*" = wrap(p1 * p2, t1 * p2 + p1 * t2),
    "/" = wrap(p1 / p2, (t1 * p2 - p1 * t2) / p2^2),
    "^" = {
      if (pow_is_variable) {
        wrap(p1^p2, p2 * p1^(p2 - 1) * t1 + p1^p2 * log(p1) * t2)
      } else {
        wrap(p1^p2, p2 * p1^(p2 - 1) * t1)
      }
    },
    # Comparisons are non-differentiable; return the primal comparison for
    # control flow (if/else). Tangent information is irrelevant.
    "<" = p1 < p2,
    "<=" = p1 <= p2,
    ">" = p1 > p2,
    ">=" = p1 >= p2,
    "==" = p1 == p2,
    "!=" = p1 != p2,
    pt_ops_unsupported(generic, array_result, unary = FALSE)
  )
}

# `Ops.PrimalTangent` and `Ops.PrimalTangentArray` are assigned the same function
# object below. R's group-generic dispatch refuses to proceed (the
# "Incompatible methods" fallback) when the two operands select different Ops
# methods; sharing one object makes any mix of a scalar pair and a tangent array
# resolve to a single method. The shared implementation reads the operand shapes
# and returns the appropriately shaped result, so a scalar-scalar operation still
# yields a scalar pair.
#' @export
Ops.PrimalTangent <- function(e1, e2) {
  if (missing(e2)) {
    return(pt_ops_unary(e1, .Generic))
  }
  pt_ops_binary(e1, e2, .Generic)
}

# ---- vectorized kink rules for abs, floor, and ceiling ----------------------

# Tangent of abs(p): sign(p) * t elementwise, where the derivative is undefined
# at p == 0 and reported as NaN (R's sign(0) is 0, so the kink is special-cased
# rather than falling out of the multiplication). Vectorized so a vector primal
# or a matrix tangent array selects the branch per element.
#' @noRd
abs_tangent <- function(p, t) {
  tangent <- sign(p) * t # +t where p > 0, -t where p < 0
  tangent[p == 0] <- NaN # derivative undefined at the kink
  tangent
}

# Tangent of floor(p) or ceiling(p): 0 away from the integers and NaN on them,
# where the step has an undefined derivative. Mirrors Python's `0 * tangent` and
# `np.nan * tangent`, so a NaN or infinite incoming tangent still propagates.
# Only a finite integer primal takes the NaN (kink) branch: Python evaluates
# `nan % 1 == 0` and `inf % 1 == 0` as False, so a NaN or infinite primal takes
# the 0-tangent branch. Guarding with `is.finite(p)` reproduces that, where R's
# `NaN %% 1 == 0` would otherwise be NA and leak an NA tangent.
#' @noRd
step_tangent <- function(p, t) {
  ifelse(is.finite(p) & p %% 1 == 0, NaN, 0) * t # NaN on finite integers, 0 elsewhere
}

# ---- shared elementwise Math rules ------------------------------------------

# Compute the primal and tangent slots for a Math group member. Every supported
# rule is elementwise, so the same formula applies whether `p` and `t` are
# scalars (a PrimalTangent), a vector, or a matrix (a PrimalTangentArray). The
# scalar and array Math methods both dispatch here, which keeps the two surfaces
# in sync: a member added here is immediately available to both. Returns a list
# with `primal` and `tangent`, or `NULL` for a member that is not supported so
# the caller can raise its own shape-appropriate abort.
#' @noRd
pt_math_rule <- function(generic, p, t) {
  switch(
    generic,
    "exp" = list(primal = exp(p), tangent = exp(p) * t),
    "log" = list(primal = log(p), tangent = t / p),
    "log2" = list(primal = log2(p), tangent = t / (p * log(2))),
    "log10" = list(primal = log10(p), tangent = t / (p * log(10))),
    "sqrt" = list(primal = sqrt(p), tangent = t / (2 * sqrt(p))),
    "sin" = list(primal = sin(p), tangent = cos(p) * t),
    "cos" = list(primal = cos(p), tangent = -sin(p) * t),
    "tan" = list(primal = tan(p), tangent = t / cos(p)^2),
    "asin" = list(primal = asin(p), tangent = t / sqrt(1 - p^2)),
    "acos" = list(primal = acos(p), tangent = -t / sqrt(1 - p^2)),
    "atan" = list(primal = atan(p), tangent = t / (1 + p^2)),
    "sinh" = list(primal = sinh(p), tangent = cosh(p) * t),
    "cosh" = list(primal = cosh(p), tangent = sinh(p) * t),
    "tanh" = list(primal = tanh(p), tangent = t / cosh(p)^2),
    "asinh" = list(primal = asinh(p), tangent = t / sqrt(p^2 + 1)),
    "acosh" = list(primal = acosh(p), tangent = t / sqrt(p^2 - 1)),
    "atanh" = list(primal = atanh(p), tangent = t / (1 - p^2)),
    "log1p" = list(primal = log1p(p), tangent = t / (1 + p)),
    "expm1" = list(primal = expm1(p), tangent = exp(p) * t),
    "lgamma" = list(primal = lgamma(p), tangent = digamma(p) * t),
    "digamma" = list(primal = digamma(p), tangent = trigamma(p) * t),
    "trigamma" = list(
      primal = trigamma(p),
      tangent = psigamma(p, deriv = 2) * t
    ),
    # cumsum is linear, so both slots run the same cumulative sum
    "cumsum" = list(primal = cumsum(p), tangent = cumsum(t)),
    "abs" = list(primal = abs(p), tangent = abs_tangent(p, t)),
    "floor" = list(primal = floor(p), tangent = step_tangent(p, t)),
    "ceiling" = list(primal = ceiling(p), tangent = step_tangent(p, t)),
    # sign is piecewise constant, so its derivative is zero away from the
    # jump at zero and treated as zero throughout, mirroring numpy's np.sign
    # on an object array of primal-tangent pairs (it collapses to a plain
    # integer sign, carrying no tangent). The penalized regression equations
    # rely on this to differentiate the bridge penalty exactly.
    "sign" = list(primal = sign(p), tangent = 0 * t),
    NULL
  )
}

# Resolve a Math group member to its primal/tangent rule, honoring the optional
# base that R threads through `...` for `log(x, base)`. Every other member is
# elementwise and delegates to `pt_math_rule`; `log` with a supplied base needs
# its own rule because the natural-log tangent `t / p` must be scaled by
# `1 / log(base)` (the primal likewise uses the base). Without this, `log(x, 10)`
# silently differentiates as the natural log, dropping the base entirely.
#' @noRd
pt_math_apply <- function(generic, p, t, dots) {
  if (generic == "log" && length(dots) > 0) {
    base <- dots[[1]]
    # The base is treated as a constant: the tangent scales the natural-log rule
    # by 1 / log(base) but carries no term for the base itself. A tangent-carrying
    # base would need that extra term, so rather than silently return the
    # constant-base derivative, abort. Differentiating with respect to a variable
    # log base is unsupported.
    if (is_tangent_container(base)) {
      cli::cli_abort(
        "{.fn log} with a tangent-carrying base is not supported."
      )
    }
    return(list(primal = log(p, base), tangent = t / (p * log(base))))
  }
  pt_math_rule(generic, p, t)
}

# ---- Math group generic (exp, log, sqrt, trig, abs, etc.) -------------------

#' @export
Math.PrimalTangent <- function(x, ...) {
  rule <- pt_math_apply(.Generic, x$primal, x$tangent, list(...))
  if (is.null(rule)) {
    # Bind to a local first: cli rejects an interpolated expression that begins
    # with a dot, so `{.Generic}` cannot be interpolated directly.
    fname <- .Generic
    cli::cli_abort(
      "Math function {.val {fname}} is not supported for a {.cls PrimalTangent}."
    )
  }
  primal_tangent(rule$primal, rule$tangent)
}

# ---- Summary group generic (sum, prod, max, min) ----------------------------

# Collect every argument of a Summary call into parallel primal/tangent vectors.
# A tangent-carrying argument (a scalar pair, a PrimalTangentVector, or a
# PrimalTangentArray) contributes its flattened primal and tangent, with a
# scalar broadcast tangent recycled to the primal length; a numeric constant
# contributes a zero tangent. Flattening every container to the same
# representation lets sum, prod, max, and min share one implementation across all
# three tangent surfaces.
#' @noRd
pt_summary_parts <- function(args) {
  all_p <- numeric(0)
  all_t <- numeric(0)
  for (a in args) {
    if (is_tangent_container(a)) {
      parts <- pt_flatten(a)
      all_p <- c(all_p, parts$primal)
      all_t <- c(all_t, parts$tangent)
    } else {
      all_p <- c(all_p, a)
      all_t <- c(all_t, rep(0, length(a)))
    }
  }
  list(primal = all_p, tangent = all_t)
}

# Apply a Summary member to the pooled primal/tangent vectors of its arguments.
# Every reduction returns a scalar pair: sum follows the linear rule, prod the
# product rule, and max and min carry the primal and tangent of the selected
# element.
#' @noRd
pt_summary <- function(generic, args, na.rm) {
  parts <- pt_summary_parts(args)
  ap <- parts$primal
  at <- parts$tangent
  # Honor na.rm for every reduction, not only sum. Base Summary decides which
  # elements to drop from the primal, so the same positions are dropped from the
  # tangent to keep the two slots aligned. Dropping by primal NA also fixes sum,
  # where the previous independent sum(at, na.rm) could keep a tangent whose
  # primal had been removed.
  if (na.rm) {
    keep <- !is.na(ap)
    ap <- ap[keep]
    at <- at[keep]
  }
  switch(
    generic,
    "sum" = primal_tangent(sum(ap), sum(at)),
    "prod" = {
      # prod(ap[-i]) is the product with the i-th factor left out
      partials <- vapply(
        seq_along(ap),
        function(i) prod(ap[-i]) * at[i],
        numeric(1)
      )
      primal_tangent(prod(ap), sum(partials))
    },
    "max" = {
      idx <- which.max(ap)
      primal_tangent(ap[idx], at[idx])
    },
    "min" = {
      idx <- which.min(ap)
      primal_tangent(ap[idx], at[idx])
    },
    cli::cli_abort(
      "Summary function {.val {generic}} is not supported for a {.cls PrimalTangent}."
    )
  )
}

#' @export
Summary.PrimalTangent <- function(..., na.rm = FALSE) {
  pt_summary(.Generic, list(...), na.rm)
}

# The same function object as `Summary.PrimalTangent`. A Summary reduction on the
# whole parameter vector (`sum(theta)`) or on a tangent-carrying array
# (`sum(X %*% theta)`) must differentiate too, and `pt_summary` already flattens
# every container to a common representation, so one implementation serves all
# three surfaces. `.Generic` is bound at dispatch time, so sharing the object
# still resolves sum, prod, max, and min correctly.
#' @export
Summary.PrimalTangentArray <- Summary.PrimalTangent

#' @export
Summary.PrimalTangentVector <- Summary.PrimalTangent

# ---- Boolean coercion (for if/else) -----------------------------------------

#' @export
as.logical.PrimalTangent <- function(x, ...) as.logical(x$primal)

#' @export
as.double.PrimalTangent <- function(x, ...) as.double(x$primal)

#' @export
as.numeric.PrimalTangent <- function(x, ...) as.numeric(x$primal)

# ---- PrimalTangentArray: primal/tangent parallel numeric arrays -------------
# A tangent-carrying vector or matrix within a single forward pass. Both slots
# hold plain numeric objects of identical shape: `primal` is the value and
# `tangent` is its directional derivative for the current pass. Because the
# tangent is scalar-per-element, every linear operation (matrix product,
# transpose, binding, reshaping, indexing) applies the same operation to both
# slots, and elementwise arithmetic applies the usual differentiation rules.

#' @noRd
primal_tangent_array <- function(primal, tangent) {
  structure(
    list(primal = primal, tangent = tangent),
    class = "PrimalTangentArray"
  )
}

#' @noRd
is_pt_array <- function(x) inherits(x, "PrimalTangentArray")

# ---- extracting primal/tangent slots from any operand -----------------------

# Return the primal and tangent slots of an operand as shape-matched arrays.
# A numeric constant contributes a zero tangent of the same shape.
#' @noRd
pt_arrays <- function(x) {
  if (is_pt_array(x)) {
    return(list(primal = x$primal, tangent = x$tangent))
  }
  if (inherits(x, "PrimalTangentVector")) {
    # Collapse the list of scalar pairs into two numeric vectors
    primal <- vapply(x$elements, function(e) e$primal, numeric(1))
    tangent <- vapply(x$elements, function(e) e$tangent, numeric(1))
    return(list(primal = primal, tangent = tangent))
  }
  if (is_pt(x)) {
    return(list(primal = x$primal, tangent = x$tangent))
  }
  if (is.list(x)) {
    # A plain list of scalar pairs (produced by `c()` on PrimalTangent objects).
    # pt_flatten already collapses this shape into parallel numeric vectors, so
    # a masked binder of such a list keeps its tangents instead of hitting the
    # numeric-constant fallback below and zeroing them.
    return(pt_flatten(x))
  }
  # Numeric constant: zero tangent, matching the constant's shape
  zero <- x
  zero[] <- 0
  list(primal = x, tangent = zero)
}

# Recycle a scalar-broadcast tangent up to the primal's length. A scalar pair
# built from a vector primal and a constant carries a length-1 tangent that
# broadcasts under elementwise arithmetic (R recycles it), but flattening and
# binding reshape the primal and tangent slots independently, so the tangent
# must first be expanded to match the primal element for element.
#' @noRd
pt_recycle_tangent <- function(primal, tangent) {
  if (length(tangent) < length(primal)) {
    tangent <- rep_len(tangent, length(primal))
    dim(tangent) <- dim(primal)
  }
  tangent
}

# Flatten a tangent-carrying container into parallel numeric vectors. Used when
# reshaping into a matrix. Accepts a PrimalTangent (a scalar pair, possibly with
# a vector primal and a scalar broadcast tangent), a PrimalTangentVector, a
# PrimalTangentArray, or a plain list of scalar pairs (produced by `c()` on
# PrimalTangent objects).
#' @noRd
pt_flatten <- function(x) {
  if (is_pt_array(x)) {
    return(list(primal = as.vector(x$primal), tangent = as.vector(x$tangent)))
  }
  if (is_pt(x)) {
    # A single scalar pair whose primal may be a length-n vector. Recycle a
    # scalar broadcast tangent to the primal length before flattening.
    primal <- as.vector(x$primal)
    tangent <- as.vector(pt_recycle_tangent(x$primal, x$tangent))
    return(list(primal = primal, tangent = tangent))
  }
  if (inherits(x, "PrimalTangentVector")) {
    x <- x$elements
  }
  # x is now a plain list of scalar pairs (or numeric constants)
  primal <- vapply(
    x,
    function(e) if (is_pt(e)) e$primal else as.numeric(e),
    numeric(1)
  )
  tangent <- vapply(x, function(e) if (is_pt(e)) e$tangent else 0, numeric(1))
  list(primal = primal, tangent = tangent)
}

# Coerce a matrix-shaped linear predictor (typically `X %*% theta`, an n-by-1
# result) to a length-n vector. A plain numeric reduces to `as.numeric()`
# exactly, so the finite-difference and undifferentiated passes are byte
# unchanged. A tangent-carrying result flattens both slots and returns a
# PrimalTangentArray, so the exact pass keeps its tangents rather than stripping
# them the way `as.numeric()` on a tangent object would.
#' @noRd
pt_as_vector <- function(x) {
  if (is_tangent_container(x)) {
    parts <- pt_flatten(x)
    return(primal_tangent_array(parts$primal, parts$tangent))
  }
  as.numeric(x)
}

# Coerce a regression design matrix argument, preserving tangents when the
# matrix is theta-derived. The second-stage design of `ee_2sls` stacks predicted
# treatment (a tangent-carrying column) with the exogenous covariates, so under
# exact mode `X` reaches `ee_regression` as a PrimalTangentArray. A plain numeric
# input reduces to `as.matrix()` exactly, so the numeric and finite-difference
# passes are byte unchanged; a tangent-carrying input passes through so the
# downstream matrix product, residual, and score keep their tangents.
#' @noRd
coerce_design <- function(X) {
  # A plain numeric matrix is the common case on a numeric-mode fit, where
  # coerce_design runs on every solver and Jacobian evaluation. A matrix is
  # never a tangent container (those are classed lists with no dim attribute),
  # so returning it directly skips both the tangent check and the as.matrix
  # redispatch, each of which would hand back the same object. Dimnames are
  # dropped so the returned score matrix is unnamed regardless of whether the
  # caller passed a plain matrix, a named matrix, or a data frame; the Python
  # reference operates on unnamed arrays, and the score dimnames otherwise leak
  # a data frame's column names into the estimating-function output.
  if (is.matrix(X)) {
    if (!is.null(dimnames(X))) {
      dimnames(X) <- NULL
    }
    return(X)
  }
  if (is_tangent_container(X)) {
    return(X)
  }
  X <- as.matrix(X)
  dimnames(X) <- NULL
  X
}

# Coerce a regression outcome argument, preserving tangents when the response is
# theta-derived. The efficient `ee_gestimation_snmm` fits an outcome model whose
# response h(phi) = Y - (V*A) phi carries tangents under exact mode, so `y`
# reaches `ee_regression` as a tangent container. A plain numeric input reduces
# to `as.numeric()` exactly; a tangent-carrying input flattens to a
# PrimalTangentArray vector so the residual keeps its tangent.
#' @noRd
coerce_outcome <- function(y) {
  # A plain, attribute-free double vector is the common numeric-mode outcome and
  # needs no coercion: as.numeric() would return it unchanged. Skipping the
  # tangent check and the redispatch avoids repeating that work on every
  # evaluation. A tangent-carrying response still flattens to a
  # PrimalTangentArray; any other numeric type still routes through as.numeric().
  if (is.double(y) && is.null(attributes(y))) {
    return(y)
  }
  if (is_tangent_container(y)) {
    return(pt_as_vector(y))
  }
  as.numeric(y)
}

# Is `x` something a tangent-aware matrix() should reshape?
#' @noRd
is_tangent_container <- function(x) {
  is_pt(x) ||
    inherits(x, "PrimalTangentVector") ||
    is_pt_array(x) ||
    (is.list(x) && length(x) > 0 && any(vapply(x, is_pt, logical(1))))
}

# ---- matrix products --------------------------------------------------------

# Differentiate a matrix product of two operands, either of which may carry
# tangents. The primal is the ordinary product; the tangent follows the matrix
# product rule d(AB) = (dA)B + A(dB), where a constant operand has zero tangent.
#' @noRd
pt_matmul <- function(x, y) {
  a <- pt_arrays(x)
  b <- pt_arrays(y)
  primal <- a$primal %*% b$primal # ordinary product
  tangent <- a$tangent %*% b$primal + a$primal %*% b$tangent
  primal_tangent_array(primal, tangent)
}

#' @export
`%*%.PrimalTangentArray` <- function(x, y) pt_matmul(x, y)

#' @export
`%*%.PrimalTangentVector` <- function(x, y) pt_matmul(x, y)

# A length-1 parameter subset (a single-column structural mean model in
# `ee_gestimation_snmm`, or the no-exogenous-covariate second stage in
# `ee_2sls`) indexes a PrimalTangentVector down to a lone scalar pair rather than
# a length-1 vector, so a matrix product `M %*% phi` reaches this method with a
# scalar PrimalTangent operand. Delegating to `pt_matmul` treats the scalar as a
# 1-by-1 factor: base `%*%` conforms it against `M`, and the tangent follows the
# same product rule, so the result matches the multi-parameter path's shape and
# tangents exactly.
#' @export
`%*%.PrimalTangent` <- function(x, y) pt_matmul(x, y)

# ---- transpose --------------------------------------------------------------

#' @export
t.PrimalTangentArray <- function(x) {
  primal_tangent_array(t(x$primal), t(x$tangent))
}

#' @export
t.PrimalTangentVector <- function(x) {
  a <- pt_arrays(x)
  primal_tangent_array(t(a$primal), t(a$tangent))
}

# ---- rbind / cbind ----------------------------------------------------------
# `rbind()` and `cbind()` are masked within the package namespace (rather than
# registered as S3 methods) so that a differentiated function reaches them under
# both a tangent-carrying argument (the exact pass) and a plain numeric argument
# (the finite-difference pass). This keeps the two passes shape- and
# name-consistent. For numeric arguments they forward to the base binder with
# `deparse.level = 0`, which suppresses the symbol-derived row and column labels
# that would otherwise leak into a Jacobian while preserving any explicit labels.

# Bind operands by applying the base binder to the primal slots and, in
# parallel, to the tangent slots. A scalar-broadcast tangent is first recycled
# to its primal's shape so the two slots bind to identical layouts (a bare
# scalar pair with a vector primal and a scalar tangent would otherwise bind a
# length-n primal against a length-1 tangent).
#' @noRd
pt_bind <- function(args, bind_fn) {
  slots <- lapply(args, function(a) {
    parts <- pt_arrays(a)
    parts$tangent <- pt_recycle_tangent(parts$primal, parts$tangent)
    parts
  })
  primal <- do.call(bind_fn, lapply(slots, `[[`, "primal"))
  tangent <- do.call(bind_fn, lapply(slots, `[[`, "tangent"))
  primal_tangent_array(primal, tangent)
}

#' @noRd
rbind <- function(..., deparse.level = 1) {
  args <- list(...)
  if (any(vapply(args, is_tangent_container, logical(1)))) {
    return(pt_bind(args, base::rbind))
  }
  base::rbind(..., deparse.level = 0)
}

#' @noRd
cbind <- function(..., deparse.level = 1) {
  args <- list(...)
  if (any(vapply(args, is_tangent_container, logical(1)))) {
    return(pt_bind(args, base::cbind))
  }
  base::cbind(..., deparse.level = 0)
}

# ---- indexing ---------------------------------------------------------------

#' @export
`[.PrimalTangentArray` <- function(x, i, j, ...) {
  # Distinguish M[i] from M[i, j] by argument count, the base-style
  # matrix-method idiom: `[`(x, i) has two arguments, while M[i, j], M[i, ],
  # and M[, j] all have three. A single index selects column-major, matching
  # base R; the two-dimensional form lets a missing i or j propagate through to
  # select a whole column or row.
  if (nargs() <= 2) {
    primal <- x$primal[i]
    tangent <- x$tangent[i]
  } else {
    primal <- x$primal[i, j, drop = FALSE]
    tangent <- x$tangent[i, j, drop = FALSE]
  }
  if (length(primal) == 1) {
    return(primal_tangent(as.vector(primal), as.vector(tangent)))
  }
  primal_tangent_array(primal, tangent)
}

#' @export
`[[.PrimalTangentArray` <- function(x, i) {
  primal_tangent(x$primal[[i]], x$tangent[[i]])
}

#' @export
length.PrimalTangentArray <- function(x) length(x$primal)

#' @export
dim.PrimalTangentArray <- function(x) dim(x$primal)

# ---- Ops group generic (elementwise arithmetic on tangent arrays) -----------

# The same function object as `Ops.PrimalTangent` (see the note there): sharing
# one method makes a tangent array and a scalar pair compatible operands under
# R's group-generic dispatch. The shared implementation returns a
# PrimalTangentArray whenever either operand carries a vector or matrix of
# tangents, which covers `(X %*% beta) / theta[k]` and the reverse order.
#' @export
Ops.PrimalTangentArray <- Ops.PrimalTangent

# The same function object again (see the note on `Ops.PrimalTangent`). A
# PrimalTangentVector is the whole parameter vector during an exact pass, and
# some estimating equations operate on it directly (the penalized regression
# equations form `theta - center` before taking the penalty). Sharing the one
# method lets those whole-vector operations differentiate exactly and stay
# compatible when a vector meets a scalar pair or a tangent array.
#' @export
Ops.PrimalTangentVector <- Ops.PrimalTangent

# ---- Math group generic (elementwise math on tangent arrays) ----------------

# Every Math member is elementwise, so a tangent array reuses the shared scalar
# rules (`pt_math_rule`) applied to its vector or matrix slots. Sharing the rule
# table keeps the scalar and array surfaces identical: a member wired for one is
# available to the other, and a member absent from both aborts here with a
# rendered cli message.
#' @export
Math.PrimalTangentArray <- function(x, ...) {
  rule <- pt_math_apply(.Generic, x$primal, x$tangent, list(...))
  if (is.null(rule)) {
    # Bind to a local first: cli rejects an interpolated expression that begins
    # with a dot, so `{.Generic}` cannot be interpolated directly.
    fname <- .Generic
    cli::cli_abort(
      "Math function {.val {fname}} is not supported for tangent arrays."
    )
  }
  primal_tangent_array(rule$primal, rule$tangent)
}

# ---- tangent-aware matrix() construction ------------------------------------

# `matrix()` is masked within the package namespace (it is not an S3 generic).
# Reshaping a tangent-carrying vector produces a PrimalTangentArray, laid out
# column-major exactly as base::matrix would; numeric inputs forward to
# base::matrix unchanged. Masking (rather than injecting into the differentiated
# function's scope) means a differentiated function reaches this under both the
# exact and the finite-difference pass.
#' @noRd
matrix <- function(data = NA, ...) {
  if (!is_tangent_container(data)) {
    return(base::matrix(data, ...))
  }
  parts <- pt_flatten(data)
  # Reshape the primal and tangent slots identically to keep column-major order
  primal <- base::matrix(parts$primal, ...)
  tangent <- base::matrix(parts$tangent, ...)
  primal_tangent_array(primal, tangent)
}

# ---- tangent-aware replication and column sums ------------------------------
# `rep()` and `colSums()` are masked within the package namespace (like matrix,
# rbind, and cbind) so a reshaping function reaches them under both the exact
# pass (a tangent-carrying argument) and the finite-difference pass (a plain
# numeric argument). Numeric arguments forward to the base functions
# bit-identically, so the undifferentiated and finite-difference paths are
# unchanged. Both operations are linear, so a tangent-carrying argument applies
# the same operation to the primal and tangent slots in parallel.

# Replicate a tangent-carrying vector, mirroring numpy's np.tile over the primal
# values. The pooled logistic equation stacks the covariate log-odds into a
# K-by-n matrix with `rep(log_odds_w, each = K)`; under exact mode log_odds_w
# carries tangents, so the primal and tangent slots are flattened and repeated
# with identical `each`/`times`/`length.out` arguments to keep them aligned.
#' @noRd
rep <- function(x, ...) {
  if (!is_tangent_container(x)) {
    return(base::rep(x, ...))
  }
  parts <- pt_flatten(x)
  primal <- base::rep(parts$primal, ...)
  tangent <- base::rep(parts$tangent, ...)
  primal_tangent_array(primal, tangent)
}

# Sum a tangent-carrying matrix down its columns, mirroring numpy's
# `np.dot(ones, residual_matrix)` reduction across the time intervals in the
# pooled logistic equation. Column summation is linear, so the primal and
# tangent slots are summed identically.
#' @noRd
colSums <- function(x, ...) {
  if (!is_tangent_container(x)) {
    return(base::colSums(x, ...))
  }
  primal_tangent_array(
    base::colSums(x$primal, ...),
    base::colSums(x$tangent, ...)
  )
}

# ---- tangent-aware parallel minimum and maximum -----------------------------
# `pmin()` and `pmax()` are masked within the package namespace (like matrix,
# rbind, and cbind) so a differentiated function reaches them under both the
# exact pass (a tangent-carrying argument) and the finite-difference pass (a
# plain numeric argument). All-numeric arguments forward to the base parallel
# extremes, so the undifferentiated and finite-difference paths are bit
# unchanged. A tangent-carrying argument selects, per element, the primal and
# tangent on the chosen side, mirroring numpy's np.maximum and np.minimum on an
# object array of primal-tangent pairs: the tangent follows the selected
# element, and at a tie the first argument is taken (numpy's object maximum
# compares with `>=` and its object minimum with `<=`).

# Combine two operands elementwise, keeping the primal and tangent of the side
# selected by `choose_max`. A numeric operand contributes a zero tangent, and a
# scalar operand broadcasts up to the longer length.
#' @noRd
pt_parallel_extreme <- function(e1, e2, choose_max) {
  a <- pt_arrays(e1)
  b <- pt_arrays(e2)
  n <- max(length(a$primal), length(b$primal))
  p1 <- rep_len(a$primal, n)
  t1 <- rep_len(a$tangent, n)
  p2 <- rep_len(b$primal, n)
  t2 <- rep_len(b$tangent, n)
  take_first <- if (choose_max) p1 >= p2 else p1 <= p2
  primal <- base::ifelse(take_first, p1, p2)
  tangent <- base::ifelse(take_first, t1, t2)
  primal_tangent_array(primal, tangent)
}

#' @noRd
pmax <- function(..., na.rm = FALSE) {
  args <- list(...)
  if (!any(vapply(args, is_tangent_container, logical(1)))) {
    return(base::pmax(..., na.rm = na.rm))
  }
  Reduce(function(e1, e2) pt_parallel_extreme(e1, e2, choose_max = TRUE), args)
}

#' @noRd
pmin <- function(..., na.rm = FALSE) {
  args <- list(...)
  if (!any(vapply(args, is_tangent_container, logical(1)))) {
    return(base::pmin(..., na.rm = na.rm))
  }
  Reduce(function(e1, e2) pt_parallel_extreme(e1, e2, choose_max = FALSE), args)
}

# ---- tangent-aware conditional selection ------------------------------------
# `ifelse()` is masked within the package namespace (like matrix, rbind, cbind,
# pmin, and pmax) so a branch-selecting function reaches it under both passes.
# With no tangent-carrying branch it forwards to base::ifelse, so the
# undifferentiated and finite-difference paths are bit unchanged. When a branch
# carries tangents, the primal and tangent of the branch the test selects are
# taken per element, mirroring numpy's np.where on an object array of
# primal-tangent pairs: both branches are evaluated and the tangent follows the
# selected branch.

#' @noRd
ifelse <- function(test, yes, no) {
  if (
    !is_tangent_container(test) &&
      !is_tangent_container(yes) &&
      !is_tangent_container(no)
  ) {
    return(base::ifelse(test, yes, no))
  }
  pt_where(test, yes, no)
}

# Select between the `yes` and `no` branches elementwise, carrying the primal
# and tangent of the chosen branch. A numeric branch contributes a zero
# tangent, and a scalar branch broadcasts up to the length of `test`.
#' @noRd
pt_where <- function(test, yes, no) {
  if (is_tangent_container(test)) {
    test <- pt_arrays(test)$primal
  }
  test <- as.logical(test)
  n <- length(test)
  y <- pt_arrays(yes)
  z <- pt_arrays(no)
  p_yes <- rep_len(y$primal, n)
  t_yes <- rep_len(y$tangent, n)
  p_no <- rep_len(z$primal, n)
  t_no <- rep_len(z$tangent, n)
  primal <- base::ifelse(test, p_yes, p_no)
  tangent <- base::ifelse(test, t_yes, t_no)
  primal_tangent_array(primal, tangent)
}

# ---- auto_differentiation ---------------------------------------------------

#' Forward-mode automatic differentiation
#'
#' Computes the exact Jacobian of a function using forward-mode automatic
#' differentiation via primal-tangent pairs. Uses one forward pass per
#' parameter direction with scalar tangents, which correctly handles
#' interactions between parameters and data vectors.
#'
#' @param theta Numeric vector of parameter values at which to evaluate the
#'   Jacobian.
#' @param f A function that takes a numeric vector and returns a numeric
#'   vector (or scalar). The function may use standard arithmetic operators and
#'   math functions (`+`, `-`, `*`, `/`, `^`, `exp`, `log`, `sqrt`, `sin`,
#'   `cos`, etc.). The matrix product `%*%`, the transpose `t()`, and
#'   two-dimensional indexing differentiate exactly in any function. The
#'   reshaping, binding, and reduction operations `matrix()`, `rbind()`,
#'   `cbind()`, `rep()`, and `colSums()` differentiate exactly for functions
#'   evaluated within the package (such as the built-in estimating equations),
#'   because those masked forms are only in
#'   scope there; a user-defined function that calls them from the global
#'   environment reaches the base versions, and the differentiation aborts
#'   rather than returning a silent approximation. `log(x, base)` differentiates
#'   with respect to `x` only; the `base` argument is treated as a constant, and
#'   a `base` that itself carries a tangent (a value derived from `theta`) is not
#'   supported and aborts rather than dropping the base's contribution.
#'
#' @returns A matrix where element `[i, j]` is the partial derivative of
#'   the `i`-th output with respect to the `j`-th parameter.
#'
#' @keywords internal
auto_differentiation <- function(theta, f) {
  theta <- as.numeric(theta)
  p <- length(theta)

  jacobian_cols <- vector("list", p)

  for (j in seq_len(p)) {
    # Create PrimalTangent input for direction j (scalar tangents)
    elements <- vector("list", p)
    for (i in seq_len(p)) {
      elements[[i]] <- primal_tangent(theta[i], if (i == j) 1 else 0)
    }
    pt_input <- primal_tangent_vector(elements)

    # Evaluate the function with PrimalTangent inputs and read off direction j
    result <- f(pt_input)
    jacobian_cols[[j]] <- extract_tangent_column(result)
  }

  do.call(cbind, jacobian_cols)
}

# Read the tangent components of an evaluated function result into one Jacobian
# column, covering every shape a differentiated function can return.
#' @noRd
extract_tangent_column <- function(result) {
  if (is_pt(result)) {
    # Single output (scalar, or a vector primal with a broadcast tangent)
    return(result$tangent)
  }
  if (is_pt_array(result)) {
    # Matrix or vector output, flattened column-major to match c()/indexing
    return(as.vector(result$tangent))
  }
  if (inherits(result, "PrimalTangentVector")) {
    return(vapply(
      result$elements,
      function(e) if (is_pt(e)) e$tangent else 0,
      numeric(1)
    ))
  }
  if (is.list(result)) {
    # A list result normally comes from c.PrimalTangent and holds at least one
    # tangent-carrying pair. A list (or list matrix) with no PrimalTangent
    # elements instead signals that a base reshaping helper such as rbind() or
    # cbind() was reached from outside the package namespace and stripped the
    # tangents. Aborting keeps the documented safety property rather than
    # returning a silent all-zero column.
    if (!any(vapply(result, is_pt, logical(1)))) {
      cli::cli_abort(
        c(
          "Automatic differentiation received a result that carries no tangents.",
          "i" = "A reshaping or binding function such as {.fn rbind} or
                 {.fn cbind} was likely reached from outside the package
                 namespace, which strips the derivative information.",
          "i" = "Combine tangent-carrying values with {.fn c} instead."
        )
      )
    }
    # Multiple outputs from c.PrimalTangent
    return(vapply(
      result,
      function(e) if (is_pt(e)) e$tangent else 0,
      numeric(1)
    ))
  }
  # Plain numeric: constant output, all derivatives are 0
  rep(0, length(result))
}
