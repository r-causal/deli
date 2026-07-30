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

# ---- exact-mode failure conditions ------------------------------------------
# The aborts raised while differentiating carry condition classes, so a test or
# a calling function can match on the class rather than on the prose. These are
# the first custom condition classes in the package; every other abort in deli
# is classless and matched by message.
#
#   deli_exact_tangent_lost
#     Derivative information is gone or would be. Either a tangent-carrying
#     value was asked to become a plain double, or a result arrived with no
#     tangent to read and evidence that one had been dropped on the way.
#
#   deli_exact_unsupported_function
#     A function reached with a tangent-carrying argument has no rule under
#     exact differentiation, either because it hands its argument to compiled
#     code without dispatching or because deli declines to differentiate it.
#
#   deli_exact_unsupported_shape
#     A result still carries its tangents but arrives in a container the
#     summing step has no rule for. Nothing was lost; the shape is the problem.

# The remedy every exact-mode abort ends on. Finite differences carry no
# tangents, so switching the derivative method resolves all of them.
#' @noRd
pt_capprox_hint <- function() {
  c(
    "i" = "Or pass {.val capprox} to {.arg deriv_method} to differentiate with
           finite differences, which need no tangent support."
  )
}

# Why deli declines to differentiate a selection among order statistics. Shared
# by the `median()` and `quantile()` methods and by the trimmed `mean()`.
#' @noRd
pt_order_statistic_hint <- function() {
  c(
    "i" = "Selecting among the order statistics of a tangent-carrying value
           differentiates whichever order statistic the current values happen to
           select. That derivative does not estimate the derivative of the
           population quantity, so a sandwich or delta-method variance built
           from it would not be valid. It is the same reason
           {.fn ee_percentile} has no usable bread."
  )
}

#' @noRd
pt_order_statistic_abort <- function(fname) {
  cli::cli_abort(
    c(
      "{.fn {fname}} cannot be differentiated when {.arg deriv_method} is
       {.val exact}.",
      pt_order_statistic_hint(),
      pt_capprox_hint()
    ),
    class = "deli_exact_unsupported_function",
    call = NULL
  )
}

# Raised when a result arrives with no tangents at all. The hint names only
# operations that carry tangents from any environment, which means the
# registered S3 methods. It deliberately does not recommend matrix(), rbind(),
# cbind(), or ifelse(), which are masked inside the namespace and so are
# unavailable to a user working in the global environment.
#' @noRd
pt_tangent_lost_abort <- function() {
  cli::cli_abort(
    c(
      "Automatic differentiation received a result that carries no tangents.",
      "i" = "A reshaping, binding, or selection function such as {.fn rbind},
             {.fn cbind}, {.fn matrix}, or {.fn ifelse} was likely reached from
             outside the package namespace, which strips the derivative
             information.",
      "i" = "Arithmetic, {.code %*%}, {.fn t}, and {.code [} are S3 methods, so
             they carry tangents from any environment. Build a p-by-n return as
             {.code t(X * resid)}.",
      "i" = "For a conditional, write {.code ind * yes + (1 - ind) * no} rather
             than using {.fn ifelse}.",
      pt_capprox_hint()
    ),
    class = "deli_exact_tangent_lost",
    call = NULL
  )
}

# Raised when a result keeps its tangents but arrives in a container the summing
# step has no reduction for. The derivatives are intact, so the remedy is to
# reshape the return rather than to give up on exact differentiation, which is
# why this is a different condition from a tangent loss.
#' @noRd
pt_unsupported_shape_abort <- function() {
  cli::cli_abort(
    c(
      "Automatic differentiation received tangents in a container shape it
       cannot sum across observations.",
      "i" = "A list holding one tangent-carrying value per equation, such as an
             {.fn lapply} over several equations, keeps its derivatives but has
             no supported reduction.",
      "i" = "Return a single tangent-carrying value instead, either built with
             {.code t()} and arithmetic or concatenated with
             {.code do.call(c, ...)}.",
      pt_capprox_hint()
    ),
    class = "deli_exact_unsupported_shape",
    call = NULL
  )
}

# Raised when a plain numeric result holds an NA. A function that does not
# recognize a tangent-carrying argument returns NA for it, and the NA then
# flows through the plain numeric arithmetic that follows, so an NA is the
# evidence separating a dropped tangent from a genuinely constant output.
#' @noRd
pt_na_tangent_abort <- function() {
  cli::cli_abort(
    c(
      "Automatic differentiation received a plain numeric result that contains
       {.code NA} and carries no tangents.",
      "i" = "A function that does not recognize a tangent-carrying argument
             returned {.code NA} for it, so the derivative was lost before the
             result was assembled.",
      "i" = "The operations deli differentiates are listed in
             {.code vignette(\"getting-started\")}.",
      pt_capprox_hint()
    ),
    class = "deli_exact_tangent_lost",
    call = NULL
  )
}

# base R functions that hand a tangent-carrying argument to compiled code,
# paired with the deli function computing the same quantity through the group
# generics. Every distribution function in stats fails the same way, so this is
# a lookup rather than the whole set: a function absent from it is still named
# in the abort, with the generic remedy.
#' @noRd
pt_exact_replacements <- list(
  plogis = "inverse_logit",
  qlogis = "logit",
  pnorm = "standard_normal_cdf",
  dnorm = "standard_normal_pdf",
  psigamma = "deli_polygamma"
)

# The name of the function a failing call invoked, or NULL when the call is not
# a simple or namespace-qualified function call. `stats::plogis(x)` puts the
# `::` call in the function position, so the name is that call's third element.
#' @noRd
pt_failing_function <- function(call) {
  if (!is.call(call)) {
    return(NULL)
  }
  fn <- call[[1]]
  if (is.call(fn) && identical(fn[[1]], as.name("::"))) {
    fn <- fn[[3]]
  }
  if (!is.name(fn)) {
    return(NULL)
  }
  as.character(fn)
}

# The members of the S3 `Math` group generic. They dispatch on class, so a
# tangent-carrying argument to any of them reaches `Math.PrimalTangent`, which
# either applies the rule or raises its own abort. A member that got as far as
# compiled code was therefore handed non-numeric data, not a tangent, and the
# rewrite below would misdiagnose it as a function deli cannot differentiate.
#
# Only `Math` is listed. An `Ops` member reports "non-numeric argument to binary
# operator" and a `Summary` member "invalid 'type' of argument", so neither can
# reach the message test the rewrite is gated on. The names are written out
# rather than read from `methods::getGroupMembers()`, which describes the S4
# groups and puts `round()` and `signif()` in `Math2` even though both dispatch
# S3 `Math` methods. `psigamma()` and every stats distribution function fall
# outside the group and keep their rewrite.
#' @noRd
pt_math_group_members <- c(
  "abs",
  "sign",
  "sqrt",
  "floor",
  "ceiling",
  "trunc",
  "round",
  "signif",
  "exp",
  "log",
  "expm1",
  "log1p",
  "log2",
  "log10",
  "cos",
  "sin",
  "tan",
  "cospi",
  "sinpi",
  "tanpi",
  "acos",
  "asin",
  "atan",
  "cosh",
  "sinh",
  "tanh",
  "acosh",
  "asinh",
  "atanh",
  "lgamma",
  "gamma",
  "digamma",
  "trigamma",
  "cumsum",
  "cumprod",
  "cummax",
  "cummin"
)

# Rewrite the compiled-code error a base R function raises when it receives a
# tangent-carrying argument, naming the offender and its deli replacement. The
# name is read off the failing call rather than matched against a fixed list,
# because every stats distribution function fails identically and a fixed list
# would leave most of them unnamed. Any other error is re-raised unchanged, so
# an unrelated failure inside a differentiated function propagates byte for
# byte, keeping its own class and message.
#' @noRd
pt_rethrow_exact_error <- function(cnd) {
  fname <- pt_failing_function(conditionCall(cnd))
  compiled <- grepl(
    "^[Nn]on-numeric argument to mathematical function",
    conditionMessage(cnd)
  )
  if (!compiled || is.null(fname) || fname %in% pt_math_group_members) {
    stop(cnd)
  }
  replacement <- pt_exact_replacements[[fname]]
  hints <- c(
    "i" = "It hands its argument to compiled code without dispatching, so it
           never sees the tangent that carries the derivative."
  )
  if (is.null(replacement)) {
    hints <- c(
      hints,
      "i" = "deli has no drop-in replacement for it. Rebuild the expression from
             arithmetic and the {.code Math} group functions deli
             differentiates."
    )
  } else {
    hints <- c(
      hints,
      "i" = "Use {.fn {replacement}} instead, which returns the same values and
             differentiates exactly."
    )
  }
  if (identical(fname, "psigamma")) {
    # The one replacement whose arguments do not line up, so a positional
    # substitution computes a different quantity and raises no error.
    hints <- c(
      hints,
      "i" = "The two take their arguments in opposite orders:
             {.code deli_polygamma(n, x)} against
             {.code psigamma(x, deriv = n)}."
    )
  }
  cli::cli_abort(
    c(
      "{.fn {fname}} cannot be differentiated when {.arg deriv_method} is
       {.val exact}.",
      hints,
      pt_capprox_hint()
    ),
    parent = cnd,
    class = "deli_exact_unsupported_function",
    call = NULL
  )
}

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

# ---- c() for combining tangent-carrying values -------------------------------
# `c()` is a base internal generic, so a registered S3 method resolves from any
# environment. That makes `c()` the one portable way to assemble a
# tangent-carrying vector: the masked reshaping helpers (`matrix`, `rbind`,
# `cbind`) are in scope only inside the package namespace, so a transform
# defined in a user's global environment reaches the base versions and loses its
# derivatives.
#
# Every operand is normalized to parallel primal/tangent vectors and the two
# slots are concatenated independently, so the result is a PrimalTangentArray
# that dispatches to the Ops, Math, and Summary methods. Returning a bare list
# instead left the result outside every group generic, so arithmetic and
# reductions over it failed.
#
# `c()` dispatches on its first argument only, so a call that leads with a plain
# numeric reaches the internal default and unpacks each pair into its two slots.
# No mask could fix that portably, since a mask applies only to code evaluated
# inside the namespace. The resulting bare list carries no tangents, and
# `extract_tangent_column()` aborts on it rather than returning a silent
# all-zero Jacobian.

# Concatenate operands by flattening each to parallel primal/tangent vectors and
# joining the two slots separately. A scalar-broadcast tangent is recycled to its
# primal's length first, mirroring `pt_bind()`, so a scalar pair holding a vector
# primal contributes one tangent per element. A numeric constant contributes a
# zero tangent per element.
#
# `recursive` and `use.names` are declared, as they are on every base `c` method,
# so that a caller writing `c(a, b, use.names = FALSE)` matches them to these
# formals. Without them both fall into `...` and are concatenated as ordinary
# operands, lengthening the result and adding a spurious row to the Jacobian.
#
# `use.names` applies to the primal only. The tangent is a parallel derivative
# array read positionally by `extract_tangent_column()`, so names on it would
# serve nothing and would leak into the Jacobian's dimnames. `recursive` has no
# analogue here and is accepted and ignored: the one list operand that
# concatenates is a list of scalar pairs, which `pt_flatten()` reduces to parallel
# numeric vectors before either slot is unlisted, so no nesting survives for
# `recursive` to act on and both settings give the same value.
#
# Every other list operand aborts. Both slots of a tangent-carrying value are
# plain numeric objects of identical shape, and a list breaks that contract
# silently in either direction: one that carries derivatives has no reduction
# down to a pair of parallel vectors, and one that carries none would be paired
# with the wrong number of zeros. Neither is a lost tangent, so the abort names
# `do.call(c, ...)`, which reaches this method with each operand in its own right.
#
# The primal's names diverge from `base::c()` when an operand carries names of
# its own: `c(a = pair, b = c(u = 1, v = 2))` yields `"a" "b1" "b2"` where base
# gives `"a" "b.u" "b.v"`, because `pt_flatten()` and the `as.vector()` fallback
# strip an operand's inner names before the two slots are concatenated. The
# divergence is cosmetic and intended. Names on a tangent-carrying value are
# carried for display only, and the Jacobian is read positionally.
#' @noRd
pt_concat <- function(..., recursive = FALSE, use.names = TRUE) {
  slots <- lapply(list(...), function(a) {
    if (is_tangent_container(a)) {
      parts <- pt_flatten(a)
      parts$tangent <- pt_recycle_tangent(parts$primal, parts$tangent)
      return(parts)
    }
    # A list operand holding no derivatives anywhere cannot be paired with a
    # tangent slot: `as.vector()` hands a list back unchanged, so the fallback
    # would pair it with one zero per top-level element while the `unlist()` below
    # flattened its internals into the primal slot. The two slots would then have
    # different lengths, and the Jacobian read off the tangent slot would have
    # fewer rows than the value it claims to differentiate.
    if (is.list(a)) {
      pt_unsupported_shape_abort()
    }
    list(primal = as.vector(a), tangent = numeric(length(a)))
  })
  primal_tangent_array(
    unlist(lapply(slots, `[[`, "primal"), use.names = use.names),
    unlist(lapply(slots, `[[`, "tangent"), use.names = FALSE)
  )
}

#' @export
c.PrimalTangent <- pt_concat

#' @export
c.PrimalTangentArray <- pt_concat

#' @export
c.PrimalTangentVector <- pt_concat

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
      # An all-NA primal reduced under na.rm leaves nothing to select. Base R
      # returns -Inf with a warning; the constant sentinel carries a zero
      # tangent. Keep base's exact wording so callers matching on it see no
      # change.
      if (length(ap) == 0L) {
        cli::cli_warn(
          "no non-missing arguments to max; returning -Inf",
          call = NULL
        )
        return(primal_tangent(-Inf, 0))
      }
      idx <- which.max(ap)
      primal_tangent(ap[idx], at[idx])
    },
    "min" = {
      # The mirror of max: an empty primal reduces to +Inf with base's warning
      # and a zero tangent.
      if (length(ap) == 0L) {
        cli::cli_warn(
          "no non-missing arguments to min; returning Inf",
          call = NULL
        )
        return(primal_tangent(Inf, 0))
      }
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

# ---- mean, median, and quantile ---------------------------------------------
# These three are ordinary S3 generics rather than members of a group generic,
# so a tangent-carrying argument reaches them only through a registered method.
# Without one, `mean.default()` warns and returns NA, and the sample-quantile
# methods index into the pair's two slots as though they were the data. Each of
# the three tangent surfaces is given the same method object, as the Summary
# methods are.

# `mean()` is linear, so both slots take the same average. A scalar-broadcast
# tangent is recycled to the primal's length first, mirroring `pt_concat()`, so
# a scalar pair holding a vector primal averages one tangent per element rather
# than dividing a single derivative by the number of elements.
#
# `na.rm = TRUE` decides which elements to drop from the primal and applies the
# same mask to the tangent, keeping the two slots aligned. This follows
# `pt_summary()`, where dropping from each slot independently would keep a
# tangent whose primal had been removed.
#
# A nonzero `trim` is not linear: it sorts and discards a fraction of the order
# statistics, which puts it in the same class as `median()` and `quantile()`
# below, so it aborts rather than quietly differentiating the untrimmed mean.
#' @noRd
pt_mean <- function(x, trim = 0, na.rm = FALSE, ...) {
  if (!isTRUE(trim == 0)) {
    cli::cli_abort(
      c(
        "{.fn mean} with a nonzero {.arg trim} cannot be differentiated when
         {.arg deriv_method} is {.val exact}.",
        "i" = "An untrimmed {.fn mean} is linear and differentiates exactly.",
        pt_order_statistic_hint(),
        pt_capprox_hint()
      ),
      class = "deli_exact_unsupported_function",
      call = NULL
    )
  }
  parts <- pt_flatten(x)
  primal <- parts$primal
  tangent <- pt_recycle_tangent(primal, parts$tangent)
  if (na.rm) {
    keep <- !is.na(primal)
    primal <- primal[keep]
    tangent <- tangent[keep]
  }
  # The primal goes through `mean()` so that it matches base R value for value:
  # `mean()` makes a second pass to correct the rounding error `sum(x) / n`
  # leaves, which the two disagree on around the tenth significant digit for a
  # cancellation-heavy vector. The tangent divides by the primal's length
  # because the recycled tangent tracks the primal element for element.
  primal_tangent(mean(primal), sum(tangent) / length(primal))
}

#' @export
mean.PrimalTangent <- pt_mean

#' @export
mean.PrimalTangentArray <- pt_mean

#' @export
mean.PrimalTangentVector <- pt_mean

# `median()` and `quantile()` select among the order statistics of their
# argument, which `pt_order_statistic_abort()` explains deli declines to
# differentiate. Registering methods that abort replaces two outcomes that were
# worse than an error: `median()` returned an empty container, so the Jacobian
# came back with no rows at all and the delta method with a 0-by-0 matrix, and
# `quantile()` read the pair's two slots as though they were the data,
# returning the primal of one element beside the tangent of another.
#' @noRd
pt_median <- function(x, na.rm = FALSE, ...) pt_order_statistic_abort("median")

#' @noRd
pt_quantile <- function(x, ...) pt_order_statistic_abort("quantile")

#' @export
median.PrimalTangent <- pt_median

#' @export
median.PrimalTangentArray <- pt_median

#' @export
median.PrimalTangentVector <- pt_median

#' @export
quantile.PrimalTangent <- pt_quantile

#' @export
quantile.PrimalTangentArray <- pt_quantile

#' @export
quantile.PrimalTangentVector <- pt_quantile

# ---- Boolean coercion (for if/else) -----------------------------------------

# Coercion to logical stays, where coercion to double aborts. A live tangent does
# reach this method, and dropping it is the right answer: a logical coercion is a
# step function, zero everywhere except at the jump, so the derivative it carries
# is zero almost everywhere. That is the same ground `ee_percentile()` stands on,
# where an indicator built from comparisons has an identically zero bread. A
# double is different, because it is the type the tangent flows through, so
# returning one would truncate a derivative that is still doing work.
#' @export
as.logical.PrimalTangent <- function(x, ...) as.logical(x$primal)

#' @export
as.logical.PrimalTangentArray <- function(x, ...) as.logical(x$primal)

# Each surface has to state the rule for itself, because base R reads the
# classed list rather than the payload: an array answered with one value per
# slot for a length-1 payload and refused every longer one, and a parameter
# vector refused them all. A parameter vector keeps its primals one per element
# rather than in a slot of its own, so the flattener assembles them before the
# same coercion runs.
#' @export
as.logical.PrimalTangentVector <- function(x, ...) {
  as.logical(pt_flatten(x)$primal)
}

# Every tangent-carrying value aborts, a genuinely scalar pair included. The
# scalar payload was the one shape a coercion could satisfy, since its primal is
# already a plain double, and satisfying it dropped a live derivative with no
# NA and no warning for any later rule to catch.
#' @export
as.double.PrimalTangent <- function(x, ...) {
  pt_coercion_abort("numeric", c("as.numeric", "as.double"))
}

# R dispatches `as.numeric()` through the `as.double()` method table, so this
# method is never reached and the rule above governs both spellings. Changing
# the behavior here alone would have no effect.
#' @export
as.numeric.PrimalTangent <- function(x, ...) as.numeric(x$primal)

# No tangent-carrying value has a representation in a plain numeric type that
# keeps its derivative, so every coercion to one aborts rather than losing the
# tangents. `target` names the type the coercion was asked for and `fns` the
# spellings that ask for it, so each method reports the rule in the caller's own
# terms while the guidance out of it is written once.
#
# Base R would otherwise report `'list' object cannot be coerced to type
# 'double'` for an array and return `NA` with a coercion warning for a parameter
# vector, neither of which points at the cause or at the remedy.
#
# `storage.mode<-` is the one spelling of a coercion no method here reaches, and
# it is a known hole rather than an accepted behavior. It is a primitive that
# coerces the object it is handed without consulting the method table, and it
# does not delegate to `mode<-`, which is an ordinary closure that looks up
# `as.<value>` and so does reach these methods. A scalar pair assigned an
# integer storage mode therefore still comes back holding its two slots as data.
#' @noRd
pt_coercion_abort <- function(target, fns) {
  cli::cli_abort(
    c(
      "A tangent-carrying value cannot be coerced to a plain {.cls {target}}.",
      "i" = "A plain {.cls {target}} has nowhere to carry a derivative, so
             {.fn {fns}} cannot preserve one.",
      "i" = "A tangent-carrying value needs no coercion: the arithmetic
             operators and the {.code Math} group functions take it directly.",
      "i" = "Use {.fn c} to flatten a tangent-carrying vector or matrix such as
             {.code X %*% theta} down to a vector; it keeps the derivative.",
      pt_capprox_hint()
    ),
    class = "deli_exact_tangent_lost",
    call = NULL
  )
}

#' @export
as.double.PrimalTangentArray <- function(x, ...) {
  pt_coercion_abort("numeric", c("as.numeric", "as.double"))
}

#' @export
as.double.PrimalTangentVector <- function(x, ...) {
  pt_coercion_abort("numeric", c("as.numeric", "as.double"))
}

# ---- coercion to a matrix, an integer, or a vector of a numeric mode --------

# A matrix of doubles is a type the tangent flows through, so the rule that
# governs `as.double()` governs this coercion too. What base R does instead is
# decided by shapes that have nothing to do with the payload:
# `as.matrix.default()` reads `length()` and `dim()` through the methods above,
# which answer for the payload, while `is.matrix()` reads the real attribute and
# sees a classed list of two slots. The two disagree, so the coercion reports
# `dims [product n] do not match the length of object [2]` from a base frame no
# deli rule runs in. A two-element payload is worse: it is the one length the
# default accepts, so it sets `dim = c(2, 1)` on the container itself and hands
# back its two slots as the cells of a matrix, with the tangents gone and
# nothing downstream to report it.
#' @noRd
pt_as_matrix <- function(x, ...) pt_coercion_abort("matrix", "as.matrix")

#' @export
as.matrix.PrimalTangent <- pt_as_matrix

#' @export
as.matrix.PrimalTangentArray <- pt_as_matrix

#' @export
as.matrix.PrimalTangentVector <- pt_as_matrix

# An integer is a plain number, so a tangent has nowhere to go in one either.
# Base R coerced the container itself for a scalar pair, so
# `as.integer(primal_tangent(2.7, 1))` returned the truncated primal with the
# tangent appended to it as a second value, which nothing downstream can tell
# from data. The other payloads stopped already, with `'list' object cannot be
# coerced to type 'integer'`, which points at neither the cause nor the remedy.
#' @noRd
pt_as_integer <- function(x, ...) pt_coercion_abort("integer", "as.integer")

#' @export
as.integer.PrimalTangent <- pt_as_integer

#' @export
as.integer.PrimalTangentArray <- pt_as_integer

#' @export
as.integer.PrimalTangentVector <- pt_as_integer

# A character and a complex are plain values in the same sense: neither has a
# place to keep a derivative, so a tangent that reaches one is gone. Base R read
# the container rather than the payload here too, and did so without stopping.
# `as.character(primal_tangent(2.7, 1))` returned the primal with the tangent
# appended to it as a second string, and on a wider payload it deparsed each
# slot into a string of its own, so `"c(1, 2)"` stood where a value belonged.
# `as.complex()` returned the two slots as two complex numbers.
#' @noRd
pt_as_character <- function(x, ...) {
  pt_coercion_abort("character", "as.character")
}

#' @export
as.character.PrimalTangent <- pt_as_character

#' @export
as.character.PrimalTangentArray <- pt_as_character

#' @export
as.character.PrimalTangentVector <- pt_as_character

#' @noRd
pt_as_complex <- function(x, ...) pt_coercion_abort("complex", "as.complex")

#' @export
as.complex.PrimalTangent <- pt_as_complex

#' @export
as.complex.PrimalTangentArray <- pt_as_complex

#' @export
as.complex.PrimalTangentVector <- pt_as_complex

# `as.vector()` dispatches on the class of its first argument, so its `mode`
# argument is another spelling of the coercions above and carries the same rule:
# every mode that forces an atomic type aborts, and `"double"` aborts here even
# though it reached base R's coercion rather than the `as.double()` method.
#
# The default mode is `"any"`, which returns a list unchanged, so
# `as.vector(X %*% theta)` keeps its tangents and differentiates exactly. That
# is the one route a differentiated function reaches this method by on purpose,
# and it is left alone. `"list"` is left to base R, which coerces nothing: a
# container already is a list of two slots, so that mode is the one faithful
# representation of it outside its own class.
#' @noRd
pt_vector_coercion <- function(x, mode = "any") {
  if (mode %in% c("integer", "double", "numeric", "character", "complex")) {
    pt_coercion_abort(mode, "as.vector")
  }
  if (identical(mode, "any")) {
    return(x)
  }
  as.vector(unclass(x), mode)
}

#' @export
as.vector.PrimalTangent <- pt_vector_coercion

#' @export
as.vector.PrimalTangentArray <- pt_vector_coercion

#' @export
as.vector.PrimalTangentVector <- pt_vector_coercion

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
    # A plain list of scalar pairs, produced by an `lapply()` or `sapply()` over
    # tangent-carrying values inside a differentiated function. (`c()` returns a
    # PrimalTangentArray, so it no longer yields this shape.) pt_flatten already
    # collapses the list into parallel numeric vectors, so a masked binder of
    # such a list keeps its tangents instead of hitting the numeric-constant
    # fallback below and zeroing them.
    return(pt_flatten(x))
  }
  # Numeric constant: zero tangent, matching the constant's shape
  zero <- x
  zero[] <- 0
  list(primal = x, tangent = zero)
}

# Recycle a scalar-broadcast tangent up to the primal's length and give it the
# primal's shape. A scalar pair built from a vector primal and a constant carries
# a length-1 tangent that broadcasts under elementwise arithmetic (R recycles
# it), but flattening, binding, and subsetting reshape the primal and tangent
# slots independently, so the tangent must first be expanded to match the primal
# element for element.
#
# The primal's `dim` decides the shape of both slots whatever their lengths, not
# only when the tangent was just recycled. A 1-by-1 primal against a length-1
# tangent, and any equal-length pairing of a shaped primal with a flat tangent,
# both leave the tangent dimensionless otherwise, and `tangent[i, j]` then
# reports `incorrect number of dimensions` where the primal selects. Assigning a
# `dim` drops the dimnames with it, so a tangent that already carries the
# primal's shape is left as it is, which is what keeps the labels a container
# records on both slots available to a selection by name.
#' @noRd
pt_recycle_tangent <- function(primal, tangent) {
  if (length(tangent) < length(primal)) {
    tangent <- rep_len(tangent, length(primal))
  }
  dims <- dim(primal)
  reshape <- !is.null(dims) &&
    !identical(dim(tangent), dims) &&
    length(tangent) == length(primal)
  if (reshape) {
    dim(tangent) <- dims
  }
  tangent
}

# Flatten a tangent-carrying container into parallel numeric vectors. Used when
# reshaping into a matrix and when concatenating with `c()`. Accepts a
# PrimalTangent (a scalar pair, possibly with a vector primal and a scalar
# broadcast tangent), a PrimalTangentVector, a PrimalTangentArray, or a plain
# list of scalar pairs (produced by an `lapply()` or `sapply()` over
# tangent-carrying values). Any other list shape has no flattening and aborts.
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
  # x is now a plain list, which flattens only when every element contributes
  # exactly one primal and one tangent: a scalar pair, or a length-1 constant
  # whose tangent is zero. A list of tangent arrays, of parameter vectors, or of
  # vector-payload pairs still holds every derivative it was given, so its
  # container shape is the problem and it is reported as such. Reaching the
  # `vapply()` below with one of those would instead raise an unclassed
  # length error, and coercing an element with `as.numeric()` would report a lost
  # tangent, which misdiagnoses an intact one.
  flattens <- vapply(
    x,
    function(e) {
      if (is_pt(e)) {
        length(e$primal) == 1L
      } else {
        length(e) == 1L && (is.numeric(e) || is.logical(e))
      }
    },
    logical(1)
  )
  if (!all(flattens)) {
    pt_unsupported_shape_abort()
  }
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
  # redispatch, each of which would hand back the same object.
  #
  # Dimnames are dropped so the returned score matrix is unnamed regardless of
  # whether the caller passed a plain matrix, a named matrix, or a data frame.
  # The Python reference operates on unnamed arrays. The two axes of the score
  # are dropped for different reasons, and only one of them is now a channel
  # for anything.
  #
  # The score's columns are the observations, and its column names would be the
  # design's row labels, which mean nothing downstream.
  #
  # The score's rows are the parameters, and `estimate()` does read their names
  # when `init` carries none, so a `p`-by-`n` return is how an estimating
  # function names the parameters it defines. What a design matrix supplies is
  # not that. The score is built as `t(X * resid)`, so the row names would be
  # whatever the caller's data frame happened to call its columns, with nothing
  # at all for an intercept, which is exactly the incomplete vector the naming
  # rule discards. An `ee_*` function that wants to name its parameters has to
  # say so deliberately rather than inherit the caller's column headings. Every
  # design argument of every built-in estimating equation runs through here, `S`
  # and `W` and `V` and the counterfactual plans as much as `X`, so that the
  # rule holds for all of them and not only for the ones whose score happens to
  # discard the labels on its way out.
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

# Does `x` belong to one of the three tangent classes?
#' @noRd
is_tangent_value <- function(x) {
  is_pt(x) || is_pt_array(x) || inherits(x, "PrimalTangentVector")
}

# Is `x` something a tangent-aware matrix() should reshape?
#
# Each of the three tangent classes is a container in its own right, and so is a
# list holding any of them: which class an element belongs to is not the author's
# choice, because `theta[k] * x` yields a scalar pair, `c(theta[k] * x)` yields a
# tangent array, and `theta[i:j]` yields a parameter vector. All three carry
# derivatives, so an `lapply()` over equations produces a container whichever of
# them it collects and wherever in the list they fall. An operand this predicate
# does not recognize takes the numeric-constant route instead, which replaces its
# tangents with zeros.
#' @noRd
is_tangent_container <- function(x) {
  is_tangent_value(x) ||
    (is.list(x) &&
      length(x) > 0 &&
      any(vapply(x, is_tangent_value, logical(1))))
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

# ---- row labels -------------------------------------------------------------
# An estimating function labels the rows of its return to name the parameters,
# and it has to be able to do so under both passes. Unlike the binders above,
# this is reached by dispatch rather than by masking. `rownames<-` is an
# ordinary closure in base R and cannot be given a method, but what it delegates
# to, `dimnames<-`, is an internal generic, so a method on it is selected from
# whatever environment the estimating function was written in. That is the whole
# point: users write estimating functions in the global environment, and a mask
# inside the package namespace would serve only the package's own code.
#
# The labels are recorded on both slots. The primal alone would not do: the
# subset methods index the two slots with the same subscript, so a character
# subscript the primal resolves and the tangent does not raises `no 'dimnames'
# attribute for array` from the tangent, at the selection rather than at the
# assignment. A container whose primal carries a shape its tangent broadcasts
# against is expanded to that shape first, so the labels have the same axes to
# sit on in both slots.
#
# Nothing reaches the Jacobian by this route. A column of it is read off the
# tangent slot through `as.vector()`, whether the slot arrives whole or already
# reduced to row sums, and that drops the labels, so a psi that names its rows
# reports the same unnamed bread the numeric pass reports. `estimate()` reads
# parameter names off the plain numeric evaluation it makes at the solved values,
# never off a differentiated one.
#
# Recording the labels also restores base R's validation of them, which returning
# the container unchanged suppressed: a label vector that does not match the axis
# is an error on the numeric pass, and the exact pass stops in the same place
# rather than being the more permissive one.
#
# Only the two containers that answer `dim()` get a setter, for the same reason.
# Base R's `rownames<-` consults `dim()` before it reaches `dimnames<-` and stops
# on anything with fewer than one dimension, so a PrimalTangentVector, and a
# PrimalTangent whose primal is a plain vector, both raise the error that a plain
# numeric vector raises on the numeric pass. Assigning `NULL` to an unlabeled
# container returns earlier still and reaches neither setter, which is why the
# `rownames(out) <- NULL` lines in `R/ee-glm.R` survive the exact pass.
#
# The reader is what lets a second assignment keep the first: base R's
# `colnames<-` builds its replacement list from `dimnames(x)`, so labels that do
# not read back are erased by the next axis to be named.

#' @noRd
pt_dimnames <- function(x) dimnames(x$primal)

#' @export
dimnames.PrimalTangent <- pt_dimnames

#' @export
dimnames.PrimalTangentArray <- pt_dimnames

#' @noRd
pt_set_dimnames <- function(x, value) {
  x$tangent <- pt_recycle_tangent(x$primal, x$tangent)
  dimnames(x$primal) <- value
  dimnames(x$tangent) <- value
  x
}

#' @export
`dimnames<-.PrimalTangent` <- pt_set_dimnames

#' @export
`dimnames<-.PrimalTangentArray` <- pt_set_dimnames

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
#'   `cos`, etc.). The matrix product `%*%`, the transpose `t()`, the
#'   concatenation `c()`, `mean()`, and two-dimensional indexing differentiate
#'   exactly in any function, because each is a registered S3 method. The
#'   reshaping, binding, and reduction operations `matrix()`, `rbind()`,
#'   `cbind()`, `rep()`, and `colSums()` differentiate exactly for functions
#'   evaluated within the package (such as the built-in estimating equations),
#'   because those masked forms are only in
#'   scope there; a user-defined function that calls them from the global
#'   environment reaches the base versions, and the differentiation aborts
#'   rather than returning a silent approximation. Coercion cannot preserve a
#'   derivative: `as.numeric()`, `as.double()`, `as.integer()`, `as.character()`,
#'   `as.complex()`, `as.matrix()`,
#'   and `as.vector()` asked for any atomic `mode` each return plain values,
#'   which have nowhere to carry one, so each aborts on a tangent-carrying value
#'   rather than dropping the tangents silently. `as.vector()` with its default `mode = "any"` returns its
#'   argument unchanged and keeps the tangents. Use `c()` to flatten a
#'   matrix-shaped result such as `X %*% theta` to a vector instead.
#'   `as.logical()` is the exception, because a logical coercion is a step
#'   function whose derivative is zero almost everywhere, so it returns the
#'   logical its payload coerces to. `log(x, base)` differentiates
#'   with respect to `x` only; the `base` argument is treated as a constant, and
#'   a `base` that itself carries a tangent (a value derived from `theta`) is not
#'   supported and aborts rather than dropping the base's contribution.
#'   `median()`, `quantile()`, and `mean(x, trim)` select among the order
#'   statistics of their argument, whose derivative is that of whichever order
#'   statistic the current values select rather than that of the population
#'   quantity, so each aborts as well.
#'
#' @returns A matrix where element `[i, j]` is the partial derivative of
#'   the `i`-th output with respect to the `j`-th parameter.
#'
#' @section Conditions:
#' The aborts raised under exact differentiation carry condition classes, which
#' are the first in the package: `deli_exact_tangent_lost` when derivative
#' information is gone or a coercion would discard it,
#' `deli_exact_unsupported_function` when a function reached with a
#' tangent-carrying argument has no rule under exact differentiation, and
#' `deli_exact_unsupported_shape` when a result keeps its tangents but arrives
#' in a container that summing the estimating equations has no reduction for.
#'
#' @keywords internal
auto_differentiation <- function(theta, f) {
  theta <- as.numeric(theta)
  p <- length(theta)

  jacobian_cols <- vector("list", p)

  # One handler around the whole sweep rather than one per direction. The
  # rewrite depends only on which function failed, not on which direction was
  # being evaluated, and installing the handler once leaves the per-direction
  # cost unchanged. `tryCatch()` evaluates its expression in this frame, so the
  # loop still fills `jacobian_cols` here.
  tryCatch(
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
    },
    error = pt_rethrow_exact_error
  )

  do.call(cbind, jacobian_cols)
}

# Read the tangent components of an evaluated function result into one Jacobian
# column, covering every shape a differentiated function can return.
#' @noRd
extract_tangent_column <- function(result) {
  if (is_pt(result)) {
    # Single output (scalar, or a vector or matrix primal, whose tangent is either
    # shaped to match or broadcast from a single derivative). Flattened
    # column-major for the same two reasons as the tangent array below: a shaped
    # tangent handed to `cbind()` unflattened is read as several columns rather
    # than one, and any labels the function recorded on the slot would ride into
    # the Jacobian's dimnames.
    return(as.vector(result$tangent))
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
    # A list result holding at least one tangent-carrying pair comes from an
    # `lapply()` or `sapply()` over per-equation values. A list (or list matrix)
    # with no PrimalTangent elements instead signals that a base reshaping,
    # binding, or selection helper such as rbind(), cbind(), or ifelse() was
    # reached from outside the package namespace and stripped the tangents.
    # Aborting keeps the documented safety property rather than returning a
    # silent all-zero column.
    if (!any(vapply(result, is_pt, logical(1)))) {
      pt_tangent_lost_abort()
    }
    # One output per list element
    return(vapply(
      result,
      function(e) if (is_pt(e)) e$tangent else 0,
      numeric(1)
    ))
  }
  # Plain numeric: either a genuinely constant output, whose derivatives are all
  # zero, or a result whose tangents were dropped before it was assembled. An NA
  # separates the two. A function that does not recognize a tangent-carrying
  # argument returns NA for it (base R's `mean.default()` was the canonical
  # case), and that NA flows through the plain numeric arithmetic that follows.
  # A tangent-free result with no NA is a genuine constant, which both
  # `function(x) 5` and `ee_percentile()`'s indicator score rely on: comparisons
  # return primal logicals by design, so that score carries no tangent and its
  # bread is identically zero under every derivative method.
  #
  # `NaN` is excluded rather than folded in with `anyNA()`, which counts it as
  # missing. A `NaN` is what arithmetic on numbers produces, `0 / 0` and
  # `log(-1)` among them, so it is evidence of a numerical problem in the
  # differentiated function rather than of a dropped tangent. Reporting a lost
  # derivative for one would be a false diagnosis, and a constant carrying a
  # `NaN` is still a constant.
  if (any(is.na(result) & !is.nan(result))) {
    pt_na_tangent_abort()
  }
  rep(0, length(result))
}
