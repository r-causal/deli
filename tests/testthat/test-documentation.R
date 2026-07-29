# A cross-package documentation link to a package that is not in DESCRIPTION
# produces a dangling reference: pkgdown cannot resolve it and `?topic` sends
# the reader to a package they may not have. R CMD check can catch this, but
# only when `_R_CHECK_XREFS_PKGS_ARE_DECLARED_` is set, and that variable
# defaults to false, so the check lives here instead.
#
# The assertion is on the whole package rather than on any one link. Every
# `\link[pkg]{...}` and `\link[pkg:topic]{...}` across all of `man/` must name
# either a declared dependency or a package that ships with R.

# ---- Helpers ----

# The parsed Rd node for `\link[pkg:topic]{text}` carries "pkg:topic" as its
# Rd_option, and `\link` is the only tag that uses that attribute. A leading
# "=" marks `\link[=topic]{text}`, which resolves within this package and names
# no external one.
rd_link_options <- function(node) {
  option <- attr(node, "Rd_option")
  here <- if (identical(attr(node, "Rd_tag"), "\\link") && !is.null(option)) {
    as.character(option)
  } else {
    character()
  }
  children <- if (is.list(node)) {
    unlist(lapply(node, rd_link_options), use.names = FALSE)
  } else {
    character()
  }
  c(here, children)
}

# Two branches because the two environments this runs in resolve the package
# differently. Under R CMD check the tests run against an installed copy, which
# keeps no `man/` but does keep the same Rd in `help/deli.rdb`, and that is what
# `tools::Rd_db(package)` reads. Under `devtools::test()` pkgload points
# `find.package()` at the source tree, which has `man/` but no `help/`, so the
# first form reports the installed help as corrupt and the sources are parsed
# instead.
deli_rd_db <- function() {
  tryCatch(
    tools::Rd_db("deli"),
    error = function(cnd) tools::Rd_db(dir = find.package("deli"))
  )
}

declared_dependencies <- function() {
  fields <- packageDescription("deli")[c("Depends", "Imports", "Suggests")]
  declared <- unlist(strsplit(unlist(fields, use.names = FALSE), ","))
  declared <- trimws(sub("\\(.*", "", declared))
  declared[nzchar(declared) & declared != "R"]
}

# The same two environments, and the same reason for two branches, as
# `deli_rd_db()` above. `_pkgdown.yml` is `.Rbuildignore`d, so it reaches
# neither the tarball nor the installed package: under `devtools::test()` it
# sits two directories above the test files in the source tree, and under
# R CMD check there is no copy of it anywhere. `NULL` says which of the two
# this is, and the test that reads it skips rather than failing on the second.
deli_pkgdown_config <- function() {
  path <- test_path("..", "..", "_pkgdown.yml")
  if (file.exists(path)) path else NULL
}

# Rd cross-references ------------------------------------------------------

test_that("every cross-package Rd link names a declared dependency", {
  db <- deli_rd_db()
  expect_gt(length(db), 0)

  found <- unlist(lapply(db, rd_link_options), use.names = FALSE)
  linked <- sort(unique(sub(":.*$", "", found[!startsWith(found, "=")])))
  # Fails if the walk stops finding links at all, which would leave the
  # assertion below trivially true.
  expect_gt(length(linked), 0)

  allowed <- c(
    declared_dependencies(),
    rownames(installed.packages(priority = "base"))
  )
  expect_equal(setdiff(linked, allowed), character())
})

# The pkgdown reference index ----------------------------------------------
#
# A topic listed twice in `_pkgdown.yml` renders twice in the reference index,
# under both groups, which says the same function belongs to two families.
# `pkgdown::check_pkgdown()` does not report it: with a duplicate injected it
# still prints "No problems found", because it checks that every topic is
# listed and that every listed topic exists, not how many times each appears.
# A duplicate `gmm_estimate` entry survived on that basis, so the assertion
# lives here.

test_that("no topic is listed twice in the pkgdown reference index", {
  skip_on_cran()
  skip_if_not_installed("yaml")
  config <- deli_pkgdown_config()
  skip_if(is.null(config), "_pkgdown.yml is not in the built package")

  reference <- yaml::read_yaml(config)$reference
  topics <- unlist(
    lapply(reference, function(group) group$contents),
    use.names = FALSE
  )
  # Fails if the groups stop yielding topics at all, which would leave the
  # assertion below trivially true.
  expect_gt(length(topics), 0)

  expect_equal(unique(topics[duplicated(topics)]), character())
})
