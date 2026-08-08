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

# Every node of a parsed Rd file carries its tag in the `Rd_tag` attribute,
# except the plain text leaves, which carry none. The text a node stands for is
# its leaves pasted back together.
rd_tag <- function(node) {
  tag <- attr(node, "Rd_tag")
  if (is.null(tag)) "" else tag
}

rd_text <- function(node) paste(unlist(node), collapse = "")

# The names a NAMESPACE makes available to a user: everything it exports, plus
# `generic.class` for every `S3method(generic, class)` registration. Either
# part of a registration arrives quoted when the name is not syntactic, which
# is the case for `%*%` and for the package-qualified class string of an S7
# object, so both are unquoted before use. A class string is a single name and
# so never contains a comma, which makes splitting on commas safe.
namespace_provides <- function() {
  namespace <- readLines(file.path(find.package("deli"), "NAMESPACE"))
  unquote <- function(x) gsub("\"", "", x)

  exports <- grep("^export\\(", namespace, value = TRUE)
  exports <- unquote(sub("^export\\((.*)\\)$", "\\1", exports))

  registrations <- grep("^S3method\\(", namespace, value = TRUE)
  registrations <- sub("^S3method\\((.*)\\)$", "\\1", registrations)
  methods <- vapply(
    strsplit(registrations, ","),
    function(parts) paste(unquote(parts[[1]]), unquote(parts[[2]]), sep = "."),
    character(1)
  )

  c(exports, methods)
}

# The names one Rd page claims to document: its aliases other than the topic
# name, the `generic.class` function behind each `\method{generic}{class}`
# usage entry, and the function named at the start of every other usage line.
# Rd backticks a class string that is not syntactic and the registration does
# not, so the backticks come off. The open parenthesis is part of the usage
# pattern because a data set page gives its data set as a bare name, which
# documents no function.
documented_names <- function(rd) {
  tags <- vapply(rd, rd_tag, character(1))
  topic <- rd_text(rd[tags == "\\name"][[1]])
  aliases <- vapply(rd[tags == "\\alias"], rd_text, character(1))

  usage <- if (any(tags == "\\usage")) rd[tags == "\\usage"][[1]] else list()
  is_method <- vapply(usage, function(x) rd_tag(x) == "\\method", logical(1))
  methods <- vapply(
    usage[is_method],
    function(x) paste0(rd_text(x[[1]]), ".", gsub("`", "", rd_text(x[[2]]))),
    character(1)
  )

  plain <- rd_text(usage[!is_method])
  calls <- regmatches(
    plain,
    gregexpr("(?m)^[[:alnum:]._]+(?=\\()", plain, perl = TRUE)
  )[[1]]

  unique(c(setdiff(aliases, topic), methods, calls))
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

test_that("no vignette is listed twice in the pkgdown article index", {
  skip_on_cran()
  skip_if_not_installed("yaml")
  config <- deli_pkgdown_config()
  skip_if(is.null(config), "_pkgdown.yml is not in the built package")

  # The article index groups vignettes the way the reference index groups
  # topics, and a vignette listed under two of those groups is the same claim
  # about two families that the reference check above refuses.
  articles <- yaml::read_yaml(config)$articles
  vignettes <- unlist(
    lapply(articles, function(group) group$contents),
    use.names = FALSE
  )
  # Fails if the groups stop yielding vignettes at all, which would leave the
  # assertion below trivially true.
  expect_gt(length(vignettes), 0)

  expect_equal(unique(vignettes[duplicated(vignettes)]), character())
})

# Documented names against the NAMESPACE ------------------------------------
#
# CRAN's incoming review reads a page that carries examples as a promise that
# a user can call what the page documents, and it looks in the NAMESPACE for
# the proof. A method registered only when the package loads leaves no trace
# there, so a page documenting one is reported as examples for an unexported
# function however correct its content is.
#
# The assertion is on the whole package rather than on any one page: for every
# page with an examples section, every name the page documents must be either
# an `export()` or an `S3method()` registration.

test_that("every page with examples documents names the NAMESPACE provides", {
  db <- deli_rd_db()
  expect_gt(length(db), 0)

  provided <- namespace_provides()
  # Fails if the NAMESPACE stops parsing at all, which would report every
  # documented name as missing rather than none.
  expect_gt(length(provided), 0)

  with_examples <- Filter(
    function(rd) any(vapply(rd, rd_tag, character(1)) == "\\examples"),
    db
  )
  # Fails if the walk stops finding pages with examples at all, which would
  # leave the assertion below trivially true.
  expect_gt(length(with_examples), 0)

  # `sprintf()` rather than `paste0()` because a page with nothing missing has
  # to contribute nothing: `paste0()` would recycle the empty second argument
  # to "" and report every page as an offender.
  undocumented <- unlist(
    Map(
      function(page, rd) {
        sprintf("%s: %s", page, setdiff(documented_names(rd), provided))
      },
      names(with_examples),
      with_examples
    ),
    use.names = FALSE
  )

  expect_equal(undocumented, character())
})
