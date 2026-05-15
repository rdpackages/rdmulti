candidate_roots <- unique(c(
  Sys.getenv("RDMULTI_REPO_ROOT", unset = NA_character_),
  normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE),
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
))
candidate_roots <- candidate_roots[!is.na(candidate_roots)]

repo_root <- NA_character_
for (root in candidate_roots) {
  if (file.exists(file.path(root, "R", "simdata_multic.csv"))) {
    repo_root <- root
    break
  }
}

if (is.na(repo_root)) {
  message("Skipping rdmulti numerical contracts: repository illustration data not found.")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(rdmulti))

assert_allclose <- function(actual, expected, tolerance = 1e-8) {
  comparison <- all.equal(
    unname(actual),
    expected,
    tolerance = tolerance,
    check.attributes = FALSE
  )
  if (!isTRUE(comparison)) {
    stop(comparison, call. = FALSE)
  }
}

data <- read.csv(file.path(repo_root, "R", "simdata_multic.csv"))
invisible(capture.output(result <- rdmc(data$y, data$x, data$c)))

assert_allclose(
  result$Coefs,
  matrix(c(484.83090323, 297.98107945, 398.91490949, 436.40049184), nrow = 1)
)
assert_allclose(
  result$Pv,
  matrix(c(7.84407236e-48, 8.54527540e-16, 0, 7.41987105e-04), nrow = 1)
)
assert_allclose(
  result$H,
  matrix(
    c(
      14.6619828, 14.6619828,
      11.9522678, 11.9522678,
      NA_real_, NA_real_,
      13.6840172, 13.6840172
    ),
    nrow = 2
  )
)
assert_allclose(
  result$Nh,
  matrix(c(149, 140, 120, 126, 269, 266, 273, 277), nrow = 2)
)

data_cumul <- read.csv(file.path(repo_root, "R", "simdata_cumul.csv"))
cvec <- c(data_cumul$c[1], data_cumul$c[2])
invisible(capture.output(result_cumul <- rdms(data_cumul$y, data_cumul$x, cvec)))

assert_allclose(
  result_cumul$Coefs,
  matrix(c(395.4918261, 342.87249624, NA_real_), nrow = 1)
)
assert_allclose(
  result_cumul$Pv,
  matrix(c(1.59995659e-145, 3.42385007e-120, NA_real_), nrow = 1)
)
assert_allclose(
  result_cumul$H,
  matrix(c(15.10927597, 15.10927597, 12.22129226, 12.22129226, NA_real_, NA_real_), nrow = 2)
)
assert_allclose(
  result_cumul$Nh,
  matrix(c(140, 146, 142, 123, NA_real_, NA_real_), nrow = 2)
)
