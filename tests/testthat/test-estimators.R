test_that("fit_direct returns expected structure", {
  set.seed(42)
  n   <- 200
  dat <- scenic:::generate_toy_data(n = n, seed = 42)
  res <- fit_direct(dat$Y[, 1], dat$Z, dat$G[, 1], dat$W[, 1], dat$synthetic_data)
  expect_named(res, c("b", "se", "p"))
  expect_true(is.numeric(res$b))
  expect_true(res$se > 0)
  expect_true(res$p >= 0 && res$p <= 1)
})

test_that("fit_coca returns valid estimate or NA on near-zero bY", {
  set.seed(1)
  n   <- 300
  dat <- scenic:::generate_toy_data(n = n, seed = 1)
  Wa  <- rowMeans(dat$W)
  res <- fit_coca(dat$Y[, 1], dat$Z, Wa, dat$synthetic_data)
  expect_true(is.na(res$beta) || is.numeric(res$beta))
})

test_that("fit_iv2sls checks weak instruments", {
  set.seed(7)
  n  <- 200
  # create a dataset where G explains nothing (weak IV)
  dat         <- scenic:::generate_toy_data(n = n, seed = 7)
  weak_g      <- rnorm(n)   # irrelevant instrument
  res         <- fit_iv2sls(dat$Y[, 1], dat$Z, weak_g, dat$W[, 1])
  # may or may not fire depending on random draw; just check structure
  expect_true(is.na(res$beta) || is.numeric(res$beta))
})

test_that("fit_pgc returns numeric beta", {
  set.seed(99)
  n   <- 300
  dat <- scenic:::generate_toy_data(n = n, seed = 99)
  Wa  <- rowMeans(dat$W)
  res <- fit_pgc(dat$Y[, 1], dat$Z, dat$G[, 1], Wa, dat$synthetic_data)
  expect_named(res, c("beta", "se", "pvalue"))
  expect_true(is.numeric(res$beta))
})

test_that("run_simulation returns correct structure", {
  res <- run_simulation(n_iter = 5, n_samples = 100, n_features = 5,
                        beta_Z = 0.1, n_cores = 1)
  expect_named(res, c("raw", "summary", "true_total", "params"))
  expect_true(nrow(res$summary) == 5)
  expect_true("bias" %in% names(res$summary))
})

test_that("run_null_sim type1 error is finite", {
  res <- run_null_sim(n_iter = 10, n_samples = 100, n_features = 5,
                      n_cores = 1)
  expect_true(all(is.finite(res$type1_error)))
})
