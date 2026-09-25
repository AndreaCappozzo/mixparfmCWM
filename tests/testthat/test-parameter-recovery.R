# Parameter recovery on well-separated simulated data.

test_that("well-separated components are recovered", {
  sim <- simulate_mixparfm_data(seed = 31, sep = 3.5)
  set.seed(31)
  fit <- suppressWarnings(fit_default(
    data = sim$data,
    class_init = sim$true_class,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(
      itermax = 50,
      tol = 1e-5,
      n_start = 1
    )
  ))

  ari <- mclust::adjustedRandIndex(fit$class, sim$true_class)
  # The classification-likelihood criterion can misassign genuinely
  # ambiguous patients (covariates compatible with either component), so
  # we require a high but not perfect adjusted Rand index.
  expect_gte(ari, 0.75)

  # mixing weights close to 0.5
  expect_lt(max(abs(fit$parameters$tau - 0.5)), 0.2)

  # regression coefficients: truth is 0.3 and 0.5 (up to label switching)
  beta <- sort(as.numeric(fit$parameters$AFT_parameters["beta.x", ]))
  expect_lt(abs(beta[1] - 0.3), 0.5)
  expect_lt(abs(beta[2] - 0.5), 0.5)

  # frailty variances are positive (truth 0.5)
  theta <- as.numeric(fit$parameters$AFT_parameters["theta", ])
  expect_true(all(theta > 0))

  # frailty predictions: one row per group, columns per component
  expect_equal(nrow(fit$frailty_effect), nlevels(sim$data$family))
  expect_true(all(c("1", "2") %in% colnames(fit$frailty_effect)))
})
