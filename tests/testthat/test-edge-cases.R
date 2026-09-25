# Edge cases surfaced by the algorithm audit: empty components, singleton
# groups, zero-event components, a single Gaussian covariate, and the
# legacy control arguments.

test_that("empty components in class_init give a clear error", {
  sim <- simulate_mixparfm_data(seed = 51)
  bad_init <- sim$true_class
  bad_init[bad_init == 2] <- 1
  expect_error(
    fit_default(
      data = sim$data,
      class_init = bad_init,
      baseline = "weibull",
      frailty = "gamma"
    ),
    "class_init leaves component"
  )
})

test_that("a degenerate start does not abort a multi-start call", {
  sim <- simulate_mixparfm_data(seed = 52)
  d <- sim$data
  # a component whose Gaussian covariates are all identical gives a
  # degenerate covariance; it has enough events to pass validation, so
  # the start must fail at initialization but the random starts must
  # still run
  bad_init <- rep(1, nrow(d))
  idx2 <- which(d$status == 1)[1:4]
  bad_init[idx2] <- 2
  d$Z1[idx2] <- 0
  d$Z2[idx2] <- 0
  set.seed(52)
  fit <- fit_default(
    data = d,
    class_init = bad_init,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(
      itermax = 5,
      tol = 1e-4,
      n_start = 2
    )
  )
  expect_true(any(fit$start_summary$error != ""))
  expect_true(is.finite(fit$classification_loglik))
})

test_that("singleton groups are handled", {
  sim <- simulate_mixparfm_data(seed = 53, n_per_group = 4)
  d <- sim$data
  # turn one group into a single-patient group
  keep <- !(d$family == levels(d$family)[1] & seq_len(nrow(d)) > 1)
  d <- droplevels(d[keep, ])
  sim_class <- sim$true_class[keep]
  set.seed(53)
  fit <- fit_default(
    data = d,
    class_init = sim_class,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(itermax = 5, tol = 1e-4)
  )
  expect_true(is.finite(fit$classification_loglik))
  expect_equal(nrow(fit$frailty_effect), nlevels(d$family))
})

test_that("a zero-event component in class_init gives a clear error", {
  sim <- simulate_mixparfm_data(seed = 54)
  d <- sim$data
  # make all patients of the second component censored: zero-event
  # components give degenerate survival-parameter estimates, so the
  # event floor must reject them before any fitting
  d$status[sim$true_class == 2] <- 0L
  expect_error(
    fit_default(
      data = d,
      class_init = sim$true_class,
      baseline = "weibull",
      frailty = "gamma",
      control_EM_algorithm = mixparfmCWM::control_EM(itermax = 5, tol = 1e-4)
    ),
    "event"
  )
})

test_that("degeneracy flags are reported in start_summary", {
  sim <- simulate_mixparfm_data(seed = 58)
  set.seed(58)
  fit <- suppressWarnings(fit_default(
    data = sim$data,
    class_init = sim$true_class,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(itermax = 5, tol = 1e-4)
  ))
  expect_true(all(
    c("degenerate", "degeneracy_reason") %in% names(fit$start_summary)
  ))
  expect_true(all(!fit$start_summary$degenerate %in% TRUE))
  expect_equal(fit$n_degenerate_starts, 0)
})

test_that("left-truncated Surv responses raise a clear error", {
  sim <- simulate_mixparfm_data(seed = 59)
  d <- sim$data
  d$entry <- 0
  expect_error(
    mixparfmCWM::fit_mixparfm(
      formula = survival::Surv(entry, time, status) ~ x,
      G = 2,
      class_init = sim$true_class,
      grouping_variable = "family",
      data = d,
      baseline = "weibull",
      frailty = "gamma"
    ),
    "right-censored"
  )
})

test_that("a single Gaussian covariate works", {
  sim <- simulate_mixparfm_data(seed = 55)
  set.seed(55)
  fit <- suppressWarnings(mixparfmCWM::fit_mixparfm(
    formula = survival::Surv(time, status) ~ x,
    G = 2,
    class_init = sim$true_class,
    grouping_variable = "family",
    X_gaussian_variables = "Z1",
    data = sim$data,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(itermax = 5, tol = 1e-4)
  ))
  expect_true(is.finite(fit$classification_loglik))
  expect_equal(length(fit$parameters$X_gaussian_parameters$mu), 2)
})

test_that("legacy control_parfm arguments are accepted and ignored", {
  ctrl <- mixparfmCWM::control_parfm(
    Fparscale = 1,
    showtime = FALSE,
    correct = 0
  )
  expect_true(isTRUE(is.list(ctrl)))
})

test_that("stochastic updates raise an informative error", {
  expect_error(
    mixparfmCWM::control_EM(E_step_update = "stochastic"),
    "stochastic"
  )
})

test_that("the frailty-variance floor is enforced and reported", {
  sim <- simulate_mixparfm_data(seed = 56)
  # a floor above any sensible estimate forces the constraint
  fit <- expect_warning(
    fit_default(
      data = sim$data,
      class_init = sim$true_class,
      baseline = "weibull",
      frailty = "gamma",
      control_EM_algorithm = mixparfmCWM::control_EM(itermax = 3, tol = 1e-4),
      control_parfm_algorithm = mixparfmCWM::control_parfm(theta_floor = 5)
    ),
    "hit the floor"
  )
  theta <- as.numeric(fit$parameters$AFT_parameters["theta", ])
  expect_true(all(theta <= 5 * (1 + 1e-6)))
})

test_that("no-frailty fits without a grouping variable", {
  sim <- simulate_mixparfm_data(seed = 57)
  set.seed(57)
  fit <- mixparfmCWM::fit_mixparfm(
    formula = survival::Surv(time, status) ~ x,
    G = 2,
    class_init = sim$true_class,
    grouping_variable = NULL,
    data = sim$data,
    baseline = "weibull",
    frailty = "none",
    control_EM_algorithm = mixparfmCWM::control_EM(itermax = 5, tol = 1e-4)
  )
  expect_true(is.finite(fit$classification_loglik))
  expect_null(fit$frailty_effect)
})
