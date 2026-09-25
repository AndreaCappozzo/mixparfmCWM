# The classification log-likelihood must never decrease at any sub-step.

test_that("the objective trace is non-decreasing at every sub-step", {
  sim <- simulate_mixparfm_data(seed = 21)
  fit <- suppressWarnings(fit_default(
    data = sim$data,
    class_init = sim$true_class,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(
      itermax = 10,
      tol = 1e-4,
      n_start = 2
    ),
    control_parfm_algorithm = mixparfmCWM::control_parfm(maxit = 500)
  ))

  trace <- fit$cem_trace
  expect_true(all(trace$c_step_delta >= -1e-8))
  expect_true(all(trace$m_step_delta >= -1e-8))
  expect_true(all(trace$total_delta >= -1e-8))
  expect_equal(fit$loglik_vec, trace$objective_after_m)
  expect_equal(fit$algorithm, "icm_cem")
})

test_that("a second start cannot beat the selected start", {
  sim <- simulate_mixparfm_data(seed = 22)
  set.seed(22)
  fit <- suppressWarnings(fit_default(
    data = sim$data,
    class_init = sim$true_class,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(
      itermax = 10,
      tol = 1e-4,
      n_start = 3
    )
  ))

  ok_rows <- !is.na(fit$start_summary$classification_loglik)
  expect_true(any(ok_rows))
  expect_equal(
    fit$classification_loglik,
    max(fit$start_summary$classification_loglik, na.rm = TRUE)
  )
  expect_equal(
    fit$selected_start,
    which.max(
      replace(fit$start_summary$classification_loglik, !ok_rows, -Inf)
    )
  )
})
