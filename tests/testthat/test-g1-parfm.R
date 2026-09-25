# With G = 1 the mixture fit must coincide with a single direct parfm fit
# of the same shared-frailty model.

test_that("G = 1 matches a direct parfm fit", {
  skip_if_not_installed("parfm")
  sim <- simulate_mixparfm_data(seed = 41)
  d <- sim$data
  formula <- survival::Surv(time, status) ~ x

  fit1 <- mixparfmCWM::fit_mixparfm(
    formula = formula,
    G = 1,
    class_init = rep(1, nrow(d)),
    grouping_variable = "family",
    data = d,
    baseline = "weibull",
    frailty = "gamma",
    control_EM_algorithm = mixparfmCWM::control_EM(itermax = 20, tol = 1e-6),
    control_parfm_algorithm = mixparfmCWM::control_parfm(maxit = 1000)
  )

  pf <- suppressWarnings(parfm::parfm(
    formula = formula,
    cluster = "family",
    data = d,
    dist = "weibull",
    frailty = "gamma",
    maxit = 1000
  ))

  # log-likelihoods agree
  pf_loglik <- as.numeric(logLik(pf))
  expect_lt(
    abs(fit1$classification_loglik - pf_loglik) / (1 + abs(pf_loglik)),
    1e-4
  )

  # parameter estimates agree
  pf_par <- as.numeric(attributes(pf)$estim_par)
  theta_pf <- exp(pf_par[1])
  expect_lt(
    abs(fit1$parameters$AFT_parameters["theta", 1] - theta_pf) /
      (1 + theta_pf),
    1e-2
  )
  rho_pf <- exp(pf_par[2])
  expect_lt(
    abs(fit1$parameters$AFT_parameters["rho", 1] - rho_pf) / (1 + rho_pf),
    1e-2
  )
  lambda_pf <- exp(pf_par[3])
  expect_lt(
    abs(fit1$parameters$AFT_parameters["lambda", 1] - lambda_pf) /
      (1 + lambda_pf),
    1e-2
  )
})
