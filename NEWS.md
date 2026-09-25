# mixparfmCWM 0.1.0

## Breaking changes

* The package no longer depends on the modified `parfm` fork at runtime
  (`parfm` moves from `Imports` to `Suggests`, where it is only used by the
  tests that check likelihood equivalence). The component-specific
  parametric frailty models are now fitted by an internal optimizer.
* The estimation algorithm is now ICM-CEM: a classification EM (CEM)
  algorithm whose classification step is carried out by iterated
  conditional modes (ICM), updating one label at a time while retaining
  the joint group-by-component shared-frailty contribution. The previous
  E-step used per-patient marginal likelihoods, which is inconsistent
  with a frailty shared at the group level and could decrease the
  classification log-likelihood (observed decreases of up to 15.1
  log-units on real data).
* `E_step_update = "stochastic"` (SEM) is no longer available;
  `control_EM()` raises an informative error for it.
* The license changed from MIT to GPL-2 because the baseline hazard and
  frailty likelihood routines are adapted from `parfm` (GPL-2), with
  attribution in `DESCRIPTION` and in the source files.
* Baseline hazards `logskewnormal` and frailties `possta`, `lognormal` and
  `loglogistic` are no longer supported (the latter was never actually
  fitted by `parfm`); supported baselines are exponential, Weibull,
  inverse Weibull / Frechet, Gompertz, lognormal and loglogistic,
  supported frailties are gamma and inverse Gaussian (plus `none`).
* In the returned object, `z` is now the hard cluster-membership indicator
  matrix (previously a matrix of posterior probabilities), and
  `fit_parfm` is a list of internal fit summaries (previously a list of
  `parfm` objects). The returned `algorithm` element is `"icm_cem"`
  (previously `"conditional_cem"`). The `bic` element is a classification BIC computed
  from the classification (not the observed-data marginal) log-likelihood,
  as before, and this is now documented explicitly.
* `control_EM()` drops the vestigial `err` argument and gains
  `objective_tolerance`, `move_tolerance`, `max_c_sweeps`,
  `minimum_component_size` and `n_start`. `control_parfm()` drops the
  parfm-specific `Fparscale`, `showtime` and `correct` arguments (they
  are accepted and ignored for backward compatibility) and gains
  `theta_floor`. The `plyr` dependency was dropped.
* Stratified models (`strata` argument) now raise an error instead of
  being passed to `parfm`.
* Left-truncated (counting-process) `Surv` responses now raise an error;
  only right-censored responses are supported. The left-truncation
  likelihood previously used (risk based on `H(t) - H(t0)`) is
  inconsistent for shared-frailty models, which require conditioning on
  cluster survival to entry, so the feature was removed.
* The default `n_start` in `control_EM()` is 5 (multiple random starts
  with best-objective selection); the classification-step tolerance
  `move_tolerance` is now relative (`1e-9` by default, scaled by
  `1 + |score|`) instead of absolute.

## New features

* `n_start` in `control_EM()` (default 5): multiple starts (the
  user-supplied partition plus random partitions) with
  best-classification-log-likelihood selection; per-start results are
  reported in `start_summary`.
* The classification step, the M-step optimizer and the reported objective
  all evaluate a single likelihood implementation, so the fast/oracle
  divergence of the prototype cannot occur by construction.
* The M-step refits are warm-started from the current component estimates
  and are protected by a per-component acceptance rule: a candidate that
  does not improve a component's survival log-likelihood is rejected and
  the previous estimates are retained, instead of aborting the fit.
* A floor (`theta_floor`, default 1e-4) is enforced on the estimated
  frailty variance, with a warning whenever the final estimate is
  constrained by the floor.
* Robustness fixes from the algorithm audit: empty components in
  `class_init` are reported with a clear error before any fitting;
  components are kept non-empty during the classification sweeps
  (`minimum_component_size`); a single cluster-weighted Gaussian covariate
  no longer crashes `mclust::cdens`; optimizer failures inside a start do
  not abort the call when `n_start > 1`; the design matrix is built once
  from the full data so component-level and full-data evaluations can
  never disagree.
* Robustness fixes from the second algorithm audit (audit 07):
  - the inverse-Gaussian frailty term no longer underflows for small
    frailty variances (the exponentially scaled Bessel function and a
    cancellation-free form of the root term replace the `parfm`
    formulation, which turned into a `-Inf` cliff near
    `theta = 1.4e-3` and trapped the M-step optimizer);
  - the ICM classification step no longer crashes with NaN
    scores when a component's parameters run away: non-finite scores are
    recomputed from exactly rebuilt cell sums, and the cell sums are
    rebuilt from scratch at the start of every sweep;
  - degenerate solutions are detected and flagged: `control_EM()` gains
    `minimum_component_events` (default: the number of survival
    parameters per component, enforced like `minimum_component_size`)
    and `runaway_parameter_limit` (default 30 on the optimizer scale).
    Flagged starts are excluded from the multi-start selection and
    reported in `start_summary` and in the new `n_degenerate_starts`
    element; previously they were silently returned as converged fits
    with absurd parameters;
  - the M-step acceptance rule is now a strict non-decrease
    (`candidate$loglik >= old_value`) instead of allowing decreases of
    up to `objective_tolerance` per component.
* The returned object gains `classification_loglik`,
  `classification_bic`, `classification_loglik_vec`, `cem_trace`
  (objective before the classification step, after it, and after the
  M-step at every iteration), `converged`, `stopping_reason`, `n_iter`,
  `reached_itermax`, `final_relative_change`, `n_start`,
  `selected_start`, `start_summary`, `n_m_step_rejections` and
  `theta_floor`.
* The empirical-Bayes frailty predictions (`frailty_effect`,
  `frailty_var_effect`) are computed with the same Laplace-transform
  routines used by the likelihood.
