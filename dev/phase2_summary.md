# Phase 2 implementation summary

Date: 2026-09-24. Follows the audit in `dev/icm_audit.md` and the decisions
recorded in its Section 10. Version bumped to 0.1.0.

## What changed

**Estimation algorithm (`R/fit_mixparfm.R`, rewritten).** Conditional
classification EM: the allocation step is the audited ICM update (one label
at a time, group-by-component frailty terms recomputed for source and
destination cells); the M-step uses closed-form updates for tau, Gaussian and
multinomial covariate parameters, and warm-started numerical refits for each
component's survival block with a per-component acceptance rule (a candidate
that does not improve a component's contribution is rejected, the previous
estimates are retained; the partition is never rolled back and the algorithm
does not halt on a rejection). Internal multi-start (`n_start`) with
best-classification-log-likelihood selection; per-start failures are
recorded in `start_summary` and do not abort the call.

**Single-path likelihood (`R/likelihood.R`, new).** Baseline hazards
(exponential, Weibull, inverse Weibull / Frechet, Gompertz, lognormal,
loglogistic), frailty Laplace-transform terms (none, gamma, inverse
Gaussian — adopting parfm's exact evaluation formulas), per-component terms,
the classification log-likelihood, natural-scale parameter transformation,
and empirical-Bayes frailty prediction (ratios of Laplace-transform
derivatives, replacing `parfm::predict.parfm`). One implementation is used
by the C-step, the M-step optimizer, the trace, the BIC and multi-start
selection — the prototype's fast/oracle divergence cannot occur by
construction. Adapted code carries GPL-2 attribution (file header and
DESCRIPTION); the package license changed from MIT to GPL-2.

**Internal component fitter (`R/component_fit.R`, new).** Replaces
`parfm::parfm()`. nlminb (default) or optim (Nelder-Mead / BFGS);
deterministic data-driven default starts; `inip`/`iniFpar` semantics
preserved; box-constrained frailty variance with `theta_floor` (default
1e-4) and a warning when the floor binds in the final fit; non-finite
objective evaluations are penalized instead of crashing; failures are
reported, not thrown.

**Audit fixes applied.** Empty components in `class_init` validated before
any fitting; components kept non-empty during sweeps
(`minimum_component_size`); the single-Gaussian-covariate `mclust::cdens`
crash fixed (matrices used consistently); design matrix built once from the
full data (full-data and component-level evaluations can never disagree);
per-start `tryCatch` so a degenerate start cannot kill a multi-start call;
`strata` raises a clear error.

**Controls (`R/control.R`).** `control_EM()` gains `objective_tolerance`,
`move_tolerance`, `max_c_sweeps`, `minimum_component_size`, `n_start`;
drops the vestigial `err`; `E_step_update` retained but only
`"classification"` accepted. `control_parfm()` gains `theta_floor`;
`Fparscale`, `showtime`, `correct` accepted and ignored for backward
compatibility.

**Breaking changes (also in NEWS.md).** No runtime dependency on the
modified parfm (moved to Suggests, tests only); `plyr` dropped (base
`merge` for the frailty tables); SEM removed; `logskewnormal` baseline and
`possta`/`lognormal`/`loglogistic` frailties removed from the supported
set; `z` returned as hard indicators (was posterior probabilities);
`fit_parfm` is a list of internal fit summaries (was a list of `parfm`
objects); AFT row names use a single `beta.` prefix (was parfm's
`beta.beta.`); `frailty_effect` rows are sorted by group (base merge, was
plyr join order); stratified models error out (previously passed through
to parfm).

## What was tested (all passing; `devtools::test()` 132/132, 0 fails/skips)

- **Likelihood equivalence with the old `loglikelihood_i()`** to < 1e-8
  across 6 baselines x 3 frailties x 3 parameter draws (ran against the
  installed parfm fork; the test skips gracefully if the installed parfm
  lacks `loglikelihood_i()`), plus group-level equivalence with
  `parfm::Mloglikelihood()` for G = 1.
- **Non-decreasing objective trace**: `cem_trace` sub-step deltas all
  >= -1e-8 in every test fit; multi-start selection returns the best start.
- **Parameter recovery** on well-separated simulated data (ARI >= 0.75 —
  the classification criterion genuinely misassigns covariate-ambiguous
  patients; a stricter threshold fights the estimator, not the code).
- **G = 1 vs a direct `parfm::parfm()` fit**: log-likelihood and
  theta/rho/lambda agreement within 1e-4 relative.
- **Edge cases**: empty `class_init` component (clear error), degenerate
  start does not abort multi-start, singleton group, zero-event component,
  single Gaussian covariate, theta-floor warning, no-frailty without
  grouping variable, legacy `control_parfm` arguments, SEM rejection.
- **All 21 baseline x frailty combinations** fit without error and stay
  monotone.
- **Replay benchmark** (not in the test suite; run manually): the new
  implementation reproduces the corrected prototype's solution on
  `conditional_cstep_share/data/replay_case.rds` exactly — same 6
  iterations, final loglik -8790.94111074466 (diff 1.5e-8), 100% class
  agreement, identical parameters, min total delta +0.068 — in 1.2 s
  versus the prototype's 48.5 s.
- **`devtools::check()`**: 0 errors, 0 warnings, 0 notes. All examples run
  (0.4 s).

## Still unresolved / deliberately left open

1. **Manuscript**: handled separately per your decision (the E-step,
   C-step, SEM and BIC sections still describe the old algorithm as of
   `wileyNJDv5_AMA.tex`).
2. **The unreproducible +341 fast/oracle case** from the audit is moot by
   construction (single likelihood path), but no root cause was found.
3. **Classification-BIC caveat** remains by design: the criterion is an
   ICL-type one and can favour over-separation/under-selection under
   overlap; documented in `fit_mixparfm()` details and README, but the
   numerical consequences were not re-quantified beyond the audit's study.
4. **`minimum_component_size = 1`** still allows singleton components when
   there are no Gaussian covariates; with Gaussian covariates a singleton
   start is rejected as degenerate (test covers this). A general
   covariance-eigenvalue floor (bounded-likelihood estimation) was not
   implemented — the classification likelihood is formally unbounded
   (demonstrated in the audit).
5. **`err` argument of `control_EM`** removed rather than ignored; old
   scripts passing it will error.
6. **ORCID warning**: `0000-0002-9200-3194` in DESCRIPTION is flagged as
   invalid by `person()` (pre-existing metadata; not guessed at). It does
   not surface in `R CMD check`.
7. **`AFT_parameters` row naming** change (`beta.beta.x` -> `beta.x`) may
   break downstream analysis code that indexes by row name.
