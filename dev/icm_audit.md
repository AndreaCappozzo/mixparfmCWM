# Audit of the conditional C-step (ICM) algorithm for mixparfmCWM

Date: 2026-09-24
Audited artifacts:

- Manuscript: `useful_for_implementation/wileyNJDv5_AMA.tex` (read-only)
- Prototype: `useful_for_implementation/conditional_cstep_share/packages/mixparfmCWMcorrected/`
  (incl. `R/shared_objective.R`, the ICM note `note/conditional_cstep_note.tex`, results in `results/`)
- Current package: `R/` at the repository root (v0.0.1, flawed EM/CEM E-step)
- parfm fork: `useful_for_implementation/parfm` (also installed as `parfm` 2.7.6 in the session library)

Everything reported below was either read from the code/manuscript or executed; commands and
outputs are given. `R/` was not modified.

---

## Verdict

**Sound after specific fixes.** The mathematical design of the proposed algorithm — sequential
conditional-mode (ICM) updates of patient labels that preserve the hospital-by-component
shared-frailty term, alternating with independent parfm refits per component — is correct for
the classification objective of the manuscript (Eq. `eq: compl loglik 5`, compiled as Eq. (7)),
and its monotonicity argument holds. I verified the objective implementation against the
manuscript expression and against parfm's own likelihood, reproduced the prototype's stored
results to ~1e-12, and confirmed empirically that the objective never decreases at any
allocation sweep or M-step in any run that completed.

However, the prototype is **not production-ready**. It contains one confirmed crash-level bug
(single Gaussian covariate), a runtime "oracle" assertion that aborts legitimate fits at
degenerate parameter corners (observed failure rates of 10/16 to 15/24 runs in my simulations),
no error handling around parfm refits (parfm optimizer failures propagate as cryptic errors),
no M-step warm starts, and several edge cases that terminate in obscure errors instead of
graceful behavior. The manuscript also disagrees with the prototype in substantive ways
(Section 1.4 below) and must be revised; as written, the manuscript's E-step is the flaw the
prototype was built to remove, and its monotonicity claim is false for its own algorithm.

Phase 2 should proceed only after the fixes in Section 8 are addressed.

---

## 1. The objective

### 1.1 The manuscript expression

The classification log-likelihood (`\label{eq: compl loglik 5}`, line 350 of the .tex; note that
despite the label it compiles as Equation (7), and the prototype's note and code correctly call
it Eq. (7)) is:

$$
\ell_c(\psi) \;=\; \sum_{g=1}^G \Bigg\{
\sum_{j=1}^J \Bigg[
\underbrace{\sum_{i \in R_{jg}} \delta_{ij}\big(\log h_0(y_{ij};\gamma_g) + x_{ij}^T\beta_g\big)}_{\text{event log-hazards, additive per patient}}
\;+\;
\log\Big[(-1)^{d_{jg}} \mathcal{L}^{(d_{jg})}\Big(\sum_{i \in R_{jg}} H_0(y_{ij};\gamma_g)e^{x_{ij}^T\beta_g};\ \theta_g\Big)\Big]
\Bigg]
\;+\;
\sum_{j=1}^J \sum_{i \in R_{jg}} \big(\log \tau_g + \log \phi(u_{ij};\mu_g,\Sigma_g) + \log \xi(v_{ij};\pi_g)\big)
\Bigg\}
$$

with $d_{jg} = \sum_{i\in R_{jg}} \delta_{ij}$. In code quantities (prototype's
`shared_objective.R`):

- additive per-patient term: `additive[i, g]` = `log tau_g` + Gaussian log-density + multinomial
  log-probability + `event[i] * (log h0 + x beta)` — exactly the first and third lines.
- frailty-integrated group likelihood, for gamma frailty with mean 1 and variance $\theta$:
  $F_g(d, s) = d\log\theta + \log\Gamma(d + 1/\theta) - \log\Gamma(1/\theta) - (d + 1/\theta)\log(1+\theta s)$,
  implemented in `.mixparfm_gamma_log_laplace_derivative(events, cumulative_hazard, theta)` and
  summed over all hospital-by-component cells as
  `.mixparfm_frailty_log_term(sum events in cell, sum risk in cell, component)`. I verified
  symbolically that $(-1)^d \mathcal L^{(d)}(s) = e^{F_g(d,s)}$ for the gamma Laplace transform
  $\mathcal L(u) = (1+\theta u)^{-1/\theta}$, i.e. that $e^{F_g(d,s)}$ is exactly
  $\int_0^\infty m^d e^{-ms} f_M(m;\theta)dm$, the group likelihood with the shared frailty
  integrated out. Constants are exact (no terms dropped); the expression is the exact
  classification log-likelihood, not up-to-constants.
- Parameterization matches parfm's optimizer scale: `p[1]` = log θ (gamma/ingau/lognormal),
  `exp(-exp(p[1]))` for possta, then transformed baseline parameters, then betas; consistent
  with `loglikelihood_i()`'s parsing in the parfm fork.

The prototype additionally recomputes the same objective through parfm itself
(`.mixparfm_classification_objective_parfm`: subset `obsdata_creator` + `-Mloglikelihood` per
component + the same covariate/tau terms) and asserts agreement at every sub-step. In healthy
parameter regions the two agree to ~1e-11 (replay data, n = 1500; and the package's own
`test-objective.R` shows the gamma term matches `parfm::fr.gamma` to <1e-10 over a wide grid).

### 1.2 Confirmed: the current package's E-step is the flaw described

The current `R/fit_mixparfm.R` E-step computes posteriors from `parfm::loglikelihood_i()` on the
**full** obsdata. In the parfm fork, `loglikelihood_i()` integrates the frailty **per patient**
(`logSurv_i[i] <- fr.gamma(k = obs$event[i], s = cumhaz_i[i], ...)`; `utils.R`, lines 290+),
although the model shares the frailty at the group level. The resulting per-patient "marginal"
is the likelihood under an *individual* frailty, which does not exist under the shared-frailty
model. Quantitatively, for a group of $n$ patients who all have events with cumulative-hazard
sums $s_i$ (θ = 1, $s_i = 0.5$):

| n (patients in group) | 2 | 5 | 10 | 20 |
|---|---|---|---|---|
| log-ratio of joint (correct) to product of per-patient marginals | 0.236 | 1.326 | 3.504 | 8.198 |

The error grows with group size; at the replay scale (hundreds of patients per hospital) the
per-patient product misstates the group likelihood by orders of magnitude.

Note the **old package's reported loglik is still Eq. (5)**: `log_density_parfm_g` is computed
from `-Mloglikelihood` on each component's subset (group-integrated). So the old algorithm
computes the right objective but allocates with the wrong likelihood — which is why its trace
can decrease (see Section 6.1).

### 1.3 Mismatch list (manuscript vs prototype vs current package)

1. **Manuscript E-step (eq:estep) uses per-patient marginals** — the same flaw as the current
   package. $p(y_{ij}|x_{ij};\gamma_g,\beta_g,\theta_g)$ with $(-1)^{\delta_{ij}}\mathcal L^{(\delta_{ij})}(H_0 e^{x'\beta};\theta)$
   is a per-patient integrated likelihood. The manuscript's own C-step (MAP over these
   posteriors) is therefore *not* what the prototype implements, and it is not an ascent
   algorithm for Eq. (5). The manuscript must be revised to describe the conditional/ICM update
   (the note in `conditional_cstep_share/note/` is the correct derivation).
2. **Manuscript claims CEM monotonicity and SEM.** With the manuscript's E-step, the
   monotonicity claim is false (reproduced decrease of −15.1164, Section 6.1). The prototype
   drops SEM entirely (`E_step_update = "stochastic"` is rejected); the manuscript still
   presents SEM and uses it to argue robustness against local optima.
3. **Equation numbering:** the label `eq: compl loglik 5` compiles as Equation (7); the
   prototype consistently says "Eq. (7)". Harmless but should be reconciled.
4. **Frailty `"loglogistic"`** is accepted by `fit_mixparfm`'s signature in the prototype but
   is unsupported in `shared_objective.R` (hard `stop()` at runtime). The corrected package
   silently advertises an option it cannot fit.
5. **Single Gaussian covariate crashes** (confirmed bug, Section 7.1): the corrected package
   passes a 1-column data.frame to `mclust::cdens` and dies with "'list' object cannot be
   coerced to type 'double'". The old package passed the same data as a vector and worked.
6. **BIC:** prototype computes a *classification* BIC from Eq. (5) (formula matches the
   manuscript's d, and `2ℓ − d ln N` matches eq:BIC; the Gaussian part `nMclustParams − (G−1)`
   agrees with the manuscript's $G\,p(p{+}3)/2$ for VVV). But the manuscript presents it as the
   Schwarz BIC of the model; the criterion is not on the observed-data likelihood (Section 5).
7. The prototype branched off an **older snapshot** of the package than the current `R/`
   (current `R/` has `tol_zero_var` guards, `converged`/`stopping_reason` fields, degenerate
   covariance validation, which the prototype copy lacks or has in different form). Phase 2
   must merge, not copy.

---

## 2. Allocation step (conditional C-step)

- **Allocation unit.** Manuscript: patient ($z_{ijg}$ per observation; `eq:cstep`). Prototype:
  patient (sequential sweep over rows). Patients of the same group *can* end up in different
  components.
- **Is a split group well-defined?** Yes, in the manuscript as written: $R_{jg}$ is defined per
  (group, component) pair, and the frailties $m_{jg}$ are cluster-specific and independent
  across $g$ ("a given group $j$ is associated with distinct frailty terms $m_{jg}$ across the
  $G$ clusters"). Eq. (5) sums $F_g(d_{jg}, s_{jg})$ over all (j,g) cells, so a split group
  contributes one integrated likelihood per piece. The prototype computes exactly this.
  (Statistically, splitting severs the dependence between the two pieces — each piece gets its
  own independent frailty. That is the model as specified; see open questions.)
- **Does a move recompute the group frailty likelihood in source and destination?** Yes.
  `.mixparfm_conditional_classification_step` maintains cell sums $D[j,g]$ (events) and
  $S[j,g]$ (cumulative hazard), removes patient $i$ from its current cell before scoring, and
  the score is $Q_{ijg} = a_{ijg} + F_g(D{+}\delta_i, S{+}w_i) - F_g(D, S)$ — destination gain
  relative to the common reduced configuration, so the move's exact objective change is
  $Q_{ij,g^*} - Q_{ij,\text{old}}$, which includes the frailty recomputation in both the
  source and destination cells. The per-patient marginal problem does not reappear.
- **Order dependence.** Updates are sequential in row order with immediate acceptance. The
  note acknowledges (and I confirmed numerically, replay data, fixed parameters, same starting
  partition): original visiting order → final objective −10413.86; permuted order →
  −10211.15, different final classification. The endpoint (not monotonicity) depends on the
  visiting order. Across 12 random starts at fixed parameters I obtained 12 distinct
  coordinate-wise optima, ranging from −10469.7 to −10411.9, versus −8840.4 from the k-prototypes
  start. Initialization and sweep order materially affect the result.
- **Empty components.** Prevented by the minimum-size constraint
  (`if (g != old && counts[old] <= minimum_component_size) next`); the admissible set is
  "no empty component" (the note's $\mathcal C$). This is enforced inside the C-step; the
  incoming partition is validated too.

**Conclusion: the allocation step is sound.** Each accepted move strictly increases Eq. (5)
(gain > `move_tolerance = 1e-12`; ties retain the current label), parameters are held fixed
during the sweep, and a post-sweep recomputation of the objective guards against drift.

---

## 3. Parameter step (M-step)

- **Does each parfm refit maximize the same terms the allocation step uses?** Yes. Given the
  hard partition, Eq. (5) is separable across components: the survival block of component $g$
  is exactly $-\texttt{Mloglikelihood}$ of the parfm model fitted to the rows assigned to $g$
  with `cluster = grouping_variable`, and the prototype refits exactly that (subset with
  `droplevels` on the grouping factor). $\tau$, $\mu$, $\Sigma$ (mclust `mstep` on the hard
  indicators) and $\pi$ (within-component frequencies, with absent levels handled) are closed-form
  MLEs of their blocks. The oracle check (`fast` vs `parfm` objective) verifies this
  equivalence at every sub-step of every iteration; agreement was ~1e-11 on the replay data.
- **Warm start?** No. Every refit passes `inip`/`iniFpar` from `control_parfm()` (defaults
  `NULL` → parfm's fixed, deterministic initial values). The current component estimates are
  never used as starts. This wastes work and increases the chance that a refit lands in a
  different local optimum of the parfm sub-problem, triggering safeguard rejections (and with
  them an early stop). Not a correctness bug — a performance and robustness defect.
- **What if the optimizer returns a lower objective?** The safeguard compares the candidate
  state's objective against the post-C-step objective and accepts only if
  $\ge$ (post-C-step objective − `objective_tolerance = 1e-7`). If rejected, the code rolls
  back **the complete previous state** — parameters *and* the pre-C-step partition — and
  terminates the whole algorithm with `stopping_reason = "m_step_rejected"`, `converged =
  FALSE`. Two issues: (i) rolling back the partition is unnecessarily conservative — the new
  partition at the old parameters is still an improvement and is monotone-safe; (ii) a single
  rejected M-step stops the entire fit rather than, e.g., retrying with warm starts. In the
  replay example the safeguard never fired (0 rejections), and on/off runs were bit-identical.

**Conclusion: the parameter step is sound in design** (blockwise exact M-step + explicit
numerical safeguard). The rollback-and-halt policy and the absence of warm starts should be
changed.

---

## 4. Monotonicity and convergence

Argument (matching the note, and verified in code):

1. Within a C-step sweep, every accepted move changes the objective by exactly
   $Q_{ij,g^*} - Q_{ij,\text{old}}$ > `move_tolerance` > 0, so the objective is strictly
   increasing across sweeps; a final recomputation checks the realized gain against
   `objective_tolerance` (assertion, never triggered in any of my runs).
2. Sweeps repeat until a full pass makes no moves (or `max_c_sweeps`); since the objective
   strictly increases per accepted move and the set of admissible partitions is finite, the
   C-step terminates at a coordinate-wise maximum (no single-patient move improves by more
   than ε_move).
3. The exact M-step is non-decreasing; the numerical one is enforced non-decreasing up to
   `objective_tolerance` by the safeguard.
4. The outer loop stops on relative-change tolerance (`tol`), `itermax`, or safeguard
   rejection. Termination is guaranteed.

**Numerical verification** (objective logged after every allocation sweep and every parameter
update):

- Replay case (paper_sim, replication 95, G = 3, start 2), safeguard on: per-iteration trace
  (`cem_trace`): min C-step delta = +0.0337, min M-step candidate delta = +0.0346, min total
  delta = +0.0682; **zero decreases** anywhere. Safeguard on/off identical (max trace
  difference 0, identical final classification).
- Sweep-level logging (my own copy of the C-step that recomputes the objective after every
  sweep): replay data, hard init: +150.807, +0.118, 0.0 over 3 sweeps; 12 random starts: 99
  sweeps total, **minimum per-sweep delta = 0 (never negative)**.
- n = 300 simulation study (16 runs, 6 completed): **0 decreases** in any completed corrected
  run; every recorded `total_delta` ≥ 0.
- By contrast, the current EM's trace (same objective values) decreased in **every scenario**;
  worst observed full-iteration decreases: −15.116 (replay), −6.3 (ovlp/low, n = 300), −45.1 and
  −46.4 (well-separated, n = 300), −42.0 (n = 100 study). The legacy trace in
  `results/classification_loglik_trace.csv` was reproduced to 4.6e-10, including its single
  −15.11642 decrease.

**Conclusion: monotonicity holds in theory and in every completed run.** The only source of
non-monotone behavior is the M-step's numerical optimizer, and the safeguard contains it.

---

## 5. Statistical soundness

What the estimator targets: the maximizer of the **classification (complete-data)
log-likelihood**, not the observed-data marginal likelihood of the mixture. Consequences:

- **Bias under overlap.** Classification-likelihood estimators are consistent only for
  well-separated components (Bryant & Williamson 1978; Celeux & Govaert 1992). When components
  overlap, hard assignments act like a plug-in that shrinks each component's effective sample
  to its "own" patients, biasing $\tau_g$, $\mu_g$, $\beta_g$, $\theta_g$ toward
  over-separation. My n = 300 study is consistent with this but is too small and too
  failure-prone to quantify bias (Section 6.3): ARI of the corrected method was *worse* than
  the EM's in one well-separated/high-frailty scenario (0.469 vs 0.621) despite a higher
  classification objective (+52.3) — maximizing Eq. (5) is not the same as recovering the true
  partition.
- **Local optima.** ICM is a coordinate-ascent method; the note says so, and my measurements
  confirm it: 12 random starts at fixed parameters → 12 distinct local optima (all far below
  the k-prototypes solution); visiting-order changes the endpoint; on the tiny brute-force
  dataset the full algorithm reached the global optimum only from favorable states. Multiple
  starts are mandatory; the prototype has no internal multi-start (single `class_init`
  argument), which is fine as an API but must be stated.
- **Unboundedness.** With unconstrained per-component Gaussian covariate densities and free
  $\Sigma_g$, the classification likelihood is unbounded (the classic mixture pathology). I
  demonstrated this directly: brute-force profile enumeration over all 4,070 viable partitions
  of a 12-patient dataset gives a "global maximum" of **+2.30 at a partition with a 2-patient
  component** (near-singular $\Sigma_g$), while the true partition ranks 6th at −32.07. The
  prototype does not constrain this; it survives because `validate_gaussian_mstep` (in the
  current `R/`) or the oracle assertion *crashes* on degenerate covariances rather than
  constraining them. Phase 2 should decide: eigenvalue floor, or documented reliance on
  initialization.
- **BIC / model selection.** The reported BIC is `2·ℓ_c − d·ln N` on the classification
  likelihood. This is not the Schwarz BIC of the marginal model (which would require the
  observed-data log-likelihood); it behaves like an ICL-type criterion and systematically
  favors fewer, more separated components when overlap is present. Comparing it across G or
  across baselines/frailties is internally consistent (same criterion), but the manuscript's
  presentation of it as "the" BIC (eq:BIC cites Schwarz 1978) is not accurate and should be
  renamed/re-derived. The prototype's documentation already calls it "classification BIC" —
  the manuscript should too.
- **Frailty identifiability.** θ_g needs several groups *with events* per component. In my
  simulations parfm drove θ to ~1e-8 corners in several runs (both methods); in the replay the
  estimated θ's were healthy (0.43–0.72). Zero-event components and θ→0 states are the fragile
  region (Section 7), and the θ→0 corner is exactly where the prototype's oracle assertion
  misfires.
- **Empty components / too few events / initialization.** Empty components are prevented
  (min size 1); but nothing prevents a component with all-censored patients, and that case
  crashes parfm (Section 7.3). Initialization is a single hard partition; the manuscript uses
  k-prototypes + 20 restarts — that is an analysis-level strategy, not in the package.

**Claims in the manuscript affected:** the E-step formula (must be replaced by the
conditional update); the statement that CEM "shares the theoretical guarantees of the standard
EM procedure" (true only for the corrected conditional C-step, and then for the classification
objective, not the observed likelihood); the SEM material (dropped by the prototype); the BIC
interpretation; and any claim that the reported log-likelihood is the model's (marginal)
log-likelihood.

---

## 6. Empirical checks

### 6.1 Reproduction of `conditional_cstep_share/results/`

I re-ran the three fits of `run_example.R` (same data, initialization, seed, controls), with
the prototype's own package copies, writing outputs to /tmp (the source tree was not modified).

| method | iterations | final objective | decreases | elapsed (s) |
|---|---|---|---|---|
| legacy (stored) | 11 | −8806.18469737851 | 1 (−15.1164221183899) | 140.5 |
| legacy (reproduced) | 11 | −8806.18469737851 | 1 (−15.1164221183899) | 82.5 |
| corrected, safeguard on (stored) | 6 | −8790.94111074466 | 0 | 81.7 |
| corrected, safeguard on (reproduced) | 6 | −8790.94111074466 | 0 | 48.5 |
| corrected, safeguard off (reproduced) | 6 | −8790.94111074466 | 0 | 48.6 |

Legacy trace reproduced to 4.6e-10; corrected final value to 3.6e-12; on/off traces identical
(max difference 0, identical classifications), matching `safeguard_effect.csv`. The stored
results are genuine and reproducible on this machine, and they support the prototype's claims:
the corrected fit is monotone, terminates sooner, and reaches a *higher* classification
objective (−8790.94 vs −8806.18) than the legacy fit from the same initialization.

### 6.2 Brute force on a tiny dataset (12 patients, 4 groups, G = 2)

Design: Weibull baselines (λ = 0.4/0.12, ρ = 1.3/1.7), one survival covariate, gamma frailty
(θ = 0.5), two Gaussian cluster-weighted covariates, components separated by 2.5 sd. ψ* taken
from parfm fits to the true components (θ pinned at the true 0.5, since at n = 6 per component
parfm collapses θ to ~1e-8).

- **Fixed ψ*, full enumeration of all 4,096 partitions:** global maximum of Eq. (5) is the
  true partition, −33.8175; runner-up −42.997 (a 9.2 log-unit gap).
- **C-step alone from 200 random partitions:** **100%** reach the global maximum.
- **Full ICM (C + M) from 30 random partitions:** 3/30 completed; 27/30 crashed (21 parfm
  optimizer errors, 6 oracle-assertion trips). None of the 3 survivors ended at the true
  partition (best −45.09 vs −32.07 achievable).
- **Profile enumeration (M-step per partition, 4,070 partitions with ≥2 patients per
  component):** 1,294/4,070 (32%) of M-steps fail outright (parfm error or non-finite
  objective); the best "profile" value is +2.30 at a pathological 2-patient component
  (Section 5, unboundedness).

Interpretation: the allocation logic itself is reliable when the objective is well posed and
parameters are sensible; the failure modes are all in the numerical M-step and the degenerate
regimes it can wander into.

### 6.3 Simulation study: corrected ICM vs current EM

Design: G = 2, J = 30 groups × 10 patients (n = 300) and a smaller n = 100 variant; scenarios
2×2 (Gaussian separation 2.5 vs 0.5; θ = 0.2 vs 2.0); Weibull baseline, gamma frailty, one
survival covariate, two Gaussian and one multinomial cluster-weighted covariates; 40%-noisy
true initialization, identical for both methods; 4–6 replications per scenario; both methods
run from the same partition on the same data. Metrics below are for n = 300 (n = 100 gave even
higher failure rates for both).

| scenario | corrected: completed runs | EM: completed runs | mean Δ loglik (corrected − EM) | mean ARI corr / EM | EM trace decreases (worst) |
|---|---|---|---|---|---|
| well-separated, θ = 0.2 | 1/4 | 4/4 | **+114.35** | 0.705 / 0.614 | 23 (−46.37) |
| well-separated, θ = 2.0 | 2/4 | 4/4 | **+52.26** | 0.469 / 0.621 | 10 (−45.12) |
| overlapping, θ = 0.2 | 3/4 | 4/4 | **+7.49** | 0.095 / 0.092 | 7 (−6.28) |
| overlapping, θ = 2.0 | 0/4 | 0/4 | — | — | — |

- In **every completed comparison the corrected algorithm reached a higher Eq. (5) value**
  (+7.5 to +114.4), and its traces never decreased; the EM's trace decreased in every scenario
  it completed (up to −46.4 in a single iteration). This is direct evidence that (a) the EM's
  allocation step is fighting the objective, and (b) the conditional C-step actually optimizes
  it.
- **Failure rates are the dominant practical problem:** corrected 10/16 failed (6 parfm
  optimizer crashes, 4 oracle-assertion trips); EM 4/16 failed (parfm crashes only). In the
  overlapping/high-θ scenario *both* methods failed on all replications.
- Parameter recovery (6 matched pairs): RMSE(log θ) 8.9 (corrected) vs 10.5 (EM); RMSE(β)
  0.165 vs 0.202 — both dominated by θ-corner degeneracy; too few pairs and too much
  degeneracy for a serious bias/RMSE statement. The one robust parameter-level conclusion is
  that **θ estimates hit the ~1e-8 corner frequently under both methods at these sizes**, and
  only the corrected method then crashes (via the oracle assertion).
- ARI: corrected better in two scenarios, worse in one (Section 5) — higher classification
  objective does not imply better clustering.

### 6.4 Order-dependence and initialization evidence

Summarized in Section 2: endpoints, not monotonicity, depend on visiting order and
initialization. The manuscript's k-prototypes + 20 restarts strategy is the right response and
should be retained (and documented as *required*, not optional).

---

## 7. Edge cases and bugs found (all reproduced)

1. **Single Gaussian covariate: confirmed crash.** `X_gaussian_variables = "Z"` (length 1):
   `shared_objective.R` passes a 1-column data.frame to `mclust::cdens` →
   `'list' object cannot be coerced to type 'double'`. Reproduced with the package's own
   state-construction path; the old package worked because it passed a vector. Any dataset
   with exactly one continuous cluster-weighted covariate cannot be fitted by the prototype.
2. **Oracle assertion misfires at degenerate corners.** The runtime check
   `fast ≈ parfm-oracle` (atol 1e-7, rtol 1e-8) holds to ~1e-11 in healthy regions but
   diverges to 1e-5–2e-4 when parfm returns θ ≈ 1e-8 (the fast and oracle survival terms are
   evaluated with different orderings, and the θ→0 corner involves cancellation between
   `d·log θ` and `lgamma(d + 1/θ) − lgamma(1/θ)`). Observed trip messages include
   `difference=-0.000184, allowed=5.11e-06`. In one (not deterministically reproducible)
   run the fast and oracle objectives differed by **885** units (fast = +341.26,
   oracle = −543.35) — a genuine large disagreement I could not pin down; treat as unresolved
   (Section 9). Consequence: the assertion converts numerical corner noise into hard crashes
   (6/15 of the corrected failures in the n = 100 study).
3. **Zero-event component: crash.** A partition whose component contains only censored
   patients makes parfm's optimizer fail (`non-finite finite-difference value`) — no handling,
   no informative message. Both packages share this failure.
4. **parfm optimizer failures propagate.** There is no `tryCatch` around the per-component
   `parfm::parfm()` calls in the M-step; failures surface as `non-finite finite-difference
   value [k]` or `'list' object cannot be coerced to type 'double'` (the latter from
   `estim_par <- as.numeric(res[...])` in parfm.R when optimx fails hard). These accounted
   for 21/27 crashes in the tiny-data ICM runs and 6/16 at n = 300.
5. **Empty component in `class_init`:** the intended validation
   (`contains empty components: ...`) is unreachable in practice, because the *initial* parfm
   fits on the empty subset fail first with the cryptic `no rows to aggregate`.
6. **Singleton group (one patient in a group):** works (the shared frailty of a one-patient
   group degenerates to an individual frailty); completed with warnings. Note that with many
   singleton groups θ is unidentifiable and the optimizer wanders.
7. **Singleton component (n_g = 1):** not reachable via the C-step (min-size constraint), but
   reachable via `class_init`; then mstep gives a degenerate covariance and the current `R/`'s
   `validate_gaussian_mstep` stops with a (clear) error, while the classification objective is
   formally unbounded there (Section 5).
8. **Subset design-matrix misalignment (structural, not reproduced):** `.mixparfm_component_terms`
   evaluates betas against the *full-data* design matrix, while the M-step fits parfm on the
   subset. If a factor in the survival formula loses a level in some component's subset, the
   subset fit returns fewer betas and the two disagree (crash or oracle trip). The corrected
   code `droplevels` only the grouping variable. In the replay data all levels are present in
   all components; the hazard is real for categorical survival covariates.

---

## 8. Required fixes before Phase 2

Ordered by severity:

1. **Downgrade the runtime oracle assertion to a diagnostic** (warning + recorded diagnostic,
   never a hard stop), and use the fast objective for all monotonicity decisions. Or keep the
   assertion but bound θ (and other parameters) away from degenerate corners so the two
   computations actually agree. As shipped, the assertion is the single largest source of fit
   failures.
2. **Wrap every per-component parfm refit in error handling**: on failure, keep the previous
   component parameters for that g (a partial-M step is still monotone if each accepted block
   does not decrease the objective) instead of crashing; retry with warm start before giving
   up; emit an informative message.
3. **Fix the single-Gaussian-covariate path** (pass a matrix/vector to `cdens` consistently;
   the crash is deterministic).
4. **Warm-start the M-step** from the current component estimates (pass `inip`/`iniFpar` from
   the current fit). This should also reduce (2) dramatically.
5. **On safeguard rejection, do not roll back the partition and do not halt**: keep
   `c^(k+1/2)` with the old parameters (still monotone), and either retry the M-step or stop
   only after repeated rejections.
6. **Handle the θ→0 corner explicitly**: lower bound on the frailty variance (e.g. θ ≥ 1e-4),
   or an explicit "frailty collapsed" state, so the algorithm degrades to a no-frailty model
   instead of wandering.
7. **Validate `class_init` before the first parfm call** (so empty components produce the
   intended message), and add the missing `loglogistic`-frailty objective (or remove it from
   the signature).
8. **Zero-event components:** detect before refitting and either merge into the closest
   component or stop with an informative error.
9. **Design-matrix alignment:** build the subset design with the full data's factor levels
   (or reject partitions that would drop survival-formula factor levels).
10. **Merge into the current `R/` code base** (which is ahead of the prototype's snapshot:
    `tol_zero_var` guards, `converged`/`stopping_reason`, etc.), preserving the corrected
    algorithm and its `cem_trace`.

---

## 9. Open questions for you

1. **Manuscript revision scope.** The E-step (eq:estep), the C-step (eq:cstep), the SEM
   material, the monotonicity claim, and the BIC presentation all need rewriting to match the
   audited algorithm. Do you want me to draft the replacement subsections as part of Phase 2
   (in a separate file — the .tex is read-only for me), or will the manuscript be handled
   separately?
2. **The BIC.** Keep the classification BIC but relabel it (and cite Celeux & Govaert / ICL-type
   criteria), or also compute the observed-data (marginal) log-likelihood for reporting? The
   latter requires the group-integrated observed likelihood summed over components — feasible
   (it is Σ_g τ_g-probabilistic mixture of the per-component group likelihoods) but changes
   what "loglik" means in the output.
3. **The unresolved +341 fast-vs-oracle case.** One simulation run showed a 885-unit
   disagreement between the fast and parfm objectives that I could not reproduce
   deterministically. Do you want a targeted hunt (instrumented runs over many seeds), or
   should the fast objective simply be re-derived to share the exact parfm evaluation path so
   the question cannot arise?
4. **Group splitting.** The model as written allows patients of one group in different
   components, each piece with an independent frailty. Should the C-step optionally constrain
   allocation at the group level (all-or-nothing per hospital), as some applications of shared
   frailties intend? This is a modeling decision, not a bug.
5. **θ degeneracy policy.** Floor on θ, or explicit fallback to `frailty = "none"` when the
   estimate collapses? This affects BIC comparability across frailty distributions.
6. **Warm starts / multi-start.** OK to add warm starting (change 4) and, optionally, an
   internal `n_start` with best-objective selection? The manuscript's current practice (20
   external restarts) would remain the recommended default.

---

## Reproduction notes

- Environment: R 4.6.0, parfm 2.7.6 (the modified fork, installed), mclust 6.1.3, macOS arm64.
- The prototype packages were loaded with `pkgload::load_all()` from
  `useful_for_implementation/conditional_cstep_share/packages/`; no files in
  `useful_for_implementation/` or `conditional_cstep_share/` were modified; scratch outputs
  were written to `/tmp` (`/tmp/icm_sim300.rds`, `/tmp/icm_sim_summary.rds`).
- Timing: legacy replay fit 82 s, corrected replay fit 49 s, profile enumeration 421 s
  (10 PSOCK workers), n = 300 simulation study 1555 s (10 workers).

---

## 10. Decisions on the open questions (2026-09-24)

1. **Manuscript revision:** handled separately by the authors; Phase 2 concerns code only.
2. **BIC:** keep computing the classification BIC and keep the `bic` element name; the
   documentation will state explicitly that it is a classification (not Schwarz) BIC.
3. **Fast/oracle divergence:** eliminated by construction — the rebuilt package has a single
   likelihood implementation used by the C-step, the M-step optimizer, and the safeguard;
   there is no second ("oracle") code path to disagree.
4. **Group splitting:** keep the model as written (patient-level allocation; split groups
   contribute one independent frailty per (group, component) piece).
5. **θ degeneracy:** floor on the frailty variance in the optimizer; a warning is emitted
   whenever the floor is hit.
6. **Warm starts / multi-start:** M-step refits are warm-started from the current component
   estimates, and `fit_mixparfm` gains an internal `n_start` (via `control_EM`) with
   best-objective selection across starts.
