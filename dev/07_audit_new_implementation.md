# Audit 07 — `new_implementation/mixparfmCWM` against `new_c_step/modified_cstep.pdf`

Date: 24 September 2026 · Package audited: `new_implementation/mixparfmCWM` v0.1.0 (uncommitted working tree, files dated 24 Sep 14:11–15:15) · Specification: *Note on the modified C-step* (3 Sep 2026) · Reproduction scripts: `audit/scripts/08_*.R`

Evidence labels: **Verified** = reproduced numerically here or read directly in the source · **Inference** = follows from verified facts · **Unverified** = plausible, not tested.

---

## 0. Bottom line

1. **The code does what `modified_cstep.pdf` asks, and it is a correct monotone algorithm for the classification log-likelihood (Eq. 7 of the manuscript).** I did not rely on the package's own tests. I rebuilt the objective independently from the model definition, and the package matches it to within 3×10⁻¹³. The package's classification sweep gives the same partition, patient by patient, as a brute-force sequential argmax of Eq. (7) in 32 of 32 cases. At termination no single-patient move improves the objective, and in 120 full fits the objective never decreased (smallest step −9×10⁻¹³, which is rounding error). The closed-form M-step blocks are exact MLEs.
2. **"Correct" here means a local, monotone maximiser, not a global one.** At fixed parameters the classification step reached the global maximum over all 4,094 admissible partitions from every random start in only 4 of 6 tiny datasets. The criterion also has degenerate boundary solutions (§3, N2).
3. **Before the revision's re-runs, fix four things** (§3):
   - **N1:** a NaN crash in the C-step when a component's parameters run away.
   - **N2:** degenerate components (zero-event, or tiny with runaway parameters) are reported as "converged".
   - **N3:** a numerical underflow in the inverse-Gaussian frailty term traps the M-step.
   - **N4:** the new left-truncation feature uses a likelihood that gives inconsistent estimates.

   None of them affects the gamma-frailty, right-censored, well-initialised fits that the paper reports. All four are local fixes.
4. **Name:** I recommend neither "conditional classification EM" nor plain "ICM". Call it a **CEM algorithm whose classification step is carried out by iterated conditional modes**, abbreviated e.g. **ICM-CEM** (§5).

---

## 1. Specification vs implementation

| # | `modified_cstep.pdf` requirement | Implementation | Status |
|---|---|---|---|
| R1 | Objective Eq. (1)–(4) = manuscript Eq. (7) | `.mixparfm_classification_loglik()` (`likelihood.R:331`) | **Verified**: equals an independent implementation to ≤3.4×10⁻¹³ at random states and ≤3.2×10⁻¹¹ at fitted states |
| R2 | One patient at a time, parameters fixed | `.mixparfm_conditional_cstep()` (`cstep.R:18`) | **Verified** |
| R3 | Score $Q_{ijg}$ from the "patient-removed" reference configuration, Eq. (6)–(8) | `cstep.R:50–74`. Cell sums $D_{jg}$ and $S_{jg}$ are updated in O(1) per move, fixing F-A11 of audit 06 | **Verified**: one sweep = brute-force sequential argmax of Eq. (7) in 32/32 cases |
| R4 | Each new label used immediately; sweeps repeated until a full pass makes no change | Sweeps repeat until no moves or `max_c_sweeps` (default 100) | **Verified**: all 32 test C-steps converged in ≤ 5 sweeps; endpoints are 1-opt |
| R5 | $\ell_c$ cannot decrease in the C-step | A move is accepted only if the gain is > `move_tolerance` (1e-12); ties keep the incumbent; a post-sweep check follows | **Verified** |
| R6 | Usual M-step plus a "numerical safeguard" | Closed-form τ, μ, Σ, π. The survival part is refitted per component with `nlminb`, warm-started, and a candidate is accepted only if it does not decrease that component's contribution (fixes F-A5, F-A6) | **Verified** exact for τ, μ, Σ, π (≤1.6×10⁻¹⁵). Survival refits sit at a local maximum, except in the N3 cases |
| R7 | Framed as ICM-type / safeguarded block-coordinate ascent | Code comments and docs call it "conditional classification EM (iterated conditional modes)" | See §5 |

**Deviations from the note.** All are defensible, and each needs one sentence in the paper:

- **(a) No empty components.** The argmax runs only over *admissible* labels: a patient cannot leave a component that is already at `minimum_component_size` (default 1). Monotonicity is unaffected, because the incumbent label is always admissible.
- **(b) Strict improvement.** A move is accepted only on strict improvement (> 1e-12), with a cap of 100 sweeps. This closes the termination gap noted in audit 06, F-A9.
- **(c) No `parfm` in the M-step.** The survival M-step uses the package's own optimiser on the same likelihood, not `parfm`.
- **(d) θ floor.** θ ≥ `theta_floor` (1e-4) is imposed, with a warning when it binds.
- **(e) Fewer model options.** SEM, the positive-stable and lognormal frailties, and strata were removed. This settles audit 06 F-A1, F-A2, F-A12, F-A13 and F-A15.
- **(f) Left truncation added.** This was not in the note, and see N4.

**Package tests:** every test I could run passes (71 expectations). The 3 tests that need `parfm` were skipped because `parfm` is not installed in my environment. The developer's summary reports them passing against the fork. My checks V1–V2 replace them with an independent reference.

---

## 2. What was verified (all reproducible)

| ID | Check | Result | Script |
|---|---|---|---|
| V1 | $F(d,s)=\log[(-1)^d\mathcal L^{(d)}(s)]$ for gamma and inverse Gaussian vs numerical integration, 200 grid points | max error 1.4×10⁻¹⁴; $F(0,0)=0$ for all families | `08a` |
| V2 | Package objective vs my implementation (written from the model, not from the code) | ≤3.4×10⁻¹³ | `08b` |
| V3 | One package sweep vs brute-force sequential $\arg\max_g \ell_c(c_{ij}=g, c_{-(ij)})$ | **32/32 identical partitions** (G = 2, 3; gamma, IG; Weibull, lognormal) | `08b` |
| V4 | Is the C-step endpoint a coordinate-wise (1-opt) maximum? | Yes in 32/32: the best single move changes $\ell_c$ by ≤ −3.9×10⁻⁴ | `08b` |
| V5 | Monotonicity over full fits: 120 fits (G = 2, 3; gamma, IG; 30%-noisy-truth and 2 random starts) plus 240 more in V9 | 0 decreases beyond rounding (min −9.1×10⁻¹³); 0 M-step rejections in the 120 fits | `08c`, `08f` |
| V6 | M-step block optimality at the final partition | τ, μ, Σ, π exact. Survival block: no improvement found from 5 perturbed restarts, except the N3 (IG) and N2 (degenerate) cases | `08c` |
| V7 | Fixed ψ, 12 patients, all 4,094 admissible partitions enumerated, 30 random ICM starts each | Global maximum reached by 100% of starts in 4/6 datasets; by 67% and 40% in the other two (worst gap 3.07) | `08g` |
| V8 | Application scale: N = 3,072, J = 32, G = 3 | 7.5 s, 16 iterations, monotone | — |

**V9 — one sweep per iteration vs sweeping to convergence** (`max_c_sweeps = 1` vs `100`). This re-runs audit 06 F-A3 with the real code: J = 20 × 20 patients, gamma frailty, 15 datasets per cell.

| G | separation | start | one sweep higher / equal / converged higher | median Δℓ (1 − 100) | ARI (1 sweep / 100) |
|---|---|---|---|---|---|
| 2 | 1.5 | 30% noisy truth | 1 / 13 / 1 | 0.00 | 0.688 / 0.684 |
| 2 | 1.5 | random | 6 / 3 / 6 | 0.00 | 0.266 / 0.239 |
| 2 | 3 | 30% noisy truth | 0 / 14 / 1 | 0.00 | 0.886 / 0.885 |
| 2 | 3 | random | 7 / 6 / 2 | 0.00 | 0.476 / 0.403 |
| 3 | 1.5 | 30% noisy truth | 8 / 6 / 1 | +0.09 | 0.531 / 0.529 |
| 3 | 1.5 | random | 8 / 0 / 3 | +10.35 | 0.298 / 0.230 |
| 3 | 3 | 30% noisy truth | 3 / 6 / 6 | 0.00 | 0.810 / 0.807 |
| 3 | 3 | random | 9 / 3 / 0 | +7.31 | 0.527 / 0.454 |

From good starts the two settings are practically the same. From random starts, one sweep is as good or better on the objective and better on ARI in every cell, at the same cost. **Suggestion:** make `max_c_sweeps = 1` the default, or at least report it as a sensitivity analysis. Both are monotone.

---

## 3. Findings

### N1 — NaN crash in the C-step when component parameters run away · **High** · Verified

**Evidence.** In V9, 5 of the 60 fits with G = 3 and a random start failed with `missing value where TRUE/FALSE needed`, raised at `cstep.R:77`. I reproduced both cases (`08h`). In one case component 1 had θ = 691 and per-patient risks up to 8.6×10⁵³; in the other, component 2 had θ = 842 and risks up to 1.2×10¹²⁵, some of them `Inf`. The incremental update `S[j, old] <- S[j, old] - risk[i]` (`cstep.R:51`) loses the smaller summands to cancellation. For example, (10⁷⁰ + 8.6×10⁵³) − 10⁷⁰ = 0 in double precision, so removing the next patient leaves S = −8.6×10⁵³. Then log(1 + θS) = NaN, `score[old]` = NaN, and `if (NaN >= …)` fails.

**Impact.**
- With `n_start = 1` the whole fit dies with an uninformative message.
- With `n_start > 1` the start is silently discarded. This is exactly what happens in a BIC sweep over G from random starts.
- There were no failures from the 30%-noisy starts.

**Fix.**
- (i) Treat a non-finite score as −Inf. If `score[old]` is not finite, recompute the cell sums exactly.
- (ii) Recompute `D` and `S` from scratch at the start of each sweep. This costs O(N), is negligible, and removes drift.
- (iii) Address the cause, N2.

### N2 — Degenerate components are reported as "converged" · **High** · Verified

**Evidence.** In 3 of 120 fits (all G = 3, all from random starts) a component ended in one of two degenerate states:

- **(a) Zero events.** 53 patients, all censored. The survival block's supremum is 0 and is not attained, so the parameters run to the boundary: θ = 7×10²⁰³, λ = 2×10⁻²⁰⁶, β = 206.
- **(b) Tiny component with a spuriously sharp fit.** For example 9 patients and 8 events: ρ = 568, β = 282, survival log-likelihood +2.87. A second case had 8 patients, 3 events and θ = 10⁴³.

Both returned `converged = TRUE` / `tolerance_reached`, with no warning. The θ-floor warning only fires at the lower bound. None of these fits won the best-objective selection in my runs.

**Why it matters.** This is the classification-likelihood counterpart of spurious maxima in mixtures. The objective has boundary solutions, and a monotone algorithm that heads towards one cannot turn back. It also feeds N1.

**Fix.**
- Require a minimum number of *events* per component (at least the number of survival parameters), not just a minimum size.
- Flag runaway parameters, for example |log θ| or |β| above a threshold.
- Exclude flagged starts from multi-start selection and from the BIC table, and report how many were excluded.
- In the paper, say that $\ell_c$ is maximised over non-degenerate solutions.

### N3 — Inverse-Gaussian frailty underflows for θ ≲ 1.43×10⁻³ and traps the M-step · **Medium-High** · Verified

**Evidence.**
- `besselK(z, k - 0.5)` underflows to 0 once z ≳ 700, that is θ ≲ 1/700. $F$ then becomes −∞ for every cell with k ≥ 1. This happened at 48 of 140 grid points (`08a`, `08d`).
- The optimiser's 1e10 penalty turns this into a cliff at θ ≈ 1.43×10⁻³. Warm-started `nlminb` sticks there.
- In 4 of 60 IG fits the M-step stopped at θ = 0.00143 while the true constrained maximum was θ = 0.0151 (loss 2.24) or θ at the 10⁻⁴ floor (loss 0.41).
- The floor warning did not fire, because 0.00143 > 10⁻⁴.
- Monotonicity is preserved (the acceptance rule), but the M-step is no longer a maximisation.
- The formula comes from `parfm::fr.ingau`, so the old code had the same problem.

**Fix (one line, tested).**
```r
log(besselK(z, k - 0.5, expon.scaled = TRUE)) - log(pi / (2 * z)) / 2 + base
# with base <- -2 * s / (1 + sqrt(1 + 2 * theta * s))   # = (1 - sqrt(1 + 2*theta*s)) / theta, cancellation-free
```
With this patch:
- The formula agrees with the original to 1.4×10⁻¹⁰ wherever the original is finite, and is finite on the whole grid.
- 0 of 60 fits get stuck.
- The previously stuck fits end 2.1 to 60.6 log-likelihood units higher.

**Relevance.** The frailty family is selected by BIC, so an under-maximised IG fit biases the gamma-vs-IG comparison.

### N4 — Left truncation uses an inconsistent likelihood · **Medium** · Verified

**Evidence.** For `Surv(t0, time, event)` the code uses risk $=(H(t)-H(t_0))e^{x'\beta}$ inside $F$, that is $F(d,\sum(H-H_0))$ (`likelihood.R:223`). The standard conditional likelihood for shared frailty conditions on the whole cluster surviving to entry: $F(d,\sum H)-F(0,\sum H_0)$. This is what `parfm` does (`Mloglikelihood.R`, `logSurv − logSurvT`); see Eriksson, Martinussen & Scheike, 2015, and van den Berg & Drepper, 2016.

Simulation (`08e`): gamma frailty, cluster-level selection at entry, J = 2,000 clusters × 4 patients, 8 replicates.

| | truth | implemented | conditional likelihood |
|---|---|---|---|
| λ | 0.200 | **0.094** (sd 0.007) | 0.213 |
| ρ | 1.50 | **1.43** (sd 0.04) | 1.51 |
| θ | 1.00 | 1.10 | 1.04 |

At the implemented estimate, the conditional log-likelihood is on average 172 units below its maximum.

**Impact.** Your paper does not use left truncation, but NEWS.md advertises it as a new feature.

**Fix.** Either `stop()` for counting-process responses, or implement it properly. A proper version tracks $S^0_{jg}$ per cell and adds $-[F(0,S^0+r^0_{ijg})-F(0,S^0)]$ to $Q_{ijg}$. The empirical-Bayes frailty predictions then also need the selected-frailty posterior.

### N5 — Tolerances · **Low-Medium** · Verified in code

- **M-step acceptance.** `candidate$loglik >= old_value - objective_tolerance` (`fit_mixparfm.R:720`) accepts decreases of up to 10⁻⁷ per component. The algorithm is then monotone only "up to 10⁻⁷". Use `>= old_value`.
- **C-step tolerance.** `move_tolerance = 1e-12` is absolute, on score differences of magnitude 10–100. That is at the rounding level of the incrementally maintained sums, so spurious "improvements" and cycling (cut off only by the 100-sweep cap) are possible. Something like `1e-9 * (1 + abs(score[old]))` would be safer.

### N6 — Visiting order · **Low** · Verified in code

Patients are visited in row order (`cstep.R:47`), so the endpoint depends on how the data are sorted. This is documented, and reproducible, which is good. For the paper, either say so or visit patients in a seeded random order and report the spread (audit 06, F-A10).

### N7 — A degenerate Gaussian covariance aborts the start · **Low** (by design) · Verified

4 of the same 60 G = 3 random-start fits in V9 stopped with `Degenerate Gaussian covariance`. Together with N1, 9 of those 60 fits (15%) failed. This is deliberate and documented. Still, with `n_start = 1` the fit fails, and the unbounded-likelihood problem remains (audit 06, F-A8). An eigenvalue-ratio constraint or `mclust::priorControl()` would turn failures into constrained solutions.

### Status of the audit 06 findings

- **Resolved by this implementation:** F-A1, F-A2, F-A5, F-A6, F-A7, F-A9, F-A11, F-A12 (by removal), F-A13, F-A14, F-A15.
- **Partially addressed:** F-A10 (order is documented, not randomised).
- **Still open by design:** F-A3 (see V9 above), F-A4 (a higher $\ell_c$ does not guarantee a higher ARI; seen again in V9) and F-A8 (unboundedness).

---

## 4. Is it a correct algorithm for maximising the classification log-likelihood?

**Yes, as a monotone block-coordinate ascent method with a local guarantee. It is not a global maximiser.** Precisely:

- **C-step.** Each accepted move raises $\ell_c(\psi^{(k)},c)$ by at least `move_tolerance` (V3–V5), and the set of admissible partitions is finite. So the sweep terminates at a 1-opt partition: no single reassignment improves $\ell_c$ at $\psi^{(k)}$.
- **M-step.** Given the partition, $\ell_c$ separates into independent blocks: τ; (μ_g, Σ_g); π_g; and each component's survival parameters. The first three are maximised exactly. The survival blocks are refitted from warm starts and only accepted if they do not decrease. So the M-step cannot decrease $\ell_c$.
- **Outcome.** $\ell_c(\psi^{(k)},c^{(k)})$ is non-decreasing. When it is bounded along the path it converges, to a *partial optimum*: the partition is 1-opt given ψ, and ψ is a (local) maximiser given the partition. This is the same kind of guarantee CEM gives, but weaker in one respect. In standard CEM the C-step maximises exactly over all partitions given ψ; here it does not (V7).
- **Caveats for the paper.** (i) $\ell_c$ is unbounded, and has degenerate boundary solutions (N2, N7), so "maximises" must mean over non-degenerate solutions. (ii) The endpoint depends on the starting partition and the visiting order, so multiple starts are required. (iii) A higher $\ell_c$ does not imply a better clustering (V9; audit 06, F-A4).

**Why not EM (useful for AE1).** $\ell_c$ couples the labels of patients in the same hospital, because $F_g(d_{jg},s_{jg})$ depends on all of them. The posterior $P(c_{ij}=g\mid y,\psi)$ would require summing over the $G^{n_j-1}$ labellings of the other patients in hospital $j$ (≈ 3⁹⁵ for $n_j≈96$). So neither an exact E-step nor an exact C-step is available. The old per-patient "posterior" was not the true posterior, and that is exactly the R1.1 objection. The full conditional $P(c_{ij}=g\mid c_{-(ij)},y,\psi)\propto\exp(Q_{ijg})$, on the other hand, is exact and cheap. That is the quantity ICM uses.

---

## 5. What to call it

**What the algorithm is, stated plainly.** It alternates two steps:

- **(i) Classification step.** Each label is set in turn to the mode of its full conditional distribution given all other labels and the current parameters, and is used immediately. This is iterated conditional modes (ICM; Besag, 1986) applied to $p(c\mid y,\psi)\propto\exp\ell_c(\psi,c)$.
- **(ii) M-step.** The parameters are updated given the partition.

Together this is block-coordinate ascent on the classification likelihood, the CML criterion that CEM maximises (Celeux & Govaert, 1992).

**"Conditional classification EM" (the current package name) — defensible, but I would not use it.**

- **For:** there is a real analogy with ECM (Meng & Rubin, 1993). An intractable joint step is replaced by a sequence of conditional maximisations that keeps the objective from decreasing. Here the joint C-step is replaced by one-label-at-a-time conditional C-steps.
- **Against:**
  - (i) It is not an established term. Referees will ask what "conditional" means, and the answer will be "ICM" anyway.
  - (ii) "Conditional" is ambiguous in a frailty paper: conditional likelihood? conditional on the frailty?
  - (iii) "Conditional EM", abbreviated CEM, is already the name of Jebara & Pentland's (1998) algorithm for maximising conditional likelihood. That collides with CEM = classification EM.
  - (iv) There is no E-step at all. That is also true of CEM, but it makes the name even less descriptive.

**"ICM" on its own — accurate for the classification step, too narrow for the algorithm.**

- **For:** the label update matches Besag's definition exactly (conditional mode, sequential, immediate use). There is also direct precedent in model-based clustering for exactly this structure, where integrating something out couples the labels. Côme & Latouche (2015) maximise the exact ICL of stochastic block models by moving one node at a time to the cluster with the largest increase, stop when a full pass brings no increase, and describe their algorithm as reminiscent of Besag's ICM.
- **Against:**
  - (i) ICM names only the label update at fixed parameters. Calling the whole estimator "ICM" hides the M-step, which does most of the statistical work: frailty, baseline and covariate effects.
  - (ii) Readers associate ICM with MAP image restoration under a Markov random field prior and will look for one. Here the dependence comes from the integrated shared frailty, not from a prior.
  - (iii) It breaks continuity with the CEM framing the paper and the referee exchange already use (AE1: why CEM/SEM rather than EM).

**Recommendation: "a CEM algorithm with an ICM classification step", shortened to ICM-CEM** (or CEM-ICM; either works if defined once). Every part of the name is an established, checkable term:

- **CEM:** what is maximised (the CML criterion) and the alternation of classification and M-steps.
- **ICM:** how the classification step is done, and why this is needed: the labels do not separate, so the standard MAP C-step is unavailable.

It also carries the right guarantee without over-claiming. Neither step is an exact maximisation: the C-step reaches a 1-opt partition, and the survival M-step is numerical with an acceptance rule. So the algorithm stands to CEM roughly as GEM stands to EM (Dempster, Laird & Rubin, 1977): monotone, with a local guarantee.

**Suggested wording for the Model Estimation section:**

> Because the integrated shared-frailty term $F_g(d_{jg},s_{jg})$ couples the labels of patients in the same hospital–cluster cell, the classification log-likelihood (7) does not separate across patients, and the usual C-step—assigning each unit by its marginal posterior probability—is not available. We therefore maximise (7) with a classification EM (CEM) algorithm (Celeux & Govaert, 1992) whose classification step is carried out by iterated conditional modes (ICM; Besag, 1986): each label is set in turn to the mode of its full conditional distribution given the current labels of all other patients and the current parameter estimates, and the new label is used immediately. We refer to the procedure as ICM-CEM. Every ICM update and every block of the M-step leaves (7) non-decreasing, so ICM-CEM is a monotone block-coordinate ascent algorithm that converges to a partial optimum: a partition that no single reassignment improves at the final parameters, together with parameters that locally maximise (7) given that partition. As with any CEM-type algorithm the solution is local, so we use multiple starts.

If you adopt this, align the package: `DESCRIPTION`, `NEWS.md`, the `control_EM()` title, the `fit_mixparfm()` details, and the returned `algorithm = "conditional_cem"`.

---

## 6. Manuscript items this implementation makes necessary (beyond the name)

- Replace the E-step (Eq. 8) and C-step (Eq. 9) with the ICM classification step. State the admissible set (no empty components), the strict-improvement rule and the visiting order.
- Remove the SEM material (Eq. 10 and the related text), or put back an S-step drawn from $\exp(Q_{ijg})$, a Gibbs sweep that reuses the same scores. The package no longer offers SEM.
- M-step paragraph: the survival blocks are no longer fitted by `parfm`. Mention the warm start, the acceptance rule and the θ floor.
- Frailty families: only gamma and inverse Gaussian remain. Remove any claim of "a wide range" of frailty distributions.
- Replace "CEM shares the theoretical guarantees of the standard EM … convergence to a stationary point" with the §4 statement.
- BIC: call it a classification (ICL-type) BIC, as already planned.

---

## 7. Recommended order of work

1. N3: the one-line IG fix. Then re-run any IG fits.
2. N1 + N2: NaN guard, recompute cell sums per sweep, minimum events per component, degeneracy flags in multi-start selection.
3. N4: disable left truncation, or implement the conditional likelihood.
4. N5: acceptance `>=`, relative move tolerance.
5. Decide `max_c_sweeps` (V9) and the visiting order (N6). Then freeze and re-run the simulation study and the application with k-prototypes + 20 starts.

## References

- Besag, J. (1986). On the statistical analysis of dirty pictures. *JRSS B* 48(3), 259–302.
- Celeux, G. & Govaert, G. (1992). A classification EM algorithm for clustering and two stochastic versions. *CSDA* 14, 315–332.
- Côme, E. & Latouche, P. (2015). Model selection and clustering in stochastic block models based on the exact integrated complete data likelihood. *Statistical Modelling* 15(6), 564–589.
- Dempster, A.P., Laird, N.M. & Rubin, D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *JRSS B* 39, 1–38.
- Jebara, T. & Pentland, A. (1998). Maximum conditional likelihood via bound maximization and the CEM algorithm. *NIPS* 11.
- Eriksson, F., Martinussen, T. & Scheike, T.H. (2015). Clustered survival data with left-truncation. *Scandinavian Journal of Statistics* 42(4), 1149–1166.
- Meng, X.-L. & Rubin, D.B. (1993). Maximum likelihood estimation via the ECM algorithm: a general framework. *Biometrika* 80(2), 267–278.
- van den Berg, G.J. & Drepper, B. (2016). Inference for shared-frailty survival models with left-truncated data. *Econometric Reviews* 35(6), 1075–1098.

*Scope limits.* All evidence is simulated: Weibull or lognormal baselines, gamma or IG frailty, J = 10–32, N = 120–3,072. The Lombardy data were not available. I ran R 4.3.3 with mclust 6.0.1 on Linux; the developer used R 4.6 on macOS. The `parfm`-equivalence tests were skipped; V1–V2 replace them with independent references.
