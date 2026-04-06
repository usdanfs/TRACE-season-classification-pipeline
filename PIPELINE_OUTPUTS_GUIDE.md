# Pipeline Outputs Interpretation Guide

This guide explains every diagnostic column and output file produced by the
season classification pipeline. It provides qualitative benchmarks, flags
combinations that warrant caution, and describes what to investigate when
results look unexpected.

---

## Quick-reference: result quality signals

Before reading individual outputs, check these five signals. If more than two
are unfavourable simultaneously, treat the winning season definition as
**provisional** and investigate before using it in downstream analysis.

| Signal | Favourable | Caution | Investigate |
|--------|-----------|---------|-------------|
| `p_top1` | ≥ 0.60 | 0.40 – 0.59 | < 0.40 |
| `bsa_min_std_quant` | ≥ 0.70 | 0.50 – 0.69 | < 0.50 |
| `bsa_min_ssa` | ≥ 0.60 | 0.40 – 0.59 | < 0.30 |
| `breakpoint_supported` | TRUE for ≥ 1 driver | — | FALSE for all drivers |
| `omega_sq` (winner) | ≥ 0.06 | 0.01 – 0.06 | < 0.01 |

*Note: `bsa_min_ssa` is only present in the full 4-stage pipeline. It is absent from the
climate-only (3-stage) decision table. See the climate-only section at the end of this guide.*

The pipeline also emits a consolidated console warning at runtime when
unfavourable signals co-occur. See the end of this guide for the interpretation
of that warning.

---

## decision_table_final.csv

The primary output. One row per retained candidate, sorted by `p_top1` then
`climate_score`.

### `climate_score`
Weighted composite rank-utility, range [0, 1]. Computed as a weighted average
of three tier scores (Tier 1: climate structure, Tier 2: internal robustness,
Tier 3: ecological verification). A value of 1.0 means the candidate ranked
first on all components; 0.0 means it ranked last on all.

**Interpretation:**
- The absolute value is less informative than the gap between candidates.
- A gap of < 0.05 between the top two candidates means the field is tight and
  the ranking is sensitive to data perturbations — check `p_top1` and
  `rank_IQR` before accepting the winner.
- A gap of ≥ 0.15 is a clear separation.
- When only 2–3 candidates remain after Stage 3, scores cluster near 0 and 1
  by construction; do not over-interpret absolute values in small candidate
  sets.

### `p_top1`
Proportion of bootstrap replicates (year-block resampling) where this
candidate ranked first. Range [0, 1].

**Benchmarks:**
- ≥ 0.60 — stable winner; the top rank is robust to inter-annual variability.
- 0.40 – 0.59 — moderate stability; consider whether the runner-up is
  scientifically distinguishable from the winner.
- < 0.40 — unstable; the ranking depends strongly on which years happen to be
  in the sample. Do not accept the winner without investigating whether the
  top candidates are ecologically equivalent.

**If the runner-up has comparable `p_top1`:** both definitions are plausible
given the data. Document both and make the choice on scientific grounds
(driver interpretability, literature precedent).

### `rank_IQR`
Inter-quartile range of bootstrap ranks. Smaller = more stable.

**Benchmarks:**
- IQR = 0 or 1 — essentially always the same rank; very stable.
- IQR = 2 — acceptable; one-rank fluctuations are common with realistic
  inter-annual variability.
- IQR ≥ 3 in a field of ≤ 6 candidates — the candidate has no clear rank
  position; the ranking is largely noise.

### `score_n_components`
Count of non-NA scoring components out of a maximum of 8 (full pipeline) or 6
(climate-only pipeline). Missing components reduce the effective weighting and
can change the winner.

**What each component requires:**
- Tier 1 (4 components: `u_mean_month`, `u_min_month`, `u_min_bin`,
  `u_switch`): always present if the candidate passed Stage 1 screening.
- Tier 2 (2 components: `u_sq_bsa`, `u_sq_ce`): requires both a `std` and a
  `quantile` candidate for the same driver × k combination to have survived
  Stage 3. If one method was dropped, these are NA for both survivors.
- Tier 3 (2 components: `u_s2_bsa`, `u_ssa_ce`): requires Stage 2 to have
  produced a valid breakpoint for the same driver × k. If
  `breakpoint_supported = FALSE` and the breakpoint estimate is non-finite,
  these are NA.

**Action:** If `score_n_components < 8`, check `ssa_reason` (for Tier 3) and
whether both threshold methods survived Stage 3 (for Tier 2). A winner scored
on 4 or fewer components should not be treated as definitive.

### `bsa_min_std_quant`
Conservative lower bound on the agreement between the `std` (literature
threshold) and `quantile` (data-driven threshold) season classifications over
the full climate record. Computed as `1 − upper Wilson CI of discord
proportion`. Range approximately [−0.05, 1].

**Interpretation:** This measures whether the result is sensitive to the
choice of threshold source. High agreement means both methods identify
essentially the same seasonal partition.

**Benchmarks:**
- ≥ 0.70 — robust; the classification is threshold-method invariant.
- 0.50 – 0.69 — acceptable; some sensitivity to threshold choice, but the
  overall pattern is consistent.
- < 0.50 — the classification changes substantially depending on how the
  threshold is set. Inspect `thresholds.csv` to understand the numerical
  difference between the two methods for this driver.

### `nce_std_quant`
Normalised conditional entropy (NCE) between the `std` and `quantile` season
label series. Range [0, 1]. Lower = more structurally similar.

- 0.00 — the two label series are identical.
- < 0.10 — near-identical structural partition; minor label swaps at
  boundary months only.
- 0.10 – 0.30 — moderate structural divergence; the two methods agree on
  the broad partition but differ in the placement of boundary months.
- > 0.30 — substantial divergence; the two methods are not producing the
  same classification of the climate record.

NCE is complementary to `bsa_min_std_quant`: BSA measures label agreement
directly; NCE measures whether the cut-points produce equivalent partitions
regardless of label identity.

### `bsa_min_ssa`
Conservative lower bound on agreement between Stage 1 (climate threshold)
and Stage 2 (ecological breakpoint) season classifications over the ecological
measurement window. Only present in the full 4-stage pipeline.

**Benchmarks:**
- ≥ 0.60 — strong cross-validation; the climate-derived and ecologically-
  derived boundaries agree. The season definition is unlikely to be an
  artefact of the threshold method.
- 0.40 – 0.59 — partial agreement. The broad seasonal pattern is shared but
  boundaries differ in some years.
- < 0.30 — investigate. The ecological breakpoint does not align with the
  climate threshold. This is not necessarily fatal (Stage 2 uses the response
  variable; if it is noisy, agreement will be low even for a correct climate-
  derived definition), but warrants investigation of `segmentation_results.csv`
  and the response variable quality.

The `poor_eco_agreement` quality flag fires at `bsa_min_ssa < 0.30`.

**If `stage2_near_constant = TRUE`:** the Stage 2 ecological labels are
dominated by one season (> `S4_NEAR_CONSTANT_THRESHOLD`). Agreement metrics
computed against a near-constant series are inflated and uninformative. This
flag does NOT indicate a problem with the Stage 1 classification; it indicates
a problem with Stage 2 validation for this driver × k combination. Possible
causes:
  1. The ecological record falls almost entirely within one climate regime.
  2. The Stage 2 breakpoint is near the edge of the driver range.
  3. The ecological record is too short to span multiple regimes.
  4. The response variable has low sensitivity to this driver.
  In all cases, treat `bsa_min_ssa` and `nce_ssa` for this candidate as
  unreliable and rely instead on Tier 1 and Tier 2 scores.

### `nce_ssa`
NCE between Stage 1 threshold-bins and Stage 2 breakpoint-bins on the driver
axis. Measures structural boundary alignment independent of label identity.
Interpretation scale is the same as `nce_std_quant` above.

When `ssa_reason != "ok"`, `nce_ssa` is NA. Common reasons:
- `no_stage2_breaks` — Stage 2 produced no finite breakpoint for this driver
  × k; check `segmentation_results.csv`.
- `degenerate_stage1_thresholds` / `degenerate_stage2_breaks` — one or both
  threshold sets are non-finite or inverted; check `thresholds.csv`.

---

## bootstrap_rank_summary.csv

### `decision_entropy`
Shannon entropy of the bootstrap winner distribution:
`−Σ p_i × log(p_i)` over all candidates that ever ranked first.

- **0.00** — one candidate always wins; maximum stability.
- **log(2) ≈ 0.69** — two candidates split the top rank 50/50.
- **log(n)** — all n candidates win equally often; the ranking is random.

**Benchmarks:**
- < 0.50 — one candidate dominates; result is stable.
- 0.50 – 1.00 — two candidates contend; review whether they are
  scientifically distinguishable.
- > 1.00 — three or more candidates contend; the ranking has no clear winner.
  Do not accept any single definition without additional evidence.

### `top_probability` and `runnerup_probability`
The gap `top_probability − runnerup_probability` is the most practical
stability indicator. A gap < 0.15 means the top two candidates are nearly
interchangeable under resampling.

---

## weight_sensitivity.csv

Shows the winning candidate under each combination of tier weights in the
sensitivity grid.

**Key question: does `top_candidate` vary across rows?**

- **No variation** — the winner is robust to the specific weight values chosen.
  This is the desired outcome.
- **Variation in ≤ 10% of combinations** — isolated sensitivity, typically at
  extreme weight values. Note which weight regions produce a different winner
  and whether those weights are scientifically defensible for your system.
- **Variation in > 25% of combinations** — the result is weight-sensitive. The
  winning candidate depends on how much you weight climate structure vs.
  ecological verification. This is a substantive finding: it means no single
  definition dominates across reasonable evidence weightings. Document all
  contending candidates.

`gap_to_second` is the score margin between first and second place. Small gaps
(< 0.05) at any weight combination indicate that the ranking is not robust at
those weights.

---

## segmentation_results.csv

### `breakpoint_supported`
TRUE if ΔAIC ≥ `MIN_DELTA_AIC` OR Davies test p < `DAVIES_ALPHA`.

**If FALSE for all drivers**, investigate in this order:

1. **Check record length.** Fewer than `MIN_MONTHS_FOR_SEG` overlapping months
   between `RESPONSE_CSV` and `CLIMATE_CSV` forces all fits to return null
   results. Check `overlap_n` in the Stage 2 console output.

2. **Check the response variable.** Plot `RESPONSE_COL` against each driver
   column. If the scatter shows no monotone or piecewise-linear relationship,
   Stage 2 will not detect a breakpoint regardless of sample size. This is
   informative: it means the selected drivers do not segment the ecological
   response, not that the pipeline failed.

3. **Check for high inter-annual variability.** If the response has strong
   year-to-year variability unrelated to the drivers, the within-season signal
   is obscured. Consider whether a transformation (log, square root) or
   aggregation (seasonal means instead of monthly) would reduce noise.

4. **Check the initial breakpoint guesses** (`psi1`, `psi2` in `SEG_DRIVERS`).
   If these are far from the true breakpoint, `segmented()` may converge to a
   degenerate solution. Set them to plausible values based on domain knowledge
   or the marginal driver distribution.

**Note:** `breakpoint_supported = FALSE` does NOT prevent Stage 3 and Stage 4
from running. The ecological candidate labels are still generated (using the
best-estimate breakpoint even if unsupported) and used in cross-validation.
It means the ecological evidence for a breakpoint is weak, which reduces
confidence in Tier 3 scores but does not invalidate the climate-derived
classification.

### `rmse_cv`
Leave-one-year-out cross-validation RMSE. Lower = better predictive accuracy
of the segmented model on held-out years. There is no universal benchmark
(depends on the response variable's scale and variance), but use it
comparatively: if `rmse_cv` for k=3 is not substantially lower than for k=2,
the additional breakpoint does not improve prediction.

### Bootstrap breakpoint CI (`b1_lo`, `b1_hi`, `b2_lo`, `b2_hi`)
2.5th and 97.5th percentile of the breakpoint estimate across B = 300
case-resampling bootstrap replicates. Quantifies uncertainty in the breakpoint
location without relying on asymptotic normality (which is violated near
boundary values of the driver range).

**Interpretation:**
- Narrow CI (width < 0.5 × IQR of the driver) — the breakpoint location is
  well-determined by the data.
- Wide CI spanning a large fraction of the driver range — the breakpoint is
  poorly constrained; the segmented model fits many plausible locations equally
  well. This reduces confidence in the threshold–breakpoint alignment diagnostic
  (`diff1_iqr`/`diff2_iqr`) even if the breakpoint is nominally supported.
- `NA` — fewer than `BOOT_CI_MIN_REPS` (default 10) bootstrap replicates
  produced a valid breakpoint estimate; CI is not reported.

---

## filter_results.csv

Drop flags: `fail_assignment`, `fail_imbalance`, `fail_calendar`,
`fail_block_collapse`, `fail_block_imbalance`.

### If many candidates fail the same screen:

**`fail_assignment`** (< `S3_MIN_PCT_ASSIGNED` months labelled in ecological
window): The ecological window has many NA-valued months. Check that
`RESPONSE_CSV` dates overlap the driver record and that the driver column has
finite values in the response period.

**`fail_imbalance`** (smallest season bin < `S3_MIN_SEASON_PROP` of months):
The threshold places almost all ecological-window months in one season. For
`std` candidates: the literature threshold may be far outside the range of
your site's driver values — inspect `thresholds.csv` and compare to the actual
driver distribution. For `quantile` candidates: the baseline period may be
climatologically unrepresentative of the ecological window.

Note: Stage-1 also applies an **absolute count screen** (`min_bin_n`) before
candidates reach Stage-3. Defaults: < **24 months** for k=2; < **18 months**
for k=3 (config params `S1_MIN_BIN_N_2S` / `S1_MIN_BIN_N_3S`). Candidates
that pass the Stage-1 absolute screen but fail the Stage-3 proportional screen
(`fail_imbalance`) typically have adequate absolute counts but unbalanced
partitions within the shorter ecological window.

**`fail_calendar`** (mean monthly consistency < `S3_MEAN_MONTH_CONS`): Season
labels flip between years for the same calendar month, meaning the threshold
does not align with the annual cycle at this site. This is most common for
transitional months. It is not always a problem with the driver — it may
indicate that the chosen driver does not have a stable seasonal cycle at this
site.

**`fail_block_collapse`** / **`fail_block_imbalance`**: At least one
multi-year block within the ecological window is dominated by one season or
loses a season level entirely. This can occur when the ecological window is
short (< 8 years) and one season is genuinely absent in some years. Consider
whether `S3_BLOCK_YEARS` should be increased (larger blocks are more tolerant
of single-year anomalies) or whether `S3_MIN_HEALTHY_BLOCKS` should be relaxed
for short records.

**If ALL candidates fail all screens:** The pipeline will stop with an error.
The most productive diagnostic path:
1. Reduce `S3_MIN_PCT_ASSIGNED` if assignment is the issue.
2. Reduce `S3_MIN_SEASON_PROP` if imbalance is the issue.
3. Check that the ecological window (response record dates) overlaps the
   climate record for at least 5 years.

---

## Stage-3 diagnostic flags (anova_summary.csv / retained_candidates.csv)

The following columns are written to `anova_summary.csv` and carried forward into
`retained_candidates.csv`. They are informational only — none triggers a drop.

### `flag_kappa_low`
TRUE when Cohen's κ (Stage-1 threshold labels vs. Stage-2 breakpoint labels) is below
`S3_FLAG_KAPPA_LOW` (default **0.10**). This corresponds to the "slight agreement"
range on the Landis & Koch scale and indicates that climate-derived and ecologically-
derived labels have negligible correspondence beyond chance.

`flag_kappa_low = TRUE` does not eliminate a candidate, but combined with
`bsa_min_ssa < 0.30` and `flag_align_far = TRUE` it is a strong signal that the
climate threshold and ecological breakpoint describe different boundaries on the
driver axis. The Stage-4 `kappa_ssa` column carries this same κ into the final
decision table for cross-reference.

### `flag_align_far`
TRUE when the IQR-normalised distance between the Stage-1 threshold and the Stage-2
breakpoint exceeds `S3_FLAG_ALIGN_IQR` (default **1.5 IQR**). This is the standard
Tukey outer fence applied to positional discrepancy: a gap larger than 1.5 × IQR of
the driver distribution is flagged as unusually large.

`diff1_iqr` (t1 vs. b1) and `diff2_iqr` (t2 vs. b2 for k=3) store the raw values.
`flag_align_far` is TRUE if either exceeds the threshold. NA values (from zero-IQR
constant drivers) are treated as non-flagged.

### `kw_p`
Kruskal-Wallis test p-value for differences in the ecological response across season
levels. A distribution-free complement to the parametric ANOVA, robust to non-normality
and unequal variances. Reported alongside `omega_sq` for each candidate.

- p < 0.05 — evidence of a rank-based difference across seasons.
- p ≥ 0.05 with low `omega_sq` — the ecological response does not differentiate seasons
  either parametrically or non-parametrically; interpret with caution.

`kw_p` is DIAGNOSTIC only and does not enter scoring or trigger flags.

### Tukey HSD pairwise comparisons (`posthoc_tbl` / `anova_summary.csv`)
For k=3 candidates where TukeyHSD can be computed, pairwise season comparisons are
written with columns `comparison`, `diff`, `lwr`, `upr`, `p_adj`. These show the
direction and magnitude of between-season differences in the ecological response.

- `diff` — mean response difference (season A minus season B).
- `lwr` / `upr` — 95% simultaneous confidence interval.
- `p_adj` — Tukey-adjusted p-value.

Use these to confirm that all three seasons are ecologically distinct (all pairwise
`p_adj < 0.05`), or to identify which pairs drive the omnibus ANOVA signal. Not
computed for k=2 (only one pair; use `anova_p` directly).

---

## anova_summary.csv

### `omega_sq`
Omega-squared: proportion of variance in the ecological response explained by
the season classification. Range [0, 1].

**Context-dependent benchmarks:**
- ≥ 0.14 — large effect (Cohen's convention); season classification explains
  substantial ecological variance.
- 0.06 – 0.13 — medium effect; moderate ecological differentiation.
- 0.01 – 0.05 — small effect; season classification explains some ecological
  variance but the signal is weak. This is common in systems with high
  inter-annual variability or where the chosen driver is only one of several
  controls on the response.
- < 0.01 — negligible effect; the season classification does not differentiate
  the ecological response at the monthly level.

**Important:** A low `omega_sq` does NOT necessarily mean the season
classification is wrong. In weakly seasonal or highly variable ecosystems,
monthly ecological data may show low between-season differentiation even when
the classification is climatologically correct. Use `omega_sq` as a
corroborating signal, not a pass/fail criterion. The pipeline flags values
below `S3_FLAG_OMEGA_SQ_LOW` (default 0.01) for attention, not for removal.

**If `omega_sq` is low for all candidates:** This may be a property of the
ecosystem (weak coupling between climate seasonality and the measured response
at monthly timescales) rather than a failure of the classification. Consider:
- Whether the response variable is measured at the right temporal scale.
- Whether a lag between climate forcing and ecological response should be
  modelled.
- Whether the chosen response variable is the most sensitive indicator of
  seasonal regime at this site.

---

## Overall result quality: when to proceed with caution

The pipeline always produces a winner. A combination of the following signals
indicates that the result should be treated as preliminary:

| Combination | Interpretation | Recommended action |
|-------------|---------------|-------------------|
| `p_top1` < 0.40 AND `rank_IQR` ≥ 3 | Ranking is essentially random under resampling | Run with more bootstrap replicates; if instability persists, the candidate field is not sufficiently differentiated |
| `bsa_min_ssa` < 0.30 AND `breakpoint_supported = FALSE` | Neither ecological label agreement nor a detectable breakpoint supports the climate-derived boundary | The classification may be climatologically coherent but lacks ecological cross-validation; treat as exploratory |
| All `breakpoint_supported = FALSE` AND `omega_sq` < 0.01 | No ecological evidence for any season boundary in any driver | Investigate whether the response variable, record length, or driver selection are appropriate |
| `stage2_near_constant = TRUE` for winner's driver × k | Ecological cross-validation is unreliable | Rely on Tier 1 and Tier 2 scores only; cross-validation metrics are uninformative |
| `score_n_components` < 6 | Multiple scoring components are missing | Identify which components are NA and why; a winner scored on 4 or fewer components lacks a full evidence base |
| `decision_entropy` > 1.0 | Three or more candidates contend for top rank | Report the top two or three candidates as a set; the data do not support a single winner |

The pipeline emits a consolidated runtime warning when unfavourable signals
co-occur. This warning appears at the end of Stage 4 (or Stage 3 for the
climate-only pipeline) console output and summarises which conditions were
triggered.

---

## Additional columns in decision_table_final.csv

### `kappa_std_quant` and `kappa_ssa`
Cohen's κ between the two label series (std vs. quantile, or Stage 1 vs. Stage 2).
κ corrects raw agreement for the proportion expected by chance.

- κ > 0.80 — near-perfect agreement.
- κ 0.60 – 0.80 — substantial agreement.
- κ 0.40 – 0.59 — moderate agreement.
- κ < 0.40 — fair or poor agreement.

κ and BSA_min are complementary: BSA_min gives a statistically conservative bound on
agreement; κ gives a chance-corrected measure. When they disagree, the BSA_min is more
conservative because it is a lower bound, not a point estimate.

### `stage2_pmax`
Highest single-class proportion in the Stage 2 (ecological breakpoint) season labels
over the ecological window. A value close to 1.0 means nearly all ecological-window
months were assigned to one season by the Stage 2 breakpoint.

When `stage2_pmax > S4_NEAR_CONSTANT_THRESHOLD` (default 0.95), `stage2_near_constant`
is set to TRUE and label-agreement metrics against Stage 2 become uninformative.

### `entropy_norm`
Normalised Shannon entropy of the season class distribution over the full climate record.
Range [0, 1]. A value of 1.0 means all k seasons are perfectly equally represented.
A value near 0 means the distribution is highly skewed.

This is a supplementary balance diagnostic alongside `min_bin_prop`. Unlike `min_bin_prop`
(which captures only the smallest bin), `entropy_norm` summarises the overall spread.
No hard threshold is applied — use it to compare candidate balance qualitatively.

### `med_run`
**File:** `output_STAGE_1/tables/validation_tbl.csv`

Median length (in consecutive months) of same-season runs in the full climate record.
Higher values mean longer unbroken seasonal stretches, which is biologically more
plausible. Very short median runs (≤ 2 months) indicate the threshold lies in a region
of high driver density where labels switch rapidly from month to month.

`med_run` is a Stage-1 DIAGNOSTIC and does not enter screening or scoring. Use it
alongside `mean_switch_per_year` to assess temporal coherence: high switching rate and
low median run length together indicate a fragmented, threshold-sensitive candidate.

---

## block_stability.csv

One row per candidate × year block. Shows whether each season level was represented
in every temporal block used for the structural stress test.

| Column | Meaning |
|--------|---------|
| `test_years` | Year range of this block (e.g. "2019-2020") |
| `n_block` | Total months in this block |
| `n_assigned_block` | Months with a non-NA season label |
| `n_levels_block` | Number of distinct season levels present |
| `min_bin_prop_block` | Smallest bin proportion within this block |
| `min_bin_n_block` | Count of months in the smallest bin |

A block where `n_levels_block < k` means one season level is entirely absent from
that period — this contributes to the `fail_block_collapse` drop criterion.

---

## retained_candidates.csv additional columns

### `prop_healthy`
Proportion of year blocks where all k season levels were represented (none collapsed).
`prop_healthy = 1.0` means every block contained all seasons; lower values indicate
periods where a season level disappeared.

### `n_collapsed`
Number of year blocks where at least one season level was absent. High values for a
short record (e.g., n_collapsed = 2 out of 3 blocks) indicate the season definition
is fragile over time.

---

## Climate-only (3-stage) pipeline differences

The climate-only pipeline produces the same `decision_table_final.csv` column layout
as the full pipeline **except** that the following columns are absent because there is
no Stage 2 ecological segmentation:

- `stage2_n`, `bsa_min_ssa`, `nce_ssa`, `kappa_ssa`
- `stage2_pmax`, `stage2_near_constant`, `ssa_reason`

The maximum `score_n_components` for the climate-only pipeline is **6** (4 Tier-1 + 2 Tier-2),
not 8. The `incomplete_scoring` quality flag uses 6 as the reference maximum.

The weight sensitivity grid in the climate-only pipeline sweeps over `w_climate` and
`w_robust` only (no `w_verify`), so `weight_sensitivity.csv` has two weight columns
instead of three.

There is no `filter_results.csv` with ecological-window columns; the climate-only
Stage 2 produces a structurally identical file but with `_val` suffixes on the
consistency and balance columns, reflecting the validation window rather than an
ecological measurement window.
