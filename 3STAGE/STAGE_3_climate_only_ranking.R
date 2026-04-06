# =============================================================================
# STAGE 3 — Climate-Only Decision Ranking and Bootstrap Stability
# =============================================================================
# Rank retained candidates using climate structure and internal robustness only.
# This is the climate-only analogue of the full-pipeline Stage 4.
#
# Inputs:
#   Stage 1: screened_tbl.rds, season_long.rds, threshold_tbl.rds
#   Stage 2: stage2_stage1_candidates_retained.rds
#
# Outputs (output_dir/tables/):
#   - decision_table_final.csv
#   - bootstrap_rank_summary.csv
#   - weight_sensitivity.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)
})

# Default path uses the 3STAGE/ prefix so the script works when run from the
# project root directory. When run from inside 3STAGE/, set SEASON_CONFIG
# explicitly or the file will still be found via the relative path.
CONFIG_FILE <- Sys.getenv("SEASON_CONFIG", unset = "3STAGE/config_climate_only.R")
source(CONFIG_FILE)
set.seed(GLOBAL_SEED)

output_dir <- stage_dir(3)
tab_dir    <- file.path(output_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

screened_tbl  <- readRDS(file.path(stage_dir(1), "screened_tbl.rds"))
season_long   <- readRDS(file.path(stage_dir(1), "season_long.rds"))
threshold_tbl <- readRDS(file.path(stage_dir(1), "threshold_tbl.rds"))
retained      <- readRDS(file.path(stage_dir(2), "stage2_stage1_candidates_retained.rds"))

base_tbl <- screened_tbl %>%
  semi_join(retained %>% dplyr::select(candidate_id), by = "candidate_id") %>%
  mutate(k            = as.integer(k),
         candidate_id = as.character(candidate_id),
         driver       = as.character(driver),
         method       = as.character(method))

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

kappa_safe <- function(tab) {
  tab <- as.matrix(tab)
  n <- sum(tab)
  if (!is.finite(n) || n <= 1) return(NA_real_)
  all_lvls <- union(rownames(tab), colnames(tab))
  sq <- matrix(0, length(all_lvls), length(all_lvls),
               dimnames = list(all_lvls, all_lvls))
  sq[rownames(tab), colnames(tab)] <- tab
  rs <- rowSums(sq); cs <- colSums(sq)
  if (sum(rs > 0) <= 1 || sum(cs > 0) <= 1) return(NA_real_)
  po <- sum(diag(sq)) / n
  pe <- sum(rs * cs) / n^2
  if (!is.finite(pe) || (1 - pe) <= 0) return(NA_real_)
  (po - pe) / (1 - pe)
}

WILSON_Z <- 1.96  # z for 95% two-sided Wilson CI (qnorm(0.975))

wilson_ci <- function(m, n, z = WILSON_Z) {
  if (!is.finite(n) || n <= 0) return(c(lo = NA_real_, hi = NA_real_))
  p <- m / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half   <- (z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)) / denom
  c(lo = center - half, hi = center + half)
}

label_agreement_bsa <- function(a, b, w = NULL) {
  ok <- !is.na(a) & !is.na(b)
  a <- as.character(a[ok]); b <- as.character(b[ok])
  if (!is.null(w)) w <- as.numeric(w[ok])
  if (length(a) <= 1)
    return(tibble(n = length(a), pct_agree = NA_real_, discord = NA_real_,
                  discord_ci_hi = NA_real_, bsa_min = NA_real_, kappa = NA_real_))
  tab <- if (is.null(w)) table(a, b)
         else xtabs(w ~ a + b, data = data.frame(a = a, b = b, w = w))
  n <- sum(tab); agree <- sum(diag(tab)); m <- n - agree
  ci <- wilson_ci(m, n)
  tibble(n = n, pct_agree = agree / n, discord = m / n,
         discord_ci_hi = ci["hi"], bsa_min = 1 - ci["hi"],
         kappa = kappa_safe(tab))
}

rank01 <- function(x, higher_better = TRUE) {
  if (!higher_better) x <- -x
  r <- rank(x, ties.method = "average", na.last = "keep")
  n_ok <- sum(!is.na(r))
  if (n_ok <= 1) return(rep(NA_real_, length(x)))
  (r - 1) / (n_ok - 1)
}

info_metrics_from_tab <- function(tab) {
  tab <- as.matrix(tab); n <- sum(tab)
  if (!is.finite(n) || n <= 1)
    return(tibble(MI = NA_real_,
                  H1_given_2 = NA_real_,
                  H2_given_1 = NA_real_,
                  nH1_given_2 = NA_real_,
                  nH2_given_1 = NA_real_))
  Pij <- tab / n; Pi <- rowSums(Pij); Pj <- colSums(Pij)
  H <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  Hi <- H(Pi); Hj <- H(Pj)
  MI <- 0
  for (i in seq_len(nrow(Pij)))
    for (j in seq_len(ncol(Pij)))
      if (Pij[i, j] > 0) MI <- MI + Pij[i, j] * log(Pij[i, j] / (Pi[i] * Pj[j]))
  H1g2 <- Hi - MI
  H2g1 <- Hj - MI
  tibble(MI = MI,
         H1_given_2 = H1g2,
         H2_given_1 = H2g1,
         nH1_given_2 = if (Hi > 0) H1g2 / Hi else NA_real_,
         nH2_given_1 = if (Hj > 0) H2g1 / Hj else NA_real_)
}

safe_tier_mean <- function(...) {
  v <- mean(c(...), na.rm = TRUE)
  if (is.nan(v)) NA_real_ else v
}

weighted_score <- function(tiers, wts) {
  ok <- !is.na(tiers)
  if (sum(ok) == 0) return(NA_real_)
  sum(tiers[ok] * wts[ok]) / sum(wts[ok])
}

std_quant_agreement <- function(sl) {
  wt_tbl <- if ("w" %in% names(sl)) {
    sl %>% distinct(DateMonth, w)
  } else {
    sl %>% distinct(DateMonth) %>% mutate(w = 1)
  }
  meta <- sl %>%
    distinct(candidate_id, driver, k, method) %>%
    filter(method %in% c("std", "quantile"))
  pairs <- meta %>%
    group_by(driver, k) %>%
    summarise(
      cid_std = candidate_id[method == "std"][1],
      cid_qtl = candidate_id[method == "quantile"][1],
      .groups = "drop"
    ) %>%
    filter(!is.na(cid_std), !is.na(cid_qtl))
  pairs %>%
    mutate(metrics = map2(cid_std, cid_qtl, function(c1, c2) {
      s1 <- sl %>% filter(candidate_id == c1) %>% dplyr::select(DateMonth, s_std = season)
      s2 <- sl %>% filter(candidate_id == c2) %>% dplyr::select(DateMonth, s_qtl = season)
      joined <- inner_join(s1, s2, by = "DateMonth") %>% left_join(wt_tbl, by = "DateMonth")
      met <- label_agreement_bsa(joined$s_std, joined$s_qtl, w = joined$w)
      lv <- union(levels(factor(joined$s_std)), levels(factor(joined$s_qtl)))
      tab <- xtabs(w ~ factor(s_std, levels = lv) + factor(s_qtl, levels = lv), data = joined)
      info <- info_metrics_from_tab(tab)
      tibble(
        stage1_n = sum(joined$w, na.rm = TRUE),
        kappa_std_quant = met$kappa,
        bsa_min_std_quant = met$bsa_min,
        pct_agree_std_quant = met$pct_agree,
        discord_ci_hi_std_quant = met$discord_ci_hi,
        nce_std_quant = mean(c(info$nH1_given_2, info$nH2_given_1), na.rm = TRUE)
      )
    })) %>%
    unnest(metrics)
}

# =============================================================================
# 3. INTERNAL ROBUSTNESS
# =============================================================================

sq_tbl <- std_quant_agreement(season_long)

# =============================================================================
# 4. ASSEMBLE DECISION SET AND SCORE
# =============================================================================

decision_set <- base_tbl %>%
  left_join(sq_tbl, by = c("driver", "k"))

decision_set <- decision_set %>%
  mutate(
    u_mean_month = rank01(mean_month_consistency, TRUE),
    u_min_month  = rank01(min_month_consistency, TRUE),
    u_min_bin    = rank01(min_bin_prop, TRUE),
    u_switch     = rank01(mean_switch_per_year, FALSE),
    u_sq_bsa     = rank01(bsa_min_std_quant, TRUE),
    u_sq_ce      = rank01(nce_std_quant, FALSE)
  ) %>%
  rowwise() %>%
  mutate(
    tier_climate = safe_tier_mean(u_mean_month, u_min_month, u_min_bin, u_switch),
    tier_robust  = safe_tier_mean(u_sq_bsa, u_sq_ce),
    climate_score = weighted_score(
      c(tier_climate, tier_robust),
      c(W_CLIMATE, W_ROBUST)
    ),
    score_n_components = sum(!is.na(c(
      u_mean_month, u_min_month, u_min_bin, u_switch,
      u_sq_bsa, u_sq_ce
    )))
  ) %>%
  ungroup()

decision_table <- decision_set %>%
  arrange(desc(climate_score), desc(bsa_min_std_quant)) %>%
  dplyr::select(
    candidate_id, driver, k, method,
    climate_score, score_n_components,
    mean_month_consistency, min_month_consistency,
    min_bin_prop, entropy_norm, mean_switch_per_year,
    stage1_n, bsa_min_std_quant, kappa_std_quant, nce_std_quant
  )

# =============================================================================
# 5. BOOTSTRAP RANK STABILITY
# =============================================================================

years_full <- sort(unique(year(season_long$DateMonth)))

resample_years <- function(years, B = BOOT_N_RANK) {
  replicate(B, sample(years, length(years), replace = TRUE), simplify = FALSE)
}

boot_sets <- resample_years(years_full, B = BOOT_N_RANK)

boot_ranks <- map_dfr(seq_along(boot_sets), function(bi) {
  yb <- boot_sets[[bi]]
  yfreq <- tibble(Year = yb) %>% count(Year, name = "w")
  sl_b <- season_long %>%
    mutate(Year = year(DateMonth)) %>%
    inner_join(yfreq, by = "Year")
  sq_b <- std_quant_agreement(sl_b)
  ds_b <- base_tbl %>%
    left_join(sq_b, by = c("driver", "k")) %>%
    mutate(
      u_mean_month = rank01(mean_month_consistency, TRUE),
      u_min_month  = rank01(min_month_consistency, TRUE),
      u_min_bin    = rank01(min_bin_prop, TRUE),
      u_switch     = rank01(mean_switch_per_year, FALSE),
      u_sq_bsa     = rank01(bsa_min_std_quant, TRUE),
      u_sq_ce      = rank01(nce_std_quant, FALSE)
    ) %>%
    rowwise() %>%
    mutate(
      tier_climate = safe_tier_mean(u_mean_month, u_min_month, u_min_bin, u_switch),
      tier_robust  = safe_tier_mean(u_sq_bsa, u_sq_ce),
      climate_score = weighted_score(c(tier_climate, tier_robust), c(W_CLIMATE, W_ROBUST))
    ) %>%
    ungroup() %>%
    arrange(desc(climate_score), desc(bsa_min_std_quant)) %>%
    mutate(rank = row_number(), boot = bi)
  ds_b %>% dplyr::select(boot, candidate_id, rank, climate_score)
})

rank_stats <- boot_ranks %>%
  group_by(candidate_id) %>%
  summarise(p_top1 = mean(rank == 1), rank_IQR = IQR(rank),
            .groups = "drop") %>%
  arrange(desc(p_top1))

top_probs <- boot_ranks %>%
  filter(rank == 1) %>%
  count(candidate_id) %>%
  mutate(prob = n / BOOT_N_RANK)

boot_summary <- tibble(
  N_BOOT               = BOOT_N_RANK,
  top_candidate        = top_probs$candidate_id[1],
  top_probability      = top_probs$prob[1],
  runnerup_probability = if (nrow(top_probs) >= 2) top_probs$prob[2] else NA_real_,
  decision_entropy     = -sum(top_probs$prob * log(top_probs$prob))
)

decision_table_final <- decision_table %>%
  left_join(rank_stats, by = "candidate_id") %>%
  arrange(desc(p_top1), desc(climate_score))

# =============================================================================
# 6. WEIGHT SENSITIVITY
# =============================================================================

weight_grid <- expand_grid(
  w_climate = seq(SENS_W_CLIMATE_RANGE[1], SENS_W_CLIMATE_RANGE[2], by = SENS_W_STEP),
  w_robust  = seq(SENS_W_ROBUST_RANGE[1],  SENS_W_ROBUST_RANGE[2],  by = SENS_W_STEP)
) %>%
  filter(abs(w_climate + w_robust - 1) < 1e-6)  # 1e-6 tolerance avoids silent

weight_sensitivity <- weight_grid %>%
  mutate(res = pmap(list(w_climate, w_robust), function(wc, wr) {
    d <- decision_set %>%
      rowwise() %>% mutate(score = weighted_score(c(tier_climate, tier_robust), c(wc, wr))) %>% ungroup()
    ord <- d %>% arrange(desc(score), desc(bsa_min_std_quant))
    tibble(top_candidate = ord$candidate_id[1],
           top_score = ord$score[1],
           gap_to_second = if (nrow(ord) > 1) ord$score[1] - ord$score[2] else NA_real_)
  })) %>%
  unnest(res)

# =============================================================================
# 7. SAVE OUTPUTS
# =============================================================================

write.csv(decision_table_final %>%
            dplyr::select(
              candidate_id, driver, n_seasons = k, method,
              climate_score, p_top1, rank_IQR, score_n_components,
              mean_month_consistency, min_month_consistency,
              min_bin_prop, mean_switch_per_year, entropy_norm,
              stage1_n, bsa_min_std_quant, nce_std_quant, kappa_std_quant),
          file.path(tab_dir, "decision_table_final.csv"), row.names = FALSE)

write.csv(weight_sensitivity,
          file.path(tab_dir, "weight_sensitivity.csv"), row.names = FALSE)

write.csv(boot_summary,
          file.path(tab_dir, "bootstrap_rank_summary.csv"), row.names = FALSE)

saveRDS(decision_set, file.path(output_dir, "decision_set.rds"))
saveRDS(boot_ranks,   file.path(output_dir, "boot_ranks.rds"))

# =============================================================================
# 8. RESULT QUALITY SYNTHESIS — consolidated runtime warning
# =============================================================================
# See PIPELINE_OUTPUTS_GUIDE.md for full interpretation of each signal.

n_candidates  <- nrow(decision_table_final)
winner_row    <- decision_table_final %>% slice(1)
quality_flags <- character(0)

# --- Single-candidate short-circuit ---
if (n_candidates == 1L) {
  message(
    "\nResult quality: SINGLE CANDIDATE — only one candidate survived all screening stages. ",
    "Bootstrap rank stability (p_top1, rank_IQR) and scoring metrics are uninformative with ",
    "a single candidate. Evaluate the season definition on scientific grounds. ",
    "Review filter_results.csv to understand why all other candidates were dropped.")
} else {

if (is.finite(winner_row$p_top1) && winner_row$p_top1 < 0.50)
  quality_flags["rank_unstable"] <- sprintf(
    "Bootstrap rank unstable: winner held rank 1 in %.0f%% of replicates (recommend >= 50%%).",
    winner_row$p_top1 * 100)

max_components <- 6L  # 4 Tier-1 + 2 Tier-2; no Tier-3 in climate-only pipeline
if (is.finite(winner_row$score_n_components) &&
    winner_row$score_n_components < max_components)
  quality_flags["incomplete_scoring"] <- sprintf(
    "Winner scored on only %d of %d components. Check whether both threshold methods survived Stage 2.",
    as.integer(winner_row$score_n_components), max_components)

if (is.finite(winner_row$bsa_min_std_quant) && winner_row$bsa_min_std_quant < 0.50)
  quality_flags["method_sensitive"] <- sprintf(
    "Low std/quantile agreement: BSA_min = %.2f (recommend >= 0.50). Classification depends strongly on threshold source.",
    winner_row$bsa_min_std_quant)

n_weight_winners <- n_distinct(weight_sensitivity$top_candidate)
if (n_weight_winners > 1) {
  pct_not_top <- mean(weight_sensitivity$top_candidate != winner_row$candidate_id) * 100
  if (pct_not_top > SENS_W_WINNER_CHANGE_PCT)
    quality_flags["weight_sensitive"] <- sprintf(
      "Winner changes in %.0f%% of weight combinations. Review weight_sensitivity.csv.",
      pct_not_top)
}

if (length(quality_flags) == 0) {
  message("\nResult quality: ACCEPTABLE — no major concerns flagged.")
} else {
  message(sprintf(
    "\nResult quality: CAUTION — %d concern(s) flagged for winner '%s':",
    length(quality_flags), winner_row$candidate_id))
  for (nm in names(quality_flags))
    message("  [", nm, "] ", quality_flags[[nm]])
  message("\nSee PIPELINE_OUTPUTS_GUIDE.md for interpretation guidance.")
}
} # end single-candidate else block

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
# =============================================================================
