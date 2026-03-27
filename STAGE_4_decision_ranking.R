# =============================================================================
# STAGE 4 — Decision Ranking and Bootstrap Stability
# =============================================================================
# Rank retained candidates using a tiered weighted composite score, then
# assess rank stability via year-block bootstrap and weight sensitivity.
#
# Scoring tiers (weighted rank-aggregation):
#   Tier 1 — Climate structure (50%): month consistency (mean + min),
#            class balance, switching rate.  Full climate record.
#   Tier 2 — Internal robustness (30%): std vs quantile threshold-method
#            agreement (BSA_min, normalized conditional entropy).
#   Tier 3 — External verification (20%): Stage 1 vs Stage 2 agreement
#            (BSA_min, normalized conditional entropy).
#
# Bootstrap: year-block resampling (B replicates). Tier 1 is fixed;
# Tiers 2 and 3 are recomputed on resampled data.
#
# Inputs (RDS from prior stages):
#   Stage 1: screened_tbl.rds, season_long.rds, threshold_tbl.rds
#   Stage 2: ecological_regime_long_seg.rds, seg_tbl.rds, monthly_response.rds
#   Stage 3: stage3_stage1_candidates_retained.rds
#
# Outputs (output_dir/tables/):
#   - decision_table_final.csv      Full decision table with bootstrap ranks
#   - bootstrap_rank_summary.csv    Rank stability summary
#   - weight_sensitivity.csv        Weight-sweep robustness check
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)
})

CONFIG_FILE <- Sys.getenv("SEASON_CONFIG", unset = "config.R")
source(CONFIG_FILE)
set.seed(GLOBAL_SEED)

output_dir <- stage_dir(4)
tab_dir    <- file.path(output_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

screened_tbl   <- readRDS(file.path(stage_dir(1), "screened_tbl.rds"))
season_long    <- readRDS(file.path(stage_dir(1), "season_long.rds"))
threshold_tbl  <- readRDS(file.path(stage_dir(1), "threshold_tbl.rds"))
stage2_long    <- readRDS(file.path(stage_dir(2), "ecological_regime_long_seg.rds")) %>%
  mutate(DateMonth = as.Date(DateMonth))
seg_tbl        <- readRDS(file.path(stage_dir(2), "seg_tbl.rds"))
monthly_response   <- readRDS(file.path(stage_dir(2), "monthly_response.rds")) %>%
  mutate(DateMonth = as.Date(DateMonth))
retained       <- readRDS(file.path(stage_dir(3), "stage3_stage1_candidates_retained.rds"))

response_months <- monthly_response %>% distinct(DateMonth)

# =============================================================================
# 2. APPLY STAGE 3 GATE
# =============================================================================

base_tbl <- screened_tbl %>%
  semi_join(retained %>% dplyr::select(candidate_id), by = "candidate_id") %>%
  mutate(k            = as.integer(k),
         candidate_id = as.character(candidate_id),
         driver       = as.character(driver),
         method       = as.character(method))

# Ecological-window level counts
n_seasons_eco_tbl <- season_long %>%
  semi_join(response_months, by = "DateMonth") %>%
  group_by(candidate_id) %>%
  summarise(n_seasons_eco = n_distinct(season[!is.na(season)]),
            .groups = "drop")

base_tbl <- base_tbl %>%
  left_join(n_seasons_eco_tbl, by = "candidate_id") %>%
  filter(n_seasons_eco >= 2)

# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

# Cohen's κ (handles non-square confusion matrices)
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

# Wilson score interval for discord proportion
wilson_ci <- function(m, n, z = 1.96) {
  if (!is.finite(n) || n <= 0) return(c(lo = NA_real_, hi = NA_real_))
  p <- m / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half   <- (z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)) / denom
  c(lo = center - half, hi = center + half)
}

# Binomial Score Agreement (BSA_min) with optional weighting
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

# Rank-to-[0,1] utility (1 = best)
rank01 <- function(x, higher_better = TRUE) {
  if (!higher_better) x <- -x
  r <- rank(x, ties.method = "average", na.last = "keep")
  n_ok <- sum(!is.na(r))
  if (n_ok <= 1) return(rep(NA_real_, length(x)))
  (r - 1) / (n_ok - 1)
}

# Bin a numeric vector by threshold cut-points
bin_from_cuts <- function(x, cuts) {
  cut(x, breaks = c(-Inf, sort(as.numeric(cuts)), Inf),
      labels = FALSE, right = TRUE)
}

# Mutual information metrics from a contingency table
info_metrics_from_tab <- function(tab) {
  tab <- as.matrix(tab); n <- sum(tab)
  if (!is.finite(n) || n <= 1)
    return(tibble(MI = NA_real_, nH1_given_2 = NA_real_,
                  nH2_given_1 = NA_real_, H1_given_2 = NA_real_, H2_given_1 = NA_real_))
  Pij <- tab / n; Pi <- rowSums(Pij); Pj <- colSums(Pij)
  H <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  Hi <- H(Pi); Hj <- H(Pj)
  MI <- 0
  for (i in seq_len(nrow(Pij)))
    for (j in seq_len(ncol(Pij)))
      if (Pij[i, j] > 0) MI <- MI + Pij[i, j] * log(Pij[i, j] / (Pi[i] * Pj[j]))
  H1g2 <- Hi - MI
  H2g1 <- Hj - MI
  tibble(MI            = MI,
         H1_given_2    = H1g2,
         H2_given_1    = H2g1,
         nH1_given_2   = if (Hi > 0) H1g2 / Hi else NA_real_,
         nH2_given_1   = if (Hj > 0) H2g1 / Hj else NA_real_)
}

# Safe tier mean: returns NA instead of NaN when all components are NA
safe_tier_mean <- function(...) {
  v <- mean(c(...), na.rm = TRUE)
  if (is.nan(v)) NA_real_ else v
}

# Proportional reweighting for missing tiers
weighted_score <- function(tiers, wts) {
  ok <- !is.na(tiers)
  if (sum(ok) == 0) return(NA_real_)
  sum(tiers[ok] * wts[ok]) / sum(wts[ok])
}

# =============================================================================
# 4. STAGE 2 LABEL AGREEMENT (per Stage 1 candidate)
# =============================================================================

stage2_best_match <- function(cid, drv, kk) {
  null_row <- tibble(
    stage2_candidate_id = NA_character_, stage2_n = 0L,
    kappa_ssa = NA_real_, bsa_min_ssa = NA_real_,
    stage2_pct_agree = NA_real_, stage2_discord_ci_hi = NA_real_,
    stage2_pmax = NA_real_, stage2_near_constant = NA,
    stage2_reason = NA_character_)
  s1 <- season_long %>%
    filter(candidate_id == cid) %>%
    semi_join(response_months, by = "DateMonth") %>%
    dplyr::select(DateMonth, s1 = season) %>% distinct()
  if (nrow(s1) == 0)
    return(null_row %>% mutate(stage2_reason = "no_stage1_overlap"))
  s2_all <- stage2_long %>%
    filter(driver == drv, k == kk) %>%
    semi_join(response_months, by = "DateMonth") %>%
    dplyr::select(stage2_candidate_id = candidate_id,
                  DateMonth, s2 = season) %>% distinct()
  if (nrow(s2_all) == 0)
    return(null_row %>% mutate(stage2_reason = "no_stage2_match"))
  scored <- s2_all %>%
    group_by(stage2_candidate_id) %>%
    group_modify(~{
      d <- inner_join(s1, .x %>% dplyr::select(DateMonth, s2), by = "DateMonth")
      met <- label_agreement_bsa(d$s1, d$s2)
      s2_tab <- table(as.character(d$s2))
      pmax2  <- if (sum(s2_tab) > 0) max(as.numeric(s2_tab) / sum(s2_tab))
                else NA_real_
      tibble(stage2_n = met$n, stage2_pct_agree = met$pct_agree,
             stage2_discord_ci_hi = met$discord_ci_hi,
             bsa_min_ssa = met$bsa_min, kappa_ssa = met$kappa,
             stage2_pmax = pmax2,
             stage2_near_constant = is.finite(pmax2) && pmax2 > 0.95)
    }) %>% ungroup() %>%
    arrange(desc(kappa_ssa), desc(bsa_min_ssa), desc(stage2_n)) %>%
    slice(1) %>%
    mutate(stage2_reason = "ok")
  scored %>%
    dplyr::select(stage2_candidate_id, stage2_n, kappa_ssa, bsa_min_ssa,
                  stage2_pct_agree, stage2_discord_ci_hi,
                  stage2_pmax, stage2_near_constant, stage2_reason)
}

stage2_match_df <- base_tbl %>%
  distinct(candidate_id, driver, k) %>%
  mutate(stage2 = pmap(list(candidate_id, driver, k), stage2_best_match)) %>%
  unnest(stage2)

# =============================================================================
# 5. THRESHOLD-METHOD ROBUSTNESS (std vs quantile agreement)
# =============================================================================

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
    summarise(cid_std = candidate_id[method == "std"][1],
      cid_qtl = candidate_id[method == "quantile"][1],
      .groups = "drop") %>% filter(!is.na(cid_std), !is.na(cid_qtl))
  pairs %>%
    mutate(metrics = map2(cid_std, cid_qtl, function(c1, c2) {
      s1 <- sl %>%
        filter(candidate_id == c1) %>% dplyr::select(DateMonth, s_std = season)
      s2 <- sl %>%
        filter(candidate_id == c2) %>% dplyr::select(DateMonth, s_qtl = season)
      joined <- inner_join(s1, s2, by = "DateMonth") %>%
        left_join(wt_tbl, by = "DateMonth")
      met <- label_agreement_bsa(joined$s_std, joined$s_qtl, w = joined$w)
      lv <- union(levels(factor(joined$s_std)), levels(factor(joined$s_qtl)))
      tab <- xtabs(w ~ factor(s_std, levels = lv) + factor(s_qtl, levels = lv),
        data = joined)
      info <- info_metrics_from_tab(tab)
      tibble(stage1_n = sum(joined$w, na.rm = TRUE),
        kappa_std_quant = met$kappa, bsa_min_std_quant = met$bsa_min,
        pct_agree_std_quant = met$pct_agree,
        discord_ci_hi_std_quant = met$discord_ci_hi,
        nce_std_quant = mean(c(info$nH1_given_2, info$nH2_given_1), na.rm = TRUE))
    })) %>%
    unnest(metrics)
}

sq_tbl <- std_quant_agreement(season_long)

# =============================================================================
# 6. STRUCTURAL SEGMENTATION ALIGNMENT (SSA)
# =============================================================================
# Bin the driver at Stage 1 thresholds and Stage 2 breakpoints independently,
# then compute normalized conditional entropy between the two binnings.
# This measures whether the classification boundaries produce structurally
# equivalent partitions of the driver, independent of season labels.

ssa_tbl <- base_tbl %>%
  distinct(candidate_id, driver, k) %>%
  pmap_dfr(function(candidate_id, driver, k) {
    cid <- candidate_id; drv <- driver; kk <- as.integer(k)
    null_ssa <- function(reason)
      tibble(candidate_id = cid, ssa_reason = reason, ssa_nce = NA_real_)
    # Stage 1 thresholds
    thr <- threshold_tbl %>% filter(candidate_id == cid) %>% slice(1)
    if (nrow(thr) == 0) return(null_ssa("no_stage1_thresholds"))
    s1_cuts <- if (kk == 2) c(thr$t1[1]) else c(thr$t1[1], thr$t2[1])
    if (any(!is.finite(s1_cuts)) || (kk == 3 && s1_cuts[1] >= s1_cuts[2]))
      return(null_ssa("degenerate_stage1_thresholds"))
    # Stage 2 breakpoints
    br <- seg_tbl %>% filter(driver == drv, k == kk) %>% slice(1)
    if (nrow(br) == 0) return(null_ssa("no_stage2_breaks"))
    s2_cuts <- if (kk == 2) c(br$b1[1]) else c(br$b1[1], br$b2[1])
    if (any(!is.finite(s2_cuts)) || (kk == 3 && s2_cuts[1] >= s2_cuts[2]))
      return(null_ssa("degenerate_stage2_breaks"))
    # Bin driver values at both cut-point sets
    d <- monthly_response %>%
      semi_join(response_months, by = "DateMonth") %>%
      distinct(DateMonth, x = .data[[drv]]) %>%
      filter(is.finite(x))
    if (nrow(d) <= 1) return(null_ssa("insufficient_data"))
    tab <- table(bin_from_cuts(d$x, s1_cuts), bin_from_cuts(d$x, s2_cuts))
    met <- info_metrics_from_tab(tab)
    nce <- mean(c(met$nH1_given_2, met$nH2_given_1), na.rm = TRUE)
    tibble(candidate_id = cid, ssa_reason = "ok",
           nce_ssa = if (is.nan(nce)) NA_real_ else nce)
  })

# =============================================================================
# 7. ASSEMBLE DECISION SET AND COMPUTE SCORE
# =============================================================================

decision_set <- base_tbl %>%
  left_join(stage2_match_df,       by = c("candidate_id", "driver", "k")) %>%
  left_join(sq_tbl,                by = c("driver", "k")) %>%
  left_join(ssa_tbl,               by = "candidate_id")

# Tiered rank-aggregation scoring
decision_set <- decision_set %>%
  mutate(
    # Tier 1: Climate structure (4 components)
    u_mean_month = rank01(mean_month_consistency, TRUE),
    u_min_month  = rank01(min_month_consistency,  TRUE),
    u_min_bin    = rank01(min_bin_prop,            TRUE),
    u_switch     = rank01(mean_switch_per_year,    FALSE),
    # Tier 2: Internal robustness (2 components)
    u_sq_bsa     = rank01(bsa_min_std_quant,       TRUE),
    u_sq_ce      = rank01(nce_std_quant,           FALSE),
    # Tier 3: External verification (2 components)
    u_s2_bsa     = rank01(bsa_min_ssa,          TRUE),
    u_ssa_ce     = rank01(nce_ssa,                 FALSE)) %>%
  rowwise() %>%
  mutate(
    tier_climate = safe_tier_mean(u_mean_month, u_min_month, u_min_bin, u_switch),
    tier_robust  = safe_tier_mean(u_sq_bsa, u_sq_ce),
    tier_verify  = safe_tier_mean(u_s2_bsa, u_ssa_ce),
    climate_score = weighted_score(
      c(tier_climate, tier_robust, tier_verify),
      c(W_CLIMATE, W_ROBUST, W_VERIFY)),
    score_n_components = sum(!is.na(c(
      u_mean_month, u_min_month, u_min_bin, u_switch,
      u_sq_bsa, u_sq_ce, u_s2_bsa, u_ssa_ce)))) %>%
  ungroup()

# =============================================================================
# 8. DECISION TABLE OUTPUT
# =============================================================================

decision_table <- decision_set %>%
  arrange(desc(climate_score), desc(bsa_min_std_quant),
          desc(bsa_min_ssa)) %>%
  dplyr::select(
    candidate_id, driver, k, method,
    # Tier 1
    climate_score, score_n_components,
    mean_month_consistency, min_month_consistency,
    min_bin_prop, entropy_norm, mean_switch_per_year,
    # Tier 2
    stage1_n, bsa_min_std_quant,
    kappa_std_quant, nce_std_quant,
    # Tier 3
    stage2_n, bsa_min_ssa,
    kappa_ssa, stage2_pmax, stage2_near_constant,
    ssa_reason, nce_ssa)

# =============================================================================
# 9. BOOTSTRAP RANK STABILITY
# =============================================================================

years_full <- sort(unique(year(season_long$DateMonth)))
years_response <- sort(unique(year(response_months$DateMonth)))

# Year-block resampler (multiplicity preserved via boot_block index)
boot_months <- function(src_tbl, years_vec) {
  yrs <- sample(years_vec, length(years_vec), replace = TRUE)
  map2_dfr(yrs, seq_along(yrs), ~{
    src_tbl %>%
      mutate(Year = year(DateMonth)) %>%
      filter(Year == .x) %>%
      mutate(.boot_block = .y)
  })
}

# Recompute resampled scoring components for one bootstrap iteration
recompute_boot_components <- function(months_full_b, months_response_b) {
  full_w <- months_full_b %>% count(DateMonth, name = "w")
  response_w <- months_response_b %>% count(DateMonth, name = "w")
  # Std vs quantile robustness (resampled full climate window)
  sq_b <- std_quant_agreement(
    season_long %>% semi_join(full_w, by = "DateMonth") %>%
      left_join(full_w, by = "DateMonth"))
  # Stage 2 label agreement (resampled ecological window)
  s2_b <- base_tbl %>%
    distinct(candidate_id, driver, k) %>%
    mutate(stage2 = pmap(list(candidate_id, driver, k), function(cid, drv, kk) {
      s1 <- season_long %>% filter(candidate_id == cid) %>%
        semi_join(response_w, by = "DateMonth") %>%
        dplyr::select(DateMonth, s1 = season) %>% distinct()
      s2_all <- stage2_long %>% filter(driver == drv, k == kk) %>%
        semi_join(response_w, by = "DateMonth") %>%
        dplyr::select(stage2_candidate_id = candidate_id,
                      DateMonth, s2 = season) %>% distinct()
      if (nrow(s1) == 0 || nrow(s2_all) == 0)
        return(tibble(stage2_candidate_id = NA_character_, stage2_n = 0L,
                      kappa_ssa = NA_real_, bsa_min_ssa = NA_real_,
                      stage2_pct_agree = NA_real_, stage2_discord_ci_hi = NA_real_,
                      stage2_pmax = NA_real_, stage2_near_constant = NA,
                      stage2_reason = "no_overlap"))
      s2_all %>%
        group_by(stage2_candidate_id) %>%
        group_modify(~{
          d <- inner_join(s1, .x %>% dplyr::select(DateMonth, s2),
                          by = "DateMonth") %>%
            left_join(response_w, by = "DateMonth")
          met <- label_agreement_bsa(d$s1, d$s2, w = d$w)
          w_by_s2 <- tapply(d$w, as.character(d$s2), sum)
          pmax2 <- if (sum(w_by_s2, na.rm = TRUE) > 0)
            max(w_by_s2 / sum(w_by_s2), na.rm = TRUE) else NA_real_
          tibble(stage2_n = met$n, stage2_pct_agree = met$pct_agree,
                 stage2_discord_ci_hi = met$discord_ci_hi,
                 bsa_min_ssa = met$bsa_min, kappa_ssa = met$kappa,
                 stage2_pmax = pmax2,
                 stage2_near_constant = is.finite(pmax2) && pmax2 > 0.95)
        }) %>% ungroup() %>%
        arrange(desc(kappa_ssa), desc(bsa_min_ssa), desc(stage2_n)) %>%
        slice(1) %>% mutate(stage2_reason = "ok") %>%
        dplyr::select(stage2_candidate_id, stage2_n, kappa_ssa,
                      bsa_min_ssa, stage2_pct_agree, stage2_discord_ci_hi,
                      stage2_pmax, stage2_near_constant, stage2_reason)
    })) %>% unnest(stage2)
  list(sq = sq_b, s2 = s2_b)
}

# Single bootstrap replicate → candidate ranks
boot_once <- function() {
  months_full_b <- boot_months(season_long, years_full)
  months_response_b <- boot_months(response_months, years_response)
  comps <- recompute_boot_components(months_full_b, months_response_b)
  decision_set %>%
    dplyr::select(candidate_id, driver, k,
                  mean_month_consistency, min_month_consistency,
                  min_bin_prop, mean_switch_per_year,
                  nce_ssa) %>%
    left_join(comps$sq, by = c("driver", "k")) %>%
    left_join(comps$s2 %>% dplyr::select(candidate_id, driver, k, bsa_min_ssa),
              by = c("candidate_id", "driver", "k")) %>%
    mutate(
      u_mean_month = rank01(mean_month_consistency, TRUE),
      u_min_month  = rank01(min_month_consistency,  TRUE),
      u_min_bin    = rank01(min_bin_prop,            TRUE),
      u_switch     = rank01(mean_switch_per_year,    FALSE),
      u_sq_bsa     = rank01(bsa_min_std_quant,       TRUE),
      u_sq_ce      = rank01(nce_std_quant,           FALSE),
      u_s2_bsa     = rank01(bsa_min_ssa,          TRUE),
      u_ssa_ce     = rank01(nce_ssa,                 FALSE)) %>%
    rowwise() %>%
    mutate(
      tier_climate  = safe_tier_mean(u_mean_month, u_min_month, u_min_bin, u_switch),
      tier_robust   = safe_tier_mean(u_sq_bsa, u_sq_ce),
      tier_verify   = safe_tier_mean(u_s2_bsa, u_ssa_ce),
      climate_score = weighted_score(
        c(tier_climate, tier_robust, tier_verify),
        c(W_CLIMATE, W_ROBUST, W_VERIFY))) %>%
    ungroup() %>%
    arrange(desc(climate_score)) %>%
    mutate(rank = row_number()) %>%
    dplyr::select(candidate_id, rank)
}

# Run bootstrap
boot_ranks <- map_dfr(seq_len(BOOT_N_RANK), ~boot_once() %>% mutate(iter = .x))

# Stability summaries
rank_stats <- boot_ranks %>%
  group_by(candidate_id) %>%
  summarise(p_top1 = mean(rank == 1), rank_IQR = IQR(rank),
            .groups = "drop") %>%
  arrange(desc(p_top1))

top_probs <- boot_ranks %>%
  filter(rank == 1) %>% count(candidate_id) %>% mutate(prob = n / BOOT_N_RANK)

boot_summary <- tibble(
  N_BOOT               = BOOT_N_RANK,
  top_candidate        = top_probs$candidate_id[1],
  top_probability      = top_probs$prob[1],
  runnerup_probability = if (nrow(top_probs) >= 2) top_probs$prob[2] else NA_real_,
  decision_entropy     = -sum(top_probs$prob * log(top_probs$prob)))

decision_table_final <- decision_table %>%
  left_join(rank_stats, by = "candidate_id") %>%
  arrange(desc(p_top1), desc(climate_score))

# =============================================================================
# 10. WEIGHT SENSITIVITY ANALYSIS
# =============================================================================

weight_grid <- expand.grid(
  w_clim = seq(0.3, 0.7, by = 0.1),
  w_rob  = seq(0.1, 0.4, by = 0.1)) %>%
  mutate(w_ver = 1 - w_clim - w_rob) %>%
  filter(w_ver >= 0.1, w_ver <= 0.4)

weight_sensitivity <- weight_grid %>%
  pmap_dfr(function(w_clim, w_rob, w_ver) {
    scored <- decision_set %>%
      rowwise() %>%
      mutate(cs = weighted_score(
        c(tier_climate, tier_robust, tier_verify),
        c(w_clim, w_rob, w_ver))) %>%
      ungroup() %>% arrange(desc(cs))
    tibble(top_candidate = scored$candidate_id[1],
           top_score     = scored$cs[1],
           w_climate     = w_clim,
           w_robust      = w_rob,
           w_verify      = w_ver,
           gap_to_second = scored$cs[1] - scored$cs[2])
  })

message("Stage 4 complete. Winner: ", boot_summary$top_candidate,
        " (P(rank 1) = ", round(boot_summary$top_probability, 3), ").",
        " Weight sensitivity: ", n_distinct(weight_sensitivity$top_candidate),
        " unique winner(s) across ", nrow(weight_grid), " weight combos.")

# =============================================================================
# 11. SAVE OUTPUTS
# =============================================================================

write.csv(decision_table_final %>%
            dplyr::select(
              candidate_id, driver, n_seasons = k, method,
              climate_score, p_top1, rank_IQR, score_n_components,
              # Tier 1 — Climate structure
              mean_month_consistency, min_month_consistency,
              min_bin_prop, mean_switch_per_year, entropy_norm,
              # Tier 2 — Internal robustness (std vs quantile)
              stage1_n, bsa_min_std_quant, nce_std_quant, kappa_std_quant,
              # Tier 3 — External verification (Stage 2 agreement)
              stage2_n, bsa_min_ssa, nce_ssa, kappa_ssa,
              stage2_pmax, stage2_near_constant, ssa_reason),
          file.path(tab_dir, "decision_table_final.csv"), row.names = FALSE)

write.csv(weight_sensitivity %>%
            dplyr::select(w_climate, w_robust, w_verify,
                          top_candidate, top_score, gap_to_second),
          file.path(tab_dir, "weight_sensitivity.csv"), row.names = FALSE)

write.csv(boot_summary,
          file.path(tab_dir, "bootstrap_rank_summary.csv"), row.names = FALSE)

# RDS objects for publication outputs
saveRDS(decision_set,  file.path(output_dir, "decision_set.rds"))
saveRDS(boot_ranks,    file.path(output_dir, "boot_ranks.rds"))
# =============================================================================
