# =============================================================================
# STAGE 3 — Season Factor Stress-Testing and Validation
# =============================================================================
# Stress-test Stage 1 candidates over the ecological measurement window and
# cross-validate against Stage 2 ecological-derived breakpoints. Only
# structural pathologies trigger drops; ecological agreement and ANOVA are
# informational flags.
#
# Inputs (RDS from prior stages):
#   Stage 1: season_long.rds, threshold_tbl.rds, screened_tbl.rds
#   Stage 2: ecological_regime_long_seg.rds, seg_tbl.rds, monthly_response.rds
#   Climate: CLIMATE_CSV
#
# Outputs (output_dir/tables/):
#   - filter_results.csv           All eco-window diagnostics + drop/flag decisions
#   - retained_candidates.csv      Survivors for Stage 4
#   - block_stability.csv          Per-block stability detail
#   - anova_summary.csv            ANOVA effect sizes per candidate
#   - anova_group_means.csv        Per-season means
#   - tukey_posthoc.csv            Pairwise comparisons (k=3 only)
#   - stage1_vs_stage2_comparison.csv  Agreement + threshold alignment
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

CONFIG_FILE <- Sys.getenv("SEASON_CONFIG", unset = "config.R")
source(CONFIG_FILE)
set.seed(GLOBAL_SEED)

output_dir <- stage_dir(3)
tab_dir    <- file.path(output_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

season_long               <- readRDS(file.path(stage_dir(1), "season_long.rds"))
threshold_tbl             <- readRDS(file.path(stage_dir(1), "threshold_tbl.rds"))
screened_tbl              <- readRDS(file.path(stage_dir(1), "screened_tbl.rds"))
ecological_regime_long_seg <- readRDS(file.path(stage_dir(2), "ecological_regime_long_seg.rds"))
seg_tbl                   <- readRDS(file.path(stage_dir(2), "seg_tbl.rds"))

monthly_clim <- read.csv(CLIMATE_CSV, stringsAsFactors = FALSE) %>%
  mutate(Year      = as.integer(Year),
         Month     = as.integer(Month),
         DateMonth = as.Date(sprintf("%d-%02d-15", Year, Month)))

monthly_response <- readRDS(file.path(stage_dir(2), "monthly_response.rds")) %>%
  mutate(DateMonth = as.Date(DateMonth))

# Ecological-month window: only months with ecological data
response_months <- monthly_response %>% distinct(DateMonth)

# =============================================================================
# 2. RESTRICT TO ECOLOGICAL WINDOW
# =============================================================================

stage1_long <- season_long %>%
  semi_join(response_months, by = "DateMonth") %>%
  dplyr::select(candidate_id, driver, k, method, DateMonth, season)

stage2_long <- ecological_regime_long_seg %>%
  mutate(DateMonth = as.Date(DateMonth)) %>%
  semi_join(response_months, by = "DateMonth") %>%
  { if (!"method" %in% names(.)) mutate(., method = "segmented") else . } %>%
  dplyr::select(candidate_id, driver, k, method, DateMonth, season)

# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

# Class balance within a season vector
balance_metrics <- function(season_vec) {
  tab <- table(droplevels(season_vec))
  if (sum(tab) == 0)
    return(tibble(min_bin_prop = NA_real_, min_bin_n = 0L,
                  n_levels_used = 0L))
  tibble(min_bin_prop = min(as.numeric(tab) / sum(tab)),
         min_bin_n = as.integer(min(tab)),
         n_levels_used = length(tab))
}

# Calendar-month consistency
month_consistency <- function(df, season_col = "season") {
  d <- df %>% filter(!is.na(.data[[season_col]]))
  if (nrow(d) == 0)
    return(tibble(mean_month_consistency = NA_real_,
                  min_month_consistency = NA_real_))
  mm <- d %>%
    count(Month, Year, .data[[season_col]], name = "n") %>%
    group_by(Month, Year) %>% slice_max(n, with_ties = FALSE) %>% ungroup() %>%
    count(Month, .data[[season_col]], name = "n") %>%
    group_by(Month) %>% mutate(prop = n / sum(n)) %>%
    summarise(max_prop = max(prop), .groups = "drop")
  tibble(mean_month_consistency = mean(mm$max_prop),
         min_month_consistency  = min(mm$max_prop))
}

# Partition dates into non-overlapping year blocks
make_year_blocks <- function(dates, block_years = 2) {
  yrs <- sort(unique(year(dates)))
  blocks <- list(); i <- 1
  while (i <= length(yrs)) {
    blocks[[length(blocks) + 1]] <- yrs[i:min(i + block_years - 1, length(yrs))]
    i <- i + block_years
  }
  blocks
}

# Cohen's κ from a contingency table
kappa_cohen <- function(tab) {
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  if (!is.finite(pe) || isTRUE(all.equal(1 - pe, 0))) return(NA_real_)
  (po - pe) / (1 - pe)
}

# =============================================================================
# 4. STRESS TESTS — BLOCK STABILITY AND CALENDAR CONSISTENCY
# =============================================================================

year_blocks <- make_year_blocks(response_months$DateMonth, block_years = S3_BLOCK_YEARS)

block_stability <- stage1_long %>%
  mutate(Year = year(DateMonth)) %>%
  group_by(candidate_id, driver, k, method) %>%
  group_modify(~{
    d <- .x
    map_dfr(year_blocks, function(yrs) {
      db  <- d %>% filter(Year %in% yrs)
      bal <- balance_metrics(db$season)
      tibble(test_years         = paste(range(yrs), collapse = "-"),
             n_block            = nrow(db),
             n_assigned_block   = sum(!is.na(db$season)),
             n_levels_block     = bal$n_levels_used,
             min_bin_prop_block = bal$min_bin_prop,
             min_bin_n_block    = bal$min_bin_n)
    })
  }) %>% ungroup()

calendar_consistency_eco <- stage1_long %>%
  mutate(Month = month(DateMonth), Year = year(DateMonth)) %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    cal = list(month_consistency(tibble(Month = Month, Year = Year, season = season))),
    .groups = "drop") %>%
  unnest_wider(cal) %>%
  rename(mean_month_consistency_eco = mean_month_consistency,
         min_month_consistency_eco  = min_month_consistency)

# =============================================================================
# 5. ECOLOGICAL DIFFERENTIATION (ANOVA)
# =============================================================================
# One-way ANOVA of monthly ecological response grouped by season level.
# Informational flag only — never triggers a drop.

ecological_anova <- stage1_long %>%
  inner_join(monthly_response, by = "DateMonth") %>%
  filter(!is.na(season), is.finite(RESPONSE_COL)) %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    anova_result = list({
      d <- tibble(RESPONSE_COL = RESPONSE_COL, season = droplevels(season))
      if (n_distinct(d$season) < 2 || any(table(d$season) < 3)) {
        tibble(n_obs = nrow(d), n_levels = n_distinct(d$season),
               p_value = NA_real_, omega_sq = NA_real_,
               kw_p = NA_real_, group_means = list(NULL))
      } else {
        fit <- aov(RESPONSE_COL ~ season, data = d)
        af  <- summary(fit)[[1]]
        SS_b <- af[["Sum Sq"]][1]; SS_t <- sum(af[["Sum Sq"]])
        df_b <- af[["Df"]][1]; MS_w <- af[["Mean Sq"]][2]
        omega_sq <- max(0, (SS_b - df_b * MS_w) / (SS_t + MS_w))
        gmeans <- d %>% group_by(season) %>%
          summarise(mean_response = mean(RESPONSE_COL), sd_response = sd(RESPONSE_COL),
                    n = n(), .groups = "drop")
        tibble(n_obs = nrow(d), n_levels = n_distinct(d$season),
               p_value = af[["Pr(>F)"]][1], omega_sq = omega_sq,
               kw_p = tryCatch(kruskal.test(RESPONSE_COL ~ season, data = d)$p.value,
                               error = function(e) NA_real_),
               group_means = list(gmeans))
      }}),
    .groups = "drop") %>%
  unnest(anova_result)

group_means_tbl <- ecological_anova %>%
  dplyr::select(candidate_id, driver, k, method, group_means) %>%
  filter(!map_lgl(group_means, is.null)) %>%
  unnest(group_means)

ecological_anova_summary <- ecological_anova %>% dplyr::select(-group_means)

# Tukey HSD post-hoc (k = 3 only)
posthoc_tbl <- stage1_long %>%
  filter(k == 3) %>%
  inner_join(monthly_response, by = "DateMonth") %>%
  filter(!is.na(season), is.finite(RESPONSE_COL)) %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    tukey = list({
      d <- tibble(RESPONSE_COL = RESPONSE_COL, season = droplevels(season))
      if (n_distinct(d$season) < 3 || any(table(d$season) < 2)) NULL
      else {
        tk <- tryCatch(TukeyHSD(aov(RESPONSE_COL ~ season, data = d))$season,
                       error = function(e) NULL)
        if (is.null(tk)) NULL
        else as.data.frame(tk) %>% rownames_to_column("comparison") %>% as_tibble()
      }
    }), .groups = "drop") %>%
  unnest(tukey)

ecological_anova_flags <- ecological_anova_summary %>%
  transmute(candidate_id, driver, k, method,
            anova_omega_sq = omega_sq, anova_p = p_value, kw_p,
            flag_low_omega_sq = is.finite(omega_sq) & omega_sq < S3_FLAG_OMEGA_SQ_LOW)

# =============================================================================
# 6. STAGE 1 vs STAGE 2 AGREEMENT (Validation Signal)
# =============================================================================

stage2_key <- stage2_long %>%
  distinct(driver, k, candidate_id) %>%
  group_by(driver, k) %>%
  summarise(stage2_candidate_id = first(candidate_id), .groups = "drop")

agreement_tbl <- stage1_long %>%
  distinct(candidate_id, driver, k, method) %>%
  left_join(stage2_key, by = c("driver", "k")) %>%
  filter(!is.na(stage2_candidate_id)) %>%
  mutate(metrics = pmap(
    list(candidate_id, stage2_candidate_id, driver, k),
    function(cid1, cid2, drv, kk) {
      d <- stage1_long %>%
        filter(candidate_id == cid1) %>%
        dplyr::select(DateMonth, s1 = season) %>%
        inner_join(
          stage2_long %>% filter(candidate_id == cid2) %>%
            dplyr::select(DateMonth, s2 = season),
          by = "DateMonth")
      lv <- union(levels(d$s1), levels(d$s2))
      tab <- table(factor(d$s1, levels = lv), factor(d$s2, levels = lv))
      tibble(n_months_overlap = sum(tab), kappa = kappa_cohen(tab))
    })) %>%
  unnest(metrics)

# =============================================================================
# 7. THRESHOLD vs BREAKPOINT ALIGNMENT (Validation Signal)
# =============================================================================

stage2_breaks <- seg_tbl %>% transmute(driver, k, b1, b2)

driver_scale <- monthly_clim %>%
  dplyr::select(DateMonth, any_of(DRIVER_META$driver)) %>%
  pivot_longer(-DateMonth, names_to = "driver", values_to = "x") %>%
  summarise(iqr = IQR(x, na.rm = TRUE), .by = driver)

alignment_tbl <- threshold_tbl %>%
  dplyr::select(candidate_id, driver, k, method, t1, t2) %>%
  inner_join(stage2_breaks, by = c("driver", "k")) %>%
  left_join(driver_scale, by = "driver") %>%
  mutate(
    abs_diff_1 = abs(b1 - t1),
    abs_diff_2 = if_else(k == 3, abs(b2 - t2), NA_real_),
    diff1_iqr  = abs_diff_1 / iqr,
    diff2_iqr  = if_else(k == 3, abs_diff_2 / iqr, NA_real_))

# =============================================================================
# 8. CONSERVATIVE FILTERING
# =============================================================================
# Drop rules target structural pathologies only.

block_flags <- block_stability %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    n_blocks    = n(),
    n_collapsed = sum(n_levels_block < k, na.rm = TRUE),
    prop_healthy = (n_blocks - n_collapsed) / n_blocks,
    min_block_min_bin_prop = {
      vals <- min_bin_prop_block[is.finite(min_bin_prop_block)]
      if (length(vals) == 0) NA_real_ else min(vals)
    },
    any_block_extreme_imbalance = any(
      is.finite(min_bin_prop_block) & min_bin_prop_block < S3_MIN_BLOCK_PROP,
      na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    min_block_min_bin_prop = if_else(is.infinite(min_block_min_bin_prop),
                                     NA_real_, min_block_min_bin_prop),
    fail_block_collapse = prop_healthy < S3_MIN_HEALTHY_BLOCKS)

# Ecological-window overview and balance
ecowin_overview <- stage1_long %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    n_months_eco     = n(),
    n_assigned_eco   = sum(!is.na(season)),
    pct_assigned_eco = n_assigned_eco / n_months_eco,
    n_levels_eco     = n_distinct(season[!is.na(season)]),
    bal = list(balance_metrics(season)),
    .groups = "drop") %>%
  unnest_wider(bal) %>%
  rename(min_bin_prop_eco = min_bin_prop,
         min_bin_n_eco    = min_bin_n)

# Assemble diagnostics and apply drop/flag rules
filters_tbl <- screened_tbl %>%
  dplyr::select(candidate_id, driver, k, method) %>%
  mutate(k = as.integer(k)) %>%
  left_join(ecowin_overview,
            by = c("candidate_id", "driver", "k", "method")) %>%
  left_join(calendar_consistency_eco,
            by = c("candidate_id", "driver", "k", "method")) %>%
  left_join(block_flags,
            by = c("candidate_id", "driver", "k", "method")) %>%
  left_join(agreement_tbl %>%
              dplyr::select(candidate_id, driver, k, method,
                            n_months_overlap, kappa),
            by = c("candidate_id", "driver", "k", "method")) %>%
  left_join(alignment_tbl %>%
              mutate(align_flag = (is.finite(diff1_iqr) & diff1_iqr > S3_FLAG_ALIGN_IQR) |
                       (k == 3 & is.finite(diff2_iqr) & diff2_iqr > S3_FLAG_ALIGN_IQR)) %>%
              dplyr::select(candidate_id, driver, k, method,
                            diff1_iqr, diff2_iqr, align_flag),
            by = c("candidate_id", "driver", "k", "method")) %>%
  left_join(ecological_anova_flags,
            by = c("candidate_id", "driver", "k", "method")) %>%
  mutate(
    # --- Drop rules ---
    fail_assignment = is.finite(pct_assigned_eco) &
      pct_assigned_eco < S3_MIN_PCT_ASSIGNED,
    fail_imbalance  = is.finite(min_bin_prop_eco) &
      min_bin_prop_eco < S3_MIN_SEASON_PROP,
    fail_calendar   = is.finite(mean_month_consistency_eco) &
      mean_month_consistency_eco < S3_MEAN_MONTH_CONS,
    fail_block_imbalance = !is.na(any_block_extreme_imbalance) &
      any_block_extreme_imbalance,
    # --- Flags ---
    flag_kappa_low = is.finite(kappa) & kappa < S3_FLAG_KAPPA_LOW,
    flag_align_far = !is.na(align_flag) & align_flag,
    # --- Final ---
    drop_candidate = fail_assignment | fail_imbalance | fail_calendar |
      (!is.na(fail_block_collapse) & fail_block_collapse) |
      fail_block_imbalance)

# =============================================================================
# 9. RETAINED CANDIDATES
# =============================================================================

retained <- filters_tbl %>%
  filter(!drop_candidate) %>%
  dplyr::select(
    candidate_id, driver, k, method,
    pct_assigned_eco, n_levels_eco, min_bin_prop_eco,
    mean_month_consistency_eco, min_month_consistency_eco,
    prop_healthy, n_collapsed, min_block_min_bin_prop,
    n_months_overlap, kappa, diff1_iqr, diff2_iqr,
    flag_kappa_low, flag_align_far,
    anova_omega_sq, anova_p, kw_p, flag_low_omega_sq)

# =============================================================================
# 10. SAVE OUTPUTS
# =============================================================================

write.csv(filters_tbl %>%
            dplyr::select(
              candidate_id, driver, n_seasons = k, method,
              pct_assigned_eco, n_levels_eco,
              min_bin_prop_eco, min_bin_n_eco,
              mean_month_consistency_eco, min_month_consistency_eco,
              prop_healthy, n_collapsed, min_block_min_bin_prop,
              fail_block_collapse, any_block_extreme_imbalance,
              n_months_overlap, kappa, diff1_iqr, diff2_iqr,
              anova_omega_sq, anova_p, kw_p,
              flag_kappa_low, flag_align_far, flag_low_omega_sq,
              fail_assignment, fail_imbalance, fail_calendar,
              fail_block_imbalance, drop_candidate),
          file.path(tab_dir, "filter_results.csv"), row.names = FALSE)

write.csv(retained %>%
            dplyr::select(
              candidate_id, driver, n_seasons = k, method,
              pct_assigned_eco, n_levels_eco, min_bin_prop_eco,
              mean_month_consistency_eco, min_month_consistency_eco,
              prop_healthy, n_collapsed, min_block_min_bin_prop,
              n_months_overlap, kappa, diff1_iqr, diff2_iqr,
              flag_kappa_low, flag_align_far,
              anova_omega_sq, anova_p, kw_p, flag_low_omega_sq),
          file.path(tab_dir, "retained_candidates.csv"), row.names = FALSE)

write.csv(block_stability %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          test_years, n_block, n_assigned_block,
                          n_levels_block, min_bin_prop_block, min_bin_n_block),
          file.path(tab_dir, "block_stability.csv"), row.names = FALSE)

write.csv(ecological_anova_summary %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          n_obs, n_levels, p_value, omega_sq, kw_p),
          file.path(tab_dir, "anova_summary.csv"), row.names = FALSE)

write.csv(posthoc_tbl %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          comparison, diff, lwr, upr, p_adj = `p adj`),
          file.path(tab_dir, "tukey_posthoc.csv"), row.names = FALSE)

write.csv(group_means_tbl %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          season, mean_response, sd_response, n),
          file.path(tab_dir, "anova_group_means.csv"), row.names = FALSE)

stage1_vs_stage2 <- agreement_tbl %>%
  left_join(alignment_tbl %>%
              dplyr::select(candidate_id, driver, k, method,
                            t1, t2, b1, b2, diff1_iqr, diff2_iqr),
            by = c("candidate_id", "driver", "k", "method"))

write.csv(stage1_vs_stage2 %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          stage2_candidate_id, n_months_overlap, kappa,
                          threshold_1 = t1, threshold_2 = t2,
                          breakpoint_1 = b1, breakpoint_2 = b2,
                          diff1_iqr, diff2_iqr),
          file.path(tab_dir, "stage1_vs_stage2_comparison.csv"),
          row.names = FALSE)

saveRDS(retained, file.path(output_dir, "stage3_stage1_candidates_retained.rds"))
# =============================================================================
