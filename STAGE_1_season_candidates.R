# =============================================================================
# STAGE 1 — Season Candidate Generation and Climate-Only Screening
# =============================================================================
#   Generate candidate season definitions from monthly climate data, screen
#   for structural quality, and Pareto-rank survivors by climatological
#   coherence.  No ecological data are used at this stage.
#
#   Each climate driver (DRIVER_META) is crossed with 2 or 3 season levels
#   and two threshold methods (literature-standard vs. baseline-quantile).
#   Candidates are evaluated on five climate-structure metrics, screened for
#   degeneracy, and Pareto-ranked.
#
# Inputs:
#   - config.R
#   - Monthly climate CSV (Year, Month, driver columns)
#
# Outputs (output_dir/tables/):
#   - season_long_data.csv         Season labels joined to climate data
#   - thresholds.csv               Threshold values per candidate
#   - validation_tbl.csv           Climate-structure metrics
#   - screened_out_candidates.csv  Removed candidates with failure reasons
#   - .rds objects for downstream stages
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)   
  library(lubridate)   
  library(zoo)         
})

CONFIG_FILE <- Sys.getenv("SEASON_CONFIG", unset = "config.R")
source(CONFIG_FILE)
set.seed(GLOBAL_SEED)

output_dir <- stage_dir(1)
tab_dir    <- file.path(output_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. CLIMATE INPUTS — load and validate the monthly driver record
# =============================================================================
# Input CSV must contain Year, Month, and all DRIVER_META$driver columns.
# Computed indices (e.g. SPEI, rolling sums, CWD) should be pre-derived.

monthly_clim <- read.csv(CLIMATE_CSV, stringsAsFactors = FALSE) %>%
  mutate(Year  = as.integer(Year), Month = as.integer(Month),
         DateMonth = as.Date(sprintf("%d-%02d-15", Year, Month)))

missing_drv <- setdiff(DRIVER_META$driver, names(monthly_clim))
if (length(missing_drv) > 0)
  stop("Drivers missing from climate CSV: ", paste(missing_drv, collapse = ", "))

# Warn on duplicate Year-Month rows: metrics are computed per calendar month and
# duplicate rows silently inflate counts and distort class-balance metrics.
dup_ym <- monthly_clim %>% count(Year, Month) %>% filter(n > 1)
if (nrow(dup_ym) > 0)
  warning(sprintf(
    "Climate CSV has %d duplicate Year-Month combination(s) (e.g. %d-%02d). ",
    nrow(dup_ym), dup_ym$Year[1], dup_ym$Month[1]),
    "Season metrics may be inflated. Deduplicate the CSV before running.")

baseline_months <- monthly_clim %>%
  filter(Year >= BASELINE_START, Year <= BASELINE_END)

# An empty baseline means all quantile thresholds will be NA, producing only
# degenerate candidates and a silent all-NA failure downstream.
if (nrow(baseline_months) == 0)
  stop(sprintf(
    "Baseline period %d--%d produced no rows. ",
    BASELINE_START, BASELINE_END),
    sprintf("Climate CSV spans %d--%d. Adjust BASELINE_START / BASELINE_END.",
            min(monthly_clim$Year, na.rm = TRUE),
            max(monthly_clim$Year, na.rm = TRUE)))

# Fewer than 24 baseline months makes quantile thresholds statistically unreliable.
n_baseline_finite <- baseline_months %>%
  dplyr::select(all_of(DRIVER_META$driver)) %>%
  summarise(across(everything(), ~sum(is.finite(.)))) %>%
  unlist()
if (any(n_baseline_finite < 24))
  warning("Baseline period has <24 finite values for driver(s): ",
          paste(names(n_baseline_finite)[n_baseline_finite < 24], collapse = ", "),
          ". Quantile thresholds for these drivers will be NA (degenerate candidates).")

# =============================================================================
# 2. SEASON-LABELLING UTILITIES — threshold assignment, quantile extraction,
#    and display rounding used throughout candidate generation
# =============================================================================

# Polarity-aware two-level labeller. lower_closed = TRUE for drivers where
# high values indicate dry conditions (e.g. CWD), so the threshold value
# itself falls in the "low" (wet) bin rather than the "high" (dry) bin.
# This keeps the boundary convention consistent with high_is_dry polarity.
# Guard: when t is non-finite (e.g. NA_real_ from an under-powered baseline),
# return all-NA rather than silently assigning every finite x to the "high"
# bin (which is what case_when does when x <= NA evaluates to NA and falls
# through to the TRUE clause).
assign_2season <- function(x, t, low = "Low", high = "High",
                           lower_closed = FALSE) {
  if (!is.finite(t)) return(rep(NA_character_, length(x)))
  if (lower_closed)
    case_when(!is.finite(x) ~ NA_character_, x <= t ~ low, TRUE ~ high)
  else
    case_when(!is.finite(x) ~ NA_character_, x < t ~ low, TRUE ~ high)
}

# Three-level generalisation of assign_2season. The same lower_closed
# convention applies at both boundaries so polarity is handled identically
# regardless of k. Returns all-NA if thresholds are non-finite or degenerate
# (t1 >= t2), which cascades to a screened-out candidate rather than silent
# mis-labelling.
assign_3season <- function(x, t1, t2, low = "Low", mid = "Mid",
                           high = "High", lower_closed = FALSE) {
  if (!is.finite(t1) || !is.finite(t2)) return(rep(NA_character_, length(x)))
  if (t1 >= t2) {
    warning(sprintf("assign_3season: t1 (%.4f) >= t2 (%.4f); degenerate.", t1, t2))
    return(rep(NA_character_, length(x)))
  }
  if (lower_closed)
    case_when(!is.finite(x) ~ NA_character_,
              x <= t1 ~ low, x <= t2 ~ mid, TRUE ~ high)
  else
    case_when(!is.finite(x) ~ NA_character_,
              x < t1 ~ low, x < t2 ~ mid, TRUE ~ high)
}

# Enforces a 24-observation minimum before computing quantile thresholds.
# Below 24 values the empirical quantile is too noisy to anchor a stable
# boundary; returning NA cascades to a degenerate candidate rather than a
# spurious-but-finite threshold that could survive screening.
get_q <- function(x, probs, type = 8) {
  x <- x[is.finite(x)]
  if (length(x) < 24) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = type))
}

# Rounds thresholds for CSV display only — not used for classification.
# Suppresses near-zero floating-point noise and adapts decimal places to
# magnitude, so the output CSV is human-readable without false precision.
clean_num <- function(x) {
  case_when(is.na(x) ~ NA_real_,
            abs(x) < 0.001 ~ 0,
            abs(x) >= 10   ~ round(x, 0),
            TRUE            ~ round(x, 3))
}

# =============================================================================
# 3. CANDIDATE GENERATION — cross driver × k × method to produce all label series
# =============================================================================
# For each driver × k × method, assign a season label to every month.
# Polarity-aware: for inverted drivers (high = Dry); or natural (high = Wet).

build_candidate <- function(df, driver, k_seasons = 2,
                            method = c("std", "quantile"),
                            baseline_df = NULL, std_list = NULL) {
  method <- match.arg(method)
  x  <- df[[driver]]
  dm <- driver_info(driver)
  labels <- if (k_seasons == 2) c(dm$label_low, dm$label_high)
  else c(dm$label_low, dm$label_mid, dm$label_high)
  meta <- list(driver = driver, k = k_seasons, method = method,
               t1 = NA_real_, t2 = NA_real_, season_levels = labels)
  if (method == "std") {
    drv_std <- std_list[[driver]]
    if (is.null(drv_std))
      return(list(season = rep(NA_character_, nrow(df)), meta = meta))
    if (k_seasons == 2) {
      meta$t1 <- suppressWarnings(as.numeric(drv_std$two$t))
      out <- assign_2season(x, t = meta$t1, low = labels[1], high = labels[2],
        lower_closed = dm$high_is_dry)
    } else {
      meta$t1 <- suppressWarnings(as.numeric(drv_std$three$t1))
      meta$t2 <- suppressWarnings(as.numeric(drv_std$three$t2))
      out <- assign_3season(x, meta$t1, meta$t2,
        low = labels[1], mid = labels[2], high = labels[3],
        lower_closed = dm$high_is_dry)
    }
    return(list(season = out, meta = meta))
  }
  if (is.null(baseline_df)) stop("baseline_df required for quantile thresholds.")
  xb <- baseline_df[[driver]]
  xb <- xb[is.finite(xb)]
  if (k_seasons == 2) {
    if (length(xb) < 24) {
      meta$t1 <- NA_real_
    } else {
      # For high_is_dry drivers (e.g. CWD), exact zeros represent "no deficit"
      # and cluster at the floor of the distribution.  Replacing 0 with a small
      # positive value (1e-6) prevents the median from collapsing to 0 when
      # >50% of baseline months are zero, while still placing the threshold at
      # effectively 0 for practical classification purposes.
      # NOTE — zero-jitter asymmetry: k=2 uses zero-jittered xb_q for the
      # quantile split, but k=3 (high_is_dry branch below) uses raw xb for t1
      # (median) and restricts t2 to strictly positive values only. This is
      # intentional: the k=3 design places t1 at the natural zero boundary
      # rather than the Q_SPLIT_2S percentile, giving a distinct biological
      # interpretation (Wet = zero-deficit; Transition/Dry = positive deficit).
      xb_q <- if (dm$high_is_dry) ifelse(xb == 0, 1e-6, xb) else xb
      meta$t1 <- suppressWarnings(as.numeric(quantile(xb_q, Q_SPLIT_2S, na.rm = TRUE)))
    }
    out <- assign_2season(x, t = meta$t1, low = labels[1], high = labels[2],
                          lower_closed = dm$high_is_dry)
  } else {
    if (dm$high_is_dry) {
      # t1 = median of all baseline values (including zeros, which represent
      # "no deficit" months and correctly anchor the lower boundary).
      # t2 = Q_HID_T2 percentile of strictly positive baseline values only; this
      # places the Dry/Transition boundary within the deficit-bearing months,
      # separating moderate from high water-stress months.  Q_HID_T2 defaults to
      # 0.66, matching the ~2/3 split used for the non-high_is_dry case (tertiles).
      # Adjust Q_HID_T2 in config with caution and document any deviation.
      xb_pos <- xb[xb > 0]
      meta$t1 <- if (length(xb) >= 24)
        suppressWarnings(as.numeric(quantile(xb, 0.50, na.rm = TRUE))) else NA_real_
      meta$t2 <- if (length(xb_pos) >= 24)
        suppressWarnings(as.numeric(quantile(xb_pos, Q_HID_T2, na.rm = TRUE))) else NA_real_
    } else {
      # t1/t2 = tertiles (33rd/67th percentiles) of baseline: equal-sized bins.
      q <- get_q(xb, probs = Q_SPLIT_3S)
      meta$t1 <- q[1]; meta$t2 <- q[2]
    }
    out <- assign_3season(x, meta$t1, meta$t2,
                          low = labels[1], mid = labels[2], high = labels[3],
                          lower_closed = dm$high_is_dry)
  }
  list(season = out, meta = meta)
}

# Enumerate all candidates (drivers × 2 k-levels × 2 methods)
candidate_grid <- expand_grid(
  driver = DRIVER_META$driver, k = c(2, 3), method = c("std", "quantile")) %>%
  mutate(candidate_id = paste0(driver, "_", k, "S_", method))

# Build each candidate and bind with time keys
candidate_results <- pmap(candidate_grid, function(driver, k, method, candidate_id) {
  res <- build_candidate(monthly_clim, driver, k, method,
                         baseline_months, STD_THRESHOLDS)
  tibble(candidate_id = candidate_id,
         Year = monthly_clim$Year, Month = monthly_clim$Month,
         season = factor(res$season, levels = res$meta$season_levels, exclude = NULL),
         season_levels  = paste(res$meta$season_levels, collapse = "|"),
         driver = res$meta$driver, k = res$meta$k, method = res$meta$method,
         t1 = res$meta$t1, t2 = res$meta$t2)
})

# Assemble long-format table with season labels joined to climate data
season_long <- bind_rows(candidate_results) %>%
  left_join(monthly_clim, by = c("Year", "Month")) %>%
  mutate(DateMonth = as.Date(sprintf("%d-%02d-15", Year, Month))) %>%
  relocate(DateMonth, Year, Month, candidate_id, driver, k, method, season)

# Threshold metadata table
threshold_tbl <- season_long %>%
  distinct(candidate_id, driver, k, method, t1, t2) %>%
  mutate(across(c(t1, t2), clean_num))

# =============================================================================
# 4. STRUCTURAL DIAGNOSTICS — temporal persistence, switching rate, calendar
#    consistency, and class balance over the full climate record
# =============================================================================
# Four complementary metrics per candidate over the full record:
#   (a) Temporal persistence  — median run length
#   (b) Switching rate        — mean seasonal transitions per year
#   (c) Calendar consistency  — fraction of years each month maps to its
#                               dominant season (mean and min across months)
#   (d) Class balance         — minimum bin proportion and normalised entropy

# Long unbroken runs indicate that the season definition does not oscillate
# on a sub-seasonal timescale — a necessary property for ecological relevance.
run_metrics <- function(season_vec) {
  s <- as.character(season_vec)
  s <- s[!is.na(s)]
  if (length(s) < 12) return(NA_real_)
  median(rle(s)$lengths)
}

# High switching per year suggests the threshold lies in a region of high
# driver density and produces noisy, rapidly alternating labels — a sign
# the boundary is not ecologically meaningful.
switch_metrics <- function(df, season_col) {
  d <- df %>% arrange(DateMonth)
  if (nrow(d) == 0) return(NA_real_)
  d %>% group_by(Year) %>%
    summarise(n_switch = sum(!is.na(.data[[season_col]]) &
                       !is.na(lag(.data[[season_col]])) &
                       .data[[season_col]] != lag(.data[[season_col]])),
      .groups = "drop") %>%
    summarise(mean_switch_per_year = mean(n_switch, na.rm = TRUE)) %>%
    pull(mean_switch_per_year)
}

# A coherent season definition should map each calendar month to the same
# label across most years. Low consistency implies the boundary is not
# aligned with the seasonal cycle, so the classification captures noise
# rather than a repeatable phenological pattern.
month_consistency <- function(df, season_col) {
  d <- df %>% filter(!is.na(.data[[season_col]]))
  if (nrow(d) == 0) {
    return(list(mean_month_consistency = NA_real_, min_month_consistency  = NA_real_))
  }
  mm <- d %>%
    count(Month, Year, .data[[season_col]]) %>%
    group_by(Month, Year) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup() %>%
    count(Month, .data[[season_col]]) %>%
    group_by(Month) %>%
    mutate(prop = n / sum(n)) %>%
    summarise(max_prop = max(prop), .groups = "drop")
  list(mean_month_consistency = mean(mm$max_prop, na.rm = TRUE),
    min_month_consistency  = min(mm$max_prop, na.rm = TRUE))
}

# Checks that no season level dominates to near-totality. A heavily
# imbalanced candidate inflates agreement metrics on the majority class
# and provides little discriminatory information for ecological analysis.
balance_metrics <- function(season_vec) {
  tab <- table(droplevels(season_vec))
  if (sum(tab) == 0) {
    return(list(min_bin_prop = NA_real_,
                n_levels_used = 0L,
                entropy_norm = NA_real_))
  }
  p <- as.numeric(tab) / sum(tab)
  H_max <- log(length(tab))
  list(min_bin_prop = min(p), n_levels_used = length(tab),
    entropy_norm = if (H_max > 0) -sum(p * log(p)) / H_max else NA_real_)
}

# Compute all metrics per candidate
validation_tbl <- season_long %>%
  group_by(candidate_id, driver, k, method) %>%
  nest() %>%
  mutate(
    n_months     = map_int(data, nrow),
    n_assigned   = map_int(data, ~sum(!is.na(.x$season))),
    pct_assigned = (n_assigned / n_months) * 100,
    med_run      = map_dbl(data, ~run_metrics(.x$season)),
    mean_switch_per_year = map_dbl(data, ~switch_metrics(.x, "season")),
    bal  = map(data, ~balance_metrics(.x$season)),
    cal  = map(data, ~month_consistency(.x, "season")),
    min_bin_prop = map_dbl(bal, "min_bin_prop"),
    n_levels_used = map_int(bal, "n_levels_used"),
    entropy_norm = map_dbl(bal, "entropy_norm"),
    mean_month_consistency = map_dbl(cal, "mean_month_consistency"),
    min_month_consistency = map_dbl(cal, "min_month_consistency")
  ) %>%
  dplyr::select(-data, -bal, -cal) %>%
  ungroup()

# =============================================================================
# 5. DEGENERACY FILTER — remove candidates that fail minimum structural quality
# =============================================================================
# Remove candidates with < 90% months assigned, fewer than k levels
# realised, or smallest bin below minimum count.

screened_tbl <- validation_tbl %>%
  mutate(
    min_bin_n     = round(min_bin_prop * n_assigned),
    min_bin_n_req = if_else(k == 2, S1_MIN_BIN_N_2S, S1_MIN_BIN_N_3S)) %>%
  filter(pct_assigned >= S1_MIN_PCT_ASSIGNED,
         n_levels_used >= k,
         min_bin_n >= min_bin_n_req) %>%
  left_join(threshold_tbl, by = c("candidate_id", "driver", "k", "method"))

# Stop if all candidates fail screening
if (nrow(screened_tbl) == 0) {
  fail_summary <- validation_tbl %>%
    mutate(
      min_bin_n     = round(min_bin_prop * n_assigned),
      min_bin_n_req = if_else(k == 2, S1_MIN_BIN_N_2S, S1_MIN_BIN_N_3S),
      fail_pct      = pct_assigned < S1_MIN_PCT_ASSIGNED,
      fail_levels   = n_levels_used < k,
      fail_bin      = min_bin_n < min_bin_n_req) %>%
    summarise(
      n_fail_pct    = sum(fail_pct),
      n_fail_levels = sum(fail_levels),
      n_fail_bin    = sum(fail_bin))
  stop(sprintf(
    "Stage 1 screening removed all %d candidate(s). ",
    nrow(validation_tbl)),
    sprintf("Failures: pct_assigned=%d, n_levels=%d, min_bin_n=%d. ",
            fail_summary$n_fail_pct,
            fail_summary$n_fail_levels,
            fail_summary$n_fail_bin),
    "Review DRIVER_META, STD_THRESHOLDS, and BASELINE_START/BASELINE_END.")
}

# Log removed candidates with failure reasons
removed_tbl <- validation_tbl %>%
  mutate(
    min_bin_n     = round(min_bin_prop * n_assigned),
    min_bin_n_req = if_else(k == 2, S1_MIN_BIN_N_2S, S1_MIN_BIN_N_3S),
    fail_pct      = pct_assigned < S1_MIN_PCT_ASSIGNED,
    fail_levels   = n_levels_used < k,
    fail_bin      = min_bin_n < min_bin_n_req) %>%
  filter(fail_pct | fail_levels | fail_bin) %>%
  dplyr::select(candidate_id, driver, k, method, min_bin_prop, min_bin_n,
         fail_pct, fail_levels, fail_bin)

message("Stage 1: ", nrow(screened_tbl), " candidates pass screening (",
        nrow(removed_tbl), " removed).")

# =============================================================================
# 6. OUTPUTS — write diagnostic tables and RDS handshakes for Stages 2–4
# =============================================================================

write.csv(season_long %>%
            dplyr::select(DateMonth, Year, Month,
                          candidate_id, driver, n_seasons = k, method, season,
                          everything()),
          file.path(tab_dir, "season_long_data.csv"), row.names = FALSE)

write.csv(threshold_tbl %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          threshold_1 = t1, threshold_2 = t2),
          file.path(tab_dir, "thresholds.csv"), row.names = FALSE)

write.csv(validation_tbl %>%
            mutate(min_bin_n = round(min_bin_prop * n_assigned)) %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          n_months, n_assigned, pct_assigned, min_bin_n,
                          n_levels_used, min_bin_prop, entropy_norm,
                          mean_month_consistency, min_month_consistency,
                          mean_switch_per_year, med_run),
          file.path(tab_dir, "validation_tbl.csv"), row.names = FALSE)

write.csv(removed_tbl %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          min_bin_prop, min_bin_n,
                          fail_pct, fail_levels, fail_bin),
          file.path(tab_dir, "screened_out_candidates.csv"), row.names = FALSE)

# RDS objects for downstream stages
saveRDS(screened_tbl,   file.path(output_dir, "screened_tbl.rds"))
saveRDS(threshold_tbl,  file.path(output_dir, "threshold_tbl.rds"))
saveRDS(season_long,    file.path(output_dir, "season_long.rds"))

# Session info — saved in every stage output directory so that each stage's
# package environment is independently recorded. Stage 1 always runs first;
# all subsequent stages also write session_info.txt to their own output_dir.
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
# =============================================================================
