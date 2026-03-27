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
# 1. LOAD CLIMATE DATA
# =============================================================================
# Input CSV must contain Year, Month, and all DRIVER_META$driver columns.
# Computed indices (e.g. SPEI, rolling sums, CWD) should be pre-derived.

monthly_clim <- read.csv(CLIMATE_CSV, stringsAsFactors = FALSE) %>%
  mutate(Year  = as.integer(Year), Month = as.integer(Month),
         DateMonth = as.Date(sprintf("%d-%02d-15", Year, Month)))

missing_drv <- setdiff(DRIVER_META$driver, names(monthly_clim))
if (length(missing_drv) > 0)
  stop("Drivers missing from climate CSV: ", paste(missing_drv, collapse = ", "))

baseline_months <- monthly_clim %>%
  filter(Year >= BASELINE_START, Year <= BASELINE_END)

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

# Assign two-level season labels from a single threshold
assign_2season <- function(x, t, low = "Low", high = "High",
                           lower_closed = FALSE) {
  if (lower_closed)
    case_when(!is.finite(x) ~ NA_character_, x <= t ~ low, TRUE ~ high)
  else
    case_when(!is.finite(x) ~ NA_character_, x < t ~ low, TRUE ~ high)
}

# Assign three-level season labels from two thresholds
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

# Robust quantiles requiring ≥ 24 finite values
get_q <- function(x, probs, type = 8) {
  x <- x[is.finite(x)]
  if (length(x) < 24) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = type))
}

# Round for display: suppress near-zero noise, adaptive precision
clean_num <- function(x) {
  case_when(is.na(x) ~ NA_real_,
            abs(x) < 0.001 ~ 0,
            abs(x) >= 10   ~ round(x, 0),
            TRUE            ~ round(x, 3))
}

# =============================================================================
# 3. BUILD SEASON CANDIDATES
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
      xb_q <- if (dm$high_is_dry) ifelse(xb == 0, 1e-6, xb) else xb
      meta$t1 <- suppressWarnings(as.numeric(quantile(xb_q, 0.5, na.rm = TRUE)))
    }
    out <- assign_2season(x, t = meta$t1, low = labels[1], high = labels[2])
  } else {
    if (dm$high_is_dry) {
      xb_pos <- xb[xb > 0]
      meta$t1 <- if (length(xb) >= 24)
        suppressWarnings(as.numeric(quantile(xb, 0.50, na.rm = TRUE))) else NA_real_
      meta$t2 <- if (length(xb_pos) >= 24)
        suppressWarnings(as.numeric(quantile(xb_pos, 0.66, na.rm = TRUE))) else NA_real_
    } else {
      q <- get_q(xb, probs = c(1/3, 2/3))
      meta$t1 <- q[1]; meta$t2 <- q[2]
    }
    out <- assign_3season(x, meta$t1, meta$t2,
                          low = labels[1], mid = labels[2], high = labels[3])
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
# 4. CLIMATE-STRUCTURE VALIDATION METRICS
# =============================================================================
# Four complementary metrics per candidate over the full record:
#   (a) Temporal persistence  — median run length
#   (b) Switching rate        — mean seasonal transitions per year
#   (c) Calendar consistency  — fraction of years each month maps to its
#                               dominant season (mean and min across months)
#   (d) Class balance         — minimum bin proportion and normalised entropy

run_metrics <- function(season_vec) {
  s <- as.character(season_vec)
  s <- s[!is.na(s)]
  if (length(s) < 12) return(NA_real_)
  median(rle(s)$lengths)
}

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
# 5. DEGENERACY SCREENING
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
  stop("Stage 1 screening removed all candidates. Review DRIVER_META, STD_THRESHOLDS, and the baseline period.")
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
# 6. SAVE OUTPUTS
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
# =============================================================================
