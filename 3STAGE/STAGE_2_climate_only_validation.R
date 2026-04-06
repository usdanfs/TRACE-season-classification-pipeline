# =============================================================================
# STAGE 2 — Climate-Only Validation
# =============================================================================
# Stress-test Stage 1 candidates for temporal stability using only structural
# criteria. This is the climate-only analogue of the full-pipeline Stage 3.
#
# Inputs:
#   Stage 1: season_long.rds, screened_tbl.rds
#   Climate: CLIMATE_CSV
#
# Outputs (output_dir/tables/):
#   - block_stability.csv
#   - filter_results.csv
#   - retained_candidates.csv
#   - .rds objects for downstream ranking
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# Default path uses the 3STAGE/ prefix so the script works when run from the
# project root directory. When run from inside 3STAGE/, set SEASON_CONFIG
# explicitly or the file will still be found via the relative path.
CONFIG_FILE <- Sys.getenv("SEASON_CONFIG", unset = "3STAGE/config_climate_only.R")
source(CONFIG_FILE)
set.seed(GLOBAL_SEED)

output_dir <- stage_dir(2)
tab_dir    <- file.path(output_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

season_long  <- readRDS(file.path(stage_dir(1), "season_long.rds"))
screened_tbl <- readRDS(file.path(stage_dir(1), "screened_tbl.rds"))

monthly_clim <- read.csv(CLIMATE_CSV, stringsAsFactors = FALSE) %>%
  mutate(Year      = as.integer(Year),
         Month     = as.integer(Month),
         DateMonth = as.Date(sprintf("%d-%02d-15", Year, Month)))

# Optional validation window: if VALIDATION_START / VALIDATION_END exist in
# config.R and are finite, restrict to that period; otherwise use all months.
use_validation_window <- exists("VALIDATION_START", inherits = TRUE) &&
  exists("VALIDATION_END", inherits = TRUE) &&
  is.finite(VALIDATION_START) && is.finite(VALIDATION_END)

validation_months <- if (use_validation_window) {
  monthly_clim %>%
    filter(Year >= VALIDATION_START, Year <= VALIDATION_END) %>%
    distinct(DateMonth)
} else {
  monthly_clim %>% distinct(DateMonth)
}

# Warn when the validation window is too short for reliable block stability testing.
# With fewer than 2 × S2_BLOCK_YEARS years, at most one block exists; a single block
# cannot distinguish temporal instability from genuine absence of a season level.
n_validation_years <- n_distinct(year(validation_months$DateMonth))
if (n_validation_years < 2L * S2_BLOCK_YEARS)
  warning(sprintf(
    "Validation window spans only %d year(s); reliable block stability requires >= %d years (%d blocks x %d yr). ",
    n_validation_years, 2L * S2_BLOCK_YEARS, 2L, S2_BLOCK_YEARS),
    "Block-collapse drop rules may be overly sensitive. Consider extending the climate ",
    "record or increasing S2_BLOCK_YEARS in config.")

# =============================================================================
# 2. RESTRICT TO VALIDATION WINDOW
# =============================================================================

stage1_long <- season_long %>%
  semi_join(validation_months, by = "DateMonth") %>%
  dplyr::select(candidate_id, driver, k, method, DateMonth, season)

# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

balance_metrics <- function(season_vec) {
  tab <- table(droplevels(season_vec))
  if (sum(tab) == 0)
    return(tibble(min_bin_prop = NA_real_, min_bin_n = 0L,
                  n_levels_used = 0L))
  tibble(min_bin_prop = min(as.numeric(tab) / sum(tab)),
         min_bin_n = as.integer(min(tab)),
         n_levels_used = length(tab))
}

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

make_year_blocks <- function(dates, block_years = 2) {
  yrs <- sort(unique(year(dates)))
  blocks <- list(); i <- 1
  while (i <= length(yrs)) {
    blocks[[length(blocks) + 1]] <- yrs[i:min(i + block_years - 1, length(yrs))]
    i <- i + block_years
  }
  blocks
}

# =============================================================================
# 4. STRUCTURAL STRESS TESTS
# =============================================================================

year_blocks <- make_year_blocks(validation_months$DateMonth, block_years = S2_BLOCK_YEARS)

overview_val <- stage1_long %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    n_months         = n(),
    n_assigned       = sum(!is.na(season)),
    pct_assigned_val = n_assigned / n_months,
    n_levels_val     = n_distinct(season[!is.na(season)]),
    .groups = "drop")

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

calendar_consistency_val <- stage1_long %>%
  mutate(Month = month(DateMonth), Year = year(DateMonth)) %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    cal = list(month_consistency(tibble(Month = Month, Year = Year, season = season))),
    .groups = "drop") %>%
  unnest_wider(cal) %>%
  rename(mean_month_consistency_val = mean_month_consistency,
         min_month_consistency_val  = min_month_consistency)

# =============================================================================
# 5. CONSERVATIVE FILTERING
# =============================================================================

block_flags <- block_stability %>%
  group_by(candidate_id, driver, k, method) %>%
  summarise(
    n_blocks    = n(),
    n_collapsed = sum(n_levels_block < first(k), na.rm = TRUE),
    prop_healthy = mean(n_levels_block >= first(k), na.rm = TRUE),
    min_block_min_bin_prop = {
      vals <- min_bin_prop_block[is.finite(min_bin_prop_block)]
      if (length(vals) == 0) NA_real_ else min(vals)
    },
    any_block_extreme_imbalance = any(min_bin_prop_block < S2_MIN_BLOCK_PROP, na.rm = TRUE),
    .groups = "drop")

filters_tbl <- overview_val %>%
  left_join(
    stage1_long %>%
      group_by(candidate_id, driver, k, method) %>%
      summarise(bal = list(balance_metrics(season)), .groups = "drop") %>%
      unnest_wider(bal) %>%
      rename(min_bin_prop_val = min_bin_prop,
             min_bin_n_val = min_bin_n,
             n_levels_used_val = n_levels_used),
    by = c("candidate_id", "driver", "k", "method")) %>%
  left_join(calendar_consistency_val, by = c("candidate_id", "driver", "k", "method")) %>%
  left_join(block_flags, by = c("candidate_id", "driver", "k", "method")) %>%
  mutate(
    fail_assignment = pct_assigned_val < S2_MIN_PCT_ASSIGNED,
    fail_imbalance  = is.finite(min_bin_prop_val) & min_bin_prop_val < S2_MIN_SEASON_PROP,
    fail_calendar   = is.finite(mean_month_consistency_val) &
      mean_month_consistency_val < S2_MEAN_MONTH_CONS,
    fail_block_collapse  = !is.na(prop_healthy) & prop_healthy < S2_MIN_HEALTHY_BLOCKS,
    fail_block_imbalance = !is.na(any_block_extreme_imbalance) & any_block_extreme_imbalance,
    drop_candidate = fail_assignment | fail_imbalance | fail_calendar |
      fail_block_collapse | fail_block_imbalance
  )

retained <- filters_tbl %>%
  filter(!drop_candidate) %>%
  dplyr::select(
    candidate_id, driver, k, method,
    pct_assigned_val, n_levels_val,
    min_bin_prop_val, min_bin_n_val,
    mean_month_consistency_val, min_month_consistency_val,
    prop_healthy, n_collapsed, min_block_min_bin_prop
  )

# =============================================================================
# 6. SAVE OUTPUTS
# =============================================================================

write.csv(block_stability %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          test_years, n_block, n_assigned_block,
                          n_levels_block, min_bin_prop_block, min_bin_n_block),
          file.path(tab_dir, "block_stability.csv"), row.names = FALSE)

write.csv(filters_tbl %>%
            dplyr::select(
              candidate_id, driver, n_seasons = k, method,
              pct_assigned_val, n_levels_val,
              min_bin_prop_val, min_bin_n_val,
              mean_month_consistency_val, min_month_consistency_val,
              prop_healthy, n_collapsed, min_block_min_bin_prop,
              any_block_extreme_imbalance,
              fail_assignment, fail_imbalance, fail_calendar,
              fail_block_collapse, fail_block_imbalance, drop_candidate),
          file.path(tab_dir, "filter_results.csv"), row.names = FALSE)

write.csv(retained %>%
            dplyr::select(
              candidate_id, driver, n_seasons = k, method,
              pct_assigned_val, n_levels_val,
              min_bin_prop_val, min_bin_n_val,
              mean_month_consistency_val, min_month_consistency_val,
              prop_healthy, n_collapsed, min_block_min_bin_prop),
          file.path(tab_dir, "retained_candidates.csv"), row.names = FALSE)

saveRDS(block_stability, file.path(output_dir, "block_stability.rds"))
saveRDS(filters_tbl,     file.path(output_dir, "filters_tbl.rds"))
saveRDS(retained,        file.path(output_dir, "stage2_stage1_candidates_retained.rds"))

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
# =============================================================================
