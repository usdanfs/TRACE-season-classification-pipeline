# =============================================================================
# STAGE 2 — Ecological Response Segmented Regression
# =============================================================================
#   Test whether climate-derived season boundaries correspond to detectable
#   regime shifts in the ecological response via piecewise linear regression.
#   This stage does NOT filter Stage 1 candidates; it generates independent
#   ecological-regime labels from breakpoints for cross-validation in Stage 3.
#
#   For each driver × k (2 or 3), fit a segmented regression of transformed
#   response vs. driver on monthly means.  Model selection via ΔAIC (linear
#   vs. segmented).  Breakpoint stability assessed by case-resampling
#   bootstrap (B iterations).  Leave-one-year-out CV-RMSE provides an
#   additional predictive diagnostic.  Breakpoints are converted to
#   ecological-regime season labels for Stage 3 agreement testing.
#
# Inputs:
#   - config.R
#   - Monthly climate CSV (CLIMATE_CSV)
#   - Pre-aggregated ecological response CSV (RESPONSE_CSV)
#
# Outputs (output_dir/tables/):
#   - segmentation_results.csv       Breakpoints, AIC, ΔAIC, CV-RMSE, bootstrap CI
#   - ecological_seg_candidates.csv  Ecological-regime season labels per month
#   - .rds objects for downstream stages
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)   
  library(lubridate)   
  library(zoo)        
  library(segmented)
})

CONFIG_FILE <- Sys.getenv("SEASON_CONFIG", unset = "config.R")
source(CONFIG_FILE)
# Seed note: set.seed() is called once here for the global sequence, and again
# inside each map() lambda below (once per driver) so that every driver's
# bootstrap starts from the same seed regardless of evaluation order or the
# number of drivers.  This makes per-driver results reproducible in isolation
# but means all drivers share the same initial RNG state for their bootstraps.
set.seed(GLOBAL_SEED)

output_dir  <- stage_dir(2)
tab_dir     <- file.path(output_dir, "tables")
seg_drivers <- as_tibble(SEG_DRIVERS)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# Internal algorithm constants — change only with care (see comments)
BOOT_CI_MIN_REPS <- 10L  # min bootstrap replicates required for a CI to be reported
CV_MIN_TEST_N    <- 3L   # min observations in a leave-one-year-out test fold

# =============================================================================
# 1. DATA INPUTS — load the ecological response and join with climate record
# =============================================================================

monthly_clim <- read.csv(CLIMATE_CSV, stringsAsFactors = FALSE) %>%
  mutate(DateMonth = as.Date(sprintf("%d-%02d-15", Year, Month)))

# Validate RESPONSE_CSV has the expected response column before attempting rename.
# Missing column gives a cryptic tidyselect error with no pointer to config.
resp_raw <- read.csv(RESPONSE_CSV, stringsAsFactors = FALSE)
if (!RESPONSE_COL %in% names(resp_raw))
  stop("RESPONSE_COL '", RESPONSE_COL, "' not found in RESPONSE_CSV. ",
       "Available columns: ", paste(names(resp_raw), collapse = ", "))

monthly_response <- resp_raw %>%
  mutate(DateMonth = as.Date(sprintf("%d-%02d-15", Year, Month))) %>%
  rename(RESPONSE_COL = all_of(RESPONSE_COL)) %>%
  left_join(monthly_clim %>%
              dplyr::select(DateMonth, Year, Month, all_of(DRIVER_META$driver)),
            by = c("DateMonth", "Year", "Month")) %>%
  arrange(DateMonth)
rm(resp_raw)

# No date overlap means all segmented fits will receive empty data and silently
# return null_results; the ecological regime candidates will be all-NA.
overlap_n <- sum(monthly_response$DateMonth %in% monthly_clim$DateMonth &
                   is.finite(monthly_response$RESPONSE_COL))
if (overlap_n == 0)
  stop("No overlapping dates with finite response values between RESPONSE_CSV and CLIMATE_CSV. ",
       sprintf("Response spans %d--%d; Climate spans %d--%d.",
               min(monthly_response$Year, na.rm = TRUE),
               max(monthly_response$Year, na.rm = TRUE),
               min(monthly_clim$Year, na.rm = TRUE),
               max(monthly_clim$Year, na.rm = TRUE)))

# Warn if total overlapping months are below the segmented-regression minimum.
if (overlap_n < MIN_MONTHS_FOR_SEG)
  warning(sprintf(
    "Only %d months overlap between RESPONSE_CSV and CLIMATE_CSV (MIN_MONTHS_FOR_SEG = %d). ",
    overlap_n, MIN_MONTHS_FOR_SEG),
    "All segmented fits will return null results and breakpoint_supported = FALSE.")

# Validate SEG_DRIVERS: all drivers must exist in DRIVER_META (so driver_info()
# works) and must be present as columns in monthly_response (so fits can run).
missing_seg_meta <- setdiff(seg_drivers$driver, DRIVER_META$driver)
if (length(missing_seg_meta) > 0)
  stop("SEG_DRIVERS contains driver(s) not in DRIVER_META: ",
       paste(missing_seg_meta, collapse = ", "))

missing_seg_col <- setdiff(seg_drivers$driver, names(monthly_response))
if (length(missing_seg_col) > 0)
  stop("SEG_DRIVERS driver(s) not found as columns in monthly_response after joining with climate data: ",
       paste(missing_seg_col, collapse = ", "),
       ". Check that these drivers are present in CLIMATE_CSV.")

# =============================================================================
# 2. SEGMENTATION UTILITIES — label assignment, jitter, bootstrap CI, and
#    leave-one-year-out CV helpers used by the fitting functions below
# =============================================================================

# Converts ecological breakpoints to season factors using the same polarity
# convention (lower_closed via high_is_dry) as Stage 1, so Stage 2 labels
# are directly comparable to Stage 1 labels in the Stage 3 agreement tests.
season_from_thresholds <- function(df, driver, k, t1, t2 = NA_real_) {
  x  <- df[[driver]]
  dm <- driver_info(driver)
  
  if (k == 2L) {
    if (dm$high_is_dry)
      out <- case_when(!is.finite(x) ~ NA_character_,
                       x <= t1 ~ dm$label_low, TRUE ~ dm$label_high)
    else
      out <- case_when(!is.finite(x) ~ NA_character_,
                       x < t1 ~ dm$label_low, TRUE ~ dm$label_high)
    return(factor(out, levels = c(dm$label_low, dm$label_high)))
  }
  
  if (k == 3L) {
    if (!is.finite(t1) || !is.finite(t2) || t1 >= t2)
      return(factor(rep(NA_character_, nrow(df))))
    if (dm$high_is_dry)
      out <- case_when(!is.finite(x) ~ NA_character_,
                       x <= t1 ~ dm$label_low, x <= t2 ~ dm$label_mid,
                       TRUE ~ dm$label_high)
    else
      out <- case_when(!is.finite(x) ~ NA_character_,
                       x < t1 ~ dm$label_low, x < t2 ~ dm$label_mid,
                       TRUE ~ dm$label_high)
    return(factor(out, levels = c(dm$label_low, dm$label_mid, dm$label_high)))
  }
  stop("k must be 2 or 3")
}

# segmented() can fail or return degenerate breakpoints when many x-values
# are identical, because tied values create flat regions in the piecewise
# likelihood surface. Minimal random jitter (< 1e-6 × SD) breaks ties
# without meaningfully shifting any breakpoint estimate.
jitter_if_tied <- function(x, frac = 1e-6) {
  x <- as.numeric(x)
  if (!all(is.finite(x))) return(x)
  if (n_distinct(x) < length(x) * 0.5) {
    sc <- max(sd(x, na.rm = TRUE), 1)
    return(x + runif(length(x), -sc * frac, sc * frac))
  }
  x
}

# Percentile-based bootstrap CI for breakpoint stability. The percentile
# method is preferred over SE-based intervals here because the bootstrap
# distribution of breakpoints is often asymmetric; min_ok = 10 ensures
# the CI is not reported when the bootstrap produced too few valid replicates.
boot_ci <- function(x, min_ok = BOOT_CI_MIN_REPS) {
  if (length(x) >= min_ok)
    c(median(x), quantile(x, 0.025), quantile(x, 0.975))
  else
    rep(NA_real_, 3)
}

# =============================================================================
# 3. FITTING FUNCTIONS — single- and double-breakpoint segmented regression
#    with bootstrap; called once per driver in section 4 and again per
#    leave-one-year-out fold in cv_seg_rmse
# =============================================================================
# fit_seg1 (k = 2) and fit_seg2 (k = 3) are kept as separate functions rather
# than a single dispatcher so their null_result structures (which differ in the
# number of bootstrap columns) can be typed explicitly without branching logic.
# Each returns: breakpoint estimate(s), AIC comparison, bootstrap CI, the
# working data frame, and the fitted segmented object.
# When B = 0 (during cross-validation folds), bootstrap is skipped.

fit_seg1 <- function(df, xvar, yvar = "RESPONSE_COL", psi_init = NULL, B = 300) {
  d0 <- df %>% filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))
  # A near-constant response has zero variance: lm() fits with 0 residual df and
  # segmented() will either crash or return a meaningless breakpoint.
  if (nrow(d0) >= 2 && sd(d0[[yvar]], na.rm = TRUE) < .Machine$double.eps^0.5) {
    warning(sprintf("fit_seg1: '%s' has near-zero variance for driver '%s'. Skipping fit.",
                    yvar, xvar))
    d0 <- d0[0L, ]   # force null_result path via nrow < MIN_MONTHS_FOR_SEG
  }
  null_result <- list(
    ok = FALSE, b1 = NA_real_,
    aic_linear = NA_real_, aic_seg = NA_real_, delta_aic = NA_real_,
    davies_p = NA_real_,   # must mirror the success-path list structure so
                           # map_dbl(res, "davies_p") does not crash on null paths
    boot_sum = tibble(b1_med = NA_real_, b1_lo = NA_real_,
                      b1_hi = NA_real_, n_boot_ok_1 = 0L),
    df0 = d0, seg_fit = NULL)
  if (nrow(d0) < MIN_MONTHS_FOR_SEG) return(null_result)
  d0 <- d0 %>% mutate(xj = jitter_if_tied(.data[[xvar]]))
  if (is.null(psi_init) || !is.finite(psi_init))
    psi_init <- as.numeric(quantile(d0$xj, 0.5, na.rm = TRUE))
  lm0  <- lm(reformulate("xj", yvar), data = d0)
  seg0 <- tryCatch(
    suppressWarnings(segmented(lm0, seg.Z = ~xj, psi = psi_init)),
    error = function(e) NULL)
  aic_linear <- AIC(lm0)
  aic_seg    <- if (!is.null(seg0)) AIC(seg0) else NA_real_
  delta_aic  <- aic_linear - aic_seg
  # Davies test for breakpoint existence (handles nuisance parameter)
  davies_p <- tryCatch({
    dt <- segmented::davies.test(lm0, seg.Z = ~xj, k = DAVIES_K)
    dt$p.value
  }, error = function(e) NA_real_)
  b1 <- NA_real_
  if (!is.null(seg0)) {
    psi <- tryCatch(summary(seg0)$psi[, "Est."], error = function(e) numeric(0))
    if (length(psi) >= 1) b1 <- as.numeric(psi[1])
  }
  if (B == 0)
    return(list(ok = TRUE, b1 = b1, aic_linear = aic_linear, aic_seg = aic_seg,
                delta_aic = delta_aic, davies_p = davies_p,
                boot_sum = tibble(b1_med = NA_real_, b1_lo = NA_real_,
                                  b1_hi = NA_real_, n_boot_ok_1 = 0L),
                df0 = d0, seg_fit = seg0))
  # Bootstrap breakpoint stability
  boot <- replicate(B, {
    db   <- d0[sample.int(nrow(d0), replace = TRUE), ]
    lm_b <- lm(reformulate("xj", yvar), data = db)
    seg_b <- tryCatch(
      suppressWarnings(segmented(lm_b, seg.Z = ~xj, psi = psi_init)),
      error = function(e) NULL)
    if (is.null(seg_b)) return(NA_real_)
    psi_b <- tryCatch(summary(seg_b)$psi[, "Est."], error = function(e) numeric(0))
    if (length(psi_b) < 1) NA_real_ else as.numeric(psi_b[1])
  })
  b1_ok <- boot[is.finite(boot)]
  ci1   <- boot_ci(b1_ok)
  boot_sum <- tibble(b1_med = ci1[1], b1_lo = ci1[2], b1_hi = ci1[3],
                     n_boot_ok_1 = length(b1_ok))
  list(ok = TRUE, b1 = b1, aic_linear = aic_linear, aic_seg = aic_seg, davies_p = davies_p,
       delta_aic = delta_aic, boot_sum = boot_sum, df0 = d0, seg_fit = seg0)
}

fit_seg2 <- function(df, xvar, yvar = "RESPONSE_COL", psi_init = NULL, B = 300) {
  d0 <- df %>% filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))
  if (nrow(d0) >= 2 && sd(d0[[yvar]], na.rm = TRUE) < .Machine$double.eps^0.5) {
    warning(sprintf("fit_seg2: '%s' has near-zero variance for driver '%s'. Skipping fit.",
                    yvar, xvar))
    d0 <- d0[0L, ]
  }
  null_result <- list(
    ok = FALSE, b1 = NA_real_, b2 = NA_real_,
    aic_linear = NA_real_, aic_seg = NA_real_, delta_aic = NA_real_,
    davies_p = NA_real_,   # must mirror the success-path list structure so
                           # map_dbl(res, "davies_p") does not crash on null paths
    boot_sum = tibble(b1_med = NA_real_, b1_lo = NA_real_, b1_hi = NA_real_,
                      b2_med = NA_real_, b2_lo = NA_real_, b2_hi = NA_real_,
                      n_boot_ok_1 = 0L, n_boot_ok_2 = 0L),
    df0 = d0, seg_fit = NULL)
  if (nrow(d0) < MIN_MONTHS_FOR_SEG) return(null_result)
  d0 <- d0 %>% mutate(xj = jitter_if_tied(.data[[xvar]]))
  if (is.null(psi_init) || length(psi_init) != 2 || any(!is.finite(psi_init)))
    psi_init <- as.numeric(quantile(d0$xj, c(1/3, 2/3), na.rm = TRUE))
  psi_init <- sort(psi_init)
  lm0  <- lm(reformulate("xj", yvar), data = d0)
  seg0 <- tryCatch(
    suppressWarnings(segmented(lm0, seg.Z = ~xj, psi = psi_init)),
    error = function(e) NULL)
  aic_linear <- AIC(lm0)
  aic_seg    <- if (!is.null(seg0)) AIC(seg0) else NA_real_
  delta_aic  <- aic_linear - aic_seg
  # Davies test for breakpoint existence (handles nuisance parameter)
  davies_p <- tryCatch({
    dt <- segmented::davies.test(lm0, seg.Z = ~xj, k = DAVIES_K)
    dt$p.value
  }, error = function(e) NA_real_)
  b1 <- NA_real_; b2 <- NA_real_
  if (!is.null(seg0)) {
    psi <- tryCatch(summary(seg0)$psi[, "Est."], error = function(e) numeric(0))
    if (length(psi) >= 1) b1 <- as.numeric(psi[1])
    if (length(psi) >= 2) b2 <- as.numeric(psi[2])
  }
  if (B == 0)
    return(list(ok = TRUE, b1 = b1, b2 = b2,  davies_p = davies_p,
                aic_linear = aic_linear, aic_seg = aic_seg, delta_aic = delta_aic,
                boot_sum = tibble(b1_med = NA_real_, b1_lo = NA_real_, b1_hi = NA_real_,
                                  b2_med = NA_real_, b2_lo = NA_real_, b2_hi = NA_real_,
                                  n_boot_ok_1 = 0L, n_boot_ok_2 = 0L),
                df0 = d0, seg_fit = seg0))
  # Bootstrap breakpoint stability
  boot <- replicate(B, {
    db   <- d0[sample.int(nrow(d0), replace = TRUE), ]
    lm_b <- lm(reformulate("xj", yvar), data = db)
    seg_b <- tryCatch(
      suppressWarnings(segmented(lm_b, seg.Z = ~xj, psi = psi_init)),
      error = function(e) NULL)
    if (is.null(seg_b)) return(c(NA_real_, NA_real_))
    psi_b <- tryCatch(summary(seg_b)$psi[, "Est."], error = function(e) numeric(0))
    if (length(psi_b) < 1) return(c(NA_real_, NA_real_))
    if (length(psi_b) == 1) return(c(as.numeric(psi_b[1]), NA_real_))
    c(as.numeric(psi_b[1]), as.numeric(psi_b[2]))
  })
  bmat  <- t(boot)
  b1_ok <- bmat[, 1][is.finite(bmat[, 1])]
  b2_ok <- bmat[, 2][is.finite(bmat[, 2])]
  ci1   <- boot_ci(b1_ok); ci2 <- boot_ci(b2_ok)
  boot_sum <- tibble(
    b1_med = ci1[1], b1_lo = ci1[2], b1_hi = ci1[3],
    b2_med = ci2[1], b2_lo = ci2[2], b2_hi = ci2[3],
    n_boot_ok_1 = length(b1_ok), n_boot_ok_2 = length(b2_ok))
  list(ok = TRUE, b1 = b1, b2 = b2, davies_p = davies_p,
       aic_linear = aic_linear, aic_seg = aic_seg, delta_aic = delta_aic,
       boot_sum = boot_sum, df0 = d0, seg_fit = seg0)
}

# Leave-one-year-out cross-validation RMSE: a predictive diagnostic separate
# from the in-sample AIC comparison, less sensitive to overfitting than ΔAIC.
cv_seg_rmse <- function(df, xvar, k_breaks) {
  years <- sort(unique(df$Year))
  if (length(years) < 5) return(NA_real_)
  rmses <- map_dbl(years, function(yr) {
    train <- df %>% filter(Year != yr)
    test  <- df %>% filter(Year == yr)
    if (nrow(test) < CV_MIN_TEST_N || nrow(train) < MIN_MONTHS_FOR_SEG) return(NA_real_)
    fit <- if (k_breaks == 1) fit_seg1(train, xvar, B = 0)
    else               fit_seg2(train, xvar, B = 0)
    breaks <- c(fit$b1, fit$b2)[is.finite(c(fit$b1, fit$b2))]
    if (length(breaks) < k_breaks) return(NA_real_)
    cuts  <- c(-Inf, sort(breaks), Inf)
    train <- train %>%
      mutate(s = factor(cut(.data[[xvar]], breaks = cuts, labels = FALSE),
                        levels = seq_len(k_breaks + 1)))
    test <- test %>%
      mutate(s = factor(cut(.data[[xvar]], breaks = cuts, labels = FALSE),
                        levels = seq_len(k_breaks + 1)))
    if (nlevels(droplevels(train$s)) < 2) return(NA_real_)
    m <- tryCatch(lm(RESPONSE_COL ~ s, data = train), error = function(e) NULL)
    if (is.null(m)) return(NA_real_)
    pred <- suppressWarnings(predict(m, newdata = test))
    ok   <- is.finite(pred) & is.finite(test$RESPONSE_COL)
    if (sum(ok) < 3) return(NA_real_)
    sqrt(mean((test$RESPONSE_COL[ok] - pred[ok])^2))
  })
  ok <- rmses[is.finite(rmses)]
  if (length(ok) < 3) NA_real_ else mean(ok)
}

# =============================================================================
# 4. REGRESSION RUNS — fit all driver × k combinations and extract breakpoints
# =============================================================================

# k = 2: single breakpoint per driver
# boot_sum here is a list-column extracted from each fit result; the local
# variable named boot_sum inside fit_seg1/fit_seg2 is its source value.
seg1_tbl <- seg_drivers %>%
  mutate(
    res      = map(driver, ~{ set.seed(GLOBAL_SEED); fit_seg1(monthly_response, .x, B = BOOT_B_SEG) }),
    b1       = map_dbl(res, "b1"),
    boot_sum = map(res, "boot_sum"),
    rmse_cv  = map_dbl(driver, ~cv_seg_rmse(monthly_response, .x, 1)),
    k        = 2L)

# k = 3: two breakpoints per driver
seg2_tbl <- seg_drivers %>%
  mutate(
    res      = pmap(list(driver, psi1, psi2), ~{
      set.seed(GLOBAL_SEED); fit_seg2(monthly_response, ..1,
                                      psi_init = c(..2, ..3), B = BOOT_B_SEG) }),
    b1       = map_dbl(res, "b1"),
    b2       = map_dbl(res, "b2"),
    boot_sum = map(res, "boot_sum"),
    rmse_cv  = map_dbl(driver, ~cv_seg_rmse(monthly_response, .x, 2)),
    k        = 3L)

# Unified table with AIC diagnostics
seg_tbl <- bind_rows(
  seg1_tbl %>% mutate(b2 = NA_real_),
  seg2_tbl) %>%
  mutate(
    aic_linear    = map_dbl(res, "aic_linear"),
    aic_seg       = map_dbl(res, "aic_seg"),
    delta_aic     = map_dbl(res, "delta_aic"),
    davies_p      = map_dbl(res, "davies_p"),
    pass_aic_gain = is.finite(delta_aic) & delta_aic >= MIN_DELTA_AIC,
    pass_davies   = is.finite(davies_p) & davies_p < DAVIES_ALPHA,
    breakpoint_supported = pass_aic_gain | pass_davies)

# =============================================================================
# 5. ECOLOGICAL-REGIME LABELS — convert breakpoints to month-level season
#    factors for cross-validation against Stage 1 labels in Stage 3
# =============================================================================
# Convert each breakpoint set into month-level season factors over the full
# climate record, for cross-validation against Stage 1 labels in Stage 3.

boot_summary <- seg_tbl %>%
  transmute(
    candidate_id = paste0("ecological_SEG_", driver, "_", k, "S"),
    driver, k, b1, b2,
    delta_aic, pass_aic_gain,
    davies_p, pass_davies, breakpoint_supported) %>%
  bind_cols(seg_tbl %>% pull(boot_sum) %>% bind_rows())

ecological_regime_candidates <- seg_tbl %>%
  transmute(
    candidate_id = paste0("ecological_SEG_", driver, "_", k, "S"),
    source = "segmented_monthly_mean",
    driver, k, t1 = b1, t2 = b2)

ecological_regime_long <- ecological_regime_candidates %>%
  mutate(
    data = pmap(list(driver, k, t1, t2), function(drv, kk, tt1, tt2) {
      d <- monthly_clim %>% dplyr::select(DateMonth, Year, Month,
                                          all_of(DRIVER_META$driver))
      d <- d[is.finite(d[[drv]]), ]
      d$season <- season_from_thresholds(d, driver = drv, k = kk,
                                         t1 = tt1, t2 = tt2)
      d
    })) %>%
  dplyr::select(candidate_id, driver, k, data) %>%
  unnest(data) %>%
  mutate(method = "segmented")

# =============================================================================
# 6. OUTPUTS — write breakpoint tables and RDS handshakes for Stage 3
# =============================================================================

seg_tbl_flat <- seg_tbl %>% dplyr::select(where(~!is.list(.)))

boot_ci_cols <- boot_summary %>%
  dplyr::select(driver, k,
                starts_with("b1_"), starts_with("b2_"), starts_with("n_boot_"))

write.csv(seg_tbl_flat %>%
            left_join(boot_ci_cols, by = c("driver", "k")) %>%
            dplyr::select(driver, n_seasons = k,
                          breakpoint_1 = b1, breakpoint_2 = b2,
                          delta_aic, pass_aic_gain,
                          davies_p, pass_davies, breakpoint_supported,
                          rmse_cv,
                          breakpoint_1_med = b1_med, breakpoint_1_lo = b1_lo,
                          breakpoint_1_hi = b1_hi, n_boot_ok_1,
                          breakpoint_2_med = b2_med, breakpoint_2_lo = b2_lo,
                          breakpoint_2_hi = b2_hi, n_boot_ok_2),
          file.path(tab_dir, "segmentation_results.csv"), row.names = FALSE)

write.csv(ecological_regime_long %>%
            dplyr::select(candidate_id, driver, n_seasons = k, method,
                          DateMonth, Year, Month, season, everything()),
          file.path(tab_dir, "ecological_seg_candidates.csv"), row.names = FALSE)

# RDS objects for downstream stages
saveRDS(ecological_regime_long, file.path(output_dir, "ecological_regime_long_seg.rds"))
saveRDS(seg_tbl,          file.path(output_dir, "seg_tbl.rds"))
saveRDS(monthly_response,     file.path(output_dir, "monthly_response.rds"))

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
# =============================================================================
