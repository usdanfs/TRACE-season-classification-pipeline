# =============================================================================
# config_climate_only.R — Pipeline Configuration (Climate-Only Version)
# =============================================================================
# Central configuration for the climate-only seasonal classification pipeline.
# This version does NOT require ecological response data.
#
# To use a different config file without editing scripts:
#   Sys.setenv(SEASON_CONFIG = "3STAGE/config_climate_only.R")
#
# Structure:
#   SECTION 1 — SITE SETTINGS
#   SECTION 2 — METHOD SETTINGS
#   SECTION 3 — ADVANCED SETTINGS
# =============================================================================

PROJECT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = TRUE),
  error = function(e) normalizePath(".", winslash = "/", mustWork = TRUE))

DIR_STAGE_1 <- "output_STAGE_1_climate_only_candidates"
DIR_STAGE_2 <- "output_STAGE_2_climate_only_validation"
DIR_STAGE_3 <- "output_STAGE_3_climate_only_ranking"

# =============================================================================
# SECTION 1 — SITE SETTINGS
# =============================================================================

CLIMATE_CSV <- file.path(PROJECT_DIR, "CLIM.csv")

DRIVER_META <- data.frame(
  driver      = c("DRIVER_1", "DRIVER_2", "DRIVER_3"),
  high_is_dry = c(TRUE,        FALSE,       TRUE),
  label_low   = c("Wet",       "Dry",       "Wet"),
  label_high  = c("Dry",       "Wet",       "Dry"),
  label_mid   = c("Transition", "Transition", "Transition"),
  stringsAsFactors = FALSE)

BASELINE_START <- 1990
BASELINE_END   <- 2020

STD_THRESHOLDS <- list(
  DRIVER_1 = list(two = list(t = 0.0), three = list(t1 = 0.0, t2 = 0.0)),
  DRIVER_2 = list(two = list(t = 0.0), three = list(t1 = 0.0, t2 = 0.0)),
  DRIVER_3 = list(two = list(t = 0.0), three = list(t1 = 0.0, t2 = 0.0)))

VALIDATION_START <- NA
VALIDATION_END   <- NA

W_CLIMATE <- 0.60
W_ROBUST  <- 0.40

# =============================================================================
# SECTION 2 — METHOD SETTINGS
# =============================================================================

Q_SPLIT_2S  <- 0.50
Q_SPLIT_3S  <- c(1/3, 2/3)
BOOT_N_RANK <- 300

S1_MIN_PCT_ASSIGNED <- 90
S1_MIN_BIN_N_2S     <- 24L
S1_MIN_BIN_N_3S     <- 18L

# =============================================================================
# SECTION 3 — ADVANCED SETTINGS
# =============================================================================

S2_MIN_PCT_ASSIGNED   <- 0.90
S2_MIN_SEASON_PROP    <- 0.10
S2_MEAN_MONTH_CONS    <- 0.55
S2_MIN_BLOCK_PROP     <- 0.05
S2_MIN_HEALTHY_BLOCKS <- 0.50
S2_BLOCK_YEARS        <- 2

SENS_W_CLIMATE_RANGE      <- c(0.40, 0.80)
SENS_W_ROBUST_RANGE       <- c(0.20, 0.60)
SENS_W_STEP               <- 0.10
SENS_W_WINNER_CHANGE_PCT  <- 25

Q_HID_T2 <- 0.66

GLOBAL_SEED <- 123

# =============================================================================
# DERIVED OBJECTS (do not edit below this line)
# =============================================================================

stage_dir <- function(stage) {
  d <- switch(as.character(stage),
              "1" = DIR_STAGE_1, "2" = DIR_STAGE_2, "3" = DIR_STAGE_3,
              stop("Unknown stage: ", stage))
  file.path(PROJECT_DIR, d)
}

driver_info <- function(drv) {
  row <- DRIVER_META[DRIVER_META$driver == drv, ]
  if (nrow(row) == 0) stop("Driver '", drv, "' not found in DRIVER_META")
  as.list(row)
}

stopifnot(abs(W_CLIMATE + W_ROBUST - 1.0) < 1e-6)

for (.drv in names(STD_THRESHOLDS)) {
  .t3 <- STD_THRESHOLDS[[.drv]][["three"]]
  if (!is.null(.t3) && is.finite(.t3$t1) && is.finite(.t3$t2) && .t3$t1 >= .t3$t2)
    stop(sprintf(
      "STD_THRESHOLDS validation: t1 (%.4f) must be < t2 (%.4f) for driver '%s' k=3.",
      .t3$t1, .t3$t2, .drv))
}
rm(.drv, .t3)
# =============================================================================
