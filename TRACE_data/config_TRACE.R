# =============================================================================
# config_TRACE.R — Pipeline Configuration for the TRACE Site
# =============================================================================
# Seasonal classification pipeline configuration for the TRACE experimental
# site, Luquillo Experimental Forest, Puerto Rico.
#
# To use this config, set the environment variable before running any stage:
#   Sys.setenv(SEASON_CONFIG = "TRAEC_data/config_TRACE.R")
#
# Structure:
#   SECTION 1 — SITE SETTINGS  (site-specific; must be reviewed for every run)
#   SECTION 2 — METHOD SETTINGS (defaults; edit only with scientific justification)
#   SECTION 3 — ADVANCED SETTINGS (internal thresholds; do not change lightly)
# =============================================================================

PROJECT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = TRUE),
  error = function(e) normalizePath(".", winslash = "/", mustWork = TRUE))

DIR_STAGE_1 <- "output_STAGE_1"
DIR_STAGE_2 <- "output_STAGE_2"
DIR_STAGE_3 <- "output_STAGE_3"
DIR_STAGE_4 <- "output_STAGE_4"

# =============================================================================
# SECTION 1 — SITE SETTINGS
# =============================================================================

CLIMATE_CSV  <- file.path(PROJECT_DIR, "CLIM_TRACE.csv")
RESPONSE_CSV <- file.path(PROJECT_DIR, "RESPONSE_TRACE.csv")
RESPONSE_COL <- "log_flux"

DRIVER_META <- data.frame(
  driver      = c("SPEI",        "CWD",        "Rain_roll"),
  high_is_dry = c(FALSE,         TRUE,          FALSE),
  label_low   = c("Dry",         "Wet",         "Dry"),
  label_high  = c("Wet",         "Dry",         "Wet"),
  label_mid   = c("Transition",  "Transition",  "Transition"),
  stringsAsFactors = FALSE)

BASELINE_START <- 1995
BASELINE_END   <- 2015

STD_THRESHOLDS <- list(
  SPEI = list(
    two   = list(t = 0),
    three = list(t1 = -1, t2 = 1)),
  CWD = list(
    two   = list(t = 0),
    three = list(t1 = 0, t2 = 20)),
  Rain_roll = list(
    two   = list(t = 300),
    three = list(t1 = 180, t2 = 300)))

SEG_DRIVERS <- data.frame(
  driver = c("Rain_roll", "CWD",  "SPEI"),
  psi1   = c(300,          10,    -1),
  psi2   = c(600,          30,     0),
  stringsAsFactors = FALSE)

W_CLIMATE <- 0.50
W_ROBUST  <- 0.30
W_VERIFY  <- 0.20

# =============================================================================
# SECTION 2 — METHOD SETTINGS
# =============================================================================

Q_SPLIT_2S   <- 0.50
Q_SPLIT_3S   <- c(1/3, 2/3)
Q_HID_T2     <- 0.66
DAVIES_ALPHA <- 0.05

MIN_MONTHS_FOR_SEG <- 45
BOOT_B_SEG         <- 300
MIN_DELTA_AIC      <- 2
BOOT_N_RANK        <- 300

S1_MIN_PCT_ASSIGNED <- 90
S1_MIN_BIN_N_2S     <- 24L
S1_MIN_BIN_N_3S     <- 18L

# =============================================================================
# SECTION 3 — ADVANCED SETTINGS
# =============================================================================

DAVIES_K <- 10

S3_MIN_PCT_ASSIGNED   <- 0.90
S3_MIN_SEASON_PROP    <- 0.10
S3_MEAN_MONTH_CONS    <- 0.55
S3_MIN_BLOCK_PROP     <- 0.05
S3_MIN_HEALTHY_BLOCKS <- 0.50
S3_BLOCK_YEARS        <- 2

S3_FLAG_KAPPA_LOW    <- 0.10
S3_FLAG_ALIGN_IQR    <- 1.5
S3_FLAG_OMEGA_SQ_LOW <- 0.01

S4_NEAR_CONSTANT_THRESHOLD <- 0.95

SENS_W_CLIMATE_RANGE      <- c(0.30, 0.70)
SENS_W_ROBUST_RANGE       <- c(0.10, 0.40)
SENS_W_VERIFY_RANGE       <- c(0.10, 0.40)
SENS_W_STEP               <- 0.10
SENS_W_WINNER_CHANGE_PCT  <- 25

GLOBAL_SEED <- 123

# =============================================================================
# DERIVED OBJECTS (do not edit below this line)
# =============================================================================

stage_dir <- function(stage) {
  d <- switch(as.character(stage),
              "1" = DIR_STAGE_1, "2" = DIR_STAGE_2,
              "3" = DIR_STAGE_3, "4" = DIR_STAGE_4,
              stop("Unknown stage: ", stage))
  file.path(PROJECT_DIR, d)
}

driver_info <- function(drv) {
  row <- DRIVER_META[DRIVER_META$driver == drv, ]
  if (nrow(row) == 0) stop("Driver '", drv, "' not found in DRIVER_META")
  as.list(row)
}

stopifnot(abs(W_CLIMATE + W_ROBUST + W_VERIFY - 1.0) < 1e-6)

for (.drv in names(STD_THRESHOLDS)) {
  .t3 <- STD_THRESHOLDS[[.drv]][["three"]]
  if (!is.null(.t3) && is.finite(.t3$t1) && is.finite(.t3$t2) && .t3$t1 >= .t3$t2)
    stop(sprintf(
      "STD_THRESHOLDS validation: t1 (%.4f) must be < t2 (%.4f) for driver '%s' k=3.",
      .t3$t1, .t3$t2, .drv))
}
rm(.drv, .t3)
# =============================================================================
