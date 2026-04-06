# =============================================================================
# config.R — Pipeline Configuration (Full 4-Stage Version)
# =============================================================================
# Central configuration for the seasonal classification pipeline.
# Edit this file for your site and dataset. All stage scripts source it.
#
# To use a different config file without editing scripts, set the environment
# variable SEASON_CONFIG before running any stage, e.g.:
#   Sys.setenv(SEASON_CONFIG = "config_TRACE.R")
#
# Structure:
#   SECTION 1 — SITE SETTINGS  (must be reviewed and edited for every new site)
#   SECTION 2 — METHOD SETTINGS (well-tested defaults; edit only if you have a
#               specific scientific reason to deviate from them)
#   SECTION 3 — ADVANCED SETTINGS (internal thresholds; do not change unless
#               you understand the downstream effects)
# =============================================================================

PROJECT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = TRUE),
  error = function(e) normalizePath(".", winslash = "/", mustWork = TRUE)
)

DIR_STAGE_1 <- "output_STAGE_1"
DIR_STAGE_2 <- "output_STAGE_2"
DIR_STAGE_3 <- "output_STAGE_3"
DIR_STAGE_4 <- "output_STAGE_4"

# =============================================================================
# SECTION 1 — SITE SETTINGS
# Review and edit all items in this section before running the pipeline at a
# new site. Every item here is site-specific and has no universal default.
# =============================================================================

# ---- Input files ------------------------------------------------------------
# Monthly climate CSV: must contain Year, Month, and the driver columns below.
CLIMATE_CSV  <- file.path(PROJECT_DIR, "CLIM.csv")
# Monthly ecological response CSV: must contain Year, Month, and RESPONSE_COL.
RESPONSE_CSV <- file.path(PROJECT_DIR, "RESPONSE.csv")
# Name of the response variable column in RESPONSE_CSV.
RESPONSE_COL <- "response_variable"

# ---- Climate driver metadata ------------------------------------------------
# One row per candidate driver. Controls polarity and season labelling
# throughout all stages.
#
# driver      : column name in CLIMATE_CSV (case-sensitive).
# high_is_dry : TRUE  if high values indicate drier conditions (e.g. CWD, VPD).
#               FALSE if high values indicate wetter conditions (e.g. rainfall).
# label_low   : season label assigned to low driver values.
# label_high  : season label assigned to high driver values.
# label_mid   : season label for the intermediate bin (k = 3 only).
#
# POLARITY RULE — must be consistent within each row:
#   high_is_dry = TRUE  → label_low  should be the wet/high-resource season
#                        → label_high should be the dry/low-resource season
#   high_is_dry = FALSE → label_low  should be the dry/low-resource season
#                        → label_high should be the wet/high-resource season
# Inconsistent polarity (e.g. high_is_dry = TRUE but label_high = "Wet") will
# not cause an error but will produce mislabelled output seasons.
#
# Add or remove rows to match your drivers. Names must match CLIMATE_CSV,
# STD_THRESHOLDS, and SEG_DRIVERS exactly.

DRIVER_META <- data.frame(
  driver      = c("DRIVER_1", "DRIVER_2", "DRIVER_3"),
  high_is_dry = c(TRUE,        FALSE,       TRUE),
  label_low   = c("Wet",       "Dry",       "Wet"),
  label_high  = c("Dry",       "Wet",       "Dry"),
  label_mid   = c("Transition", "Transition", "Transition"),
  stringsAsFactors = FALSE)

# ---- Baseline period --------------------------------------------------------
# Years used to compute climatological thresholds (quantile method).
# Choose a period long enough to be climatologically representative (>=20 yr).
BASELINE_START <- 1990
BASELINE_END   <- 2020

# ---- Standard thresholds ----------------------------------------------------
# Literature- or expert-based thresholds for each driver and season count.
# These define the "std" candidate set; they are independent of the data.
# two$t    : single threshold separating k = 2 seasons.
# three$t1 : lower threshold for k = 3 (boundary between low and mid season).
# three$t2 : upper threshold for k = 3 (boundary between mid and high season).

STD_THRESHOLDS <- list(
  DRIVER_1 = list(
    two   = list(t = 0.0),
    three = list(t1 = 0.0, t2 = 0.0)),
  DRIVER_2 = list(
    two   = list(t = 0.0),
    three = list(t1 = 0.0, t2 = 0.0)),
  DRIVER_3 = list(
    two   = list(t = 0.0),
    three = list(t1 = 0.0, t2 = 0.0)))

# ---- Segmented regression initial breakpoints (Stage 2) ---------------------
# One row per driver x k combination used in segmented regression.
# psi1/psi2 are initial guesses; they must be plausible but need not be exact.
# The algorithm searches locally around these values.
# k = 2: only psi1 is used.  k = 3: both psi1 and psi2 are used.

SEG_DRIVERS <- data.frame(
  driver = c("DRIVER_1", "DRIVER_1", "DRIVER_2"),
  psi1   = c(0.0,         0.0,        0.0),
  psi2   = c(0.0,         0.0,        0.0),
  stringsAsFactors = FALSE)

# ---- Tier weights (Stage 4 ranking) ----------------------------------------
# Controls the relative influence of each evidence tier on the composite score.
# Tier 1 (climate structure) should dominate to avoid circular reasoning
# when the response variable is also used in downstream analysis.
# Must sum to exactly 1.0.

W_CLIMATE <- 0.50    # Tier 1: Climate structure
W_ROBUST  <- 0.30    # Tier 2: Internal robustness (std vs quantile agreement)
W_VERIFY  <- 0.20    # Tier 3: Ecological verification (Stage 1 vs Stage 2)

# =============================================================================
# SECTION 2 — METHOD SETTINGS
# These parameters have well-tested defaults. Edit only if you have a specific
# scientific reason (e.g., a very short record requires lower S1_MIN_BIN_N_*).
# =============================================================================

# ---- Quantile split points --------------------------------------------------
Q_SPLIT_2S <- 0.50          # Median split for k = 2
Q_SPLIT_3S <- c(1/3, 2/3)  # Tertile split for k = 3

# Upper percentile for the high_is_dry k = 3 quantile threshold (t2).
# Applied to strictly positive baseline values only, separating moderate from
# high water-stress months. 0.66 matches the ~2/3 split used for the
# non-high_is_dry case (tertiles). Change with caution and document if adjusted.
Q_HID_T2 <- 0.66

# ---- Davies test significance level -----------------------------------------
DAVIES_ALPHA <- 0.05

# ---- Segmented regression controls -----------------------------------------
MIN_MONTHS_FOR_SEG <- 45
BOOT_B_SEG         <- 300
MIN_DELTA_AIC      <- 2

# ---- Rank bootstrap ---------------------------------------------------------
BOOT_N_RANK <- 300

# ---- Stage 1 screening thresholds ------------------------------------------
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

# Validate that every k=3 threshold pair has t1 < t2
for (.drv in names(STD_THRESHOLDS)) {
  .t3 <- STD_THRESHOLDS[[.drv]][["three"]]
  if (!is.null(.t3) && is.finite(.t3$t1) && is.finite(.t3$t2) && .t3$t1 >= .t3$t2)
    stop(sprintf(
      "STD_THRESHOLDS validation: t1 (%.4f) must be < t2 (%.4f) for driver '%s' k=3.",
      .t3$t1, .t3$t2, .drv))
}
rm(.drv, .t3)
# =============================================================================
