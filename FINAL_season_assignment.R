# =============================================================================
# FINAL — Season Assignment from Pipeline Winner
# =============================================================================
#   Read the pipeline winner from the bootstrap rank summary and produce a
#   clean, publication-ready CSV mapping every month to its assigned season
#   label.  This is the file you join to your data for downstream analysis.
#
# Compatible with both pipelines:
#   Full pipeline  → source("config.R"),  winner from Stage 4
#   Climate-only   → source("config_climate_only.R"), winner from Stage 3
#
# Output:
#   - season_assignment_final.csv
#     Columns: Year, Month, season, candidate_id
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})

# ---------------------------------------------------------------------------
# 0.  CONFIGURATION
# ---------------------------------------------------------------------------
# No editing is needed here for standard use.
# The config file is read from the SEASON_CONFIG environment variable
# (default: "config.R" for the full pipeline, or set it to
# "3STAGE/config_climate_only.R" for the climate-only pipeline).
# The pipeline version is detected automatically from the config: if DIR_STAGE_4
# is defined, the Stage 4 ranking output is used; otherwise Stage 3 is used.
CONFIG_FILE <- Sys.getenv("SEASON_CONFIG", unset = "config.R")
source(CONFIG_FILE)
set.seed(GLOBAL_SEED)

# Detect final ranking stage automatically
ranking_stage <- if (exists("DIR_STAGE_4")) 4L else 3L

# =============================================================================
# 1. IDENTIFY WINNER
# =============================================================================

# Detect which ranking stage actually exists
boot_csv_4 <- if (exists("DIR_STAGE_4")) {
  file.path(stage_dir(4), "tables", "bootstrap_rank_summary.csv")
} else {
  NA_character_
}

boot_csv_3 <- file.path(stage_dir(3), "tables", "bootstrap_rank_summary.csv")

if (!is.na(boot_csv_4) && file.exists(boot_csv_4)) {
  ranking_stage <- 4L
  boot_csv <- boot_csv_4
} else if (file.exists(boot_csv_3)) {
  ranking_stage <- 3L
  boot_csv <- boot_csv_3
} else {
  stop("No bootstrap_rank_summary.csv found in Stage 3 or Stage 4 output.")
}

boot_summary <- read.csv(boot_csv, stringsAsFactors = FALSE)

if (nrow(boot_summary) == 0 || is.na(boot_summary$top_candidate[1]))
  stop("bootstrap_rank_summary.csv contains no winner. ",
       "The ranking stage may have run with no surviving candidates. ",
       "Re-run the ranking stage and check for upstream errors.")

winner_id <- boot_summary$top_candidate[1]
message("Pipeline winner: ", winner_id)

# =============================================================================
# 2. LOAD SEASON LABELS AND METADATA
# =============================================================================

season_long   <- readRDS(file.path(stage_dir(1), "season_long.rds"))
threshold_tbl <- readRDS(file.path(stage_dir(1), "threshold_tbl.rds"))

winner_meta <- threshold_tbl %>%
  filter(candidate_id == winner_id) %>%
  slice(1)

if (nrow(winner_meta) == 0)
  stop("Winner '", winner_id, "' not found in threshold_tbl.rds from Stage 1. ",
       "Ensure Stage 1 was run with the same config and that its output files are intact.")

message("  Driver : ", winner_meta$driver,
        "\n  k      : ", winner_meta$k,
        "\n  Method : ", winner_meta$method)

# =============================================================================
# 3. BUILD FINAL ASSIGNMENT TABLE
# =============================================================================

season_final <- season_long %>%
  filter(candidate_id == winner_id, !is.na(season)) %>%
  transmute(
    Year         = year(DateMonth),
    Month        = month(DateMonth),
    season       = as.character(season),
    candidate_id = winner_id) %>%
  arrange(Year, Month)

if (nrow(season_final) == 0)
  stop("No labelled months found for winner '", winner_id, "' in season_long.rds. ",
       "The candidate may have been screened out or the season labels are all NA. ",
       "Re-run Stage 1 and check screened_out_candidates.csv.")

# =============================================================================
# 4. SUMMARY
# =============================================================================

season_summary <- season_final %>%
  count(season, name = "n_months") %>%
  mutate(pct = round(n_months / sum(n_months) * 100, 1))

message("\nSeason distribution:")
for (i in seq_len(nrow(season_summary)))
  message("  ", season_summary$season[i], ": ",
          season_summary$n_months[i], " months (",
          season_summary$pct[i], "%)")

message("\nCoverage: ",
        min(season_final$Year), "-", max(season_final$Year),
        " (", nrow(season_final), " months)")

# =============================================================================
# 5. SAVE
# =============================================================================

write.csv(season_final,
          file.path(PROJECT_DIR, "season_assignment_final.csv"),
          row.names = FALSE)

message("Saved: ", file.path(PROJECT_DIR, "season_assignment_final.csv"))
# =============================================================================
