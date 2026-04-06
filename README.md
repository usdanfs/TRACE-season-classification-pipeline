# Season classification pipeline

Reproducible pipeline for climate-based season classification in weakly seasonal ecosystems, with optional ecological breakpoint verification.

This repository contains:

- a **4-stage pipeline** for climate-based season classification with independent ecological verification
- a **3-stage climate-only pipeline** for study systems where no ecological response time series is available

## Software

The pipeline is written in **R**. Package requirements are loaded within the scripts.

## Overview

The pipeline selects a season definition from a predefined candidate set using climate structure, internal robustness, and, when available, independent ecological corroboration.

### Full 4-stage pipeline
1. **Stage 1**: Generate candidate season definitions from climate data and screen them for structural stability
2. **Stage 2**: Fit segmented regressions to an ecological response variable and evaluate breakpoint support
3. **Stage 3**: Stress-test surviving candidates within the ecological measurement window
4. **Stage 4**: Rank candidates using a tiered composite score, bootstrap rank stability, and weight sensitivity analysis

### Climate-only 3-stage pipeline
1. **Stage 1**: Generate candidate season definitions from climate data and screen them for structural stability
2. **Stage 2**: Stress-test surviving candidates within a validation window using structural criteria only
3. **Stage 3**: Rank candidates using climate structure, internal robustness, bootstrap rank stability, and weight sensitivity analysis

## Core repository files

### Full 4-stage pipeline (run in this order)
- `config.R`
- `STAGE_1_season_candidates.R`
- `STAGE_2_ecological_segmentation.R`
- `STAGE_3_season_validation.R`
- `STAGE_4_decision_ranking.R`
- `FINAL_season_assignment.R`

### Climate-only 3-stage pipeline (run in this order)
Scripts for this pipeline are in the `3STAGE/` subfolder.
- `3STAGE/config_climate_only.R`
- `3STAGE/STAGE_1_climate_only_candidates.R`
- `3STAGE/STAGE_2_climate_only_validation.R`
- `3STAGE/STAGE_3_climate_only_ranking.R`
- `FINAL_season_assignment.R`

### General documentation
- `SOP_Pipeline.docx` — step-by-step procedural guide: how to set up your data, fill in `config.R`, run the stages in order, and check outputs. Start here for a new study.
- `PIPELINE_OUTPUTS_GUIDE.md` — reference guide for interpreting every output column and CSV file (benchmarks, flag meanings, climate-only differences). Read this after your first run to understand what the numbers mean.

## Input requirements

### Full 4-stage pipeline
Two input tables are required:

1. **Climate data**
   - one row per `Year`–`Month`
   - monthly time step
   - required columns defined in `config.R`
   - must include the selected climate drivers or the variables needed to derive them

2. **Ecological response data**
   - one row per `Year`–`Month`
   - must contain `Year`, `Month`, and the response column specified in `RESPONSE_COL`
   - must already be aggregated to monthly scale before running the pipeline

### Climate-only 3-stage pipeline
One monthly climate table is required:
- one row per `Year`–`Month`
- required columns defined in `3STAGE/config_climate_only.R`

## General configuration

All site-specific settings are edited in the config file (not in the stage scripts). Open the config file for your pipeline version and review all parameters before the first run:
- `config.R` for the full 4-stage pipeline
- `3STAGE/config_climate_only.R` for the climate-only pipeline

At minimum, review and set:
- `CLIMATE_CSV` and `RESPONSE_CSV` / `RESPONSE_COL` (full pipeline only)
- `DRIVER_META` — driver names, polarity, and season labels
- `STD_THRESHOLDS` — literature-based classification thresholds
- `SEG_DRIVERS` and initial breakpoint guesses (full pipeline only)
- `BASELINE_START`, `BASELINE_END`
- Tier weights (`W_CLIMATE`, `W_ROBUST`, `W_VERIFY`)

Driver names must match exactly across:
- the climate CSV column headers
- `DRIVER_META`
- `STD_THRESHOLDS`
- `SEG_DRIVERS` (full pipeline only)

### Switching config files without editing scripts

Set the `SEASON_CONFIG` environment variable before sourcing any stage script to point to a different config file:

```r
Sys.setenv(SEASON_CONFIG = "TRAEC_data/config_TRACE.R")
source("STAGE_1_season_candidates.R")
```

If `SEASON_CONFIG` is not set, each script defaults to `config.R` (full pipeline) or `3STAGE/config_climate_only.R` (climate-only pipeline; path includes the subfolder prefix so scripts can be run from the project root).

## Pipeline outputs

The final deliverable is `season_assignment_final.csv` written to the project root directory by `FINAL_season_assignment.R`.

| Column | Description |
|--------|-------------|
| `Year` | Calendar year |
| `Month` | Calendar month (1–12) |
| `season` | Season label assigned to this month (e.g. "Wet", "Dry", "Transition") |
| `candidate_id` | Identifier of the winning season definition |

All intermediate diagnostic outputs are written to `output_STAGE_1/` through `output_STAGE_4/` (or `output_STAGE_1_climate_only_candidates/` etc. for the 3-stage pipeline). See `PIPELINE_OUTPUTS_GUIDE.md` for a full description of every output column and recommended interpretation benchmarks.

## Reproducibility notes

- All season definitions are derived from a predefined candidate set specified in the configuration file.
- The ecological tier is used only as an independent verification line and is weighted below climatological structure to avoid circularity in downstream analysis.
- The climate-only version is intended for cases where no ecological response time series is available for Stage 2 segmentation.
- The pipeline is applicable to any weakly seasonal ecosystem.

## Test data

The folder `test_data/` contains synthetic test data and a dedicated SOP describing the exact config adjustments required to run that dataset:

- `test_data/CLIM_test_alt_1990_2025.csv`
- `test_data/RESPONSE_test_alt_2019_2025.csv`
- `test_data/SOP_TestData_ConfigAdjustment.docx`
- `test_data/config_testdata.R`

## Code availability

This repository contains the code accompanying the manuscript:

**A Four-Stage Fraemwork for Methodical Season Classification in Humid Tropical Forests: Case study in the Luquillo Experimental Forest, Puerto Rico**
Daniel Minikaev, Debjani Sihi, Sasha C. Reed, Tana E. Wood
Submitted to *Agricultural and Forest Meteorology*

## Manuscript data

This repository contains the code used to classify season definitions at the TRACE site in the Luquillo Experimental Forest, Puerto Rico.

The folder `TRACE_data/` contains the TRACE climate and ecological data used in the manuscript, along with a dedicated SOP and config file:

- `TRACE_data/SOP_ManuscriptData_ConfigAdjustment.docx`
- `TRACE_data/CLIM_TRACE.csv`
- `TRACE_data/RESPONSE_TRACE.csv`
- `TRACE_data/config_TRACE.R`

## Citation

If you use this code, please cite this repository release.

Repository citation:
> Minikaev, D. 2026. Season Classification Pipeline (Version 1.0.0) [R 4.5.1]. GitHub.
> https://github.com/usdanfs/TRACE-season-classification-pipeline.git.

## License

This repository is released under the MIT License. See `LICENSE` for details.
