# Season classification pipeline

Reproducible pipeline for climate-based season classification in weakly seasonal ecosystems, with optional ecological breakpoint verification.

This repository contains:

- a **4-stage pipeline** for climate-based season classification with independent ecological verification
- a **3-stage climate-only pipeline** for study systems where no ecological response time series is available

## Software

The pipeline was developed in **R**. Package requirements are loaded within the scripts.

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

### Shared configuration
- `config.R`

### Full 4-stage pipeline (run in this order)
- `config.R`
- `STAGE_1_season_candidates.R`
- `STAGE_2_ecological_segmentation.R`
- `STAGE_3_season_validation.R`
- `STAGE_4_decision_ranking.R`
- `FINAL_season_assignment.R`

### Climate-only 3-stage pipeline (run in this order)
- `config_climate_only.R`
- `STAGE_1_climate_only_candidates.R`
- `STAGE_2_climate_only_validation.R`
- `STAGE_3_climate_only_ranking.R
- `FINAL_season_assignment.R`

### General documentation
- `SOP_Season_Pipeline.docx`

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
- required columns defined in `config_climate_only.R`

## General configuration

All site-specific settings are edited in config.
- project directory and file paths
- climate drivers and polarity
- standard thresholds
- baseline period
- ecological response column and related settings
- scoring weights and bootstrap settings

At minimum, users must review:
- `CLIMATE_CSV`
- `RESPONSE_CSV` (full pipeline only)
- `RESPONSE_COL` (full pipeline only)
- `DRIVER_META`
- `STD_THRESHOLDS`
- `SEG_DRIVERS` (full pipeline only)
- `BASELINE_START`, `BASELINE_END`
- stage thresholds and tier weights

Driver names must match exactly across:
- the climate CSV
- `DRIVER_META`
- `STD_THRESHOLDS`
- `SEG_DRIVERS` (full pipeline only)

## Reproducibility notes

- All season definitions are derived from a predefined candidate set specified in the configuration file.
- The ecological tier is used only as an independent verification line and is weighted below climatological structure to avoid circularity in downstream analysis.
- The climate-only version is intended for cases where no ecological response time series is available for Stage 2 segmentation.
- The pipeline is applicable to any weakly seasonal ecosystem.

## Test data

The folder `test_data/` contains synthetic test data and one separate SOP with the exact config adjustments required to run that dataset successfully:

- `test_data/CLIM_test_alt_1990_2025.csv`
- `test_data/RESPONSE_test_alt_2019_2025.csv`
- `test_data/SOP_TestData_ConfigAdjustment.docx`
- `config_testdata.R`

## Code availability

This repository contains the code accompanying the manuscript:

**A Four-Stage Framework for Methodical Season Classification in Humid Tropical Forests: Case study in the Luquillo Experimental Forest, Puerto Rico**  
Daniel Minikaev, Debjani Sihi, Sasha C. Reed, Tana E. Wood  
Submitted to *Agricultural and Forest Meteorology*

## Manuscript data

This repository contains the code used to classify season definitions at the TRACE site in the Luquillo Experimental Forest, Puerto Rico.

The folder `TRACE_data/` contains the TRACE climate and ecological data used to classify season definitions at the site in the Luquillo Experimental Forest, PR.
An SOP with the exact config adjustments required to run that dataset successfully, and the dedicated config are provided.

- `SOP_ManuscriptData_ConfigAdjustment.docx`
- `test_data/CLIM_TRACE.csv`
- `test_data/RESPONSE_TRACE.csv`
- `config_TRACE.R`

## Citation

If you use this code, please cite this repository release.

Repository citation:
> Minikaev, D. 2026. Season Classification Pipeline (Version 1.0.0) [R 4.5.1]. GitHub.
> https://github.com/usdanfs/TRACE-season-classification-pipeline.git.

## License

This repository is released under the MIT License. See `LICENSE` for details.
