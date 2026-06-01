# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`kinfitr` is an R package for kinetic modelling of PET (Positron Emission Tomography) Time Activity Curves (TACs). It implements reference tissue models, arterial input models, blood data processing, and BIDS-compliant data handling.

## Common Commands

```r
# Load package for interactive development
devtools::load_all()

# Regenerate documentation from roxygen2 comments
devtools::document()

# Run full R CMD check
devtools::check()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-refreg.R")

# Run tests matching a pattern
testthat::test_file("tests/testthat/test-refreg.R", filter = "srtm")
```

## Architecture

### Model Categories

**Reference Tissue Models** (no blood input required): SRTM, SRTM2, SRTMV, FRTM, MRTM1, MRTM2, refLogan, refPatlak, refmlLogan

**Arterial Input Models** (require blood input function): 1TCM, 2TCM variants, Logan, mlLogan, Patlak, MA1, MA2, lin2TCM

**Blood Data** (`kinfitr_blooddata.R`, `kinfitr_bloodfuncs.R`, `kinfitr_bloodmodels.R`): BIDS-compliant blood data structures, interpolation, dispersion correction, metabolite correction

**BIDS Support** (`kinfitr_bids.R`): PET BIDS data parsing and extraction

### Consistent Model API

All model functions follow a uniform pattern:

**Inputs:**
- `t_tac`: Time vector in minutes
- `tac`/`reftac`/`roitac`: Activity in kBq/ml
- `input`: Blood input function object (arterial models only)
- `weights`: Optional frame weights
- `frameStartEnd`: Frame index subsetting (or `timeStartEnd` for time-based)
- `Param.start`, `Param.lower`, `Param.upper`: Bounds for each fitted parameter
- `multstart_iter`: Number of random starts for optimization (default 1)

**Output** — all models return a named list with:
- `par`: Data frame of fitted parameters
- `par.se`: Standard errors as fractions (0–1 scale, where 1 = 100%)
- `fit`: Raw fit object (nls or similar)
- `weights`: Applied weights
- `tacs`: Data frame with columns `Times`, `Measured`, `Fitted`
- `model`: Model name string (used for S3 dispatch)

### Input Tidying

Two validation functions standardize inputs before fitting:
- `tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)` — reference tissue models
- `tidyinput_art(t_tac, tac, input, weights, frameStartEnd)` — arterial input models

Both auto-add a zero time point if missing and handle weight defaults.

### Plotting

S3 dispatch: `plot(fit_object)` → `plot.kinfit()` → `plot_kinfit()` → `plot_modelname()`. All plot functions return ggplot2 objects.

### Optimization

Uses `nls.multstart` (from the `nls.multstart` package) for robust nonlinear regression with optional multi-start. `fix_multstartpars()` handles multi-start parameter setup.

### Blood Data Structure

```r
blooddata$Data$Blood$Discrete$Values   # tibble: time, activity
blooddata$Data$Blood$Continuous$Values # tibble: time, activity
blooddata$Data$Plasma$Values           # tibble: time, activity
blooddata$Data$Metabolite$Values       # tibble: time, parentFraction
blooddata$TimeShift                    # numeric
```

Key functions: `create_blooddata_components()`, `bids_create_blooddata()`, `bd_getdata()`, `bd_extract()`, `bd_addfit()`.

### Weights

Nine weight calculation methods in `kinfitr_weights.R` via `weights_create()`. Some models (Logan, Patlak, refLogan, refPatlak) apply weight transformations to account for variable transformations on the dependent variable.

## Testing Conventions

Tests validate that fitted parameters fall within plausible ranges (e.g., SRTM BP between 1.5–2.5). Tests cover: basic fitting, with/without weights, `frameStartEnd` subsetting, multi-start, and plot output type.

## Key Dependencies

`nls.multstart`, `minpack.lm`, `ggplot2`, `dplyr`, `tidyr`, `purrr`, `tibble`, `broom`, `pracma`, `mgcv`, `jsonlite`
