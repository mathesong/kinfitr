# Implementation Plan: Nested Multi-Region Kinetic Models

## Overview

Implement nested multi-region kinetic model fitting where an outer `optim()` loop
estimates shared parameters across regions while an inner `nlsLM`/`nls_multstart`
loop fits per-region kinetic parameters. Input data is in long format (matching the
`longdata_multiTAC_models` branch conventions).

## Models to Implement

| # | Function | File | Outer param(s) | Inner model function | Inner params (per region) |
|---|---|---|---|---|---|
| 1 | `nested_1tcm_delay` | `R/kinfitr_nested_1tcm.R` | `inpshift` | `onetcm_model` | K1, k2 |
| 2 | `nested_2tcm_delay` | `R/kinfitr_nested_2tcm.R` | `inpshift` | `twotcm_model` | K1, k2, k3, k4 |
| 3 | `nested_2tcm` | `R/kinfitr_nested_2tcm.R` | Vnd and/or k4 | `twotcm_macro_model` | K1, BPp (+ whichever of Vnd/k4 not outer) |
| 4 | `nested_srtm` | `R/kinfitr_nested_srtm.R` | `k2prime` | `srtm2_model` (k2prime fixed) | R1, bp |

Tests go in `tests/testthat/test-nested.R`.

---

## Shared Architecture

All nested functions follow the same algorithmic pattern. Each function is
self-contained (no generic dispatcher), but they share an internal helper for
per-region fitting.

### Common Algorithm

```
1. Validate & prepare input data
2. Build outer objective function for optim()
3. Run optim() to find optimal shared parameters
4. Final pass: re-fit all regions at optimal outer params
5. Compute derived parameters
6. Assemble and return output
```

### Internal Helper: `.fit_region()`

Defined once (in `kinfitr_nested_1tcm.R`, the first file created), used by all
nested functions. Not exported.

```r
#' Fit a single region using nlsLM or nls_multstart
#'
#' @param formula_str Character string of the nlsLM formula
#' @param modeldata Named list with t_tac, tac, weights, input (or reftac)
#' @param start Named numeric vector of starting values
#' @param lower Named numeric vector of lower bounds
#' @param upper Named numeric vector of upper bounds
#' @param multstart_iter Numeric. If > 1, use nls_multstart
#' @param multstart_lower Named numeric vector of multstart lower bounds
#' @param multstart_upper Named numeric vector of multstart upper bounds
#'
#' @return The nls fit object, or NULL if fitting fails
.fit_region <- function(formula_str, modeldata, start, lower, upper,
                        multstart_iter = 1,
                        multstart_lower = NULL, multstart_upper = NULL) {
  tryCatch({
    if (prod(multstart_iter) == 1) {
      minpack.lm::nlsLM(
        as.formula(formula_str),
        data = modeldata,
        start = start, lower = lower, upper = upper,
        weights = weights,
        control = minpack.lm::nls.lm.control(maxiter = 200)
      )
    } else {
      nls.multstart::nls_multstart(
        as.formula(formula_str),
        data = modeldata,
        supp_errors = "Y",
        start_lower = multstart_lower,
        start_upper = multstart_upper,
        iter = multstart_iter, convergence_count = FALSE,
        lower = lower, upper = upper,
        modelweights = weights
      )
    }
  }, error = function(e) NULL)
}
```

### Common Output Structure

All nested functions return a list with class `c("<model_name>", "kinfit")`:

```r
list(
  par        = data.frame(),   # One row per region. Columns: Region, fitted params,
                               #   outer params (repeated), fixed params, derived params
  par.se     = data.frame(),   # One row per region. Columns: Region, <param>.se for
                               #   inner params only (outer param SEs not available from
                               #   inner nlsLM fits)
  outer_par  = data.frame(),   # Single row with optimal outer parameter value(s)
  optim      = optim_result,   # Raw optim() output for diagnostics (includes
                               #   convergence, counts, value)
  tacs       = data.frame(),   # Long-format: Time, Region, Radioactivity, Fitted
  input      = input,          # Blood input after applying optimal shift (or original
                               #   if no shift). For SRTM: not included
  weights    = weights_per_region, # Per-region weight vector
  roiweights = roiweights,     # Named vector of ROI weights
  model      = "<model_name>"  # String for plot_kinfit() dispatch
)
```

### Common Input Interface

Following the long-format convention from `SIME()` on this branch:

```r
# For arterial input models (1tcm, 2tcm):
nested_Xtcm_delay(
  t_tac,             # Numeric vector of times (minutes), repeated per region
  tac,               # Numeric vector of radioactivity, stacked per region
  region,            # Character vector identifying regions
  input,             # interpblood object from blood_interp()
  ...
)

# For reference tissue models (srtm):
nested_srtm(
  t_tac,             # Numeric vector of times (minutes), repeated per region
  roitac,            # Numeric vector of target ROI radioactivity, stacked per region
  reftac,            # Numeric vector of reference tissue radioactivity (same for all
                     #   regions, but repeated per region in long format)
  region,            # Character vector identifying regions
  ...
)
```

### Common Data Preparation Steps

All functions perform these steps at the start:

```r
# 1. Convert timeStartEnd to frameStartEnd
if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
  t_first <- t_tac[region == regions[1]]
  frameStartEnd <- c(which(t_first >= timeStartEnd[1])[1],
                     tail(which(t_first <= timeStartEnd[2]), 1))
}

# 2. Tidy input (adds zero frame, handles weights, validates lengths)
tidied <- tidyinput_long(t_tac, tac, region, weights, frameStartEnd)
t_tac   <- tidied$t_tac
tac     <- tidied$tac
region  <- tidied$region
weights <- tidied$weights
regions <- unique(region)

# 3. Handle ROI weights
if (is.null(roiweights)) {
  roiweights <- rep(1, length(regions))
  names(roiweights) <- regions
} else if (!is.null(names(roiweights))) {
  roiweights <- roiweights[regions]
} else {
  names(roiweights) <- regions
}
roiweights <- roiweights / max(roiweights)

# 4. Extract per-region data
weights_per_region <- weights[region == regions[1]]
n_frames <- sum(region == regions[1])
```

### Common Outer Objective Pattern

```r
outer_objective <- function(outer_vals) {
  # Apply outer parameters (model-specific: shift input, fix k2prime, etc.)
  # ...

  total_rss <- 0
  for (r in regions) {
    region_tac <- tac_shifted[region_shifted == r]
    region_t   <- t_shifted[region_shifted == r]

    # Build modeldata list for nlsLM
    modeldata <- list(
      tac = region_tac,
      t_tac = region_t,
      weights = weights_per_region,
      input = shifted_input       # or reftac for reference models
    )

    fit <- .fit_region(formula_str, modeldata, start, lower, upper,
                       multstart_iter, multstart_lower, multstart_upper)

    if (is.null(fit)) {
      total_rss <- total_rss + 1e10 * roiweights[r]
    } else {
      region_rss <- sum(weights(fit) * residuals(fit)^2)
      total_rss <- total_rss + region_rss * roiweights[r]
    }
  }
  return(total_rss)
}
```

### Common Final Pass Pattern

After `optim()` converges, re-run inner fits at optimal outer params to extract
full results:

```r
# Apply optimal outer params
# ...

par_list <- list()
par_se_list <- list()
fitted_tacs <- list()

for (r in regions) {
  # Build modeldata, fit with .fit_region()
  fit <- .fit_region(...)

  # Extract parameters
  region_par <- as.data.frame(as.list(coef(fit)))
  region_par$Region <- r

  # Add outer params
  region_par$<outer_param> <- optimal_value

  # Add fixed params
  region_par$vB <- vB  # (or whichever is fixed)

  # Compute derived params (model-specific)
  # ...

  # Extract SEs
  region_se <- region_par
  region_se[1, ] <- purrr::map_dbl(names(coef(fit)), ~ get_se(fit, .x))
  names(region_se) <- paste0(names(region_se), ".se")
  region_se$Region <- r

  # Store fitted values
  fitted_tacs[[r]] <- data.frame(
    Time = region_t,
    Region = r,
    Radioactivity = region_tac,
    Fitted = as.numeric(fitted(fit)),
    stringsAsFactors = FALSE
  )

  par_list[[r]] <- region_par
  par_se_list[[r]] <- region_se
}

par <- do.call(rbind, par_list)
par.se <- do.call(rbind, par_se_list)
tacs <- do.call(rbind, fitted_tacs)
```

---

## Model 1: `nested_1tcm_delay`

**File:** `R/kinfitr_nested_1tcm.R`

### Function Signature

```r
#' Nested 1TCM with Delay Estimation
#'
#' Fits the One Tissue Compartment Model to multiple regions simultaneously,
#' estimating a shared input delay (inpshift) across all regions using optim()
#' in the outer loop, while fitting K1 and k2 per region using nlsLM in the
#' inner loop.
#'
#' @param t_tac Numeric vector of times for each frame in minutes (repeated for
#'   each region). If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations, stacked across
#'   regions.
#' @param region Character vector identifying the region for each observation.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time. Generated using \code{blood_interp}.
#' @param vB The blood volume fraction. Default is 0.05.
#' @param inpshift.start Starting value for the delay parameter. Default is 0.
#' @param inpshift.lower Lower bound for the delay parameter. Default is -0.5.
#' @param inpshift.upper Upper bound for the delay parameter. Default is 0.5.
#' @param K1.start Starting parameter for K1. Default is 0.1.
#' @param K1.lower Lower bound for K1. Default is 0.0001.
#' @param K1.upper Upper bound for K1. Default is 1.
#' @param k2.start Starting parameter for k2. Default is 0.1.
#' @param k2.lower Lower bound for k2. Default is 0.0001.
#' @param k2.upper Upper bound for k2. Default is 0.5.
#' @param weights Optional. Numeric vector of weights. Same conventions as SIME.
#' @param roiweights Optional. Named numeric vector of ROI weights.
#' @param optim_method Optimization method for optim(). Default is "L-BFGS-B".
#' @param optim_control List of control parameters for optim().
#' @param multstart_iter Number of multistart iterations for inner fits.
#'   Default is 1.
#' @param multstart_lower Optional. Lower bounds for multistart starting params.
#' @param multstart_upper Optional. Upper bounds for multistart starting params.
#' @param frameStartEnd Optional. Frame range c(start, end). Applied per region.
#' @param timeStartEnd Optional. Time range for frame selection.
#'
#' @return A list with class c("nested_1tcm_delay", "kinfit") containing:
#'   \code{$par} (per-region parameters), \code{$par.se} (standard errors),
#'   \code{$outer_par} (optimal inpshift), \code{$optim} (raw optim result),
#'   \code{$tacs} (long-format TACs with fitted values), \code{$input},
#'   \code{$weights}, \code{$roiweights}, \code{$model}.
#'
#' @export
nested_1tcm_delay <- function(
    t_tac, tac, region, input,
    vB = 0.05,
    inpshift.start = 0, inpshift.lower = -0.5, inpshift.upper = 0.5,
    K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
    k2.start = 0.1, k2.lower = 0.0001, k2.upper = 0.5,
    weights = NULL, roiweights = NULL,
    optim_method = "L-BFGS-B",
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL)
```

### Outer Objective

```r
outer_objective <- function(outer_vals) {
  inpshift_val <- outer_vals[["inpshift"]]

  # Shift timings for all regions
  shifted <- shift_timings_long(t_tac, tac, region, input, inpshift_val)

  total_rss <- 0
  for (r in regions) {
    region_mask <- shifted$region == r
    modeldata <- list(
      tac = shifted$tac[region_mask],
      t_tac = shifted$t_tac[region_mask],
      weights = weights_per_region,
      input = shifted$input
    )

    formula_str <- paste0("tac ~ onetcm_model(t_tac, input, K1, k2, vB=", vB, ")")

    fit <- .fit_region(formula_str, modeldata, start, lower, upper,
                       multstart_iter, multstart_lower, multstart_upper)

    if (is.null(fit)) {
      total_rss <- total_rss + 1e10 * roiweights[r]
    } else {
      total_rss <- total_rss + sum(weights(fit) * residuals(fit)^2) * roiweights[r]
    }
  }
  return(total_rss)
}
```

Where:
- `start <- c(K1 = K1.start, k2 = k2.start)`
- `lower <- c(K1 = K1.lower, k2 = k2.lower)`
- `upper <- c(K1 = K1.upper, k2 = k2.upper)`

### optim() Call

```r
optim_result <- optim(
  par = c(inpshift = inpshift.start),
  fn = outer_objective,
  method = optim_method,
  lower = c(inpshift = inpshift.lower),
  upper = c(inpshift = inpshift.upper),
  control = optim_control
)
```

### Derived Parameters

Per region:
- `Vt = K1 / k2`

SE for Vt: `get_se(fit, "K1/k2")`

### Plotting: `plot_nested_1tcm_delayfit()`

**Single combined plot** (NOT faceted) with all regions overlaid in different
colours. Points are measured TAC values, lines are fitted values, coloured by
region. The input function (AIF) is shown as a shared line. This works because
the delay is shared — all regions are on the same time axis after shifting.

```r
#' @export
plot_nested_1tcm_delayfit <- function(nested_out, roiname = NULL) {
  tacs <- nested_out$tacs

  # Measured points: colour by Region
  # Fitted lines: colour by Region (matching)
  # AIF: separate line

  measured <- data.frame(
    Time = tacs$Time,
    Radioactivity = tacs$Radioactivity,
    Region = tacs$Region,
    stringsAsFactors = FALSE
  )

  fitted_df <- data.frame(
    Time = tacs$Time,
    Radioactivity = tacs$Fitted,
    Region = tacs$Region,
    stringsAsFactors = FALSE
  )

  aif_df <- data.frame(
    Time = nested_out$input$Time,
    Radioactivity = nested_out$input$AIF,
    Region = "AIF"
  )

  outplot <- ggplot(measured, aes(x = Time, y = Radioactivity, colour = Region)) +
    geom_point() +
    geom_line(data = fitted_df) +
    geom_line(data = aif_df, colour = "black", linetype = "dashed") +
    labs(title = paste0("Nested 1TCM Delay (inpshift = ",
                        round(nested_out$outer_par$inpshift, 4), ")"))

  return(outplot)
}
```

This dispatches via `plot_kinfit()` because `model = "nested_1tcm_delay"` →
`plot_nested_1tcm_delayfit()`.

---

## Model 2: `nested_2tcm_delay`

**File:** `R/kinfitr_nested_2tcm.R`

### Function Signature

```r
#' Nested 2TCM with Delay Estimation
#'
#' Fits the Two Tissue Compartment Model to multiple regions simultaneously,
#' estimating a shared input delay (inpshift) across all regions using optim()
#' in the outer loop, while fitting K1, k2, k3, k4 per region using nlsLM in
#' the inner loop.
#'
#' @inheritParams nested_1tcm_delay
#' @param k3.start Starting parameter for k3. Default is 0.1.
#' @param k3.lower Lower bound for k3. Default is 0.0001.
#' @param k3.upper Upper bound for k3. Default is 0.5.
#' @param k4.start Starting parameter for k4. Default is 0.1.
#' @param k4.lower Lower bound for k4. Default is 0.0001.
#' @param k4.upper Upper bound for k4. Default is 0.5.
#'
#' @return A list with class c("nested_2tcm_delay", "kinfit").
#'
#' @export
nested_2tcm_delay <- function(
    t_tac, tac, region, input,
    vB = 0.05,
    inpshift.start = 0, inpshift.lower = -0.5, inpshift.upper = 0.5,
    K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
    k2.start = 0.1, k2.lower = 0.0001, k2.upper = 0.5,
    k3.start = 0.1, k3.lower = 0.0001, k3.upper = 0.5,
    k4.start = 0.1, k4.lower = 0.0001, k4.upper = 0.5,
    weights = NULL, roiweights = NULL,
    optim_method = "L-BFGS-B",
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL)
```

### Differences from nested_1tcm_delay

- Inner formula: `tac ~ twotcm_model(t_tac, input, K1, k2, k3, k4, vB=<value>)`
- Inner params: `start <- c(K1 = K1.start, k2 = k2.start, k3 = k3.start, k4 = k4.start)`
- Derived parameters per region:
  - `Vt = (K1/k2) * (1 + k3/k4)` — SE: `get_se(fit, "(K1/k2) * (1 + k3/k4)")`
  - `Vnd = K1/k2` — SE: `get_se(fit, "K1/k2")`
  - `BPnd = k3/k4` — SE: `get_se(fit, "k3/k4")`
  - `BPp = (K1/k2) * (k3/k4)` — SE: `get_se(fit, "(K1/k2) * (k3/k4)")`

### Plotting: `plot_nested_2tcm_delayfit()`

Same single combined plot as `plot_nested_1tcm_delayfit()` (NOT faceted) — all
regions overlaid with coloured points and lines, plus shared AIF. Title reflects
2TCM.

---

## Model 3: `nested_2tcm`

**File:** `R/kinfitr_nested_2tcm.R` (same file as nested_2tcm_delay)

This model uses the **macro parameterisation** (`twotcm_macro_model`) and
optimizes Vnd and/or k4 in the outer loop. The inner loop fits K1 and BPp
(and whichever of Vnd/k4 is not in the outer loop) per region.

### Function Signature

```r
#' Nested 2TCM with Shared Macroparameters
#'
#' Fits the Two Tissue Compartment Model (macro parameterisation) to multiple
#' regions simultaneously. Optimizes shared parameters (Vnd and/or k4) across
#' all regions using optim() in the outer loop, while fitting per-region
#' parameters (K1, BPp, and any non-shared parameters) using nlsLM in the
#' inner loop.
#'
#' @param t_tac Numeric vector of times (minutes), repeated per region.
#' @param tac Numeric vector of radioactivity, stacked per region.
#' @param region Character vector identifying regions.
#' @param input Data frame from blood_interp().
#' @param inpshift The delay in minutes. Must be pre-fitted (not estimated here).
#'   Default is 0.
#' @param vB Blood volume fraction. Default is 0.05.
#' @param Vnd Optional. If NULL, Vnd is estimated in the outer loop (shared
#'   across regions). If a numeric value, Vnd is fixed.
#' @param k4 Optional. If NULL, k4 is estimated in the outer loop (shared
#'   across regions). If a numeric value, k4 is fixed.
#' @param Vnd.start Starting value for Vnd (if fitted). Default is 1.
#' @param Vnd.lower Lower bound for Vnd. Default is 0.0001.
#' @param Vnd.upper Upper bound for Vnd. Default is 10.
#' @param k4.start Starting value for k4 (if fitted). Default is 0.1.
#' @param k4.lower Lower bound for k4. Default is 0.0001.
#' @param k4.upper Upper bound for k4. Default is 0.5.
#' @param K1.start Starting parameter for K1. Default is 0.1.
#' @param K1.lower Lower bound for K1. Default is 0.0001.
#' @param K1.upper Upper bound for K1. Default is 1.
#' @param BPp.start Starting parameter for BPp. Default is 1.
#' @param BPp.lower Lower bound for BPp. Default is 0.0001.
#' @param BPp.upper Upper bound for BPp. Default is 50.
#' @param weights Optional. Numeric vector of weights.
#' @param roiweights Optional. Named numeric vector of ROI weights.
#' @param optim_method Optimization method. Default is "L-BFGS-B".
#' @param optim_control Control parameters for optim().
#' @param multstart_iter Multistart iterations for inner fits. Default is 1.
#' @param multstart_lower Optional. Lower bounds for multistart.
#' @param multstart_upper Optional. Upper bounds for multistart.
#' @param frameStartEnd Optional. Frame range c(start, end).
#' @param timeStartEnd Optional. Time range for frame selection.
#'
#' @return A list with class c("nested_2tcm", "kinfit").
#'
#' @export
nested_2tcm <- function(
    t_tac, tac, region, input,
    inpshift = 0, vB = 0.05,
    Vnd = NULL, k4 = NULL,
    Vnd.start = 1, Vnd.lower = 0.0001, Vnd.upper = 10,
    k4.start = 0.1, k4.lower = 0.0001, k4.upper = 0.5,
    K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
    BPp.start = 1, BPp.lower = 0.0001, BPp.upper = 50,
    weights = NULL, roiweights = NULL,
    optim_method = "L-BFGS-B",
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL)
```

### Key Design Details

**Determining what goes in the outer loop:**

At least one of `Vnd` or `k4` must be NULL (fitted in outer loop). If both are
NULL, both are optimized jointly. If both are provided, there's nothing to nest
— error with a helpful message.

```r
# Determine outer parameters
outer_start <- c()
outer_lower <- c()
outer_upper <- c()

if (is.null(Vnd)) {
  outer_start <- c(outer_start, Vnd = Vnd.start)
  outer_lower <- c(outer_lower, Vnd = Vnd.lower)
  outer_upper <- c(outer_upper, Vnd = Vnd.upper)
}
if (is.null(k4)) {
  outer_start <- c(outer_start, k4 = k4.start)
  outer_lower <- c(outer_lower, k4 = k4.lower)
  outer_upper <- c(outer_upper, k4 = k4.upper)
}

if (length(outer_start) == 0) {
  stop("At least one of Vnd or k4 must be NULL to be estimated in the outer loop.")
}
```

**Inner formula construction:**

The formula for `twotcm_macro_model(t_tac, input, K1, Vnd, BPp, k4, vB)` needs
to fix whichever of Vnd/k4 comes from the outer loop:

```r
# Inside outer_objective:
outer_objective <- function(outer_vals) {
  Vnd_val <- if (is.null(Vnd)) outer_vals[["Vnd"]] else Vnd
  k4_val  <- if (is.null(k4))  outer_vals[["k4"]]  else k4

  # Shift timings (inpshift is pre-fitted for this model)
  shifted <- shift_timings_long(t_tac, tac, region, input, inpshift)

  total_rss <- 0
  for (r in regions) {
    region_mask <- shifted$region == r
    modeldata <- list(
      tac = shifted$tac[region_mask],
      t_tac = shifted$t_tac[region_mask],
      weights = weights_per_region,
      input = shifted$input
    )

    # Build formula with fixed outer params
    formula_str <- paste0(
      "tac ~ twotcm_macro_model(t_tac, input, K1",
      ", Vnd=", Vnd_val,
      ", BPp",
      ", k4=", k4_val,
      ", vB=", vB, ")"
    )

    # Inner params are K1 and BPp only (Vnd and k4 fixed from outer)
    fit <- .fit_region(formula_str, modeldata,
                       start = c(K1 = K1.start, BPp = BPp.start),
                       lower = c(K1 = K1.lower, BPp = BPp.lower),
                       upper = c(K1 = K1.upper, BPp = BPp.upper),
                       multstart_iter, multstart_lower, multstart_upper)

    # ... aggregate RSS as usual
  }
  return(total_rss)
}
```

**Note:** If only one of Vnd/k4 is in the outer loop, the other goes to the
inner loop. The formula and inner start/lower/upper vectors must be adjusted
accordingly. Specifically:

- If only Vnd is outer: inner params are K1, BPp, k4
- If only k4 is outer: inner params are K1, Vnd, BPp
- If both Vnd and k4 are outer: inner params are K1, BPp

### Derived Parameters

Per region (using the optimal outer values + per-region inner values):
- `Vt = Vnd + BPp` — or computed from whichever params are available
- `BPnd = BPp / Vnd`
- `k2 = K1 / Vnd`
- `k3 = (BPp / Vnd) * k4`

### Plotting: `plot_nested_2tcmfit()`

**Faceted plot** with one panel per region (unlike the delay models which use a
single combined plot). Each panel shows measured points and fitted line. Faceting
is appropriate here because different regions may have different parameter values
for Vnd/k4 from the outer loop vs inner loop, and the shared parameter story is
less visual than delay.

---

## Model 4: `nested_srtm`

**File:** `R/kinfitr_nested_srtm.R`

### Function Signature

```r
#' Nested SRTM with Shared k2prime Estimation
#'
#' Fits the Simplified Reference Tissue Model 2 (SRTM2) to multiple regions
#' simultaneously, estimating a shared k2prime across all regions using optim()
#' in the outer loop, while fitting R1 and bp per region using nlsLM in the
#' inner loop. This is the proper way to estimate k2prime for SRTM2 — by
#' leveraging information across multiple regions.
#'
#' @param t_tac Numeric vector of times (minutes), repeated per region.
#' @param roitac Numeric vector of target ROI radioactivity, stacked per region.
#' @param reftac Numeric vector of reference tissue radioactivity. Should be the
#'   same reference region values repeated for each target region (i.e., same
#'   length as roitac). Alternatively, can be a single vector of per-frame
#'   values that will be repeated for each region.
#' @param region Character vector identifying the region for each observation.
#' @param k2prime.start Starting value for k2prime. Default is 0.1.
#' @param k2prime.lower Lower bound for k2prime. Default is 0.001.
#' @param k2prime.upper Upper bound for k2prime. Default is 1.
#' @param R1.start Starting parameter for R1. Default is 1.
#' @param R1.lower Lower bound for R1. Default is 0.
#' @param R1.upper Upper bound for R1. Default is 10.
#' @param bp.start Starting parameter for bp. Default is 1.5.
#' @param bp.lower Lower bound for bp. Default is 0.
#' @param bp.upper Upper bound for bp. Default is 15.
#' @param weights Optional. Numeric vector of weights.
#' @param roiweights Optional. Named numeric vector of ROI weights.
#' @param optim_method Optimization method. Default is "L-BFGS-B".
#' @param optim_control Control parameters for optim().
#' @param multstart_iter Multistart iterations for inner fits. Default is 1.
#' @param multstart_lower Optional. Lower bounds for multistart.
#' @param multstart_upper Optional. Upper bounds for multistart.
#' @param frameStartEnd Optional. Frame range c(start, end).
#' @param timeStartEnd Optional. Time range for frame selection.
#'
#' @return A list with class c("nested_srtm", "kinfit").
#'
#' @export
nested_srtm <- function(
    t_tac, roitac, reftac, region,
    k2prime.start = 0.1, k2prime.lower = 0.001, k2prime.upper = 1,
    R1.start = 1, R1.lower = 0, R1.upper = 10,
    bp.start = 1.5, bp.lower = 0, bp.upper = 15,
    weights = NULL, roiweights = NULL,
    optim_method = "L-BFGS-B",
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL)
```

### Key Design Details

**Reference tissue handling:**

The reference TAC is the same for all target regions. The `reftac` vector should
be provided in long format (repeated for each region). If the user provides a
single vector (length = frames per region), it should be repeated:

```r
n_per_region <- length(roitac) / length(unique(region))
if (length(reftac) == n_per_region) {
  reftac <- rep(reftac, times = length(unique(region)))
}
```

**Data tidying:**

Uses `tidyinput_long()` on the `roitac` data. The `reftac` is handled separately
(tidied to match the same frames). We could also use `tidyinput_ref()` per
region in the inner loop.

**Outer objective:**

```r
outer_objective <- function(outer_vals) {
  k2prime_val <- outer_vals[["k2prime"]]

  total_rss <- 0
  for (r in regions) {
    region_mask <- region == r
    region_roitac <- roitac[region_mask]
    region_reftac <- reftac[region_mask]
    region_t_tac  <- t_tac[region_mask]

    modeldata <- list(
      roitac = region_roitac,
      reftac = region_reftac,
      t_tac = region_t_tac,
      weights = weights_per_region
    )

    formula_str <- paste0("roitac ~ srtm2_model(t_tac, reftac, R1, k2prime=",
                          k2prime_val, ", bp)")

    fit <- .fit_region(formula_str, modeldata,
                       start = c(R1 = R1.start, bp = bp.start),
                       lower = c(R1 = R1.lower, bp = bp.lower),
                       upper = c(R1 = R1.upper, bp = bp.upper),
                       multstart_iter, multstart_lower, multstart_upper)

    if (is.null(fit)) {
      total_rss <- total_rss + 1e10 * roiweights[r]
    } else {
      total_rss <- total_rss + sum(weights(fit) * residuals(fit)^2) * roiweights[r]
    }
  }
  return(total_rss)
}
```

**Note on srtm2_model:** The existing `srtm2_model(t_tac, reftac, R1, k2prime, bp)`
accepts k2prime as a regular argument. When we fix it in the formula string as
`k2prime=<value>`, it becomes a fixed value in the nlsLM call and R1, bp are
estimated.

### Derived Parameters

Per region:
- `k2a = (R1 * k2prime) / (bp + 1)` — SE: `get_se(fit, paste0("(R1 * ", k2prime_val, ") / (bp + 1)"))`

### Output Structure

Same as other nested models but:
- No `$input` field (reference tissue model)
- Has `$reftac` field instead (the reference tissue TAC)
- `$tacs` includes columns: Time, Region, Radioactivity (roitac), Reference
  (reftac), Fitted

### Plotting: `plot_nested_srtmfit()`

**Faceted plot** with one panel per region. Each panel shows measured target ROI
points, fitted line, and the shared reference TAC line.

---

## Test File: `tests/testthat/test-nested.R`

```r
context("test-nested.R")

data("pbr28")

meas <- 2

# Prepare long-format data from pbr28 (selecting a few regions)
tac_wide <- pbr28$tacs[[meas]]
t_tac_single <- tac_wide$Times / 60
weights_single <- tac_wide$Weights
input <- pbr28$input[[meas]]

# Select 3-4 regions for testing
selected_regions <- c("FC", "TC", "TH")

long_data <- do.call(rbind, lapply(selected_regions, function(r) {
  data.frame(
    t_tac = t_tac_single,
    tac = tac_wide[[r]],
    region = r,
    weights = weights_single,
    stringsAsFactors = FALSE
  )
}))

set.seed(42)

# --- nested_1tcm_delay tests ---

test_that("nested_1tcm_delay works", {
  out <- nested_1tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    weights = long_data$weights,
    vB = 0.05,
    K1.start = 0.1, K1.lower = 0.05, K1.upper = 0.5,
    k2.start = 0.1, k2.lower = 0.01, k2.upper = 0.3
  )

  # Check structure
  expect_true("nested_1tcm_delay" %in% class(out))
  expect_true("kinfit" %in% class(out))
  expect_true(all(c("par", "par.se", "outer_par", "optim", "tacs") %in% names(out)))

  # Check inpshift is within bounds
  expect_gte(out$outer_par$inpshift, -0.5)
  expect_lte(out$outer_par$inpshift, 0.5)

  # Check per-region results
  expect_equal(nrow(out$par), length(selected_regions))
  expect_true(all(out$par$Vt > 0))
  expect_true(all(out$par$Vt < 20))

  # Check plotting
  expect_true(any(class(plot_kinfit(out)) == "ggplot"))
})

# --- nested_2tcm_delay tests ---

test_that("nested_2tcm_delay works", {
  out <- nested_2tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    weights = long_data$weights,
    vB = 0.05,
    K1.start = 0.1, K1.lower = 0.05, K1.upper = 0.5,
    k2.start = 0.1, k2.lower = 0.01, k2.upper = 0.3,
    k3.start = 0.05, k3.lower = 0.001, k3.upper = 0.3,
    k4.start = 0.05, k4.lower = 0.001, k4.upper = 0.3
  )

  expect_true("nested_2tcm_delay" %in% class(out))
  expect_equal(nrow(out$par), length(selected_regions))
  expect_true(all(out$par$Vt > 0))
  expect_true(any(class(plot_kinfit(out)) == "ggplot"))
})

# --- nested_2tcm tests ---

test_that("nested_2tcm with Vnd in outer loop works", {
  out <- nested_2tcm(
    long_data$t_tac, long_data$tac, long_data$region, input,
    weights = long_data$weights,
    inpshift = 0.14, vB = 0.05,
    Vnd = NULL, k4 = 0.05  # Vnd in outer, k4 fixed
  )

  expect_true("nested_2tcm" %in% class(out))
  expect_equal(nrow(out$par), length(selected_regions))
  expect_true(!is.null(out$outer_par$Vnd))
})

test_that("nested_2tcm with both Vnd and k4 in outer loop works", {
  out <- nested_2tcm(
    long_data$t_tac, long_data$tac, long_data$region, input,
    weights = long_data$weights,
    inpshift = 0.14, vB = 0.05,
    Vnd = NULL, k4 = NULL  # Both in outer
  )

  expect_true(!is.null(out$outer_par$Vnd))
  expect_true(!is.null(out$outer_par$k4))
})

# --- nested_srtm tests ---
# (Uses simref data instead of pbr28)

test_that("nested_srtm works", {
  data("simref")

  # Build long-format data from simref
  sim_t <- simref$tacs[[2]]$Times
  sim_ref <- simref$tacs[[2]]$Reference
  sim_weights <- simref$tacs[[2]]$Weights

  sim_regions <- c("ROI1", "ROI2", "ROI3")
  sim_long <- do.call(rbind, lapply(sim_regions, function(r) {
    data.frame(
      t_tac = sim_t,
      roitac = simref$tacs[[2]][[r]],
      reftac = sim_ref,
      region = r,
      weights = sim_weights,
      stringsAsFactors = FALSE
    )
  }))

  out <- nested_srtm(
    sim_long$t_tac, sim_long$roitac, sim_long$reftac, sim_long$region,
    weights = sim_long$weights
  )

  expect_true("nested_srtm" %in% class(out))
  expect_equal(nrow(out$par), length(sim_regions))
  expect_true(out$outer_par$k2prime > 0)
  expect_true(out$outer_par$k2prime < 1)
  expect_true(any(class(plot_kinfit(out)) == "ggplot"))
})
```

---

## Implementation Order

1. **`R/kinfitr_nested_1tcm.R`** — Contains `nested_1tcm_delay()`,
   `plot_nested_1tcm_delayfit()`, and the shared `.fit_region()` helper
2. **`R/kinfitr_nested_2tcm.R`** — Contains `nested_2tcm_delay()`,
   `nested_2tcm()`, `plot_nested_2tcm_delayfit()`, `plot_nested_2tcmfit()`
3. **`R/kinfitr_nested_srtm.R`** — Contains `nested_srtm()`,
   `plot_nested_srtmfit()`
4. **`tests/testthat/test-nested.R`** — Tests for all models

## Existing Functions to Reuse

| Function | File | Used for |
|---|---|---|
| `tidyinput_long()` | `R/kinfitr_miscfuncs.R` | Data validation, zero-frame, weights |
| `shift_timings_long()` | `R/kinfitr_bloodfuncs.R` | Shifting input by inpshift |
| `onetcm_model()` | `R/kinfitr_1tcm.R` | 1TCM predicted values |
| `twotcm_model()` | `R/kinfitr_2tcm.R` | 2TCM predicted values |
| `twotcm_macro_model()` | `R/kinfitr_2tcm_macro.R` | 2TCM macro predicted values |
| `srtm2_model()` | `R/kinfitr_srtm2.R` | SRTM2 predicted values |
| `get_se()` | `R/kinfitr_miscfuncs.R` | Delta method SEs |
| `fix_multstartpars()` | `R/kinfitr_miscfuncs.R` | Multistart parameter handling |
| `plot_kinfit()` | `R/kinfitr_miscfuncs.R` | Generic plot dispatch |

## Verification

1. `devtools::document()` — regenerate NAMESPACE with new exports
2. `devtools::test(filter = "nested")` — run nested model tests
3. Manual: compare `nested_1tcm_delay` inpshift to single-region `onetcm()` fits
4. Manual: compare `nested_srtm` k2prime to single-region `srtm2()` fits
5. `devtools::check()` — full package check
