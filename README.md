
[![R-CMD-check](https://github.com/mathesong/kinfitr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mathesong/kinfitr/actions/workflows/R-CMD-check.yaml)
[![Coverage
Status](https://img.shields.io/codecov/c/github/mathesong/kinfitr/master.svg)](https://codecov.io/github/mathesong/kinfitr?branch=master)
[![Launch Rstudio
Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mathesong/kinfitr_vignettes/master?urlpath=rstudio)

# kinfitr : PET Kinetic Modelling using R

## Overview

kinfitr is a package for PET Kinetic Modelling Using R. The main goal of
this package is to equip PET modellers with great flexibility, while
simultaneously making it easier to produce, present and share their
results in a highly transparent manner using R and its ecosystem of
tools for computational reproducibility.

## Installation

This package is currently only available on GitHub. It can be downloaded
as follows:

``` r
remotes::install_github("mathesong/kinfitr")
```

Additionally, this package can be cloned with git and setup with renv:

```bash
git clone https://github.com/mathesong/kinfitr.git
cd kinfitr
Rscript -e "renv::install()"
```

There is also a docker container if you would prefer not to have to
download everything. If you download Docker, you can pull the container
and start it with the following commands. The rstudio session can be
opened in a browser pointing to `http://localhost:8787`.

    docker pull mathesong/kinfitr_docker
    docker run -e PASSWORD=XYZ --name rstudio -p 8787:8787 mathesong/kinfitr_docker

Alternatively, if you want to mess around with the kinfitr package and
its vignettes in order to get a feel for it before installing it, click
the “launch binder” button above to open the vignettes in a binder
instance.

## Citation

At present, there are two preprints out about *kinfitr*. Please cite at
least one of these if you use *kinfitr* in your study.

An introduction to the package:

> Matheson, G. J. (2019). *Kinfitr: Reproducible PET Pharmacokinetic
> Modelling in R*. bioRxiv: 755751. <https://doi.org/10.1101/755751>

A validation study compared against commercial software:

> Tjerkaski, J., Cervenka, S., Farde, L., & Matheson, G. J. (2020).
> *Kinfitr – an open source tool for reproducible PET modelling:
> Validation and evaluation of test-retest reliability*. bioRxiv:
> 2020.02.20.957738. <https://doi.org/10.1101/2020.02.20.957738>

## Background and Usage

I’ve written up a series of blog posts describing some PET modelling
theory and demonstrating how the package can be used. These will
eventually be included as vignettes in the package. You can launch them
in a cloud instance to play around with them yourself by clicking the
“launch binder” button above.

  - [Part 1: PET Modelling
    Theory](https://www.granvillematheson.com/post/pharmacokinetic-modelling-of-pet-data-in-r-using-kinfitr-part-1-theory)

  - [Part 2: Basics and Iteration using
    kinfitr](https://www.granvillematheson.com/post/pharmacokinetic-modelling-of-pet-data-in-r-using-kinfitr-part-2-basics-and-iteration)

  - [Part 3: Finding tstar using
    kinfitr](https://www.granvillematheson.com/post/pharmacokinetic-modelling-of-pet-data-in-r-using-kinfitr-part-3-finding-tstar)

  - [Part 4: Blood Processing using
    kinfitr](https://www.granvillematheson.com/post/pharmacokinetic-modelling-of-pet-data-in-r-using-kinfitr-part-4-blood-processing)

<!-- ## Example Usage -->

<!-- ### Data Structure -->

<!-- The optimal data structure for _kinfitr_ is that of a nested tibble, with rows representing the desired level of chunking, e.g. whether modelling be across individuals, or across ROIs within individuals). The package contains two datasets: `pbr28` containing [$^{11}$C]PBR28 TACs for testing models involving arterial input function, and `simref` containing simulated TACs of a tracer with a reference region, for testing reference tissue models. -->

<!-- ```{r, message=F} -->

<!-- library(tidyverse) -->

<!-- library(kinfitr) -->

<!-- library(knitr) -->

<!-- data(simref) -->

<!-- ``` -->

<!-- Thus the data looks as follows: -->

<!-- ```{r} -->

<!-- head(simref) -->

<!-- ``` -->

<!-- ...and inside each nested tibble is the following: -->

<!-- ```{r} -->

<!-- head(simref$tacs[[1]]) -->

<!-- ``` -->

<!-- ### Fitting a Model for a single TAC -->

<!-- As a conscious decision, almost all input arguments of times or radioactivity concentrations are as numeric vectors. This allows the functions to be used at any stage of an analysis, and does not require complicated data structures created in earlier steps. So let's create vectors and run a model. -->

<!-- ```{r srtmfit} -->

<!-- times <- simref$tacs[[1]]$Times -->

<!-- tac <- simref$tacs[[1]]$ROI1 -->

<!-- reference <- simref$tacs[[1]]$Reference -->

<!-- weights <- simref$tacs[[1]]$Weights -->

<!-- srtmfit <- srtm(t_tac = times, reftac = reference, -->

<!--       roitac = tac,weights = weights) -->

<!-- plot_kinfit(srtmfit) -->

<!-- ``` -->

<!-- ### Fitting a Model to Many TACs -->

<!-- I recommend iterating through multiple TACs using the `purrr` package to go through the nested tibble.  Let's also examine how we might decide to chunk the data, using MRTM1 and MRTM2. -->

<!-- What we want to do from here is to model the data using MRTM1 and MRTM2.  Our plan is as follows: -->

<!-- * MRTM1 fits BP~ND~ and k2' -->

<!-- * MRTM2 fits BP~ND~ (using k2' from MRTM1 from a high-binding region) -->

<!-- So what we'll do: -->

<!-- * Fit MRTM1 to one region of each PET Measurement -->

<!-- * Fit MRTM2 to all regions of each PET Measurement -->

<!-- #### Fitting k2prime using MRTM1 using purrr::map -->

<!-- ```{r mrtm2fit} -->

<!-- simref <- simref %>% -->

<!--   group_by(PET) %>%     # Group by each PET measurement -->

<!--   mutate(mrtm1_fit = map(tacs, ~mrtm1(t_tac = .x$Times, reftac = .x$Reference,      # Add MRTM1 fit column -->

<!--                                       roitac = .x$ROI1, weights = .x$Weights))) %>% -->

<!--   mutate(k2prime = map_dbl(mrtm1_fit, c('par', 'k2prime')))     # Extract k2prime from the fit output -->

<!-- plot_kinfit(simref$mrtm1_fit[[1]])     # Plot the first TAC -->

<!-- ``` -->

<!-- Now we have a k2' value for each individual from a single region (ROI1). -->

<!-- ```{r} -->

<!-- head(simref) -->

<!-- ``` -->

<!-- #### Tidy Data: Gathering into Long format -->

<!-- Now we want to use the k2prime from fitting MRTM1 to one region from each PET measurement.  Now, we want to chunk the data a little bit differently: we want to make the arrangement a little bit longer: each TAC which we wish to model should be a row. -->

<!-- ```{r} -->

<!-- long_simref <- simref %>% -->

<!--   select(-mrtm1_fit) %>% # Remove the fit object -->

<!--   unnest() %>%    # Unnest -->

<!--   select(-StartTime, -Duration) %>%  -->

<!--   gather(Region, TAC, -PET, -Subjname, -PETNo, -Weights,  -->

<!--          -Times, -k2prime, -Reference) %>%    # Gather the data into long-er format -->

<!--   group_by(PET, Subjname, PETNo, Region, k2prime) %>%    # Group by more than just PET -->

<!--   nest(.key=tacs) %>%       # Nest the data again %>%  -->

<!--   arrange(Subjname, Region) -->

<!-- ``` -->

<!-- This produces data which looks like this: -->

<!-- ```{r} -->

<!-- head(long_simref) -->

<!-- ``` -->

<!-- For which the TACs nested object looks like this: -->

<!-- ```{r} -->

<!-- head( long_simref$tacs[[1]] ) -->

<!-- ``` -->

<!-- #### Fitting MRTM2 using purrr::pmap -->

<!-- Now we can iterate through this list, using the fitted k2prime. We can't use `map()` any longer, as we'll be iterating through 2 lists.  For this, we could either use `map2()`, or we can use `pmap()`, which allows us to iterate through an unlimited numberof columns.  Let's go with the latter for now to show how we could do this. -->

<!-- First we define a function for the iteration: -->

<!-- ```{r} -->

<!-- dofit_mrtm2 <- function(tacs, k2prime) { -->

<!--   mrtm2(t_tac = tacs$Times, reftac = tacs$Reference,  -->

<!--             roitac = tacs$TAC, weights = tacs$Weights, -->

<!--             k2prime = k2prime) -->

<!-- } -->

<!-- ``` -->

<!-- ... and then we apply it: -->

<!-- ```{r} -->

<!-- long_simref <- long_simref %>% -->

<!--   mutate(fit_mrtm2 = pmap(list(tacs, k2prime), dofit_mrtm2)) %>% -->

<!--   mutate(bp_mrtm2 = map_dbl(fit_mrtm2, c('par', 'bp'))) -->

<!-- plot_kinfit(long_simref$fit_mrtm2[[6]]) -->

<!-- ``` -->

<!-- In this way, we can run several different models on the data, keep the fits in their own separate columns, and plot them at will.  And we can do all of this within the tidyverse paradigm. -->

## Implemented Models

**Reference Region Models**

  - Simplified Reference Tissue Model (SRTM) *(Lammertsma & Hume, 1996)*

  - Non-Invasive Logan Plot *(Logan et al., 1996)*

  - Non-Invasive Multilinear Logan Plot *(Turkheimer et al., 2003)*

  - Ichise Multilinear Reference Tissue Model 1 (MRTM1) *(Ichise et al.,
    2003)*

  - Ichise Multilinear Reference Tissue Model 2 (MRTM2) *(Ichise et al.,
    2003)*

  - Patlak Reference Tissue Model *(Patlak & Blasberg, 1985)*

**Models Requiring Arterial Input**

  - One-Tissue Compartment Model

  - Two-Tissue Compartment Model

  - Logan Plot *(Logan et al., 1990)*

  - Multilinear Logan Plot *(Turkheimer et al., 2003)*

  - Ichise Multilinear Analysis 1 (MA1) *(Ichise et al., 2002)*

  - Ichise Multilinear Analysis 2 (MA2) *(Ichise et al., 2002)*

  - Patlak Plot *(Patlak et al., 1983)*

**Other Models**

  - Simultaneous Estimation of Non-Displaceable Binding (SIME) *(Ogden
    et al., 2015)*

## To-Do

**General**

  - Vignettes

**Additions**

**Improvements**

  - Tidy up
    
      - Documentation still quite rough and several inconsistencies
    
      - T-star finders have lots of code duplication: should be more
        generic

  - Revise vignette

  - Update plyr functions to dplyr functions

  - Get the package CRAN-ready
