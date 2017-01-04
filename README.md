
kinfitr : PET Kinetic Modelling using R
=======================================

Overview
--------

kinfitr is a a new package for PET Kinetic Modelling Using R. Aside from allowing me to use R more in my day-to-day work, the main goal of this package is to make it easier to perform more reproducible research in PET research. Furthermore, both kinfitr and the language it is written in, [R](https://cran.r-project.org), are free to download. Furthermore, the primary IDE for R, [RStudio](https://www.rstudio.com/), now comes with [R Notebooks](http://rmarkdown.rstudio.com/r_notebooks.html), which makes it particularly easy to generate reproducible reports.

Status
------

Please note that this package is currently under development. Many of the models have not been comprehensively tested to make sure that they produce the correct output in all situations. I plan to validate the performance of the kinfitr package against other kinetic modelling software at some point in the future. Please refer to the 'To Do' section below to check on what has not yet been done and is planned to be completed as soon as I have time.

Installation
------------

This package is currently only available on GitHub. It can be downloaded as follows:

``` r
# install.packages("devtools")  # If you do not already have devtools
devtools::install_github("mathesong/kinfitr")
```

Usage
-----

I have aimed to make the package as intuitive as possible to use. Almost all input arguments of times or radioactivity concentrations are as numeric vectors.

``` r
library(kinfitr)

srtmout <- srtm(times, reftac, roitac, weights=weights)
plot_srtmfit(srtmout)
```

Blood data is interpolated into an `input` object using the `blood_interp` command to make it easier to work with.

``` r
input <- blood_interp(t_blood = blooddata$Time.sec./60, 
                      blood = blooddata$Cbl.nCi.cc.,
                      t_plasma = plasmadata$Time.sec./60, 
                      plasma = plasmadata$Cpl.nCi.cc.,
                      t_parentfrac = parentdata$Times/60,
                      parentfrac = parentdata$Fraction)

onetcmout <- onetcm(times, tac, input, weights=weights)
plot_1tcmfit(onetcmout)
```

All the parameter output values are in data frames, and are thus easily manipulated with tidyverse tools. Tidy data is data where

1.  Each variable forms a column.
2.  Each observation forms a row.
3.  Each type of observational unit forms a table.

Thus, I can create a tidy dataframe from a list of TAC dataframes for each individual

``` r
library(tidyverse)
library(stringr)

tidytacs <- map(tacdata, c('tacdf')) %>%  # Extract the TAC dataframe from each list element
  bind_rows(.id = "id") %>%   # Add an id column of the list element name (= Acro_PETNo_results)
  mutate(Acronym = str_extract(id, "(^[a-z]*)")) %>% # Extract and add an acronym column
  mutate(PETNo = str_extract(id, "\\d")) %>% # Extract and add a PETNo column
  select(-id) %>% # Remove the id column
  gather(Region, TAC, -Acronym, -PETNo, -weights, -RefCBL, -times) %>% # Make into long format
  arrange(Acronym, PETNo, Region, times) # Order by Acronym, PET number, Region and then time

head(tidytacs)
```

|    RefCBL|   weights|      times| Acronym | PETNo | Region |       TAC|
|---------:|---------:|----------:|:--------|:------|:-------|---------:|
|    0.0000|  0.000000|  0.0000000| Subj01  | 1     | DLPFC  |    0.0000|
|  246.2982|  0.730071|  0.5002667| Subj01  | 1     | DLPFC  |  223.0276|
|  467.8367|  0.721462|  1.5002667| Subj01  | 1     | DLPFC  |  424.3305|
|  473.6503|  0.720655|  2.5002667| Subj01  | 1     | DLPFC  |  437.5260|
|  444.4454|  0.793644|  4.5002667| Subj01  | 1     | DLPFC  |  438.9326|
|  368.2536|  0.789719|  7.5002667| Subj01  | 1     | DLPFC  |  404.4193|

From here, we can model it, and extract the parameters:

``` r
# Fit SRTM to each TAC of each individual
srtmout <- tidytacs %>%
  group_by(Acronym, PETNo, Region) %>%
  do(modeloutput = srtm(t_tac = .$times, reftac=.$RefCBL, roitac=.$TAC, weights=.$weights))

# Extract a dataframe of all the BP values
srtm_pars <- srtmout$modeloutput %>%
  map('par') %>% # Extract the parameters
  do.call(rbind, .) %>% # Bind them together
  as.data.frame() %>% # Make them a data frame
  cbind( select(srtmout, -modeloutput) ) # Remove the modeloutput column

head(srtm_pars)
```

|         R1|         k2|         bp| Acronym | PETNo | Region   |
|----------:|----------:|----------:|:--------|:------|:---------|
|  0.8184882|  0.1073673|  0.2352065| Subj01  | 1     | DLPFC    |
|  0.8818911|  0.1164704|  1.5536115| Subj01  | 1     | gmfslSTR |
|  0.7996433|  0.0959518|  0.2558649| Subj01  | 2     | DLPFC    |
|  0.8972514|  0.1241690|  1.5615239| Subj01  | 2     | gmfslSTR |
|  0.8292853|  0.0997430|  0.3185725| Subj02  | 1     | DLPFC    |
|  0.8867102|  0.1088630|  1.9806780| Subj02  | 1     | gmfslSTR |

Implemented Models
------------------

**Reference Region Models**

-   Simplified Reference Tissue Model (SRTM) *(Lammertsma & Hume, 1996)*

-   Non-Invasive Logan Plot *(Logan et al., 1996)*

-   Non-Invasive Multilinear Logan Plot *(Turkheimer et al., 2003)*

-   Ichise Multilinear Reference Tissue Model 1 (MRTM1) *(Ichise et al., 2003)*

-   Ichise Multilinear Reference Tissue Model 2 (MRTM2) *(Ichise et al., 2003)*

-   Patlak Reference Tissue Model *(Patlak & Blasberg, 1985)*

**Models Requiring Arterial Input**

-   One-Tissue Compartment Model

-   Two-Tissue Compartment Model

-   Logan Plot *(Logan et al., 1990)*

-   Multilinear Logan Plot *(Turkheimer et al., 2003)*

-   Ichise Multilinear Analysis 1 (MA1) *(Ichise et al., 2002)*

-   Ichise Multilinear Analysis 2 (MA2) *(Ichise et al., 2002)*

-   Patlak Plot *(Patlak et al., 1983)*

**Other Models**

-   Simultaneous Estimation of Non-Displaceable Binding (SIME) *(Ogden et al., 2015)*

To-Do
-----

**General**

-   Validate model output against other software

    -   Reference models output very similar to our group's in-house MATLAB tools

    -   Arterial models produce reasonable values, but not yet fully validated

    -   Irreversible methods completely unvalidated (no data to try these out on)

-   Write up a tidyverse workflow for the README for models with arterial input

-   Add some sample data for testing and for the vignette

**Additions**

-   Add vB correction into the remaining linearised arterial models (only implemented in Logan Plot thus far)

-   Add function for creating weights

-   Add functions for processing blood data

    -   Combination of automatic and manual blood samples

    -   Dispersion Correction

-   Add more models

    -   More kinetic models

    -   Models of arterial input function

    -   Models of plasma parent fraction

-   Add code tests

**Improvements**

-   Steamline 1TCM, 2TCM and SIME models: currently quite slow

    -   SIME should be parallelised
-   Tidy up

    -   Functions a little messy

    -   Documentation still quite rough and several inconsistencies

    -   T-star finders have lots of code duplication: should be more generic

-   Revise vignette

-   Update plyr functions to dplyr functions
