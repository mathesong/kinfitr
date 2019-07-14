library(tidyverse)

load("data/pbr28.RData")

set.seed(123)

tacs <- pbr28$tacs[[1]]
input <- pbr28$input[[1]]

noisify_tac <- function(t_start, t_end, tac, maxprop) {

  t_tac <- t_start + (t_end - t_start)/2
  weights <- weights_create(t_start, t_end, tac, "C11")
  stdev <- sqrt(1/weights)
  stdev[1] <- 0
  stdev <- stdev / max(stdev)
  stdev <- maxprop * max(tac) * stdev
  noise <- rnorm(length(t_tac), mean=tac, sd = stdev)
  tac <- tac + noise
  tac[1] <- 0

  return(tac)

}

sim_fake_reftacs <- function(tacs, input) {

  times <- tacs$Times/60
  t_start <- tacs$StartTime/60
  t_end <- tacs$StartTime + tacs$Duration/60
  weights <- tacs$Weights
  onetcm_fit <- onetcm(times, tacs$FC, input)

  reftac <- onetcm_model(onetcm_fit$tacs$Time, onetcm_fit$input, K1=0.5, k2=0.2, vB=0.05)

  targ1 <- srtm_model(times, reftac,
                      R1 = rnorm(1, 1.1, 0.1),
                      k2 = rnorm(1, 0.1, 0.001),
                      bp = rnorm(1, 1.5, 0.15))

  targ2 <- srtm_model(t_tac = times, reftac = reftac,
                      R1 = rnorm(1, 1.1, 0.1),
                      k2 = rnorm(1, 0.1, 0.001),
                      bp = rnorm(1, 0.8, 0.08))

  targ3 <- srtm_model(times, reftac,
                      R1 = rnorm(1, 1.1, 0.1),
                      k2 = rnorm(1, 0.1, 0.001),
                      bp = rnorm(1, 0.3, 0.03))

  weights_create(t_start, t_end, targ1, "C11")

  reftac_noisy <- abs(noisify_tac(t_start, t_end, reftac, 0.02))
  targ1_noisy <- abs(noisify_tac(t_start, t_end, targ1, 0.02))
  targ2_noisy <- abs(noisify_tac(t_start, t_end, targ2, 0.02))
  targ3_noisy <- abs(noisify_tac(t_start, t_end, targ3, 0.02))

  weights <- weights_create(t_start, t_end, targ1_noisy, radioisotope = "C11")

  out <- tibble::tibble(
    Times = times,
    Reference = reftac_noisy,
    ROI1 = targ1_noisy,
    ROI2 = targ2_noisy,
    ROI3 = targ3_noisy,
    Weights = weights,
    StartTime = t_start,
    Duration = tacs$Duration/60
  )

  return(out)

}

subjnames <- stringi::stri_rand_strings( 100, 4, pattern = "[a-z]")
subjnames <- subjnames[!duplicated(subjnames)][1:nrow(pbr28)]

sim_refdata <- pbr28 %>%
  mutate(simtacs = map2(tacs, input, ~sim_fake_reftacs(.x, .y)))

sim_refdata <- sim_refdata %>%
  mutate(Subjname = subjnames, PET = paste0(Subjname, "_1"),
         PETNo = 1) %>%
  select(Subjname, PETNo, PET, tacs = simtacs)

simref <- sim_refdata

save(simref, file = "data/simref.RData", compress = "xz")
