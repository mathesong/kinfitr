context("test-blooddata")

data(pbr28)
data(oldbids_json)

suppressMessages(
  blooddata_old <- create_blooddata_bids(oldbids_json) )

blooddata <- pbr28$blooddata[[1]]

blooddata <- bd_blood_dispcor(blooddata)


test_that("plotting blooddata works", {
  bdplot <- plot(blooddata)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("creating blooddata from vectors works", {

  blooddata2 <- create_blooddata_components(
     Blood.Discrete.Values.time =
       blooddata$Data$Blood$Discrete$Values$time,
     Blood.Discrete.Values.activity =
       blooddata$Data$Blood$Discrete$Values$activity,
     Plasma.Values.time =
       blooddata$Data$Plasma$Values$time,
     Plasma.Values.activity =
       blooddata$Data$Plasma$Values$activity,
     Metabolite.Values.time =
       blooddata$Data$Metabolite$Values$time,
     Metabolite.Values.parentFraction =
       blooddata$Data$Metabolite$Values$parentFraction,
     Blood.Continuous.Values.time =
       blooddata$Data$Blood$Continuous$Values$time,
     Blood.Continuous.Values.activity =
       blooddata$Data$Blood$Continuous$Values$activity,
     Blood.Continuous.DispersionConstant =
       blooddata$Data$Blood$Continuous$DispersionConstant,
     Blood.Continuous.DispersionCorrected = FALSE,
     TimeShift = 0)

  expect_true(class(blooddata2) == "blooddata")

  bdplot <- plot(blooddata2)
  expect_true(any(class(bdplot) == "ggplot"))
})


test_that("blooddata from vectors with no continuous works", {

  blood_discrete <- tibble::tibble(
    time = c(blooddata$Data$Blood$Discrete$Values$time,  # d & c as discrete
             blooddata$Data$Blood$Continuous$Values$time),
    activity = c(blooddata$Data$Blood$Discrete$Values$activity,  # d & c as discrete
                 blooddata$Data$Blood$Continuous$Values$activity)
  ) %>%
    dplyr::distinct(time, .keep_all = TRUE)

  blooddata2 <- create_blooddata_components(
    Blood.Discrete.Values.time =
      blood_discrete$time,
    Blood.Discrete.Values.activity =
      blood_discrete$activity,
    Plasma.Values.time =
      blooddata$Data$Plasma$Values$time,
    Plasma.Values.activity =
      blooddata$Data$Plasma$Values$activity,
    Metabolite.Values.time =
      blooddata$Data$Metabolite$Values$time,
    Metabolite.Values.parentFraction =
      blooddata$Data$Metabolite$Values$parentFraction,
    #Blood.Continuous.Values.time =
    #  blooddata$Data$Blood$Continuous$Values$time,
    #Blood.Continuous.Values.activity =
    #  blooddata$Data$Blood$Continuous$Values$activity,
    #Blood.Continuous.DispersionConstant =
    #  blooddata$Data$Blood$Continuous$DispersionConstant,
    #Blood.Continuous.DispersionCorrected = FALSE,
    TimeShift = 0)

  expect_true(class(blooddata2) == "blooddata")

  bdplot <- plot(blooddata2)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("blooddata with missing plasma works", {

  blooddata2 <- create_blooddata_components(
    Blood.Discrete.Values.time =
      blooddata$Data$Blood$Discrete$Values$time,
    Blood.Discrete.Values.activity =
      blooddata$Data$Blood$Discrete$Values$activity,
    #Plasma.Values.time =
    #  blooddata$Data$Plasma$Values$time,
    #Plasma.Values.activity =
    #  blooddata$Data$Plasma$Values$activity,
    Metabolite.Values.time =
      blooddata$Data$Metabolite$Values$time,
    Metabolite.Values.parentFraction =
      blooddata$Data$Metabolite$Values$parentFraction,
    Blood.Continuous.Values.time =
      blooddata$Data$Blood$Continuous$Values$time,
    Blood.Continuous.Values.activity =
      blooddata$Data$Blood$Continuous$Values$activity,
    Blood.Continuous.DispersionConstant =
      blooddata$Data$Blood$Continuous$DispersionConstant,
    Blood.Continuous.DispersionCorrected = FALSE,
    TimeShift = 0)

  expect_true(class(blooddata2) == "blooddata")

  bdplot <- plot(blooddata2)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("blooddata with missing metabolite works", {

  blooddata2 <- create_blooddata_components(
    Blood.Discrete.Values.time =
      blooddata$Data$Blood$Discrete$Values$time,
    Blood.Discrete.Values.activity =
      blooddata$Data$Blood$Discrete$Values$activity,
    Plasma.Values.time =
     blooddata$Data$Plasma$Values$time,
    Plasma.Values.activity =
     blooddata$Data$Plasma$Values$activity,
    # Metabolite.Values.time =
    #   blooddata$Data$Metabolite$Values$time,
    # Metabolite.Values.parentFraction =
    #   blooddata$Data$Metabolite$Values$parentFraction,
    Blood.Continuous.Values.time =
      blooddata$Data$Blood$Continuous$Values$time,
    Blood.Continuous.Values.activity =
      blooddata$Data$Blood$Continuous$Values$activity,
    Blood.Continuous.DispersionConstant =
      blooddata$Data$Blood$Continuous$DispersionConstant,
    Blood.Continuous.DispersionCorrected = FALSE,
    TimeShift = 0)

  expect_true(class(blooddata2) == "blooddata")

  bdplot <- plot(blooddata2)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("blooddata with missing WB works", {

  blooddata2 <- create_blooddata_components(
    #Blood.Discrete.Values.time =
    #  blooddata$Data$Blood$Discrete$Values$time,
    #Blood.Discrete.Values.activity =
    #  blooddata$Data$Blood$Discrete$Values$activity,
    Plasma.Values.time =
      blooddata$Data$Plasma$Values$time,
    Plasma.Values.activity =
      blooddata$Data$Plasma$Values$activity,
    Metabolite.Values.time =
      blooddata$Data$Metabolite$Values$time,
    Metabolite.Values.parentFraction =
      blooddata$Data$Metabolite$Values$parentFraction,
    # Blood.Continuous.Values.time =
    #   blooddata$Data$Blood$Continuous$Values$time,
    # Blood.Continuous.Values.activity =
    #   blooddata$Data$Blood$Continuous$Values$activity,
    # Blood.Continuous.DispersionConstant =
    #   blooddata$Data$Blood$Continuous$DispersionConstant,
    # Blood.Continuous.DispersionCorrected = FALSE,
    TimeShift = 0)

  expect_true(class(blooddata2) == "blooddata")

  bdplot <- plot(blooddata2)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("blooddata with missing WB and metabolite works", {

  blooddata2 <- create_blooddata_components(
    #Blood.Discrete.Values.time =
    #  blooddata$Data$Blood$Discrete$Values$time,
    #Blood.Discrete.Values.activity =
    #  blooddata$Data$Blood$Discrete$Values$activity,
    Plasma.Values.time =
      blooddata$Data$Plasma$Values$time,
    Plasma.Values.activity =
      blooddata$Data$Plasma$Values$activity,
    # Metabolite.Values.time =
    #   blooddata$Data$Metabolite$Values$time,
    # Metabolite.Values.parentFraction =
    #   blooddata$Data$Metabolite$Values$parentFraction,
    # Blood.Continuous.Values.time =
    #   blooddata$Data$Blood$Continuous$Values$time,
    # Blood.Continuous.Values.activity =
    #   blooddata$Data$Blood$Continuous$Values$activity,
    # Blood.Continuous.DispersionConstant =
    #   blooddata$Data$Blood$Continuous$DispersionConstant,
    # Blood.Continuous.DispersionCorrected = FALSE,
    TimeShift = 0)

  expect_true(class(blooddata2) == "blooddata")

  bdplot <- plot(blooddata2)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("updating blooddata works", {

  blooddata2 <- update_blooddata(blooddata_old)

  expect_true(class(blooddata2) == "blooddata")

  bdplot <- plot(blooddata2)
  expect_true(any(class(bdplot) == "ggplot"))
})


test_that("getting data from blooddata works", {
  blood <- bd_extract(blooddata, output = "Blood")
  expect_true(any(class(blood) == "tbl"))

  bpr <- bd_extract(blooddata, output = "BPR")
  expect_true(any(class(bpr) == "tbl"))

  pf <- bd_extract(blooddata, output = "parentFraction")
  expect_true(any(class(pf) == "tbl"))

  aif <- bd_extract(blooddata, output = "AIF")
  expect_true(any(class(aif) == "tbl"))
})

test_that("getting input data from blooddata works", {
  input <- bd_create_input(blooddata)
  expect_true(any(class(input) == "interpblood"))
})




test_that("addfit works", {
  pf <- bd_extract(blooddata, output = "parentFraction")
  pf_fit <- metab_sigmoid(pf$time, pf$parentFraction)
  blooddata <- bd_addfit(blooddata, fit = pf_fit, modeltype = "parentFraction")

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitted works", {
  pf <- bd_extract(blooddata, output = "parentFraction")
  pf_fit <- metab_sigmoid(pf$time, pf$parentFraction)

  fitted <- tibble::tibble(
    time = seq(min(pf$time), max(pf$time), length.out = 100)
  )
  fitted$pred <- predict(pf_fit, newdata = list(time = fitted$time))

  blooddata <- bd_addfitted(blooddata,
    time = fitted$time,
    predicted = fitted$pred,
    modeltype = "parentFraction"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitpars works", {
  pf <- bd_extract(blooddata, output = "parentFraction")
  pf_fit <- metab_sigmoid(pf$time, pf$parentFraction)

  fitpars <- as.list(coef(pf_fit))

  blooddata <- bd_addfitpars(blooddata,
    modelname = "metab_sigmoid_model", fitpars = fitpars,
    modeltype = "parentFraction"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitpars works for the AIF", {
  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 1)

  bd_pars <- bd_addfitpars(blooddata,
    modelname = "blmod_triexp_model", fitpars = as.list(blood_fit$par),
    modeltype = "AIF"
  )

  # The interpolated AIF is named the same whichever Method produced it
  expect_true("aif" %in% names(bd_extract(bd_pars, output = "AIF",
                                          what = "interp")))

  # Parameters and the fit itself describe the same curve
  bd_fit <- bd_addfit(blooddata, fit = blood_fit, modeltype = "AIF")
  expect_equal(bd_create_input(bd_pars)$AIF, bd_create_input(bd_fit)$AIF)

  bdplot <- plot(bd_pars)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitpars works for the BPR", {
  # metab_sigmoid_model is used here only for its shape: kinfitr has no BPR
  # model function of its own, and what is being checked is the plumbing
  fitpars <- list(a = 1, b = 0.001, c = 1, ppf0 = 1, delay = 0)

  bd_pars <- bd_addfitpars(blooddata,
    modelname = "metab_sigmoid_model", fitpars = fitpars,
    modeltype = "BPR"
  )

  i_bpr <- bd_extract(bd_pars, output = "BPR", what = "interp")

  expect_true("bpr" %in% names(i_bpr))

  input <- bd_create_input(bd_pars)

  expect_true(all(is.finite(input$Plasma)))
  expect_equal(input$Plasma, input$Blood / i_bpr$bpr)
})


# Blood/AIF Model Tests ----------------------------------------------

test_that("the blood and metabolite models require sorted times", {
  # These models describe the curve in segments either side of t0, and return
  # the segments in time order, so unsorted times would silently be given the
  # wrong values
  tt <- c(0, 5, 20, 64, 300, 1200)
  unsorted <- tt[c(5, 2, 6, 1, 4, 3)]

  expect_error(blmod_triexp_model(unsorted, 2, 64, 126, 80, 0.08, 20, 0.005,
                                  2, 2e-4), "ascending")
  expect_error(blmod_feng_model(unsorted, 2, 80, 0.08, 20, 0.005, 2, 2e-4),
               "ascending")
  expect_error(blmod_fengconv_model(unsorted, 2, 80, 0.08, 20, 0.005, 2, 2e-4,
                                    45), "ascending")
  expect_error(blmod_fengconvplus_model(unsorted, 2, 80, 0.08, 20, 0.005, 2,
                                        2e-4, 45, 30, 0.001), "ascending")
  expect_error(metab_sigmoid_model(unsorted, 1, 0.001, 1), "ascending")
  expect_error(metab_power_model(unsorted, 1, 0.001, 1), "ascending")
  expect_error(metab_exponential_model(unsorted, 0.02, 0, 0.001), "ascending")
  expect_error(metab_gamma_model(unsorted, 1, 0.001, 1, 1), "ascending")
  expect_error(metab_invgamma_model(unsorted, 1, 0.001, 1, 1), "ascending")

  # Sorted times are unaffected
  expect_length(blmod_triexp_model(tt, 2, 64, 126, 80, 0.08, 20, 0.005,
                                   2, 2e-4), length(tt))

  # metab_hill_model handles either, so it is not restricted
  expect_equal(metab_hill_model(unsorted, 1, 0.001, 1),
               metab_hill_model(tt, 1, 0.001, 1)[c(5, 2, 6, 1, 4, 3)])
})

test_that("multstart bounds may name every starting parameter", {
  aif <- bd_extract(blooddata, output = "AIF")

  startpars <- blmod_exp_startpars(aif$time, aif$aif,
                                   fit_exp3 = TRUE,
                                   expdecay_props = c(1 / 60, 0.1))

  # peaktime and peakval are not fitted by default, so bounds naming all nine
  # starting parameters have to be pruned to those which are
  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 25,
                             multstart_lower = purrr::map(startpars,
                                                          ~(.x - abs(.x * 0.8))),
                             multstart_upper = purrr::map(startpars,
                                                          ~(.x + abs(.x * 2))))

  expect_true(all(c("A", "alpha", "B", "beta") %in% names(blood_fit$par)))

  # But a bound which is genuinely absent is still an error
  expect_error(blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 25,
                             multstart_lower = purrr::map(
                               startpars, ~(.x - abs(.x * 0.8)))[c("A", "alpha")],
                             multstart_upper = purrr::map(
                               startpars, ~(.x + abs(.x * 2)))),
               "should include a value")
})

test_that("the default bounds never cross, and crossed bounds are refused", {
  aif <- bd_extract(blooddata, output = "AIF")

  startpars <- blmod_exp_startpars(aif$time, aif$aif,
                                   fit_exp3 = TRUE,
                                   expdecay_props = c(1 / 60, 0.1))

  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 1)

  # A is bounded below at half the peak, so its ceiling has to clear that even
  # when three times the peeled starting value does not
  expect_gt(blood_fit$upper$A, blood_fit$lower$A)
  expect_gte(blood_fit$upper$A, startpars$peakval)

  expect_true(all(unlist(blood_fit$upper) >= unlist(blood_fit$lower)))

  # Bounds which cross are pinned rather than fitted, so say so
  crossed <- purrr::map(startpars, ~(.x - abs(.x * 0.8)))
  crossed$A <- 1e6

  expect_error(blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 1,
                             lower = crossed),
               "below their lower bounds")
})

test_that("a single inftime fixes ti rather than fitting it", {
  aif <- bd_extract(blooddata, output = "AIF")

  # A single value is a known infusion time, so it is substituted into the
  # model; two values are the limits within which ti should be fitted
  fixed <- blmod_fengconv(aif$time,
                          aif$aif,
                          Method = aif$Method,
                          multstart_iter = 25,
                          inftime = 45)

  expect_equal(fixed$par$ti, 45)
  expect_false("ti" %in% names(fixed$start))

  ranged <- blmod_fengconv(aif$time,
                           aif$aif,
                           Method = aif$Method,
                           multstart_iter = 25,
                           inftime = c(30, 60))

  expect_true("ti" %in% names(ranged$start))
  expect_gte(ranged$par$ti, 30)
  expect_lte(ranged$par$ti, 60)

  # Fixing ti costs one parameter fewer than fitting it
  expect_equal(ncol(ranged$start), ncol(fixed$start) + 1)

  expect_error(blmod_fengconv(aif$time,
                              aif$aif,
                              Method = aif$Method,
                              multstart_iter = 1,
                              inftime = c(1, 2, 3)),
               "single value")

  # The same holds for fengconvplus
  fixedplus <- blmod_fengconvplus(aif$time,
                                  aif$aif,
                                  Method = aif$Method,
                                  multstart_iter = 25,
                                  inftime = 45)

  expect_equal(fixedplus$par$ti, 45)
  expect_false("ti" %in% names(fixedplus$start))
})

test_that("fengconvplus bounds reach the right parameters", {
  aif <- bd_extract(blooddata, output = "AIF")

  # asymptote, slope and ti appear in the model formula in a different order
  # from the one the bounds are built in, so each was applied to whichever
  # parameter shared its position. The bounds for slope are relative to the
  # length of the measurement, so a differently scaled time axis is what
  # brings the mismatch out.
  set.seed(7)
  blood_fit <- blmod_fengconvplus(aif$time / 60,
                                  aif$aif,
                                  Method = aif$Method,
                                  multstart_iter = 50)

  pars <- names(blood_fit$start)

  expect_true(all(unlist(blood_fit$par[pars]) >=
                    unlist(blood_fit$lower[pars]) - 1e-9))
  expect_true(all(unlist(blood_fit$par[pars]) <=
                    unlist(blood_fit$upper[pars]) + 1e-9))
})

test_that("the fengconvplus rise bounds are relative to the measurement", {
  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_fengconvplus(aif$time,
                                  aif$aif,
                                  Method = aif$Method,
                                  multstart_iter = 1)

  # risefunc reaches half of asymptote log(3)/slope after t0, and the ceiling
  # is set at a hundredth of the measurement
  risedur <- max(aif$time) - blood_fit$par$t0

  expect_equal(log(3) / blood_fit$upper$slope, 0.01 * risedur,
               tolerance = 0.05)
  expect_lt(blood_fit$start$slope, blood_fit$upper$slope)

  # A starting value given by the user is not overwritten
  given <- as.list(blood_fit$start)
  given$slope <- 0.5 * blood_fit$start$slope

  expect_equal(blmod_fengconvplus(aif$time,
                                  aif$aif,
                                  Method = aif$Method,
                                  start = given,
                                  check_startpars = TRUE)$start$slope,
               given$slope)
})

test_that("bd_getdata is defunct", {
  expect_error(bd_getdata(blooddata), "defunct")
})



test_that("bloodsplines works", {
  blood <- bd_extract(blooddata, output = "Blood")
  blood_fit <- blmod_splines(blood$time,
    blood$activity,
    Method = blood$Method
  )

  blooddata <- bd_addfit(blooddata,
    fit = blood_fit,
    modeltype = "Blood"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("starting parameters for expontial when AIF contains zeros works", {

  aif <- bd_extract(blooddata, output = "AIF")

  aif$aif[100] <- 0
  aif$aif[500] <- 0
  aif$aif[length(aif$aif)] <- -1
  aif$aif[length(aif$aif)-1] <- 0

  start <- blmod_exp_startpars(aif$time,
                               aif$aif,
                               fit_exp3 = T,
                               expdecay_props = c(1/60, 0.1))



  expect_true(any(class(start) == "list"))

})

test_that("exponential works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 1)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("exponential 2exp works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                                 aif$aif,
                                 Method = aif$Method,
                                 multstart_iter = 1, fit_exp3 = F)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})


test_that("exponential peakfitting works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                                 aif$aif,
                                 Method = aif$Method,
                                 multstart_iter = 1,
                                 fit_peaktime = T, fit_peakval = T)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("multstart bounds stay matched to their parameters", {
  # Defining peakval as A+B+C puts A, B and C ahead of the rate constants in
  # the model formula, which is not the order in which the bounds are built.
  # nls_multstart takes its parameters from the formula and pairs them with
  # the bounds positionally, so the bounds have to be reordered to match: a
  # bound landing on the wrong parameter shows up as a fit which breaches the
  # bounds it was given.
  aif <- bd_extract(blooddata, output = "AIF")

  set.seed(99)
  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 25,
                             fit_t0 = FALSE,
                             peakval_set = FALSE)

  pars <- c("A", "alpha", "B", "beta", "C", "gamma")

  expect_true(all(unlist(blood_fit$par[pars]) >=
                    unlist(blood_fit$lower[pars]) - 1e-10))
  expect_true(all(unlist(blood_fit$par[pars]) <=
                    unlist(blood_fit$upper[pars]) + 1e-10))

  # If the bounds cannot be matched to parameters by name, they cannot be
  # reordered, so say so rather than fitting against the wrong ones
  unnamed <- as.list(unname(unlist(blood_fit$lower[pars])))

  expect_error(blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 1,
                             fit_t0 = FALSE,
                             peakval_set = FALSE,
                             lower = unnamed),
               "must be named")

  # But unnamed bounds are fine when no reordering is called for
  expect_no_error(blmod_exp(aif$time,
                                aif$aif,
                                Method = aif$Method,
                                multstart_iter = 1,
                                fit_t0 = FALSE,
                                lower = unnamed))
})

test_that("exponential with interpolated rise works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 1,
                             rise = "interp")

  # The nudged peak sample is used for fitting, but excluded from predictions
  expect_equal(length(predict(blood_fit)), nrow(blood_fit$blood))
  expect_lt(stats::nobs(blood_fit$fit), nrow(blood_fit$blood))

  expect_equal(blood_fit$fit_details$rise, "interp")
  expect_equal(blood_fit$par$t0, 0)
  expect_equal(blood_fit$par$peaktime, aif$time[which.max(aif$aif)])
  expect_setequal(names(blood_fit$par),
                  c("A", "alpha", "B", "beta", "C", "gamma",
                    "peaktime", "peakval", "t0"))

  # The rise passes exactly through the measured samples before the peak
  expect_equal(predict(blood_fit,
                       newdata = list(time = blood_fit$rise_samples$time)),
               blood_fit$rise_samples$activity)

  # The nudged sample anchors A+B+C close to the measured peak
  expect_lt(with(blood_fit$par, abs(A + B + C - peakval)),
            0.1 * blood_fit$par$peakval)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitpars refuses a model it cannot fully describe", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             Method = aif$Method,
                             multstart_iter = 1,
                             rise = "interp")

  # An interpolated rise is data, not parameters, so it cannot be smuggled in
  # alongside them
  expect_error(bd_addfitpars(blooddata,
    modelname = "blmod_triexp_model",
    fitpars = c(as.list(blood_fit$par),
                list(risetime = blood_fit$rise_samples$time,
                     riseval = blood_fit$rise_samples$activity)),
    modeltype = "AIF"
  ), "single value")

  # bd_addfit is the route for these fits
  blooddata <- bd_addfit(blooddata, fit = blood_fit, modeltype = "AIF")

  expect_true(any(class(plot(blooddata)) == "ggplot"))
})

test_that("interpolated rise 2exp and without method works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                             aif$aif,
                             multstart_iter = 1,
                             rise = "interp", fit_exp3 = F)

  expect_equal(ncol(blood_fit$start), 4)
  expect_equal(blood_fit$par$C, 0)
  expect_equal(blood_fit$par$gamma, 0)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("interpolated rise cannot fit the peak", {

  aif <- bd_extract(blooddata, output = "AIF")

  expect_warning(
    blood_fit <- blmod_exp(aif$time,
                               aif$aif,
                               Method = aif$Method,
                               multstart_iter = 1,
                               rise = "interp",
                               fit_peaktime = T))

  expect_false(blood_fit$fit_details$fit_peaktime)
  expect_false(blood_fit$fit_details$fit_t0)
})

test_that("interpolated rise needs samples after the peak", {

  aif <- bd_extract(blooddata, output = "AIF")

  beforepeak <- aif$time <= aif$time[which.max(aif$aif)]

  expect_error(blmod_exp(aif$time[beforepeak],
                             aif$aif[beforepeak],
                             multstart_iter = 1,
                             rise = "interp"),
               "no samples after the peak")
})

test_that("the rise argument does not change the linear rise", {

  aif <- bd_extract(blooddata, output = "AIF")

  default <- blmod_exp(aif$time,
                           aif$aif,
                           Method = aif$Method,
                           multstart_iter = 1)

  linear <- blmod_exp(aif$time,
                          aif$aif,
                          Method = aif$Method,
                          multstart_iter = 1,
                          rise = "linear")

  expect_equal(default$par, linear$par)
  expect_equal(predict(default), predict(linear))
  expect_null(default$rise_samples)
  expect_equal(length(predict(default)), nrow(default$blood))
})

test_that("exponential without method works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_exp(aif$time,
                                 aif$aif,
                                 multstart_iter = 1)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("exponential with start parameters works", {

  aif <- bd_extract(blooddata, output = "AIF")

  startpars <- blmod_exp_startpars(aif$time,
                                   aif$aif)

  blood_fit <- blmod_exp(aif$time,
                         aif$aif,
                         Method = aif$Method,
                         multstart_iter = 1,
                         start = startpars)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("Feng works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_feng(aif$time,
                          aif$aif,
                          Method = aif$Method,
                          multstart_iter = 1)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("Feng works with startpars", {

  aif <- bd_extract(blooddata, output = "AIF")

  startpars <- blmod_feng_startpars(aif$time,
                                   aif$aif)

  blood_fit <- blmod_feng(aif$time,
                         aif$aif,
                         Method = aif$Method,
                         multstart_iter = 1,
                         start = startpars)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("Fengconv works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_fengconv(aif$time,
                          aif$aif,
                          Method = aif$Method,
                          multstart_iter = 1)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("Fengconv works with startpars", {

  aif <- bd_extract(blooddata, output = "AIF")

  startpars <- blmod_feng_startpars(aif$time,
                                    aif$aif)

  blood_fit <- blmod_fengconv(aif$time,
                          aif$aif,
                          Method = aif$Method,
                          multstart_iter = 1,
                          start = startpars)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("Fengconvplus works", {

  aif <- bd_extract(blooddata, output = "AIF")

  blood_fit <- blmod_fengconvplus(aif$time,
                              aif$aif,
                              Method = aif$Method,
                              inftime = 20,
                              multstart_iter = 1)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("Fengconv works with startpars", {

  aif <- bd_extract(blooddata, output = "AIF")

  startpars <- blmod_feng_startpars(aif$time,
                                    aif$aif)

  blood_fit <- blmod_fengconvplus(aif$time,
                              aif$aif,
                              Method = aif$Method,
                              inftime = 20,
                              multstart_iter = 1,
                              start = startpars)

  blooddata <- bd_addfit(blooddata,
                         fit = blood_fit,
                         modeltype = "AIF"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})



# Parent Fraction Model Tests ----------------------------------------------

pf <- bd_extract(blooddata, "parentFraction")
set.seed(12345)

test_that("hill function works", {

  fit <- metab_hill(pf$time, pf$parentFraction)

  expect_true(any(class(fit) == "nls"))

})

test_that("exponential function works", {

  fit <- metab_exponential(pf$time, pf$parentFraction)

  expect_true(any(class(fit) == "nls"))

})

test_that("power function works", {

  fit <- metab_power(pf$time, pf$parentFraction)

  expect_true(any(class(fit) == "nls"))

})

test_that("sigmoid function works", {

  fit <- metab_sigmoid(pf$time, pf$parentFraction)

  expect_true(any(class(fit) == "nls"))

})

test_that("gamma function works", {

  fit <- metab_gamma(pf$time, pf$parentFraction)

  expect_true(any(class(fit) == "nls"))

})

test_that("invgamma function works", {

  fit <- metab_invgamma(pf$time, pf$parentFraction)

  expect_true(any(class(fit) == "nls"))

})

# test_that("exp_sep works", {
#
#   aif <- bd_extract(blooddata, output = "AIF")
#
#   blood_fit <- blmod_exp_sep(aif$time,
#                                  aif$aif,
#                                  Method = aif$Method,
#                                  multstart_iter = 1)
#
#   blooddata <- bd_addfit(blooddata,
#                          fit = blood_fit,
#                          modeltype = "AIF"
#   )
#
#   bdplot <- plot(blooddata)
#
#   expect_true(any(class(bdplot) == "ggplot"))
# })


# test_that("exp_sep without method works", {
#
#   aif <- bd_extract(blooddata, output = "AIF")
#
#   blood_fit <- blmod_exp_sep(aif$time,
#                                  aif$aif,
#                                  multstart_iter = 1)
#
#   blooddata <- bd_addfit(blooddata,
#                          fit = blood_fit,
#                          modeltype = "AIF"
#   )
#
#   bdplot <- plot(blooddata)
#
#   expect_true(any(class(bdplot) == "ggplot"))
# })

# test_that("exp_sep with start parameters works", {
#
#   aif <- bd_extract(blooddata, output = "AIF")
#
#   startpars <- blmod_exp_startpars(aif$time,
#                                    aif$aif)
#
#   blood_fit <- blmod_exp_sep(aif$time,
#                          aif$aif,
#                          Method = aif$Method,
#                          multstart_iter = 1,
#                          start = startpars)
#
#   blooddata <- bd_addfit(blooddata,
#                          fit = blood_fit,
#                          modeltype = "AIF"
#   )
#
#   bdplot <- plot(blooddata)
#
#   expect_true(any(class(bdplot) == "ggplot"))
# })


test_that("dispcor with different intervals works", {

  time <- 1:20
  activity <- rnorm(20)
  tau <- 2.5

  time <- time[-c(15, 17, 19)]
  activity <- activity[-c(15, 17, 19)]

  out <- blood_dispcor(time, activity, tau, keep_interpolated = T)

  expect_true(nrow(out)==20)
})

test_that("dispcor with different intervals works with orig times", {

  time <- 1:20
  activity <- rnorm(20)
  tau <- 2.5

  time <- time[-c(15, 17, 19)]
  activity <- activity[-c(15, 17, 19)]

  out <- blood_dispcor(time, activity, tau, keep_interpolated = F)

  expect_true(nrow(out)==17)
})

# Missing blood measurements and BPR--------------------------------------------

blooddata_pl <- blooddata

blooddata_pl$Data$Blood$Discrete$Values <- blooddata_pl$Data$Blood$Discrete$Values[sample(1:18, 5),]

test_that("plotting blooddata with missing blood samples works", {
  bdplot <- plot(blooddata_pl)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfit bpr with missing blood samples works", {
  bpr <- bd_extract(blooddata_pl, output = "BPR")
  bpr_fit <- lm(bpr ~ time, data=bpr)
  blooddata_pl <- bd_addfit(blooddata_pl, fit = bpr_fit, modeltype = "BPR")

  bdplot <- plot(blooddata_pl)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitted bpr with missing blood samples works", {
  bpr <- bd_extract(blooddata_pl, output = "BPR")
  bpr_fit <- lm(bpr ~ time, data=bpr)
  fitted <- tibble::tibble(
    time = seq(min(bpr$time), max(bpr$time), length.out = 100)
  )
  fitted$pred <- predict(bpr_fit, newdata = list(time = fitted$time))

  blooddata_pl <- bd_addfitted(blooddata_pl,
                            time = fitted$time,
                            predicted = fitted$pred,
                            modeltype = "BPR"
  )

  bdplot <- plot(blooddata_pl)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("bloodstream_import_inputfunctions pads explicit AIF when first sample is after t=0", {
  tmp <- tempfile("blstream-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  base <- file.path(tmp, "sub-01_ses-baseline_inputfunction")

  # First sample at 25.2 s (no t=0 row), mimicking PMOD .km -> BIDS output.
  utils::write.table(
    data.frame(
      time = c(25.2, 60, 120, 300, 600),
      whole_blood_radioactivity = c(10, 8, 6, 4, 2),
      plasma_radioactivity = c(12, 10, 7, 5, 2),
      metabolite_parent_fraction = c(1, 0.95, 0.85, 0.7, 0.5),
      AIF = c(12, 9.5, 5.95, 3.5, 1)
    ),
    file = paste0(base, ".tsv"), sep = "\t",
    row.names = FALSE, quote = FALSE
  )

  jsonlite::write_json(
    list(
      time = list(Units = "s"),
      whole_blood_radioactivity = list(Units = "kBq"),
      plasma_radioactivity = list(Units = "kBq"),
      metabolite_parent_fraction = list(Units = "arbitrary"),
      AIF = list(Units = "kBq")
    ),
    paste0(base, ".json"),
    auto_unbox = TRUE, pretty = TRUE
  )

  result <- suppressWarnings(bloodstream_import_inputfunctions(tmp))

  expect_equal(nrow(result), 1)
  interp <- result$input[[1]]
  expect_equal(interp$Time[1], 0)
  expect_equal(interp$AIF[1], 0)
  expect_false(any(is.na(interp$AIF)))
})
