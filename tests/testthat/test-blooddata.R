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


# Blood/AIF Model Tests ----------------------------------------------


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

