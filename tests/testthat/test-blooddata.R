context("test-blooddata")

data(pbr28)

blooddata <- create_blooddata_bids(pbr28$jsondata[[1]])

blooddata_nocont <- blooddata
blooddata_nocont$Data$Blood$Continuous <- NULL

test_that("creating blooddata from BIDS works", {
  blooddata <- create_blooddata_bids(pbr28$jsondata[[1]])
  expect_true(class(blooddata) == "blooddata")
})


test_that("creating blooddata from vectors works", {
  blooddata2 <- create_blooddata_components(
    Blood.Discrete.Data.Values.sampleStartTime =
      blooddata$Data$Blood$Discrete$Data$Values$sampleStartTime,
    Blood.Discrete.Data.Values.sampleDuration =
      blooddata$Data$Blood$Discrete$Data$Values$sampleDuration,
    Blood.Discrete.Data.Values.activity =
      blooddata$Data$Blood$Discrete$Data$Values$activity,
    Plasma.Data.Values.sampleStartTime =
      blooddata$Data$Plasma$Data$Values$sampleStartTime,
    Plasma.Data.Values.sampleDuration =
      blooddata$Data$Plasma$Data$Values$sampleDuration,
    Plasma.Data.Values.activity =
      blooddata$Data$Plasma$Data$Values$activity,
    Metabolite.Data.Values.sampleStartTime =
      blooddata$Data$Metabolite$Data$Values$sampleStartTime,
    Metabolite.Data.Values.sampleDuration =
      blooddata$Data$Metabolite$Data$Values$sampleDuration,
    Metabolite.Data.Values.parentFraction =
      blooddata$Data$Metabolite$Data$Values$parentFraction,
    Blood.Continuous.Data.Values.time =
      blooddata$Data$Blood$Continuous$Data$Values$time,
    Blood.Continuous.Data.Values.activity =
      blooddata$Data$Blood$Continuous$Data$Values$activity,
    Blood.Continuous.WithdrawalRate =
      blooddata$Data$Blood$Continuous$WithdrawalRate,
    Blood.Continuous.DispersionConstant =
      blooddata$Data$Blood$Continuous$DispersionConstant,
    Blood.Continuous.DispersionConstantUnits =
      blooddata$Data$Blood$Continuous$DispersionConstantUnits,
    Blood.Continuous.DispersionCorrected = FALSE,
    TimeShift = 0
  )
  expect_true(class(blooddata2) == "blooddata")
})

test_that("plotting blooddata works", {
  bdplot <- plot(blooddata)
  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("getting data from blooddata works", {
  blood <- blooddata_getdata(blooddata, output = "Blood")
  expect_true(any(class(blood) == "tbl"))

  bpr <- blooddata_getdata(blooddata, output = "BPR")
  expect_true(any(class(bpr) == "tbl"))

  pf <- blooddata_getdata(blooddata, output = "parentFraction")
  expect_true(any(class(pf) == "tbl"))

  aif <- blooddata_getdata(blooddata, output = "AIF")
  expect_true(any(class(aif) == "tbl"))

  input <- blooddata_getdata(blooddata)
  expect_true(any(class(aif) == "tbl"))
})


test_that("addfit works", {
  pf <- blooddata_getdata(blooddata, output = "parentFraction")
  pf_fit <- metab_hillguo(pf$time, pf$parentFraction)
  blooddata <- blood_addfit(blooddata, fit = pf_fit, modeltype = "parentFraction")

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitted works", {
  pf <- blooddata_getdata(blooddata, output = "parentFraction")
  pf_fit <- metab_hillguo(pf$time, pf$parentFraction)

  fitted <- tibble::tibble(
    time = seq(min(pf$time), max(pf$time), length.out = 100)
  )
  fitted$pred <- predict(pf_fit, newdata = list(time = fitted$time))

  blooddata <- blood_addfitted(blooddata,
    time = fitted$time,
    predicted = fitted$pred,
    modeltype = "parentFraction"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("addfitpars works", {
  pf <- blooddata_getdata(blooddata, output = "parentFraction")
  pf_fit <- metab_hillguo(pf$time, pf$parentFraction)

  fitpars <- as.list(coef(pf_fit))

  blooddata <- blood_addfitpars(blooddata,
    modelname = "metab_hillguo_model", fitpars = fitpars,
    modeltype = "parentFraction"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})

test_that("bloodsplines works", {
  blood <- blooddata_getdata(blooddata, output = "Blood")
  blood_fit <- blmod_splines(blood$time,
    blood$activity,
    Method = blood$Method
  )

  blooddata <- blood_addfit(blooddata,
    fit = blood_fit,
    modeltype = "Blood"
  )

  bdplot <- plot(blooddata)

  expect_true(any(class(bdplot) == "ggplot"))
})
