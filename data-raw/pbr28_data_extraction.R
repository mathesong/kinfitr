# Read in the data and save to a nested data.frame

library(tidyverse)
library(readr)
library(kinfitr)

demog <- read_csv("data-raw/pbr28_demographics.csv") %>%
  select(
    Subjname, Genotype,
    MBq_PET1, MBq_PET2
  ) %>%
  gather(PETNo, injRad, contains("MBq")) %>%
  mutate(PETNo = as.numeric(str_match(PETNo, "\\d"))) %>%
  mutate(PET = paste(Subjname, PETNo, sep = "_"))


tacdat <- read_csv("data-raw/pbr28_tacdata.csv") %>%
  mutate_at(
    .vars = vars(FC:CBL),
    .funs = ~ . * 0.037
  ) %>%
  group_by(PET) %>%
  nest(tacs = -PET)

procblood <- read_csv("data-raw/pbr28_blooddata.csv") %>%
  mutate(
    Cbl_dispcorr = ifelse(Cbl_dispcorr < 0, 0, Cbl_dispcorr),
    Cpl_metabcorr = ifelse(Cpl_metabcorr < 0, 0, Cpl_metabcorr)
  ) %>%
  mutate(
    Cbl_dispcorr = Cbl_dispcorr * 0.037,
    Cpl_metabcorr = Cpl_metabcorr * 0.037
  ) %>%
  group_by(PET) %>%
  nest(procblood = -PET) %>%
  mutate(input = map(procblood,
                     ~blood_interp(t_blood = .x$Time / 60,
                                   blood = .x$Cbl_dispcorr,
                                   t_plasma = .x$Time / 60,
                                   plasma = .x$Cpl_metabcorr,
                                   t_parentfrac = 1, parentfrac = 1
  )))

petdat <- inner_join(tacdat, procblood) %>%
  separate(PET, c("Subjname", "PETNo"), sep = "_", remove = F, convert = T)






jsondat <- readRDS("data-raw/pbr28_jsondata.rds")

jsondat2petinfo <- function(j) {
  p <- list()

  p$Manufacturer <- "Siemens"
  p$ManufacturersModelName <- "High-Resolution Research Tomograph (HRRT, CTI/Siemens)"
  p$BodyPart <- "Brain"
  p$Unit <- "kBq/ml"
  p$TracerName <- "PBR28"
  p$TracerRadionuclide <- "C11"
  p$TracerMolecularWeight <- 346.41
  p$TracerMolecularWeightUnit <- "g/mol"
  p$InjectedRadioactivity <- j$Radiochem$InjectedRadioactivity
  p$InjectedRadioactivityUnit <- j$Radiochem$InjectedRadioactivityUnits
  p$InjectedMass <- NA
  p$InjectedMassUnit <- j$Radiochem$InjectedMassUnits
  p$SpecificRadioactivity <- NA
  p$SpecificRadioactivityUnit <- "GBq/ug"
  p$ModeOfAdministration <- "bolus"
  p$MolarActivity <- NA
  p$MolarActivityUnit <- "Bq/mol"
  p$MolarActivityMeasTime <- NA
  p$TimeZero <- "12:00:00"
  p$ScanStart <- 0
  p$InjectionStart <- 0
  p$FrameTimesStart <- j$Time$FrameTimes$Values[,1]
  p$FrameDuration <- j$Time$FrameTimes$Values[,2] - j$Time$FrameTimes$Values[,1]
  p$AcquisitionMode <- "list mode"
  p$ImageDecayCorrected <- TRUE
  p$ImageDecayCorrectionTime <- 0
  p$ReconMatrixSize <- NA
  p$ImageVoxelSize <- NA
  p$ReconMethodName <- "3D-OSEM-PSF"
  p$ReconMethodParameterLabels <- c("subsets", "iterations")
  p$ReconMethodParameterUnit <- c("none", "none")
  p$ReconMethodParameterValues <- c(10, 16)
  p$ReconFilterType <- "none"
  p$ReconFilterSize <- 0
  p$AttenuationCorrection <- "[137Cs]transmission scan-based"
  p$PlasmaAvail <- TRUE
  p$MetaboliteAvail <- TRUE
  p$MetaboliteMethod <- "HPLC"
  p$MetaboliteRecoveryCorrectionApplied <- TRUE
  p$ContinuousBloodAvail <- TRUE
  p$ContinuousBloodDispersionCorrected <- FALSE
  p$DiscreteBloodAvail <- TRUE

  return(p)

}

jsondat2bd <- function(j) {

  b <- create_blooddata_components(
    Blood.Discrete.Values.time = j$Blood$Discrete$Data$Values[,1],
    Blood.Discrete.Values.activity = j$Blood$Discrete$Data$Values[,3],
    Plasma.Values.time = j$Plasma$Data$Values[,1],
    Plasma.Values.activity = j$Plasma$Data$Values[,3],
    Metabolite.Values.time = j$Metabolite$Data$Values[,1],
    Metabolite.Values.parentFraction = j$Metabolite$Data$Values[,3],
    Blood.Continuous.Values.time = j$Blood$Continuous$Data$Values[,1],
    Blood.Continuous.Values.activity = j$Blood$Continuous$Data$Values[,2],
    Blood.Continuous.DispersionConstant = j$Blood$Continuous$DispersionConstant,
    Blood.Continuous.DispersionCorrected = FALSE,
    TimeShift = 0)


  b$Data$Blood$Continuous$activity$Description <- "Radioactivity in uncorrected whole blood samples from Allogg autosampler."
  b$Data$Blood$Continuous$time$Description <- "Time in relation to time zero defined by the _pet.json"

  b$Data$Blood$Discrete$activity$Description <- "Radioactivity in whole blood samples."
  b$Data$Blood$Discrete$time$Description <- "Time in relation to time zero defined by the _pet.json"

  b$Data$Plasma$activity$Description <- "Radioactivity in plasma samples."
  b$Data$Plasma$activity$Units <- "Time in relation to time zero defined by the _pet.json"

  b$Data$Metabolite$parentFraction$Description <- "Parent fraction of the radiotracer."
  b$Data$Metabolite$time$Description <- "Time in relation to time zero defined by the _pet.json"
  b$Data$Metabolite$Method <- "HPLC"
  b$Data$Metabolite$RecoveryCorrectionApplied <- TRUE

  return(b)

}

jsondat2tactimes <- function(j) {

  t <- tibble::tibble(
    start = j$Time$FrameTimes$Values[,1],
    dur = j$Time$FrameTimes$Values[,2] - j$Time$FrameTimes$Values[,1],
    time = j$Time$FrameTimes$Values[,1] + 0.5*dur
  )

  return(t)

}

blooddata2input <- function(bd) {
  bd <- bd_blood_dispcor(bd)
  input <- bd_create_input(bd)
  return(input)
}

bidsdat <- jsondat %>%
  group_by(Subjname, PETNo) %>%
  mutate(petinfo = map(jsondata,   jsondat2petinfo),
         blooddata = map(jsondata, jsondat2bd),
         tactimes = map(jsondata,  jsondat2tactimes)) %>%
  select(-jsondata)


pbr28 <- petdat %>%
  inner_join(demog) %>%
  inner_join(bidsdat) %>%
  arrange(PET)

save(pbr28, file = "data/pbr28.RData", compress = "xz")



oldbids_json <- jsondat$jsondata[[1]]

save(oldbids_json, file = "data/oldbids_json.RData", compress = "xz")
