# Read in the data and save to a nested data.frame

library(tidyverse)
library(readr)

demog <- read_csv("data-raw/pbr28_demographics.csv") %>%
  select(
    Subjname, Genotype,
    MBq_PET1, MBq_PET2) %>%
  gather(PETNo, injRad, contains("MBq")) %>%
  mutate(PETNo = as.numeric(str_match(PETNo, "\\d"))) %>%
  mutate(PET = paste(Subjname, PETNo, sep = "_"))


tacdat <- read_csv("data-raw/pbr28_tacdata.csv") %>%
  mutate_at(.vars = vars(FC:CBL),
            .funs = ~ . * 0.037 ) %>%
  group_by(PET) %>%
  nest(.key = "tacs")

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
  nest(.key = "procblood") %>%
  mutate(input = map(procblood, ~blood_interp(
    t_blood = .x$Time / 60, blood = .x$Cbl_dispcorr,
    t_plasma = .x$Time / 60, plasma = .x$Cpl_metabcorr,
    t_parentfrac = 1, parentfrac = 1
  )))

petdat <- inner_join(tacdat, procblood) %>%
  separate(PET, c("Subjname", "PETNo"), sep = "_", remove = F, convert = T)

jsondat <- readRDS('data-raw/pbr28_jsondata.rds')

pbr28 <- petdat %>%
  inner_join(demog) %>%
  inner_join(jsondat) %>%
  arrange(PET)

save(pbr28, file="data/pbr28.RData", compress='xz')
