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
  group_by(PET) %>%
  nest(.key = "tacs")

blooddat <- read_csv("data-raw/pbr28_blooddata.csv") %>%
  mutate(
    Cbl_dispcorr = ifelse(Cbl_dispcorr < 0, 0, Cbl_dispcorr),
    Cpl_metabcorr = ifelse(Cpl_metabcorr < 0, 0, Cpl_metabcorr)
  ) %>%
  group_by(PET) %>%
  nest(.key = "blooddata") %>%
  mutate(input = map(blooddata, ~blood_interp(
    t_blood = .x$Time / 60, blood = .x$Cbl_dispcorr,
    t_plasma = .x$Time / 60, plasma = .x$Cpl_metabcorr,
    t_parentfrac = 1, parentfrac = 1
  )))

petdat <- inner_join(tacdat, blooddat) %>%
  separate(PET, c("Subjname", "PETNo"), sep = "_", remove = F, convert = T)

pbr28 <- petdat %>%
  inner_join(demog) %>%
  arrange(PET)

save(pbr28, file="data/pbr28.RData", compress='xz')
