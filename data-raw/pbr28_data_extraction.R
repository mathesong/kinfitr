# Read in the data and save to a nested data.frame

library(tidyverse)
library(readr)

demog <- read_csv2('TrT_chemistry_demographics.csv') %>%
  select(Subjname=Akronym, Genotype, Comment,
         `MBq PET1`, `MBq PET2`, bodyMass=Weight_kg) %>%
  gather(PETNo, injRad, contains('MBq')) %>%
  mutate(PETNo = as.numeric(str_match(PETNo, '\\d'))) %>%
  mutate(PET = paste(Subjname, PETNo, sep='_'))


tacdat <- read_csv('tacdata.csv') %>%
  group_by(PET) %>%
  nest(.key = 'tacs')

blooddat <- read_csv('blooddata.csv') %>%
  mutate(Cbl.disp.corr = ifelse(Cbl.disp.corr < 0, 0, Cbl.disp.corr),
         Cpl..metabcorr. = ifelse(Cpl..metabcorr. < 0, 0, Cpl..metabcorr.)) %>%
  group_by(PET) %>%
  nest(.key='blooddata') %>%
  mutate(input = map(blooddata, ~blood_interp(
    t_blood = .x$ABSS.sec/60, blood=.x$Cbl.disp.corr,
    t_plasma=.x$ABSS.sec/60, plasma=.x$Cpl..metabcorr.,
    t_parentfrac = 1, parentfrac=1 ) ))

petdat <- inner_join(tacdat, blooddat) %>%
  separate(PET, c("Subjname", "PETNo"), sep='_', remove = F, convert=T)

pbr28 <- petdat %>%
  inner_join(demog) %>%
  arrange(PET)

