library(dplyr)
library(tidyverse)

# Aplha/ Omega

AO <- readRDS("data_new/processed_data/Alpha_Omega_splicing_data.RDS")

AO_splice <- AO %>%
  drop_na()


AO_pre <- readRDS("data_new/Pre_Exercise/Alpha_Omega_PreExc_splicing_data.RDS") %>%
  drop_na()
# Conteratrain

Ct <- readRDS("data_new/processed_data/contratrain_splicing_data.RDS")

Ct_splice <- Ct %>%
  drop_na()

ct_pre <- readRDS("data_new/Pre_Exercise/ct_PreExc_splicing_data.RDS") %>%
  drop_na()


# COPD
COPD <- readRDS("data_new/processed_data/copd_splicing_data.RDS")

COPD_splice <- COPD %>%
  drop_na()


COPD_pre <- readRDS("data_new/Pre_Exercise/copd_preExc_splicing_data.RDS") %>%
  drop_na()

# Relief

Relief <- readRDS("data_new/processed_data/Relief_splicing_data.RDS")
Relief_splice <- Relief %>%
  drop_na()

ReliefPre <- readRDS("data_new/Pre_Exercise/Relief_PreExc_splicing_data.RDS") %>%
  drop_na()




# SRP102542
SRP102542 <- readRDS("data_new/processed_data/SRP102542_splicing_data.RDS")

SRP102542_splice <- SRP102542 %>%
  drop_na()

SRP102542_pre <- readRDS("data_new/Pre_Exercise/SRP102542_preExc_splicing_data.RDS") %>%
  drop_na()


# LOAD vOLUME DATA
Volume <- readRDS("data_new/processed_data/volume_splicing_data.RDS")

Volume_spliceq <- Volume %>%
  drop_na()

Volume_PreExc_spliceq <- readRDS("data_new/Pre_Exercise/vol_preExc_splicing_data.RDS") %>%
  drop_na()





