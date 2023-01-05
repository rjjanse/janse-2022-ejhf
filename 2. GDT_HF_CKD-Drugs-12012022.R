############################### General information ###################################
# Guideline directed therapy in HF and CKD
# Roemer Janse - 12/01/2022
# Code for drugs (determining initiators and non-initiators)
### -------------------------------------------------------------------------------- ###

rm(list = ls())

pacman::p_load("dplyr", "tidyverse")

memory.limit(size = 60000)
setwd("~/Research/[] Nephrology/5. GDT_HF_CKD/Codes/Dataframes/")
load("cohort.Rdata")
load("~/Datasets/SwedeHF2019/lmsel.Rdata")

##### 1. Prepare medication data #####
# load("~/Datasets/SwedeHF2019/lmsel.RData")
# drugs <- lmsel %>% rename(lopnr = LopNr) %>% filter(lopnr %in% cohort$lopnr) %>% rename(edatum = EDATUM, antal = ANTAL, atc = ATC) %>%
#     dplyr::select(lopnr, atc, edatum, forpddd, antnum, antal) %>%
#     filter(grepl("^C09[A-D]", atc) | # RASi
#                grepl("^C03DA", atc) | # MRA
#                grepl("^C07", atc) | # BB
#                grepl("^C09DX04", atc)) %>% # ARNI
#     mutate(drug = ifelse(grepl("^C09[A-D]", atc), "RASi",
#                          ifelse(grepl("^C03DA", atc), "MRA",
#                                 ifelse(grepl("^C07", atc), "BB", NA))),
#            drug = ifelse(grepl("^C09DX04", atc), "ARNi", drug)) %>%
#     # Correct wrong prescriptions (code from Faizan)
#     mutate(tot_forddd = forpddd * antal,
#            tot_antnum = antnum * antal) %>% arrange(lopnr, atc, edatum) %>% group_by(lopnr, atc, edatum) %>%
#     mutate(purchase_ddd = sum(tot_forddd),
#            purchase_antal = sum(antal),
#            purchase_antnum = sum(tot_antnum)) %>% # value = zero means that the whole purchase on that date has been cancelled
#     filter(purchase_antal > 0 & purchase_antnum > 0) %>% # here we would normally choose purchase_antnum > 0, but for this study we are only
#     # interested in the amount of packages, so we only keep those with at least a package left with at least more than 0 pills
#     slice(1L) %>% ungroup() %>%  # keep only one indication of purchase per date
#     mutate(antal = purchase_antal, antnum = (purchase_antnum / antal)) %>%
#     dplyr::select(-tot_forddd, -tot_antnum, -purchase_ddd, -purchase_antal, -purchase_antnum) %>%
#     # Remove dispensations after death
#     left_join(cohort %>% dplyr::select(lopnr, death_dt), "lopnr") %>%
#     filter(edatum <= death_dt) %>% dplyr::select(-death_dt)
# 
# drugs.rasi_arni <- drugs %>% filter(drug == "RASi" | drug == "ARNi") %>% mutate(drug = "RASi/ARNi")
# 
# drugs <- rbind(drugs, drugs.rasi_arni)
# 
# save(drugs, file = "drugs.Rdata")

load("drugs.Rdata")

##### 2. Determine initiators and non-initiators #####
### Determine how long a dispensation lasts
# Keep first 3 dispensations per person per ATC
# Calculate total amount of days and packages over these three dispensations
# Calculate how long it takes to finish a package (days per package)
# Calculate average amount of days over all people per ATC
# atc_dur <- dplyr::select(cohort, lopnr, index_dt) %>% left_join(drugs, "lopnr") %>% filter(edatum >= index_dt) %>%
#     arrange(lopnr, atc, edatum) %>% group_by(lopnr, atc) %>% mutate(n = row_number()) %>% filter(n >= 1 & n <= 3) %>%
#     mutate(days = as.numeric(lead(edatum) - edatum),
#            tdays = sum(days, na.rm = TRUE),
#            npack = sum(antal),
#            dpp = tdays / npack,
#            dpp = ifelse(dpp == 0, NA, dpp)) %>% slice(1L) %>% ungroup() %>%
#     arrange(atc) %>% group_by(atc) %>% mutate(dpp_avg = mean(dpp, na.rm = TRUE)) %>% slice(1L) %>% ungroup() %>%
#     mutate(dpp = ifelse(dpp_avg < 45, 30,
#                         ifelse(dpp_avg >= 45 & dpp_avg < 75, 60,
#                                ifelse(dpp_avg >= 75, 90, NA)))) %>%
#     dplyr::select(atc, drug, dpp, dpp_avg)
# 
# save(atc_dur, file = "atc_dur.Rdata")
# 
# ### Add a row per drug to cohort
# lopnrs <- unique(cohort$lopnr)
# 
# for(i in lopnrs){
#     new <- cbind(rep(i, 6), c("RASi", "BB", "MRA", "ARNi", "TT", "RASi/ARNi"))
# 
#     if(i == lopnrs[[1]]){
#         dat <- new
#     }
# 
#     else{
#         dat <- rbind(dat, new)
#     }
# }
# 
# dat <- as.data.frame(dat) %>% rename(lopnr = V1, drug = V2) %>% mutate(lopnr = as.numeric(lopnr))
# save(dat, file = "dat.Rdata")

load("dat.Rdata")
load("atc_dur.Rdata")

# Join data
sample <- dat %>% left_join(cohort, "lopnr") 

### Add drugs to sample and determine end date of each dispensation with dpp_avg
## First only keep last dispensation before index and all dispensations after
# Last dispensation before index date
# Keep end date of the last dispensation per drug (without grace period)
ldisp <- dplyr::select(sample, lopnr, adm_dt) %>% left_join(drugs, "lopnr") %>% filter(edatum < adm_dt) %>%
    arrange(lopnr, drug, desc(edatum)) %>% group_by(lopnr, drug) %>% slice(1L) %>% ungroup() %>% left_join(atc_dur, c("atc", "drug")) %>%
    mutate(end_dt_ldisp = as.Date(edatum + (antal * dpp), origin = "1970-01-01")) %>% dplyr::select(lopnr, drug, end_dt_ldisp)

# All dispensations after index
ndisp <- dplyr::select(sample, lopnr, index_dt) %>% left_join(drugs, "lopnr") %>% filter(edatum >= index_dt) %>% dplyr::select(-index_dt)

# Add drugs to sample
# Determine whether a drug is still used after index date (out-of-medication date: oom_dt)
sample <- sample %>% left_join(ldisp, c("lopnr", "drug")) %>% left_join(ndisp, c("lopnr", "drug")) %>%
    mutate(oom_dt = as.Date(ifelse(end_dt_ldisp <= index_dt | is.na(end_dt_ldisp), index_dt, end_dt_ldisp), origin = "1970-01-01"),
           edatum = as.Date(ifelse(edatum < oom_dt, NA, edatum), origin = "1970-01-01"))

# Placeholder for number of drugs
ph <- sample %>% arrange(lopnr, drug) %>% group_by(lopnr, drug) %>% slice(1L) %>% ungroup()
table(ph$drug)

# Exclude oom_dt after censor_dt
sample <- sample %>% filter(!(oom_dt > censor_dt))
ph <- sample %>% arrange(lopnr, drug) %>% group_by(lopnr, drug) %>% slice(1L) %>% ungroup()
table(ph$drug)
length(unique(sample$lopnr))

ph <- sample %>% arrange(lopnr) %>% group_by(lopnr) %>% slice(1L) %>% ungroup()
table(ph$cov_ef)
# n = 43,738
# rEF = 22,631
# mEF = 21,107

# Determine whether a dispensation is collected within 90 days after the out-of-medication date (post-discharge: init_pd)
sample <- sample %>%
    arrange(lopnr, drug, edatum) %>% group_by(lopnr, drug) %>% slice(1L) %>% ungroup() %>%
    mutate(prev = ifelse(oom_dt == index_dt, 0, 1),
           init_pd = ifelse(prev == 0 & edatum <= (oom_dt + 90), 1,
                            ifelse(prev == 1 & edatum <= (oom_dt + 60), 1, 0)),
           init_pd = ifelse(is.na(init_pd) & drug != "TT", 0, init_pd),
           init_dt = as.Date(ifelse(init_pd == 1, edatum, NA), origin = "1970-01-01"))

tts <- sample %>% filter(drug %in% c("RASi/ARNi", "MRA", "BB", "TT")) %>% arrange(lopnr, drug) %>% group_by(lopnr) %>%
    mutate(concdrug = ifelse(drug == "TT", sum(init_pd, na.rm = TRUE), NA), 
           init_pd = ifelse(drug == "TT" & concdrug >= 3, 1,
                            ifelse(drug == "TT" & concdrug < 3, 0, init_pd))) %>% filter(drug == "TT" & init_pd == 1)

sample <- sample %>% arrange(lopnr, drug) %>% group_by(lopnr) %>% 
    mutate(init_pd = ifelse(drug == "TT" & lopnr %in% tts[["lopnr"]], 1, 
                            ifelse(drug == "TT" & !(lopnr %in% tts[["lopnr"]]), 0, init_pd)),
           init_dt = as.Date(ifelse(drug == "TT" & init_pd == 1, max(init_dt, na.rm = TRUE), init_dt), origin = "1970-01-01")) %>%
    ungroup() %>% dplyr::select(-end_dt_ldisp)

# Change medication prescription from wide to long
sample <- sample %>% mutate(# Prescription for RASi
    pres_pd = ifelse(drug == "RASi" & med_rasi == 1, 1, NA),
    # Prescription for BBs
    pres_pd = ifelse(drug == "BB" & med_bb == 1, 1, pres_pd),
    # Prescription for MRAs
    pres_pd = ifelse(drug == "MRA" & med_mra == 1, 1, pres_pd),
    # Prescription for ARNis
    pres_pd = ifelse(drug == "ARNi" & med_arni == 1, 1, pres_pd),
    # Prescription for TT
    pres_pd = ifelse(drug == "TT" & med_tt == 1, 1, pres_pd),
    # Prescription for RASi/ARNi
    pres_pd = ifelse(drug == "RASi/ARNi" & (med_rasi == 1 | med_arni == 1), 1, pres_pd)) %>%
    tidyr::replace_na(list(pres_pd = 0)) %>% dplyr::select(-(med_bb:med_tt), -(atc:antal))

# Determine outcomes
sample <- sample %>% mutate(time2mort = as.numeric(death_dt - oom_dt),
                            # For time to heart failure hospitalisation, extract days from index to oom_dt from time2hosphf
                            time2hosphf = time2hosphf - as.numeric(oom_dt - index_dt),
                            time2cvmort = as.numeric(death_dt - oom_dt),
                            # Composite outcome is cv death or heart failure hospitalisation
                            oc_comp = ifelse(oc_hosphf == "Yes" | oc_cvmort == "Yes", 1, 0),
                            time2comp = ifelse(oc_comp == 1, pmin(time2hosphf, time2cvmort, na.rm = TRUE),
                                               as.numeric(censor_dt - oom_dt)))

save(sample, file = "sample_interim.Rdata")
load("sample_interim.Rdata")

### Only C03DA01 and C03DA04 in lakemedelsregistret
# ## Check diuretics: we have a variable ldiu and diu, respectively loop diuretics and diuretics.
# # C03C: loop diuretics
# # Other C03: other diuretics
# check <- sample %>% dplyr::select(lopnr, index_dt, cov_diu, cov_ldiu) %>% arrange(lopnr) %>% group_by(lopnr) %>% slice(1L) %>% ungroup()
# 
# check2 <- check %>% left_join(lmsel %>% rename(lopnr = 1, edatum = 2, atc = 8) %>% dplyr::select(lopnr, edatum, atc) %>%
#                                   filter(grepl("^C03", atc)), "lopnr") %>%
#     mutate(atc = ifelse(edatum < (index_dt - 365.25) | edatum > (index_dt + 30), NA, atc),
#            edatum = as.Date(ifelse(is.na(atc), NA, edatum), origin = "1970-01-01"),
#            lmed_ldiu = ifelse(grepl("^C03C", atc), "Yes", "No"), # Loop diuretics in LMED
#            lmed_diu = ifelse(!is.na(atc), "Yes", "No"), # Any diuretic in LMED
#            lmed_odiu = ifelse(!grepl("^C03C", atc) & !is.na(atc), "Yes", "No")) # Other diuretics in LMED
# 
# table(check2$cov_ldiu, check2$lmed_ldiu, useNA = "always")
# table(check2$cov_diu, check2$lmed_diu, useNA = "always")
# table(check2$cov_diu, check2$lmed_odiu, useNA = "always")

save(sample, file = "sample.Rdata")

##### 3. Sensitivity analysis with 180d grace period #####
# load("dat.Rdata")
# load("atc_dur.Rdata")
# 
# # Join data
# sample <- dat %>% left_join(cohort, "lopnr") 
# 
# ### Add drugs to sample and determine end date of each dispensation with dpp_avg
# ## First only keep last dispensation before index and all dispensations after
# # Last dispensation before index date
# # Keep end date of the last dispensation per drug (without grace period)
# ldisp <- dplyr::select(sample, lopnr, adm_dt) %>% left_join(drugs, "lopnr") %>% filter(edatum < adm_dt) %>%
#     arrange(lopnr, drug, desc(edatum)) %>% group_by(lopnr, drug) %>% slice(1L) %>% ungroup() %>% left_join(atc_dur, c("atc", "drug")) %>%
#     mutate(end_dt_ldisp = as.Date(edatum + (antal * dpp), origin = "1970-01-01")) %>% dplyr::select(lopnr, drug, end_dt_ldisp)
# 
# # All dispensations after index
# ndisp <- dplyr::select(sample, lopnr, index_dt) %>% left_join(drugs, "lopnr") %>% filter(edatum >= index_dt) %>% dplyr::select(-index_dt)
# 
# # Add drugs to sample
# # Determine whether a drug is still used after index date (out-of-medication date: oom_dt)
# sample <- sample %>% left_join(ldisp, c("lopnr", "drug")) %>% left_join(ndisp, c("lopnr", "drug")) %>%
#     mutate(oom_dt = as.Date(ifelse(end_dt_ldisp <= index_dt | is.na(end_dt_ldisp), index_dt, end_dt_ldisp), origin = "1970-01-01"),
#            edatum = as.Date(ifelse(edatum < oom_dt, NA, edatum), origin = "1970-01-01"))
# 
# # Exclude death within 90 days (this is also censoring within 90 days)
# sample <- sample %>% filter(!(death_dt < (oom_dt + 90)))
# 
# length(unique(sample$lopnr))
# # n = 71,096
# 
# # Exclude oom_dt after censor_dt
# sample <- sample %>% filter(!(oom_dt > censor_dt))
# 
# length(unique(sample$lopnr))
# # n = 71,096
# 
# # Determine whether a dispensation is collected within 90 days after the out-of-medication date (post-discharge: init_pd)
# sample <- sample %>%
#     arrange(lopnr, drug, edatum) %>% group_by(lopnr, drug) %>% slice(1L) %>% ungroup() %>%
#     mutate(prev = ifelse(oom_dt == index_dt, 0, 1),
#            init_pd = ifelse(prev == 0 & edatum <= (oom_dt + 90), 1,
#                             ifelse(prev == 1 & edatum <= (oom_dt + 180), 1, 0)),
#            init_pd = ifelse(is.na(init_pd) & drug != "TT", 0, init_pd),
#            init_dt = as.Date(ifelse(init_pd == 1, edatum, NA), origin = "1970-01-01")) %>%
#     # Determine TT with init_pd for arni/rasi, bb & mra by taking the sum of init_pd. If this is 3 or 4, all drugs have been given
#     arrange(lopnr, drug) %>% group_by(lopnr) %>% mutate(concdrug = ifelse(drug == "TT", sum(init_pd, na.rm = TRUE), NA),
#                                                         init_pd = ifelse(drug == "TT" & concdrug >= 3, 1,
#                                                                          ifelse(drug == "TT" & concdrug < 3, 0, init_pd)),
#                                                         init_dt = as.Date(ifelse(drug == "TT" & init_pd == 1, max(init_dt, na.rm = TRUE), 
#                                                                                  init_dt), origin = "1970-01-01")) %>% ungroup() %>%
#     dplyr::select(-concdrug, -end_dt_ldisp)
# 
# # Change medication prescription from wide to long
# sample <- sample %>% mutate(# Prescription for RASi
#     pres_pd = ifelse(drug == "RASi" & med_rasi == 1, 1, NA),
#     # Prescription for BBs
#     pres_pd = ifelse(drug == "BB" & med_bb == 1, 1, pres_pd),
#     # Prescription for MRAs
#     pres_pd = ifelse(drug == "MRA" & med_mra == 1, 1, pres_pd),
#     # Prescription for ARNis
#     pres_pd = ifelse(drug == "ARNi" & med_arni == 1, 1, pres_pd),
#     # Prescription for TT
#     pres_pd = ifelse(drug == "TT" & med_tt == 1, 1, pres_pd)) %>%
#     tidyr::replace_na(list(pres_pd = 0)) %>% dplyr::select(-(med_bb:med_tt), -(atc:antal))
# 
# # Determine outcomes
# sample <- sample %>% mutate(time2mort = as.numeric(death_dt - oom_dt),
#                             # For time to heart failure hospitalisation, extract days from index to oom_dt from time2hosphf
#                             time2hosphf = time2hosphf - as.numeric(oom_dt - index_dt),
#                             time2cvmort = as.numeric(death_dt - oom_dt),
#                             # Composite outcome is cv death or heart failure hospitalisation
#                             oc_comp = ifelse(oc_hosphf == "Yes" | oc_cvmort == "Yes", 1, 0),
#                             time2comp = ifelse(oc_comp == 1, pmin(time2hosphf, time2cvmort, na.rm = TRUE),
#                                                as.numeric(censor_dt - oom_dt)))
# 
# # Exclude hospitalization for heart failure within 90 days 
# sample_sa180gp <- sample %>% filter(!(time2hosphf < 90))
# 
# length(unique(sample_sa180gp$lopnr))
# # n = 64,135
# 
# save(sample_sa180gp, file = "sample_sa180gp.Rdata")
# 
# rm(list = ls())
