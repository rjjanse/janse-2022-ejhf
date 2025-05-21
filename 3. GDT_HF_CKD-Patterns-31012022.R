############################### General information ###################################
# Guideline directed therapy in HF and CKD
# Roemer Janse - 31/01/2022
# Code for adherence, persistence, and late initiation 
### -------------------------------------------------------------------------------- ###

rm(list = ls()); gc()

pacman::p_load("dplyr", "tidyverse")

setwd("C:/Users/rjjanse.lumcnet/onedrive - lumc/research/projects/7. gdt_hf_ckd/codes/dataframes")
load("sample.Rdata")
load("sample_arni.Rdata")
load("sample_rasi.Rdata")
load("drugs_new.Rdata")
load("drugs_arni_new.Rdata")
load("drugs_rasi_new.Rdata")
load("atc_dur.Rdata")
load("atc_dur_arni.Rdata")
load("atc_dur_rasi.Rdata")

##### 1. Adherence and persistence among initiators #####
### Function for adherence ###
adh <- function(df, drug_df, yr, df_atc_dur){
    # Prepare drugs data
    drugs_df <- df %>% mutate(yr = yr,
                              year_dt = as.Date(ifelse(yr == 0, censor_dt, (init_dt + (yr * 365.25))), origin = "1970-01-01"),
                              year_dt = as.Date(ifelse(year_dt > censor_dt, censor_dt, year_dt), origin = "1970-01-01")) %>% 
        left_join(drug_df, "lopnr") %>% filter(edatum >= init_dt & edatum <= year_dt) %>% dplyr::select(-(drug:year_dt)) %>% 
        arrange(lopnr, edatum) %>% left_join(df_atc_dur %>% dplyr::select(atc, dpp), "atc")
    
    # Prepare data for the loop
    adh_drug <- df %>% mutate(yr = yr,
                              year_dt = as.Date(ifelse(yr == 0, censor_dt, (init_dt + (yr * 365.25))), origin = "1970-01-01"),
                              year_dt = as.Date(ifelse(year_dt > censor_dt, censor_dt, year_dt), origin = "1970-01-01"),
                              totday = as.numeric(year_dt - init_dt)) %>% left_join(drugs_df, "lopnr") %>%
        mutate(days = antal * dpp, 
               end_dt = as.Date(edatum + days),
               # Days without pills
               missed_days = 0,
               # Out of pills date
               out_of_pills_dt = end_dt) %>% 
        arrange(lopnr) %>% group_by(lopnr) %>% mutate(rownr = row_number(),
                                                      nrow = max(rownr),
                                                      out_of_pills_dt = as.Date(ifelse(rownr == 1, out_of_pills_dt, NA), 
                                                                                origin = "1970-01-01")) %>% ungroup()
    
    # Loop to also take into account stockpiling
    for(i in 2:max(adh_drug$rownr)){
        adh_drug <- adh_drug %>% arrange(lopnr) %>% group_by(lopnr) %>% 
            mutate(lag_out_of_pills_dt = lag(out_of_pills_dt),
                   # If new edatum is later than out of pills date (previous observation), calculate missed days, otherwise 0
                   missed_days = ifelse(rownr == i & edatum > lag_out_of_pills_dt, as.numeric(edatum - lag_out_of_pills_dt), 
                                        ifelse(rownr == i & edatum <= lag_out_of_pills_dt, 0, missed_days)),
                   # If no days missed, out of pills date is lag out of pills date + days
                   # If days were missed, out of pills date starts anew from edatum
                   out_of_pills_dt = as.Date(ifelse(rownr == i & missed_days == 0, lag_out_of_pills_dt + days,
                                                    ifelse(rownr == i & missed_days != 0, edatum + days, out_of_pills_dt)), 
                                             origin = "1970-01-01")) %>%
            ungroup()
    }
    
    # Determine adherence                                                  
    adh_drug <- adh_drug %>% arrange(lopnr, desc(rownr)) %>% group_by(lopnr) %>% 
        mutate(out_of_pills_dt = as.Date(ifelse(out_of_pills_dt > year_dt, year_dt, out_of_pills_dt), origin = "1970-01-01"),
               md = sum(missed_days),
               md = md + as.numeric(year_dt - out_of_pills_dt)) %>% 
        slice(1L) %>% ungroup() %>%
        mutate(dc = totday - md,
               pdc = round(dc / totday * 100),
               # Patients who died on the initiation date
               pdc = ifelse(is.na(pdc) & init_dt == death_dt, 100, pdc),
               pdc_ind = ifelse(pdc >= 80, 0, 1)) %>% dplyr::select(lopnr, drug, pdc, pdc_ind)
    
    return(adh_drug)
}

### Function for persistence ###
per <- function(df, drug_df, yr, gp, df_atc_dur){
    # To get persistence (60 day gap)
    # First determine when a prescription ends
    # Then calculate the gap between each end date and subsequent start date
    drugs_df <- df %>% mutate(yr = yr,
                              year_dt = as.Date(ifelse(yr == 0, censor_dt, (init_dt + (yr * 365.25))), origin = "1970-01-01"),
                              year_dt = as.Date(ifelse(year_dt > censor_dt, censor_dt, year_dt), origin = "1970-01-01")) %>% 
        left_join(drug_df, "lopnr") %>% filter(edatum >= init_dt & edatum <= year_dt) %>% dplyr::select(-(drug:year_dt)) %>% 
        arrange(lopnr, edatum) %>% left_join(df_atc_dur %>% dplyr::select(atc, dpp), "atc")
    
    per_drug <- df %>% mutate(yr = yr,
                              year_dt = as.Date(ifelse(yr == 0, censor_dt, (init_dt + (yr * 365.25))), origin = "1970-01-01"),
                              year_dt = as.Date(ifelse(year_dt > censor_dt, censor_dt, year_dt), origin = "1970-01-01")) %>% 
        left_join(drugs_df, "lopnr") %>%
        mutate(days = antal * dpp,
               end_dt = as.Date(edatum + days)) %>%
        arrange(lopnr, edatum) %>% group_by(lopnr) %>% mutate(rownr = row_number(),
                                                              out_of_pills_dt = as.Date(ifelse(rownr == 1, end_dt, NA), 
                                                                                        origin = "1970-01-01"),
                                                              gap = NA) %>% ungroup()
    
    # Take into account stockpiling
    for(i in 1:max(per_drug$rownr)){
        if(i == 1){
            per_drug <- per_drug %>% arrange(lopnr) %>% group_by(lopnr) %>%
                mutate(lead_edatum = lead(edatum),
                       lead_edatum = as.Date(ifelse(is.na(lead_edatum), year_dt, lead_edatum), origin = "1970-01-01"),
                       # If new edatum is later than out of pills date, calculate gap
                       gap = ifelse(rownr == 1 & lead_edatum > out_of_pills_dt, as.numeric(lead_edatum - out_of_pills_dt), 
                                    ifelse(rownr == 1 & lead_edatum <= out_of_pills_dt, 0, gap))) %>% ungroup()
        }
        
        else{
            per_drug <- per_drug %>% arrange(lopnr) %>% group_by(lopnr) %>%
                mutate(lag_gap = lag(gap),
                       lag_out_of_pills_dt = lag(out_of_pills_dt),
                       # If there was a gap, no matter how long, out of pills date is reset to end date of current edatum, otherwise add
                       # days of current dispensation to out of pills date
                       out_of_pills_dt = as.Date(ifelse(rownr == i & lag_gap == 0, lag_out_of_pills_dt + days, 
                                                        ifelse(rownr == i & lag_gap != 0, end_dt, out_of_pills_dt)), origin = "1970-01-01"),
                       lead_edatum = lead(edatum),
                       lead_edatum = as.Date(ifelse(is.na(lead_edatum), year_dt, lead_edatum), origin = "1970-01-01"),
                       # If new edatum is later than out of pills date, calculate gap, otherwise 0
                       gap = ifelse(rownr == i & lead_edatum > out_of_pills_dt, as.numeric(lead_edatum - out_of_pills_dt), 
                                    ifelse(rownr == i & lead_edatum <= out_of_pills_dt, 0, gap))) %>% ungroup()
        }
    }
    
    per_drug <- per_drug %>%
        mutate(gap_ind = ifelse(gap >= gp, 1, 0),
               # Gap date for people with a gap is end date plus grace period, otherwise end date
               gap_dt = as.Date(ifelse(gap_ind == 1, (out_of_pills_dt + gp), year_dt), origin = "1970-01-01"),
               # Gap indicators for people of whom the grace period ends after the year date are 0 (no gap)
               gap_ind = ifelse(gap_dt > year_dt & gap_ind == 1, 0, gap_ind),
               # End date for end dates past year date is year date
               # This also sets end dates with grace period after one year for people who went from gap (1) to no gap (0) back to right date
               gap_dt = as.Date(ifelse(gap_dt > year_dt, year_dt, gap_dt), origin = "1970-01-01")) %>%
        arrange(lopnr, desc(gap_ind), gap_dt) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
        mutate(time2gap = ifelse(gap_ind == 1, as.numeric(gap_dt - init_dt), as.numeric(year_dt - init_dt))) %>%
        dplyr::select(lopnr, drug, gap_ind, time2gap, gap_dt)
    
    return(per_drug)
}

### Run functions
inits <- sample %>% filter(init_pd == 1)
inits_arni <- sample_arni %>% filter(init_pd == 1)
inits_rasi <- sample_rasi %>% filter(init_pd == 1)

adh_bb <- adh(inits %>% filter(drug == "BB"), drugs %>% filter(drug == "BB") %>% dplyr::select(-drug), yr = 1, atc_dur)
adh_rasi <- adh(inits_rasi %>% filter(drug == "RASi"), drugs_rasi %>% filter(drug == "RASi") %>% dplyr::select(-drug), yr = 1, atc_dur_rasi)
adh_arni <- adh(inits_arni %>% filter(drug == "ARNi"), drugs_arni %>% filter(drug == "ARNi") %>% dplyr::select(-drug), yr = 1, atc_dur_arni)
adh_mra <- adh(inits %>% filter(drug == "MRA"), drugs %>% filter(drug == "MRA") %>% dplyr::select(-drug), yr = 1, atc_dur)
adh_rasi_arni <- adh(inits %>% filter(drug == "RASi/ARNi"), drugs %>% filter(drug == "RASi/ARNi") %>% dplyr::select(-drug), yr = 1, atc_dur)
adherence <- rbind(adh_bb, adh_rasi, adh_arni, adh_mra, adh_rasi_arni)
save(adherence, file = "adherence.Rdata")
load("adherence.Rdata")

# Save individual adherence files
save(adh_bb, file = "adh_bb.Rdata")
save(adh_rasi, file = "adh_rasi.Rdata")
save(adh_arni, file = "adh_arni.Rdata")
save(adh_mra, file = "adh_mra.Rdata")
save(adh_rasi_arni, file = "adh_rasi_arni.Rdata")

per_bb <- per(inits %>% filter(drug == "BB"), drugs %>% filter(drug == "BB") %>% dplyr::select(-drug), yr = 1, gp = 60, atc_dur)
per_rasi <- per(inits_rasi %>% filter(drug == "RASi"), drugs_rasi %>% filter(drug == "RASi") %>% dplyr::select(-drug), yr = 1, gp = 60, atc_dur_rasi)
per_arni <- per(inits_arni %>% filter(drug == "ARNi"), drugs_arni %>% filter(drug == "ARNi") %>% dplyr::select(-drug), yr = 1, gp = 60, atc_dur_arni)
per_mra <- per(inits %>% filter(drug == "MRA"), drugs %>% filter(drug == "MRA") %>% dplyr::select(-drug), yr = 1, gp = 60, atc_dur)
per_rasi_arni <- per(inits %>% filter(drug == "RASi/ARNi"), drugs %>% filter(drug == "RASi/ARNi") %>% dplyr::select(-drug), yr = 1, gp = 60, atc_dur)
persistence <- rbind(per_bb, per_rasi, per_arni, per_mra, per_rasi_arni)
save(persistence, file = "persistence.Rdata")
load("persistence.Rdata")

# Save individual persistence files
save(per_bb, file = "per_bb.Rdata")
save(per_rasi, file = "per_rasi.Rdata")
save(per_arni, file = "per_arni.Rdata")
save(per_mra, file = "per_mra.Rdata")
save(per_rasi_arni, file = "per_rasi_arni.Rdata")

load("sample_full.Rdata")

### Add data together and determine adherence and persistence for triple therapy initiators
pop <- sample_full %>% left_join(adherence, c("lopnr", "drug")) %>% left_join(persistence, c("lopnr", "drug")) %>%
    arrange(lopnr, drug) %>% group_by(lopnr) %>% 
    mutate(# PDC for triple therapy is the lowest PDC observed
           pdc = ifelse(drug == "TT" & init_pd == 1, min(pdc, na.rm = TRUE), pdc),
           # Based on the new PDC determine PDC inidcator
           pdc_ind = ifelse(drug == "TT" & init_pd == 1 & pdc >= 80, 0, ifelse(drug == "TT" & init_pd == 1 & pdc < 80, 1, pdc_ind)),
           # Gap indicator is the maximum of the others; i.e., if one drug has a gap, TT has a gap
           gap_ind = ifelse(drug == "TT" & init_pd == 1, max(gap_ind, na.rm = TRUE), gap_ind),
           # Gap date is either one year after the initiation date in the case of no gap, or the first time a gap occurs
           gap_dt = as.Date(ifelse(drug == "TT" & gap_ind == 0, init_dt + 365.25, 
                                   ifelse(drug == "TT" & gap_ind == 1, min(gap_dt, na.rm = TRUE), gap_dt)), origin = "1970-01-01"),
           # Calculate time to gap
           time2gap = ifelse(drug == "TT" & init_pd == 1, as.numeric(gap_dt - init_dt), time2gap),
           # If the time to gap is negative, these people discontinued a drug before initiating the last, so TT was not actually initiated
           init_pd = ifelse(drug == "TT" & init_pd == 1 & time2gap < 0, 0, 
                            ifelse(drug == "TT" & init_pd == 1 & time2gap >= 0, 1,
                                   ifelse(drug == "TT" & init_pd == 0, init_pd, init_pd))),
           # Reset other variables in this case
           pdc = ifelse(drug == "TT" & init_pd == 0, NA, pdc),
           pdc_ind = ifelse(drug == "TT" & init_pd == 0, NA, pdc_ind),
           gap_ind = ifelse(drug == "TT" & init_pd == 0, NA, gap_ind),
           gap_dt = as.Date(ifelse(drug == "TT" & init_pd == 0, NA, gap_dt), origin = "1970-01-01"),
           time2gap = ifelse(drug == "TT" & init_pd == 0, NA, time2gap)) %>% ungroup()

### Determine if death occurred before end of observation
# Only matters if they had a gap (which might be due to death)
pop <- pop %>% mutate(mort_ind = ifelse(death_dt <= gap_dt & gap_ind == 1, 1, 0),
                      mort_ind = ifelse(!is.na(gap_ind) & is.na(mort_ind), 0, mort_ind),
                      time2event = ifelse(mort_ind == 1, as.numeric(death_dt - init_dt), time2gap))

### Among discontinuers, determine who restarts within 3 months
rest <- pop %>% filter(gap_ind == 1) %>% left_join(drugs, c("lopnr", "drug")) %>% filter(edatum > gap_dt & edatum <= (gap_dt + 90)) %>%
    arrange(lopnr, drug, edatum) %>% group_by(lopnr, drug) %>% slice(1L) %>% ungroup() %>% mutate(time2restart = as.numeric(edatum - gap_dt),
                                                                                                  restart_ind = 1) %>% 
    dplyr::select(lopnr, drug, time2restart, restart_ind)

# Restarts for triple therapy
rest <- pop %>% dplyr::select(lopnr, drug, gap_ind) %>% filter(gap_ind == 1) %>% left_join(rest, c("lopnr", "drug")) %>%
    arrange(lopnr, drug) %>% group_by(lopnr) %>% mutate(# Calculate amount of drugs that were discontinued (gap_ind == 1) minus TT
                                                        # This will give wrong results for non-TT lopnrs but we are not interested in those 
                                                        # now anyway
                                                        ndisc = max(row_number()) - 1,
                                                        # Amount of drugs restarted
                                                        nrest = sum(restart_ind, na.rm = TRUE),
                                                        # If drugs restarted is equal to drugs discontinued, triple therapy is also restarted
                                                        # if it was discontinued in the first place
                                                        restart_ind = ifelse(drug == "TT" & ndisc == nrest, 1,
                                                                             ifelse(drug == "TT" & ndisc != nrest, 0,
                                                                                    restart_ind)),
                                                        time2restart = ifelse(drug == "TT" & restart_ind == 1, 
                                                                              min(time2restart, na.rm = TRUE), time2restart)) %>%
    ungroup() %>% dplyr::select(lopnr, drug, restart_ind, time2restart) %>% replace_na(list(restart_ind = 0))

### Among non initiators, determine who starts within 6 months
rints <- pop %>% filter(init_pd == 0) %>% left_join(drugs, c("lopnr", "drug")) %>% filter(edatum > oom_dt & edatum <= oom_dt + 180) %>%
    dplyr::select(lopnr, drug, oom_dt, edatum) %>% arrange(lopnr, edatum) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
    mutate(time2claim = as.numeric(edatum - oom_dt),
           claim_ind = 1) %>% dplyr::select(lopnr, drug, time2claim, claim_ind)

# Rints for triple therapy
rints <- pop %>% dplyr::select(lopnr, drug, init_pd) %>% filter(init_pd == 0) %>% left_join(rints, c("lopnr", "drug")) %>% 
    arrange(lopnr) %>% group_by(lopnr) %>% mutate(# Determine amount of drugs not claimed
                                                  # -2: one for TT, one for RASi/ARNi (as only one of two should be claimed)
                                                  nninit = max(row_number()) - 2,
                                                  # Amount of drugs claimed later
                                                  nlclaim = sum(claim_ind, na.rm = TRUE),
                                                  # If drugs not claimed & drugs claimed later are equal, then TT is late claimed
                                                  claim_ind = ifelse(drug == "TT" & nninit == nlclaim, 1,
                                                                     ifelse(drug == "TT" & nninit != nlclaim, 0, claim_ind)),
                                                  time2claim = ifelse(drug == "TT" & claim_ind == 1, min(time2claim, na.rm = TRUE),
                                                                      time2claim),
                                                  # Some people initiate all three drugs but discontinue one before all three are initiated. 
                                                  # They thus do not initiate triple therapy as not all three drugs are used at the same time.
                                                  # When this code runs, these people get their TT row, but none of the other drugs (except
                                                  # maybe ARNi/RASi), but the code still recognizes them as late claiming triple therapy. 
                                                  # Here, we recognize them by infinite time2claim (as it does not exist) and set them to no 
                                                  # late claim.
                                                  claim_ind = ifelse(is.infinite(time2claim), 0, claim_ind),
                                                  time2claim = ifelse(claim_ind == 0, NA, time2claim)) %>%
    ungroup() %>% dplyr::select(lopnr, drug, claim_ind, time2claim) %>% replace_na(list(claim_ind = 0))

# Also set some missing triple therapy init_pds to 0
pop <- pop %>% left_join(rints, c("lopnr", "drug")) %>% left_join(rest, c("lopnr", "drug")) %>% replace_na(list(init_pd = 0))

save(pop, file = "pop_new.Rdata")

# ##### 2. Sensitivity analysis with 180d grace period #####
# ### Run functions
# inits <- sample %>% filter(init_pd == 1)
# 
# adh_bb <- adh(inits %>% filter(drug == "BB"), drugs %>% filter(drug == "BB") %>% dplyr::select(-drug), yr = 1)
# adh_rasi <- adh(inits %>% filter(drug == "RASi"), drugs %>% filter(drug == "RASi") %>% dplyr::select(-drug), yr = 1)
# adh_arni <- adh(inits %>% filter(drug == "ARNi"), drugs %>% filter(drug == "ARNi") %>% dplyr::select(-drug), yr = 1)
# adh_mra <- adh(inits %>% filter(drug == "MRA"), drugs %>% filter(drug == "MRA") %>% dplyr::select(-drug), yr = 1)
# adherence_sa <- rbind(adh_bb, adh_rasi, adh_arni, adh_mra)
# save(adherence_sa, file = "adherence_sa.Rdata")
# load("adherence_sa.Rdata")
# 
# per_bb <- per(inits %>% filter(drug == "BB"), drugs %>% filter(drug == "BB") %>% dplyr::select(-drug), yr = 1, gp = 180)
# per_rasi <- per(inits %>% filter(drug == "RASi"), drugs %>% filter(drug == "RASi") %>% dplyr::select(-drug), yr = 1, gp = 180)
# per_arni <- per(inits %>% filter(drug == "ARNi"), drugs %>% filter(drug == "ARNi") %>% dplyr::select(-drug), yr = 1, gp = 180)
# per_mra <- per(inits %>% filter(drug == "MRA"), drugs %>% filter(drug == "MRA") %>% dplyr::select(-drug), yr = 1, gp = 180)
# persistence_sa <- rbind(per_bb, per_rasi, per_arni, per_mra)
# save(persistence_sa, file = "persistence_sa.Rdata")
# load("persistence_sa.Rdata")
# 
# ### Add data together and determine persistence for triple therapy initiators
# pop <- sample %>% left_join(adherence_sa, c("lopnr", "drug")) %>% left_join(persistence_sa, c("lopnr", "drug")) %>%
#     arrange(lopnr, drug) %>% group_by(lopnr) %>% 
#     mutate(# PDC for triple therapy is the lowest PDC observed
#         pdc = ifelse(drug == "TT" & init_pd == 1, min(pdc, na.rm = TRUE), pdc),
#         # Based on the new PDC determine PDC inidcator
#         pdc_ind = ifelse(drug == "TT" & init_pd == 1 & pdc >= 80, 0, ifelse(drug == "TT" & init_pd == 1 & pdc < 80, 1, pdc_ind)),
#         # Gap indicator is the maximum of the others; i.e., if one drug has a gap, TT has a gap
#         gap_ind = ifelse(drug == "TT" & init_pd == 1, max(gap_ind, na.rm = TRUE), gap_ind),
#         # Gap date is either one year after the initiation date in the case of no gap, or the first time a gap occurs
#         gap_dt = as.Date(ifelse(drug == "TT" & gap_ind == 0, init_dt + 365.25, 
#                                 ifelse(drug == "TT" & gap_ind == 1, min(gap_dt, na.rm = TRUE), gap_dt)), origin = "1970-01-01"),
#         # Calculate time to gap
#         time2gap = ifelse(drug == "TT" & init_pd == 1, as.numeric(gap_dt - init_dt), time2gap),
#         # If the time to gap is negative, these people discontinued a drug before initiating the last, so TT was not actually initiated
#         init_pd = ifelse(drug == "TT" & init_pd == 1 & time2gap < 0, 0, 
#                          ifelse(drug == "TT" & init_pd == 1 & time2gap >= 0, 1,
#                                 ifelse(drug == "TT" & init_pd == 0, init_pd, init_pd))),
#         # Reset other variables in this case
#         pdc = ifelse(drug == "TT" & init_pd == 0, NA, pdc),
#         pdc_ind = ifelse(drug == "TT" & init_pd == 0, NA, pdc_ind),
#         gap_ind = ifelse(drug == "TT" & init_pd == 0, NA, gap_ind),
#         gap_dt = as.Date(ifelse(drug == "TT" & init_pd == 0, NA, gap_dt), origin = "1970-01-01"),
#         time2gap = ifelse(drug == "TT" & init_pd == 0, NA, time2gap)) %>% ungroup()
# 
# ### Determine if death occurred before end of observation
# pop <- pop %>% mutate(mort_ind = ifelse(death_dt <= gap_dt & gap_ind == 1, 1, 0),
#                       mort_ind = ifelse(!is.na(gap_ind) & is.na(mort_ind), 0, mort_ind),
#                       time2event = ifelse(mort_ind == 1, as.numeric(death_dt - init_dt), time2gap))
# 
# ### Among discontinuers, determine who restarts within 3 months
# rest <- pop %>% filter(gap_ind == 1) %>% left_join(drugs, c("lopnr", "drug")) %>% filter(edatum > gap_dt & edatum <= (gap_dt + 90)) %>%
#     arrange(lopnr, edatum) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>% mutate(time2restart = as.numeric(edatum - gap_dt),
#                                                                                       restart_ind = 1) %>% 
#     dplyr::select(lopnr, drug, time2restart, restart_ind)
# 
# # Restarts for triple therapy
# rest <- pop %>% dplyr::select(lopnr, drug, gap_ind) %>% filter(gap_ind == 1) %>% left_join(rest, c("lopnr", "drug")) %>%
#     arrange(lopnr, drug) %>% group_by(lopnr) %>% mutate(# Calculate amount of drugs that were discontinued (gap_ind == 1) minus TT
#         # This will give wrong results for non-TT lopnrs but we are not interested in those 
#         # now anyway
#         ndisc = max(row_number()) - 1,
#         # Amount of drugs restarted
#         nrest = sum(restart_ind, na.rm = TRUE),
#         # If drugs restarted is equal to drugs discontinued, triple therapy is also restarted
#         # if it was discontinued in the first place
#         restart_ind = ifelse(drug == "TT" & ndisc == nrest, 1,
#                              ifelse(drug == "TT" & ndisc != nrest, 0,
#                                     restart_ind)),
#         time2restart = ifelse(drug == "TT" & restart_ind == 1, 
#                               min(time2restart, na.rm = TRUE), time2restart)) %>%
#     ungroup() %>% dplyr::select(lopnr, drug, restart_ind, time2restart) %>% replace_na(list(restart_ind = 0))
# 
# # Also set some triple therapy inits to 0
# pop_sa180gp <- pop %>% left_join(rest, c("lopnr", "drug")) %>% replace_na(list(init_pd = 0))
# 
# save(pop_sa180gp, file = "pop_sa180gp.Rdata")
# 
# rm(list = ls())