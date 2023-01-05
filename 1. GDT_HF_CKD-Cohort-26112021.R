############################### General information ###################################
# Guideline directed therapy in HF and CKD
# Roemer Janse - 25/11/2021
# Code for cohort derivation
### -------------------------------------------------------------------------------- ###

rm(list = ls())

pacman::p_load("dplyr", "tidyverse")

memory.limit(size = 60000)
setwd("~/Research/[] Nephrology/5. GDT_HF_CKD/Codes/Dataframes/")

##### 1. Exclusion criteria #####
# 1. Last registration with available eGFR & HF
# 2. Missing eGFR
# 3. Missing EF
# 4. eGFR >= 90
# 5. Maintenance dialysis
load("rs.data.Rdata")

length(unique(rs.data$LopNr))
# n = 90,383

# Before we keep only last registration, get death date from all information in SwedeHF.
death <- rs.data %>% mutate(death_dt = shf_indexdtm + sos_outtime_death) %>% arrange(LopNr, desc(sos_out_death), death_dt) %>% 
    group_by(LopNr) %>% slice(1L) %>% ungroup() %>% dplyr::select(LopNr, sos_out_death, death_dt, sos_out_deathcv, sos_deathcause) 

# Index date from 2009-01-01
rs.data <- rs.data %>% filter(shf_indexdtm >= as.Date("2009-01-01"))

length(unique(rs.data$LopNr))
# n = 71,889

# Missing eGFR
rs.data <- rs.data %>% filter(!is.na(shf_egfr))

length(unique(rs.data$LopNr))
# n = 71,042

# Missing EF
rs.data <- rs.data %>% filter(!is.na(EF))

length(unique(rs.data$LopNr))
# n = 63,782

# Death within 90 days
rs.data <- rs.data %>% dplyr::select(-(sos_out_death:sos_deathcause)) %>% left_join(death, "LopNr") %>% filter(death_dt > shf_indexdtm + 90)

length(unique(rs.data$LopNr))
# n = 60,397

# HFpEF
rs.data <- rs.data %>% filter(EF != "pEF")

length(unique(rs.data$LopNr))
ph <- rs.data %>% arrange(LopNr, desc(shf_indexdtm)) %>% group_by(LopNr) %>% slice(1L) %>% ungroup()
table(ph$EF)
# n = 47,272
# rEF = 32,004
# mEF = 15,268

# Last registration
rs.data <- rs.data %>% arrange(LopNr, desc(shf_indexdtm)) %>% group_by(LopNr) %>% slice(1L) %>% ungroup()

length(unique(rs.data$LopNr))
table(rs.data$EF)
# n = 47,272
# rEF = 32,004
# mEF = 15,268

# Residence care home
rs.data <- rs.data %>% filter(!(shf_residence == "Other home") | is.na(shf_residence))

length(unique(rs.data$LopNr))
table(rs.data$EF)
# n = 46,782
# rEF = 31,668
# mEF = 15,114

##### 2. eGFR classification #####
rs.data <- rs.data %>% mutate(cat_ef = EF,
                              cat_egfr = ifelse(shf_egfr >= 60, "eGFR >60",
                                                ifelse(shf_egfr >= 45 & shf_egfr < 60, "eGFR 45-59",
                                                       ifelse(shf_egfr >= 30 & shf_egfr < 45, "eGFR 30-44",
                                                              ifelse(shf_egfr >= 15 & shf_egfr < 30, "eGFR 15-29",
                                                                     ifelse(shf_egfr < 15, "eGFR <15", NA))))),
                              exp_egfr = ifelse(shf_egfr >= 60, "eGFR >60",
                                                ifelse(shf_egfr >= 45 & shf_egfr < 60, "eGFR 45-59",
                                                       ifelse(shf_egfr >= 30 & shf_egfr < 45, "eGFR 30-44",
                                                              ifelse(shf_egfr < 30, "eGFR <30", NA))))) %>% rename(lopnr = LopNr) 

table(rs.data$cat_ef, rs.data$exp_egfr, useNA = "ifany") 

cohort <- as_tibble(rs.data) 

##### 3. Create cohort and keep relevant variables #####
### Determine admission and censor (death or censor_dtm) date
# Index_dt is admission date for outpatient clinics and discharge date for hospitalizations (then shf_indexhosptime is not NA)
cohort <- cohort %>% rename(index_dt = shf_indexdtm) %>% 
    mutate(hosp = ifelse(is.na(shf_indexhosptime), 0, 1),
           adm_dt = as.Date(ifelse(hosp == 1, index_dt - shf_indexhosptime, index_dt), origin = "1970-01-01"),
           censor_dt = pmin(as.Date(ifelse(sos_out_death == "Yes", death_dt, censdtm), origin = "1970-01-01"),
                            as.Date("2019-12-31")))

### Rename variables
cohort <- cohort %>% rename(
    cov_age = shf_age,
    cov_loc = shf_location,
    cov_smok = shf_smoking,
    cov_cs = scb_famtype,
    cov_hfdur = shf_durationhf,
    cov_pe = shf_primaryetiology,
    cov_bpsys = shf_bpsys,
    cov_bpdia = shf_bpdia,
    cov_map = shf_map,
    cov_hb = shf_hb,
    cov_creat = shf_creatinine,
    cov_ntprobnp = shf_ntprobnp,
    med_arni = shf_arni,
    med_bb = shf_bbl,
    med_mra = shf_mra,
    cov_ld = shf_lungdisease,
    cov_edu = scb_education,
    cov_phhf = sos_prevhosphf,
    cov_egfr = shf_egfr,
    cov_an = shf_anemia, # Anemia
    cov_af = AF, # Atrial fibrillation
    cov_dm = DM, # Diabetes mellitus
    cov_revasc = revasc, # Revascularization
    cov_ihd = ischaemic, # Ischemic heart disease
    cov_ht = HT, # Hypertension
    med_rasi = ace_arb, 
    cov_dev = shf_device, 
    cov_valv = valvular_disease, # Valvular disease
    cov_cevd = stroke_tia, # Cerebrovascular disease
    cov_copd = sos_com_copd, # COPD
    cov_pad = sos_com_peripheralartery, # Peripheral artery disease
    cov_inc = income,
    cov_ca = sos_com_cancer3y, # Cancer in previous 3 years
    cov_hr = shf_heartrate,
    cov_bnp = bnp_cat,
    cov_bmi = shf_bmi,
    cov_age_cat = age.10,
    cov_nyha = nyha_cat,
    cov_ef = cat_ef,
    cov_egfr_cat = cat_egfr,
    cov_dx = shf_digoxin,
    cov_diu = shf_diuretic,
    cov_ldiu = shf_loopdiuretic,
    cov_st = shf_statin,
    cov_ac = shf_anticoagulantia,
    cov_ap = shf_asaantiplatelet,
    cov_nt = shf_nitrate,
    cov_pot = shf_potassium,
    cov_ob = bmi_cat,
    oc_mort = sos_out_death,
    oc_cvmort = sos_out_deathcv,
    oc_hosphf = sos_out_hosphf) %>%
    mutate(time2hosphf = sos_outtime_hosphf,
           cov_female = ifelse(shf_sex == "Female", 1, 0),
           cov_age_cat = ifelse(cov_age < 45, "<45 years",
                                ifelse(cov_age >= 45 & cov_age < 66, "45-65 years",
                                       ifelse(cov_age >= 66 & cov_age < 76, "66-75 years",
                                              ifelse(cov_age >= 76, ">75 years", NA)))),
           cov_dev = ifelse(cov_dev == "CRT" | cov_dev == "CRT & ICD" | cov_dev == "ICD" | cov_dev == "Pacemaker", 1,
                            ifelse(cov_dev == "No", 0, NA)),
           med_bb = ifelse(med_bb == "Yes", 1, 0),
           med_bb = ifelse(is.na(med_bb), 0, med_bb),
           med_rasi = ifelse(med_rasi == "Yes", 1, 0),
           med_rasi = ifelse(is.na(med_rasi), 0, med_rasi),
           med_mra = ifelse(med_mra == "Yes", 1, 0),
           med_mra = ifelse(is.na(med_mra), 0, med_mra),
           med_arni = ifelse(med_arni == "Yes", 1, 0),
           med_arni = ifelse(is.na(med_arni), 0, med_arni),
           med_tt = ifelse(med_bb == 1 & med_rasi == 1 & med_mra == 1, 1, 0),
           cov_cs = ifelse(cov_cs == "Cohabitating", "Cohabitating",
                           ifelse(cov_cs == "Living alone", "Living alone", NA)), 
           mod_cs = ifelse(is.na(cov_cs), "Missing", cov_cs),
           mod_cs = factor(mod_cs, levels = c("Living alone", "Cohabitating", "Missing")),
           cov_edu = ifelse(cov_edu == "Compulsory school", "Compulsory school",
                            ifelse(cov_edu == "Secondary school", "Secondary school",
                                   ifelse(cov_edu == "University", "University", NA))),
           mod_edu = ifelse(is.na(cov_edu), "Missing", cov_edu),
           mod_edu = factor(mod_edu, levels = c("Compulsory school", "Secondary school", "University", "Missing")),
           cov_smok = ifelse(cov_smok == "Current", "Current",
                             ifelse(cov_smok == "Former", "Former",
                                    ifelse(cov_smok == "Never", "Never", NA))),
           mod_smok = ifelse(is.na(cov_smok), "Missing", cov_smok),
           mod_smok = factor(mod_smok, levels = c("Never", "Former", "Current", "Missing"))) %>%
    dplyr::select(# Basic information & prescription
                  lopnr, adm_dt, index_dt, exp_egfr, censor_dt, med_bb, med_rasi, med_mra, med_arni, med_tt, 
                  # General information
                  cov_age, cov_age_cat, cov_female, cov_egfr, cov_egfr_cat, cov_loc,
                  # Heart failure information
                  cov_nyha, cov_ef, cov_hfdur, cov_pe, cov_phhf, cov_ntprobnp,
                  # Other comorbidities
                  cov_ld, cov_ob, cov_an, cov_af, cov_dm, cov_ihd, cov_ht, cov_dcm, cov_valv, cov_cevd, cov_copd, cov_pad, cov_ca,
                  # Procedures
                  cov_revasc, cov_dev,
                  # Medication
                  cov_dx, cov_ldiu, cov_diu, cov_st, cov_ac, cov_ap, cov_nt,
                  # Measurements
                  cov_hr, cov_bnp, cov_bpsys, cov_bpdia, cov_map, cov_hb, cov_bmi, cov_pot, 
                  # Other information
                  cov_smok, cov_edu, cov_cs, cov_inc,
                  # Outcomes
                  oc_mort, oc_cvmort, oc_hosphf, time2hosphf, death_dt,
                  # Model
                  mod_cs, mod_edu, mod_smok)

## Refactor variables
# Small function for binary variables
fac <- function(var){
    ph <- ifelse(var == "No", "No",
                 ifelse(var == "Yes", "Yes", NA))
    ph <- ifelse(is.na(ph), "Missing", ph)
    ph <- factor(ph, levels = c("No", "Yes", "Missing"))
    
    return(ph)
}

cohort <- cohort %>% mutate(cov_year = format(cohort$index_dt, "%Y"),
                            cov_year_cat = ifelse(index_dt <= as.Date("2013-12-31"), "2009-2013",
                                           ifelse(index_dt >= as.Date("2014-01-01"), "2014-2018", NA)),
                            cov_egfr_cat = factor(cov_egfr_cat, levels = c("eGFR >60", "eGFR 45-59", "eGFR 30-44",
                                                                           "eGFR 15-29", "eGFR <15")),
                            exp_egfr = factor(exp_egfr, levels = c("eGFR >60", "eGFR 45-59", "eGFR 30-44", "eGFR <30")),
                            cov_ef = factor(cov_ef, levels = c("pEF", "mEF", "rEF")),
                            cov_age_cat = factor(cov_age_cat, levels = c("<45 years", "45-65 years", "66-75 years", ">75 years")),
                            cov_cs = factor(cov_cs, levels = c("Cohabitating", "Living alone", "Missing")),
                            cov_edu = factor(cov_edu, levels = c("Compulsory school", "Secondary school", "University", "Missing")),
                            cov_loc = factor(cov_loc, levels = c("Out-patient", "In-patient")),
                            cov_smok = factor(cov_smok, levels = c("Never", "Former", "Current", "Missing")),
                            ## Impute missings with either median for continuous or new category for binary
                            # Continuous
                            mod_phhf = ifelse(is.na(cov_phhf), median(cov_phhf, na.rm = TRUE), cov_phhf),
                            mod_ntprobnp = ifelse(is.na(cov_ntprobnp), median(cov_ntprobnp, na.rm = TRUE), cov_ntprobnp),
                            mod_hr = ifelse(is.na(cov_hr), median(cov_hr, na.rm = TRUE), cov_hr),
                            mod_bpsys = ifelse(is.na(cov_bpsys), median(cov_bpsys, na.rm = TRUE), cov_bpsys),
                            mod_bpdia = ifelse(is.na(cov_bpdia), median(cov_bpdia, na.rm = TRUE), cov_bpdia),
                            mod_map = ifelse(is.na(cov_map), median(cov_map, na.rm = TRUE), cov_map),
                            mod_hb = ifelse(is.na(cov_hb), median(cov_hb, na.rm = TRUE), cov_hb),
                            mod_bmi = ifelse(is.na(cov_bmi), median(cov_bmi, na.rm = TRUE), cov_bmi),
                            mod_pot = ifelse(is.na(cov_pot), median(cov_pot, na.rm = TRUE), cov_pot),
                            # Categorical
                            cov_nyha = ifelse(cov_nyha == "Class I/II", "Class I/II",
                                              ifelse(cov_nyha == "Class III/IV", "Class III/IV", NA)),
                            mod_nyha = ifelse(is.na(cov_nyha), "Missing", cov_nyha),
                            mod_nyha = factor(mod_nyha, levels = c("Class I/II", "Class III/IV", "Missing")),
                            cov_hfdur = ifelse(cov_hfdur == "<6mo", "<6mo",
                                               ifelse(cov_hfdur == ">6mo", ">6mo", NA)),
                            mod_hfdur = ifelse(is.na(cov_hfdur), "Missing", cov_hfdur),
                            mod_hfdur = factor(mod_hfdur, levels = c("<6mo", ">6mo", "Missing")),
                            cov_pe = ifelse(cov_pe == "DCM", "DCM", 
                                            ifelse(cov_pe == "Heart valve disease", "Heart valve disease",
                                                   ifelse(cov_pe == "Hypertension", "Hypertension",
                                                          ifelse(cov_pe == "IHD", "IHD",
                                                                 ifelse(cov_pe == "Known alcoholic cardiomyopathy", 
                                                                        "Known alcoholic cardiomyopathy",
                                                                        ifelse(cov_pe == "Other", "Other", NA)))))),
                            mod_pe = ifelse(is.na(cov_pe), "Missing", cov_pe),
                            mod_pe = factor(mod_pe, levels = c("DCM","Heart valve disease", "Hypertension", "IHD",
                                                               "Known alcoholic cardiomyopathy", "Other", "Missing")),
                            # Binary
                            mod_ld = fac(cov_ld), mod_dx = fac(cov_dx), mod_diu = fac(cov_diu), mod_st = fac(cov_st), mod_ac = fac(cov_ac), 
                            mod_ap = fac(cov_ap), mod_nt = fac(cov_nt),
                            cov_an = ifelse(cov_an == "no", "No", ifelse(cov_an == "yes", "Yes", NA)),
                            mod_an = fac(cov_an),
                            mod_dev = ifelse(is.na(cov_dev), 3, cov_dev),
                            mod_dev = ifelse(mod_dev == 0, "No", ifelse(mod_dev == 1, "Yes", ifelse(mod_dev == 3, "Missing", NA))),
                            mod_dev = factor(mod_dev, levels = c("No", "Yes", "Missing")),
                            mod_ob = ifelse(cov_ob == "below_equal30", "below_equal30",
                                            ifelse(cov_ob == "over30", "over30", NA)),
                            mod_ob = ifelse(is.na(mod_ob), "Missing", mod_ob),
                            mod_ob = factor(mod_ob, levels = c("below_equal30", "over30", "Missing")))

### Save file
save(cohort, file = "cohort.Rdata")

rm(list = ls())
