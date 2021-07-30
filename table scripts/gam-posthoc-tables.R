rm(list=ls())

library('flextable')
library('officer')
source(here::here("0-config.R"))
source(here::here("table scripts/table-functions.R"))


# load enrollment characteristics and results
H <- readRDS(here('results/gam_results/unadjusted/posthoc_res.RDS'))
Hadj <- readRDS(here('results/gam_results/adjusted/posthoc_adj_res.RDS'))


#### TABLES ####
exposure <- c("laz_t1", "waz_t1", "whz_t1" ,"hcz_t1")
outcome <- c("TS_t2", "TS_t3")
expo_var <- c("LAZ Month 3", "WAZ Month 3", "WLZ Month 3", "HCZ Month 3")
out_var <- c("Telomere length at Year 1", "Telomere length at Year 2")

posthoc1 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H, Hadj)
posthoc1flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H, Hadj, F)


exposure <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2")
outcome <- c("TS_t3")
expo_var <- c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1")
out_var <- c("Telomere length at Year 2")

posthoc2 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H, Hadj)
posthoc2flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H, Hadj, F)


#### GAM Table 3 ####

exposure <- c("delta_laz_t1_t2", "delta_waz_t1_t2", "delta_whz_t1_t2", "delta_hcz_t1_t2",
              "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3",
              "delta_laz_t1_t3", "delta_waz_t1_t3", "delta_whz_t1_t3", "delta_hcz_t1_t3")
outcome <- c("TS_t2", "TS_t3")
expo_var <- c("Change in LAZ between Month 3 and Year 1", "Change in WAZ between Month 3 and Year 1", 
              "Change in WLZ between Month 3 and Year 1", "Change in HCZ between Month 3 and Year 1",
              "Change in LAZ between Year 1 and Year 2", "Change in WAZ between Year 1 and Year 2", 
              "Change in WLZ between Year 1 and Year 2", "Change in HCZ between Year 1 and Year 2",
              "Change in LAZ between Month 3 and Year 2", "Change in WAZ between Month 3 and Year 2", 
              "Change in WLZ between Month 3 and Year 2", "Change in HCZ between Month 3 and Year 2")
out_var <- c("Telomere length at Year 1", "Telomere length at Year 2")

posthoc3 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H, Hadj)
posthoc3flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H, Hadj, F, 1.3, 1)

#### GAM Table 4 ####

exposure <- c("len_velocity_t1_t2", "wei_velocity_t1_t2", "hc_velocity_t1_t2",
              "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3",
              "len_velocity_t1_t3", "wei_velocity_t1_t3", "hc_velocity_t1_t3")
outcome <- c("TS_t2", "TS_t3")
expo_var <- c("Length velocity between Month 3 and Year 1", "Weight velocity between Month 3 and Year 1", "Head circumference velocity between Month 3 and Year 1",
              "Length velocity between Year 1 and Year 2", "Weight velocity between Year 1 and Year 2", "Head circumference velocity between Year 1 and Year 2",
              "Length velocity between Month 3 and Year 2", "Weight velocity between Month 3 and Year 2", "Head circumference velocity between Month 3 and Year 2")
out_var <- c("Telomere length at Year 1", "Telomere length at Year 2")

posthoc4 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H, Hadj)
posthoc4flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H, Hadj, F, 1.3, 1)




#### SAVE TABLES ####
save_as_docx("Association Between Growth at Month 3 and Subsequent Telomere Length" = posthoc1flex, 
             "Association Between Growth at Year 1 and Subsequent Telomere Length" = posthoc2flex, 
             "Association Between Change in Growth and Telomere Length" = posthoc3flex, 
             "Association Between Change in Telomere Length and Growth" = posthoc4flex,
             path="C:/Users/Sophia/Documents/WASH/WASH Telomeres and Growth/telo-growth gam posthoc tables.docx",
             pr_section = sect_properties)

