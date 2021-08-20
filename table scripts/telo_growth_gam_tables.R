rm(list=ls())

library('flextable')
library('officer')
source(here::here("0-config.R"))
source(here::here("table scripts/table-functions.R"))


# load enrollment characteristics and results
H1 <- readRDS(here('results/gam_results/unadjusted/telot2_res.RDS'))
H2 <- readRDS(here('results/gam_results/unadjusted/telot3_res.RDS'))
H3 <- readRDS(here('results/gam_results/unadjusted/dtelo_res.RDS'))
H1adj <- readRDS(here('results/gam_results/adjusted/telot2_adj_res.RDS'))
H2adj <- readRDS(here('results/gam_results/adjusted/telot3_adj_res.RDS'))
H3adj <- readRDS(here('results/gam_results/adjusted/dtelo_adj_res.RDS'))


#### TABLES ####
#### GAM Table 1 ####

exposure <- c("TS_t2")
outcome <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2", "laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
             "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3",
             "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")
expo_var <- c("Telomere length at Year 1 (T/S Ratio)")
out_var <- c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1",
             "LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2", 
             "Change in LAZ between Year 1 and Year 2", "Change in WAZ between Year 1 and Year 2", 
             "Change in WLZ between Year 1 and Year 2", "Change in HCZ between Year 1 and Year 2",
             "Length velocity between Year 1 and Year 2", "Weight velocity between Year 1 and Year 2",
             "Head circumference velocity between Year 1 and Year 2")

tbl1 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl1flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H1, H1adj, T, exp_col_size = 1.3, out_col_size = 1.5)
tbl1supp <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl1flexsupp <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H1, H1adj, F, exp_col_size = 1, out_col_size = 1.3)

#### GAM Table 2 ####

exposure <- c("TS_t3")
outcome <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3")
expo_var <- c("Telomere length at Year 2 (T/S Ratio)")
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

tbl2 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl2flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H2, H2adj, T, 1.1, .7)
tbl2supp <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl2flexsupp <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H2, H2adj, F, 1.1, .7)


#### GAM Table 3 ####

exposure <- c("delta_TS")
outcome <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
             "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3", 
             "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")
expo_var <- c("Change in telomere length between Year 1 and Year 2 (T/S Ratio)")
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2", 
             "Change in LAZ between Year 1 and Year 2", "Change in WAZ between Year 1 and Year 2", 
             "Change in WLZ between Year 1 and Year 2", "Change in HCZ between Year 1 and Year 2",
             "Length velocity between Year 1 and Year 2", "Weight velocity between Year 1 and Year 2",
             "Head circumference velocity between Year 1 and Year 2")

tbl3 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H3, H3adj, T)
tbl3flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H3, H3adj, T, 1.1, 1.3)
tbl3supp <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl3flexsupp <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H3, H3adj, F, 1.1, 1.3)




#### SAVE TABLES ####

write.csv(tbl1, file=here("tables/gam/telo-growth-table2.csv"))
write.csv(tbl2, here('tables/gam/telo-growth-table3.csv'))
write.csv(tbl3, here('tables/gam/telo-growth-table4.csv'))
write.csv(tbl1supp, file=here("tables/gam/telo-growth-supptable1.csv"))
write.csv(tbl2supp, here('tables/gam/telo-growth-supptable2.csv'))
write.csv(tbl3supp, here('tables/gam/telo-growth-supptable3.csv'))

save_as_docx("Table 2: Association Between Telomere Length at Year 1 and Growth" = tbl1flex, 
             "Table 3: Association Between Telomere Length at Year 2 and Growth" = tbl2flex, 
             "Table 4: Association Between Change in Telomere Length and Growth" = tbl3flex, 
             path="C:/Users/Sophia/Documents/WASH/WASH Telomeres and Growth/telo-growth main gam tables 08.19.21.docx",
             pr_section = sect_properties)

save_as_docx("Table S1: Association Between Telomere Length at Year 1 and Growth" = tbl1flexsupp, 
             "Table S2: Association Between Telomere Length at Year 2 and Growth" = tbl2flexsupp, 
             "Table S3: Association Between Change in Telomere Length and Growth" = tbl3flexsupp, 
             path="C:/Users/Sophia/Documents/WASH/WASH Telomeres and Growth/telo-growth supplementary gam tables 08.19.21.docx",
             pr_section = sect_properties)

