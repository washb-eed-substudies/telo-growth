rm(list=ls())

library('flextable')
library('officer')
source(here::here("0-config.R"))

# load enrollment characteristics and results
load(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.RData"))
H1 <- readRDS(here('results/gam_results/unadjusted/H1_res_BH.RDS'))
H2 <- readRDS(here('results/gam_results/unadjusted/H2_res_BH.RDS'))
H3 <- readRDS(here('results/gam_results/unadjusted/H3_res_BH.RDS'))
H1adj <- readRDS(here('results/gam_results/adjusted/H1_adj_res_BH.RDS'))
H2adj <- readRDS(here('results/gam_results/adjusted/H2_adj_res_BH.RDS'))
H3adj <- readRDS(here('results/gam_results/adjusted/H3_adj_res_BH.RDS'))


#### Functions for growth tables ####
growth_tbl <- function(name, expo_var, out_var, exposure, outcome, results, results_adj){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as a csv
  
  tbl <- data.table(" " = character(), " " = character(), " " = character(), " " = character(),
                    " Outcome, 90th pctl v. 10th pctl" = character(), " " = character(), " " = character(), " " = character(), 
                    " " = character(), " " = character(), " " = character(), " " = character())
  tbl <- rbind(tbl, list(" ", " ", " ", " ", "Unadjusted", " ", " ", " ", "Fully adjusted", " ", " ", " "))
  tbl <- rbind(tbl, list(" ", "Outcome", "10th pctl Mean", "90th pctl Mean", 
                         "Predicted Outcome at 10th pctl", "Predicted Outcome at 90th pctl", "Coefficient (95% CI)", "P-value", 
                         "Predicted Outcome at 10th pctl", "Predicted Outcome at 90th pctl", "Coefficient (95% CI)", "P-value"))
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (j==1){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }else {
        tbl <- rbind(tbl, list("", out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","","","","",""))
    }
  }
  tbl
}

growth_tbl_flex <- function(name, expo_var, out_var, exposure, outcome, results, results_adj){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as an image or 
  ### directly to a word document!
  
  # build table
  tbl <- data.table(matrix(nrow=0, ncol=12))
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (j==1){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }else {
        tbl <- rbind(tbl, list(" ", out_var[j], round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$BH.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","", "","","",""))
    }
  }
  
  # format for export
  flextbl<-flextable(tbl, col_keys=names(tbl))
  flextbl <- set_header_labels(flextbl,
                               values = list("V1" = name, "V2" = "Outcome", "V3" = "10th pct1 Mean", "V4" = "90th pctl Mean",
                                             "V5" = "Predicted Outcome at 10th pctl", "V6" = "Predicted Outcome at 90th pctl", "V7" = "Coefficient (95% CI)", "V8" = "P-value",
                                             "V9" = "Predicted Outcome at 10th pctl", "V10" = "Predicted Outcome at 90th pctl", "V11" = "Coefficient (95% CI)", "V12" = "P-value"))
  flextbl <- add_header_row(flextbl, values = c("","","","", "Unadjusted", "Fully adjusted"), colwidths=c(1,1,1,1,4,4))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- add_header_row(flextbl, values = c("","","","", "Outcome, 90th pctl v. 10th pctl"), colwidths=c(1,1,1,1,8))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
  flextbl <- autofit(flextbl, part = "all")
  flextbl <- fit_to_width(flextbl, max_width=8)
  
  flextbl
}



#### TABLES ####
#### GAM Table 1 ####

exposure <- c("TS_t2", "TS_t2_Z")
outcome <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2", "laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
             "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3",
             "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")
expo_var <- c("Telomere length at Year 1", "Z-score of telomere length at Year 1")
out_var <- c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1",
             "LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2", 
             "Change in LAZ between Year 1 and Year 2", "Change in WAZ between Year 1 and Year 2", 
             "Change in WLZ between Year 1 and Year 2", "Change in HCZ between Year 1 and Year 2",
             "Length velocity between Year 1 and Year 2", "Weight velocity between Year 1 and Year 2",
             "Head circumference velocity between Year 1 and Year 2")

tbl1 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl1flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H1, H1adj)

#### GAM Table 2 ####

exposure <- c("TS_t3", "TS_t3_Z")
outcome <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3")
expo_var <- c("Telomere length at Year 2", "Z-score of telomere length at Year 2")
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

tbl2 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl2flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H2, H2adj)


#### GAM Table 3 ####

exposure <- c("delta_TS", "delta_TS_Z")
outcome <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
             "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3", 
             "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")
expo_var <- c("Change in telomere length between Year 1 and Year 2", "Z-score change in telomere length between Year 1 and Year 2")
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2", 
             "Change in LAZ between Year 1 and Year 2", "Change in WAZ between Year 1 and Year 2", 
             "Change in WLZ between Year 1 and Year 2", "Change in HCZ between Year 1 and Year 2",
             "Length velocity between Year 1 and Year 2", "Weight velocity between Year 1 and Year 2",
             "Head circumference velocity between Year 1 and Year 2")

tbl3 <- growth_tbl("Exposure", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl3flex <- growth_tbl_flex("Exposure", expo_var, out_var, exposure, outcome, H3, H3adj)




#### SAVE TABLES ####

write.csv(tbl1, file=here("tables/gam/telo-growth-table1.csv"))
write.csv(tbl2, here('tables/gam/telo-growth-table2.csv'))
write.csv(tbl3, here('tables/gam/telo-growth-table3.csv'))

save_as_docx("Table 1" = tbl1flex, "Table 2" = tbl2flex, "Table 3" = tbl3flex, path="C:/Users/Sophia/Documents/WASH/WASH Telomeres and Growth/GAM Tables v3.docx")
