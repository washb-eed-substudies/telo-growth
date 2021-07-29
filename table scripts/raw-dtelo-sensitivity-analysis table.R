rm(list=ls())

source(here::here("0-config.R"))

library(boxr)
box_auth()
d <- box_read_csv(839767614700)

library(flextable)
library(officer)

# make table comparing results
exposure <- c("delta_TS", "raw_delta_TS")
outcome <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
             "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3", 
             "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")
expo_var <- c("Corrected Change in telomere length between Year 1 and Year 2", "Uncorrected Change in telomere length between Year 1 and Year 2")
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2", 
             "Change in LAZ between Year 1 and Year 2", "Change in WAZ between Year 1 and Year 2", 
             "Change in WLZ between Year 1 and Year 2", "Change in HCZ between Year 1 and Year 2",
             "Length velocity between Year 1 and Year 2", "Weight velocity between Year 1 and Year 2",
             "Head circumference velocity between Year 1 and Year 2")

results <- rbind(readRDS(here("results/gam_results/unadjusted/raw_dtelo_res.RDS")), readRDS(here("results/gam_results/unadjusted/dtelo_res.RDS")) %>% select(!BH.Pval))
results_adj <- rbind(readRDS(here('results/gam_results/adjusted/raw_dtelo_adj_res.RDS')), readRDS(here('results/gam_results/adjusted/dtelo_adj_res.RDS')) %>% select(!BH.Pval))

tbl <- data.table(matrix(nrow=0, ncol=13))
skipped<-F
for (i in 1:length(exposure)) {
  for (j in 1:length(outcome)) {
    exp <- exposure[i]
    out <- outcome[j]
    filtered_res <- results[results$Y==out & results$X==exp,]
    filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
    unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
    adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
    if (nrow(filtered_res)==0){
      skipped<-T
      next
    }
    if(j==1|skipped==T){
      tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                             round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), 
                             round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2)))
      skipped=F
    }else {
      tbl <- rbind(tbl, list(" ", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                             round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), 
                             round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2)))
    }
  }
  if (i != length(exposure)) {
    tbl <- rbind(tbl, list("","","","","","","","", "","","","",""))
  }
}

flextbl<-flextable(tbl, col_keys=names(tbl))
flextbl <- set_header_labels(flextbl,
                             values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                           "V6" = "Predicted Outcome at 10th Percentile", "V7" = "Predicted Outcome at 90th Percentile", "V8" = "Coefficient (95% CI)", "V9" = "P-value",
                                           "V10" = "Predicted Outcome at 10th Percentile", "V11" = "Predicted Outcome at 90th Percentile", "V12" = "Coefficient (95% CI)", "V13" = "P-value"))
flextbl <- add_header_row(flextbl, values = c("","","","","", "Unadjusted", "Adjusted"), colwidths=c(1,1,1,1,1,4,4))
flextbl <- add_header_row(flextbl, values = c("Exposure", "Outcome","N","10th Percentile","90th Percentile", "Outcome, 90th Percentile v. 10th Percentile"), colwidths=c(1,1,1,1,1,8))

flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
flextbl <- align(flextbl, align = "center", part = "all")
flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
flextbl <- fontsize(flextbl, part = "all", size = 6)
flextbl <- width(flextbl, 1:13, width=c(1.1, 1.3, .3, .5, .5, .5, .5, 1, .3, .5, .5, 1, .3))

sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)
save_as_docx("Comparison of Regression-to-the-mean Corrected Change in Telomere Length and Raw Change in Telomere Length" = flextbl, 
             path="C:/Users/Sophia/Documents/WASH/WASH Telomeres and Growth/telo-growth change in telo sensitivity analysis tables.docx",
             pr_section = sect_properties)
