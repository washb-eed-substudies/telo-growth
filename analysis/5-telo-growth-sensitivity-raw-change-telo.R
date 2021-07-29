rm(list=ls())

source(here::here("0-config.R"))

library(boxr)
box_auth()
d <- box_read_csv(839767614700)

# Z-score telomere length
d <- d %>%
  mutate(TS_t2_Z = scale(TS_t2, center=T, scale=T)[,1]) %>%
  mutate(TS_t3_Z = scale(TS_t3, center=T, scale=T)[,1]) %>%
  mutate(delta_TS_Z = scale(delta_TS, center=T, scale=T)[,1])


#### Association between change in telomere length between years 1 and 2 and growth ####
# immune ratios at y1 and growth velocity outcomes between y1 and y2
Xvars <- c("raw_delta_TS")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
           "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3", 
           "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")

#Fit models
H3_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H3_models <- bind_rows(H3_models, res)
  }
}

#Get primary contrasts
H3_res <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  preds <- predict_gam_diff(fit=H3_models$fit[i][[1]], d=H3_models$dat[i][[1]], quantile_diff=c(0.10,0.90), Xvar=res$X, Yvar=res$Y)
  H3_res <-  bind_rows(H3_res , preds$res)
}

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  simul_plot <- gam_simul_CI(H3_models$fit[i][[1]], H3_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H3_models, here("results/gam_models/unadjusted/raw_dtelo_models.RDS"))

#Save results
saveRDS(H3_res, here("results/gam_results/unadjusted/raw_dtelo_res.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("results/gam_figure_data/unadjusted/raw_dtelo_unadj_spline_data.RDS"))

#clean covariates
for(i in 1:ncol(d)){
  cat(colnames(d)[i],"  ",class(d[,i]),"\n")
}
#set variables as factors/numeric
d$sex<-as.factor(d$sex)
d$sex <- factor(d$sex, labels = c("female", "male"))
d$birthord<-as.factor(d$birthord)
d$momage<-as.numeric(d$momage)
d$momheight<-as.numeric(d$momheight)
d$momedu<-as.factor(d$momedu)
d$hfiacat<-as.factor(d$hfiacat)
d$Nlt18<-as.numeric(d$Nlt18)
d$Ncomp<-as.numeric(d$Ncomp)
d$watmin<-as.numeric(d$watmin)
d$floor<-as.factor(d$floor)
d$walls<-as.factor(d$walls)
d$elec<-as.factor(d$elec)
d$asset_wardrobe<-as.factor(d$asset_wardrobe)
d$asset_table<-as.factor(d$asset_table)
d$asset_chair<-as.factor(d$asset_chair)
d$asset_clock<-as.factor(d$asset_clock)
d$asset_khat<-as.factor(d$asset_khat)
d$asset_chouki<-as.factor(d$asset_chouki)
d$asset_radio<-as.factor(d$asset_radio)
d$asset_tv<-as.factor(d$asset_tv)
d$asset_refrig<-as.factor(d$asset_refrig)
d$asset_bike<-as.factor(d$asset_bike)
d$asset_moto<-as.factor(d$asset_moto)
d$asset_sewmach<-as.factor(d$asset_sewmach)
d$asset_mobile<-as.factor(d$asset_mobile)
d$n_cattle<-as.numeric(d$n_cattle)
d$n_goat<-as.numeric(d$n_goat)
d$n_chicken<-as.numeric(d$n_chicken)

d$lenhei_med_t2<-as.numeric(d$lenhei_med_t2)
d$weight_med_t2<-as.numeric(d$weight_med_t2)

d$monsoon_ht2<-as.factor(d$monsoon_ht2)
d$monsoon_ht2<-addNA(d$monsoon_ht2)
levels(d$monsoon_ht2)[length(levels(d$monsoon_ht2))]<-"Missing"

d$monsoon_ht3<-as.factor(d$monsoon_ht3)
d$monsoon_ht3<-addNA(d$monsoon_ht3)
levels(d$monsoon_ht3)[length(levels(d$monsoon_ht3))]<-"Missing"

d$ageday_ht2<-as.numeric(d$ageday_ht2)
d$ageday_ht3<-as.numeric(d$ageday_ht3)

d$anthro_days_btwn_t2_t3<-as.numeric(d$anthro_days_btwn_t2_t3)

d$tr <- factor(d$tr,levels=c("Control","Nutrition + WSH"))

d$cesd_sum_t2<-as.numeric(d$cesd_sum_t2)
d$cesd_sum_ee_t3<-as.numeric(d$cesd_sum_ee_t3)
d$pss_sum_mom_t3<-as.numeric(d$pss_sum_mom_t3)

d$diar7d_t2<-as.factor(d$diar7d_t2)
d$diar7d_t2<-addNA(d$diar7d_t2)
levels(d$diar7d_t2)[length(levels(d$diar7d_t2))]<-"Missing"

d$diar7d_t3<-as.factor(d$diar7d_t3)
d$diar7d_t3<-addNA(d$diar7d_t3)
levels(d$diar7d_t3)[length(levels(d$diar7d_t3))]<-"Missing"

d$life_viol_any_t3<-as.factor(d$life_viol_any_t3)
d$life_viol_any_t3<-addNA(d$life_viol_any_t3)
levels(d$life_viol_any_t3)[length(levels(d$life_viol_any_t3))]<-"Missing"

#### Association between change in telomere length between years 1 and 2 and growth ####
Xvars <- c("raw_delta_TS")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
           "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3", 
           "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")

Wvars_23_23<-c("sex","birthord", "momage","momheight","momedu", 
               "hfiacat", "Nlt18", "Ncomp", "watmin", 
               "floor", "walls", "elec", "asset_wardrobe", "asset_table", "asset_chair", "asset_clock", "asset_khat", 
               "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
               "asset_bike", "asset_moto", "asset_sewmach", "asset_mobile", 
               "n_cattle", "n_goat", "n_chicken", "monsoon_ht2", "monsoon_ht3", "ageday_ht2", 
               "ageday_ht3", "tr", "cesd_sum_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", "diar7d_t2", "diar7d_t3", 
               "life_viol_any_t3", "anthro_days_btwn_t2_t3")

Wvars_23_3<-c("sex","birthord", "momage", "momheight","momedu", 
              "hfiacat", "Nlt18", "Ncomp", "watmin", 
              "floor", "walls", "elec", "asset_wardrobe", "asset_table", "asset_chair", "asset_clock", "asset_khat", 
              "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
              "asset_bike", "asset_moto", "asset_sewmach", "asset_mobile", 
              "n_cattle", "n_goat", "n_chicken", "monsoon_ht2", "monsoon_ht3", "ageday_ht2", 
              "ageday_ht3", "tr", "cesd_sum_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", "diar7d_t2", "diar7d_t3", 
              "life_viol_any_t3", "lenhei_med_t2", "weight_med_t2", "anthro_days_btwn_t2_t3")

#Fit models
H3_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    if(grepl("delta|velocity", j)){
      Wvars <- Wvars_23_23
    }else{
      Wvars <- Wvars_23_3
    }
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.10,0.90), Xvar=res$X, Yvar=res$Y)
  H3_adj_res <-  bind_rows(H3_adj_res , preds$res)
}

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H3_adj_models$fit[i][[1]], H3_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H3_adj_models, here("results/gam_models/adjusted/raw_dtelo_adj_models.RDS"))

#Save results
saveRDS(H3_adj_res, here("results/gam_results/adjusted/raw_dtelo_adj_res.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("results/gam_figure_data/adjusted/raw_dtelo_adj_spline_data.RDS"))


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

results <- rbind(H3_res, readRDS(here("results/gam_results/unadjusted/dtelo_res.RDS")) %>% select(!BH.Pval))
results_adj <- rbind(H3_adj_res, readRDS(here('results/gam_results/adjusted/dtelo_adj_res.RDS')) %>% select(!BH.Pval))

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
