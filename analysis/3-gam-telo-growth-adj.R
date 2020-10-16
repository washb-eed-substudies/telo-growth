rm(list=ls())

source(here::here("0-config.R"))
source(here::here("analysis/1-gam-functions.R"))

load(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.RData"))
# Z-score telomere length
d <- d %>%
  mutate(TS_t2_Z = scale(TS_t2, center=T, scale=T)[,1]) %>%
  mutate(TS_t3_Z = scale(TS_t3, center=T, scale=T)[,1])


#Loop over exposure-outcome pairs

#### Association between telomere length at year 1 and child growth ####
Xvars <- c("TS_t2", "TS_t2_Z")            
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2", "laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
           "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3",
           "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3") 

Wvars_2<-c("sex","birthord", "momage", "momheight","momedu", 
           "hfiacat", "Nlt18", "Ncomp", "watmin", "floor", 
           "walls", "elec", "asset_wardrobe", "asset_table", 
           "asset_chair", "asset_clock", "asset_khat", 
           "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
           "asset_bike", "asset_moto", "asset_sewmach",
           "asset_mobile", "n_cattle", "n_goat", "n_chicken",
           "monsoon_ht2", "ageday_ht2", "tr", "cesd_sum_t2", "diar7d_t2", 
           "life_viol_any_t3")

#Fit models
H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars_2)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H1_adj_res <-  bind_rows(H1_adj_res , preds$res)
}
H1_adj_res$adjusted <- 0

#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H1_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/H1_adj_models.RDS"))

#Save results
saveRDS(H1_adj_res, here("results/gam_results/adjusted/H1_adj_res.RDS"))


#Save plots
#saveRDS(H1_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H1_adj_splines.RDS"))

#Save plot data
#saveRDS(H1_adj_plot_data, paste0(dropboxDir,"results/stress-growth-models/figure-data/H1_adj_spline_data.RDS"))



#### Association between telomere length at year 2 and growth ####
# all immune outcomes at y1 v. growth at y2
Xvars <- c("TS_t3", "TS_t3_Z")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 

Wvars_H3<-c("sex","birthord", "momage", "momheight","momedu", 
            "hfiacat", "Nlt18", "Ncomp", "watmin", 
            "floor", "walls", "elec", "asset_wardrobe", "asset_table", "asset_chair", "asset_clock", "asset_khat", 
            "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
            "asset_bike", "asset_moto", "asset_sewmach", "asset_mobile", 
            "n_cattle", "n_goat", "n_chicken", "monsoon_ht2", "monsoon_ht3", "ageday_ht2", 
            "ageday_ht3", "tr", "cesd_sum_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", "diar7d_t2", "diar7d_t3", 
            "life_viol_any_t3", "lenhei_med_t2", "weight_med_t2")

#Fit models
H2_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars_H3)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

#Get primary contrasts
H2_adj_res <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_adj_res <-  bind_rows(H2_adj_res , preds$res)
}
H2_adj_res$adjusted <- 0

#Make list of plots
H2_plot_list <- NULL
H2_plot_data <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_models$fit[i][[1]], H2_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_plot_list[[i]] <-  simul_plot$p
  H2_plot_data <-  rbind(H2_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H2_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_H2_adj_models.RDS"))

#Save results
saveRDS(H2_adj_res, here("results/gam_results/adjusted/H2_adj_res.RDS"))


#Save plots
#saveRDS(H2_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H2_adj_splines.RDS"))

#Save plot data
#saveRDS(H2_plot_data, paste0(dropboxDir,"results/stress-growth-models/figure-data/H2_adj_spline_data.RDS"))



#### Association between change in telomere length between years 1 and 2 and growth ####
# immune ratios at y1 and growth velocity outcomes between y1 and y2
Xvars <- c("delta_TS")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
           "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3", 
           "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")

Wvars_2_3<-c("sex","birthord", "momage","momheight","momedu", 
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
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars_2_3)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_adj_res <-  bind_rows(H3_adj_res , preds$res)
}
H3_adj_res$adjusted <- 0

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
#saveRDS(H3_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_H3_adj_models.RDS"))

#Save results
saveRDS(H3_adj_res, here("results/gam_results/adjusted/H3_adj_res.RDS"))


#Save plots
#saveRDS(H3_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H3_adj_splines.RDS"))

#Save plot data
#saveRDS(H3_plot_data, paste0(dropboxDir,"results/stress-growth-models/figure-data/H3_adj_spline_data.RDS"))

