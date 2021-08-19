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

outliers <- read.csv("C:/Users/Sophia/Downloads/outliers.csv") %>% select(!X)
outliers %>% head()

laz <- filter(outliers, lhflag ==1)
hc <- filter(outliers, hcflag ==1)
wei <- filter(outliers, weiflag ==1)
wh <- filter(outliers, whflag ==1)

d$laz_t2 <- ifelse(d$childid %in% laz$childid, NA, d$laz_t2)
d$laz_t3 <- ifelse(d$childid %in% laz$childid, NA, d$laz_t3)
d$delta_laz_t2_t3 <- ifelse(d$childid %in% laz$childid, NA, d$delta_laz_t2_t3)
d$hcz_t2 <- ifelse(d$childid %in% hc$childid, NA, d$hcz_t2)
d$hcz_t3 <- ifelse(d$childid %in% hc$childid, NA, d$hcz_t3)
d$delta_hcz_t2_t3 <- ifelse(d$childid %in% hc$childid, NA, d$delta_hcz_t2_t3)
d$waz_t2 <- ifelse(d$childid %in% wei$childid, NA, d$waz_t2)
d$waz_t3 <- ifelse(d$childid %in% wei$childid, NA, d$waz_t3)
d$delta_waz_t2_t3 <- ifelse(d$childid %in% wei$childid, NA, d$delta_waz_t2_t3)
d$whz_t2 <- ifelse(d$childid %in% wh$childid, NA, d$whz_t2)
d$whz_t3 <- ifelse(d$childid %in% wh$childid, NA, d$whz_t3)
d$delta_whz_t2_t3 <- ifelse(d$childid %in% wh$childid, NA, d$delta_whz_t2_t3)



#Example:

#Fit GAM model with random effects for childid
#res_unadj <- fit_RE_gam(d=d, X="t3_cort_z01", Y="laz_t3",  W=NULL)

#Get predictions of differences from the 25th percentile of exposure
#preds_unadj <- predict_gam_diff(fit=res_unadj$fit, d=res_unadj$dat, quantile_diff=c(0.10,0.90), Xvar="delta_TS", Yvar="laz_t3")


#Primary parameter we are estimating: difference between 25th and 75th percentile of the exposure
#preds_unadj$res

#Plot the difference from the 25th percentile for the full range of the exposure:
#NOTE: not making these plots anymore, just using for diagnostics
#p <- plot_gam_diff(preds_unadj$plotdf)
#print(p)

#Fit spline with simultaneous confidence intervals
#simul_plot <- gam_simul_CI(res_unadj$fit, res_unadj$dat, xlab="delta_TS", ylab="laz_t3", title="example title")
#simul_plot$p


#### Loop over exposure-outcome pairs ####

#### Association between telomere length at year 1 and child growth ####
Xvars <- c("TS_t2", "TS_t2_Z")            
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2", "laz_t3", "waz_t3", "whz_t3", "hcz_t3", 
           "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3",
           "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3") 

#Fit models
H1_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}

#Get primary contrasts
H1_res <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  preds <- predict_gam_diff(fit=H1_models$fit[i][[1]], d=H1_models$dat[i][[1]], quantile_diff=c(0.10,0.90), Xvar=res$X, Yvar=res$Y)
  H1_res <-  bind_rows(H1_res , preds$res)
}

#Make list of plots
H1_plot_list <- NULL
H1_plot_data <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  simul_plot <- gam_simul_CI(H1_models$fit[i][[1]], H1_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_plot_list[[i]] <-  simul_plot$p
  H1_plot_data <-  rbind(H1_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H1_models, here("results/gam_models/unadjusted/telot2_models_nooutliers.RDS"))

#Save results
saveRDS(H1_res, here("results/gam_results/unadjusted/telot2_res_nooutliers.RDS"))

#Save plot data
saveRDS(H1_plot_data, here("results/gam_figure_data/unadjusted/telot2_unadj_spline_data_nooutliers.RDS"))



#### Association between telomere length at year 2 and growth ####
# all immune outcomes at y1 v. growth at y2
Xvars <- c("TS_t3", "TS_t3_Z")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 

#Fit models
H2_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H2_models <- bind_rows(H2_models, res)
  }
}

#Get primary contrasts
H2_res <- NULL
for(i in 1:nrow(H2_models)){
  res <- data.frame(X=H2_models$X[i], Y=H2_models$Y[i])
  preds <- predict_gam_diff(fit=H2_models$fit[i][[1]], d=H2_models$dat[i][[1]], quantile_diff=c(0.10,0.90), Xvar=res$X, Yvar=res$Y)
  H2_res <-  bind_rows(H2_res , preds$res)
}

#Make list of plots
H2_plot_list <- NULL
H2_plot_data <- NULL
for(i in 1:nrow(H2_models)){
  res <- data.frame(X=H2_models$X[i], Y=H2_models$Y[i])
  simul_plot <- gam_simul_CI(H2_models$fit[i][[1]], H2_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_plot_list[[i]] <-  simul_plot$p
  H2_plot_data <-  rbind(H2_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H2_models, here("results/gam_models/unadjusted/telot3_models_nooutliers.RDS"))

#Save results
saveRDS(H2_res, here("results/gam_results/unadjusted/telot3_res_nooutliers.RDS"))

#Save plot data
saveRDS(H2_plot_data, here("results/gam_figure_data/unadjusted/telot3_unadj_spline_data_nooutliers.RDS"))



#### Association between change in telomere length between years 1 and 2 and growth ####
# immune ratios at y1 and growth velocity outcomes between y1 and y2
Xvars <- c("delta_TS", "delta_TS_Z")            
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
saveRDS(H3_models, here("results/gam_models/unadjusted/dtelo_models_nooutliers.RDS"))

#Save results
saveRDS(H3_res, here("results/gam_results/unadjusted/dtelo_res_nooutliers.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("results/gam_figure_data/unadjusted/dtelo_unadj_spline_data_nooutliers.RDS"))
