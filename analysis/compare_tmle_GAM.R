

rm(list=ls())
source(here::here("0-config.R"))

source(here::here("analysis/0-base-quartileTMLE_functions.R"))

#Set seed to replicate with R versions less than 3.6.1
RNGkind(sample.kind = "Rounding")


#load covariates, exposures, outcomes dataset
load(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.RData"))



#check if variables are factors/numeric
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



#Hypothesis 1
#The change in telomere length between Year 1 and Year 2 is associated with the change in child length-for-age Z score (LAZ), weight-for-age Z score (WAZ), weight-for-length Z score (WLZ), and head circumference-for-age Z score from Year 1 to Year 2.
#Exposure: Quartiles of change in telomere length between Year 1 and Year 2
#Outcome: Change in child length-for-age Z score (LAZ), weight-for-age Z score (WAZ), weight-for-length Z score (WLZ), and head circumference-for-age Z score from Year 1 to Year 2

#unadjusted, with continuous outcome

#Null data.frame
h1unadj.res = NULL

h1unadj.res <- tmle_quart(dat=d, 
                          Y="delta_laz_t2_t3", 
                          W=NULL, 
                          A="delta_TS", 
                          id="block",
                          Alevels=c("Q1","Q2","Q3","Q4"),
                          outputdf = h1unadj.res,
                          family="gaussian", 
                          SLlibrary="SL.gam")




rm(list=ls())

source(here::here("0-config.R"))
source(here::here("analysis/1-gam-functions.R"))

load(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.RData"))
# Z-score telomere length
d <- d %>%
  mutate(TS_t2_Z = scale(TS_t2, center=T, scale=T)[,1]) %>%
  mutate(TS_t3_Z = scale(TS_t3, center=T, scale=T)[,1]) %>%
  mutate(delta_TS_Z = scale(delta_TS, center=T, scale=T)[,1])


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

#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}

