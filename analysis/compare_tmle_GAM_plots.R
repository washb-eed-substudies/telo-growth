rm(list=ls())

source(here::here("0-config.R"))
source(here::here("analysis/1-gam-functions.R"))
source(here::here("analysis/0-base-quartileTMLE_functions.R"))

load(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.RData"))
# Z-score telomere length
d <- d %>%
  mutate(TS_t2_Z = scale(TS_t2, center=T, scale=T)[,1]) %>%
  mutate(TS_t3_Z = scale(TS_t3, center=T, scale=T)[,1]) %>%
  mutate(delta_TS_Z = scale(delta_TS, center=T, scale=T)[,1])


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




#Loop over exposure-outcome pairs

#### Association between telomere length at year 1 and child growth ####


#### Association between change in telomere length between years 1 and 2 and growth ####
# immune ratios at y1 and growth velocity outcomes between y1 and y2
H3_adj_models <- NULL
i <- Xvars <- c("TS_t2")            
j <- Yvars <- c("laz_t2")


Wvars_2<-c("sex","birthord", "momage", "momheight","momedu", 
           "hfiacat", "Nlt18", "Ncomp", "watmin", "floor", 
           "walls", "elec", "asset_wardrobe", "asset_table", 
           "asset_chair", "asset_clock", "asset_khat", 
           "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
           "asset_bike", "asset_moto", "asset_sewmach",
           "asset_mobile", "n_cattle", "n_goat", "n_chicken",
           "monsoon_ht2", "ageday_ht2", "tr", "cesd_sum_t2", "diar7d_t2", 
           "life_viol_any_t3")


Wvars_H7<-c("sex","birthord", "momage", "momheight","momedu", 
            "hfiacat", "Nlt18", "Ncomp", "watmin", 
            "floor", "walls", "elec", "asset_wardrobe", "asset_table", "asset_chair", "asset_clock", "asset_khat", 
            "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
            "asset_bike", "asset_moto", "asset_sewmach", "asset_mobile", 
            "n_cattle", "n_goat", "n_chicken", "monsoon_ht2", "monsoon_ht3", "ageday_ht2", 
            "ageday_ht3", "tr", "cesd_sum_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", "diar7d_t2", "diar7d_t3", 
            "life_viol_any_t3", "lenhei_med_t2", "weight_med_t2")


    res_adj <- fit_RE_gam(d=d, X=Xvars, Y=Yvars,  W=Wvars_2)
    res <- data.frame(X=Xvars, Y=Yvars, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.1,0.9), Xvar=res$X, Yvar=res$Y)
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



h1adj.res=NULL
h1adj.res <- tmle_quart(dat=d, 
                        Y=Yvars, 
                        W=Wvars_2, 
                        A=Xvars, 
                        id="block",
                        Alevels=c("Q1","Q2","Q3","Q4"),
                        outputdf = h1adj.res,
                        family="gaussian", 
                        SLlibrary="SL.gam")





j <- Yvars <- c("delta_whz_t2_t3")


H3_adj_models=NULL
res_adj <- fit_RE_gam(d=d, X=Xvars, Y=Yvars,  W=Wvars_H7)
res <- data.frame(X=Xvars, Y=Yvars, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
H3_adj_models <- bind_rows(H3_adj_models, res)

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.1,0.9), Xvar=res$X, Yvar=res$Y)
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

preds$res
plot_gam_diff(preds$plotdf)


h3adj.res=NULL
h3adj.res <- tmle_quart(dat=d, 
                        Y=Yvars, 
                        W=Wvars_H7, 
                        A=Xvars, 
                        id="block",
                        Alevels=c("Q1","Q2","Q3","Q4"),
                        outputdf = h3adj.res,
                        family="gaussian", 
                        SLlibrary="SL.gam")



h3adj.res$cutpoints <- gsub("\\[","",h3adj.res$cutpoints)
h3adj.res$cutpoints <- gsub(")","",h3adj.res$cutpoints)
h3adj.res$cutpoints <- gsub("<","",h3adj.res$cutpoints)
h3adj.res$cutpoints <- gsub(">=","",h3adj.res$cutpoints)
h3adj.res$cutpoints <- str_split(h3adj.res$cutpoints, ",", simplify = T)[,1]
h3adj.res$cutpoints <- as.numeric(h3adj.res$cutpoints)

tmle_res <- h3adj.res %>% select(level, meanY, ATE, cutpoints) %>% pivot_wider(names_from=level, values_from=c(meanY, ATE, cutpoints)) %>% 
  mutate(Y="delta_whz_t2_t3", X="TS_t2")

plotdf <- left_join(preds$plotdf, tmle_res, by=c("X","Y"))

plotdf <- plotdf %>%
  mutate(
    ATE_Q2 = meanY_Q1 + ATE_Q2,
    ATE_Q3 = meanY_Q1 + ATE_Q3,
    ATE_Q4 = meanY_Q1 + ATE_Q4
  )

head(plotdf)

p<-ggplot(plotdf) + geom_ribbon(aes(x=x, ymin=pred.q1+lb.diff, ymax=pred.q1+ub.diff), alpha=0.5) + 
  geom_path(aes(x=x, y=pred.q1+lb.diff), color="blue")+
  geom_path(aes(x=x, y=pred.q1+ub.diff), color="red")+
  geom_path(aes(x=x, y=pred.q1+point.diff), color="black") + 
  geom_vline(aes(xintercept=q1)) +
  geom_vline(aes(xintercept=q3)) + 
  geom_vline(aes(xintercept=cutpoints_Q2), linetype="dashed", color="green") +
  geom_vline(aes(xintercept=cutpoints_Q3), linetype="dashed", color="green") +
  geom_vline(aes(xintercept=cutpoints_Q4), linetype="dashed", color="green") +
  geom_point(aes(x= cutpoints_Q2 - (cutpoints_Q2 - min(x))/2, y=meanY_Q1)) +
  geom_point(aes(x= cutpoints_Q3 - (cutpoints_Q3 - cutpoints_Q2)/2, y=meanY_Q2)) +
  geom_point(aes(x= cutpoints_Q4 - (cutpoints_Q4 - cutpoints_Q3)/2, y=meanY_Q3)) +
  geom_point(aes(x= cutpoints_Q4 + (max(x) - cutpoints_Q4)/2, y=meanY_Q4)) +
  
  geom_point(aes(x= cutpoints_Q3 - (cutpoints_Q3 - cutpoints_Q2)/2, y=ATE_Q2), color="green") +
  geom_point(aes(x= cutpoints_Q4 - (cutpoints_Q4 - cutpoints_Q3)/2, y=ATE_Q3), color="green") +
  geom_point(aes(x= cutpoints_Q4 + (max(x) - cutpoints_Q4)/2, y=ATE_Q4), color="green") +
  xlab("Exposure") + 
  ylab("GAM-estimated differences from 25th percentile of exposure")
p



