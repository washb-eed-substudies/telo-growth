




#---------------------------------------
# quartileTMLE_telo_growth.R
#
# audrie lin (audrielin@berkeley.edu)
#
# Summarize unadjusted and adjusted mean growth outcomes by quartiles of telomere length. 
# Calculate means within each quartile and estimate adjusted means differences 
# between quartiles using TMLE
#---------------------------------------

#---------------------------------------

# input files
# 
# ~/Dropbox/WBB-EE-analysis/Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.csv (from 3-bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.do)
# "~/Dropbox/WBB-EE-analysis/Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.RData (from ~/ee-secondary/audrie R scripts/observational/0-file_conversion.R)
#
# output files:
#	
# 
# 
#---------------------------------------




#---------------------------------------
# preamble
#---------------------------------------

rm(list=ls())
source(here::here("0-config.R"))

source(here::here("analysis/0-base-quartileTMLE_functions.R"))
source(here::here("analysis/1-gam-functions.R"))

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





#null dataframe
h1adj.res = NULL 

Wvars_2_3<-c("sex","birthord", "momage","momheight","momedu", 
             "hfiacat", "Nlt18", "Ncomp", "watmin", 
             "floor", "walls", "elec", "asset_wardrobe", "asset_table", "asset_chair", "asset_clock", "asset_khat", 
             "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
             "asset_bike", "asset_moto", "asset_sewmach", "asset_mobile", 
             "n_cattle", "n_goat", "n_chicken", "monsoon_ht2", "monsoon_ht3", "ageday_ht2", 
             "ageday_ht3", "tr", "cesd_sum_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", "diar7d_t2", "diar7d_t3", 
             "life_viol_any_t3", "lenhei_med_t2", "weight_med_t2", "anthro_days_btwn_t2_t3")


h1adj.res <- tmle_quart(dat=d, 
                        Y="delta_laz_t2_t3", 
                        W=Wvars_2_3, 
                        A="delta_TS", 
                        id="block",
                        Alevels=c("Q1","Q2","Q3","Q4"),
                        outputdf = h1adj.res,
                        family="gaussian", 
                        SLlibrary="SL.gam")

A="TS_t2"
dat=d
Acuts=as.numeric(summary( subset(dat, select=A)[,1])[c(2,3,5)])
cat("Cutpoints: ", Acuts, "\nNumbers per category:\n")
a<-subset(dat, select=A)
a[,1]<-findInterval(a[,1], Acuts)
d$quartiles<-factor(a[,1])
#if(!is.null(Alevels)){levels(a[,1])<-Alevels[as.numeric(levels(a[,1]))+1]}

res_adj <- fit_RE_gam(d=d, X=A, Y="laz_t2",  W=Wvars_2_3)
res_adj_cat <- fit_RE_gam(d=d, X="quartiles", Y="laz_t2",  W=Wvars_2_3)

library(emmeans)
# res = emmeans(res_adj$fit, "X", by=colnames(res_adj$dat)[-c(1:4)], data=res_adj$dat)
# head(res)
preddf = res_adj_cat$dat %>% mutate(dummy=0)
res = emmeans(res_adj_cat$fit, "X", data=preddf)
res
res_adj_cat$dat %>% group_by(X) %>%
  summarise(mean(Y, na.rm=T))

#res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])

simul_plot <- gam_simul_CI(res_adj$fit, res_adj$dat, xlab=A, ylab="laz_t2", title="")
simul_plot$p



m=res_adj$fit
newdata=res_adj$dat
nreps=10000
xlab=A
ylab="laz_t2"
title=""

gam_simul_CI <- function(m,newdata,nreps=10000, xlab="", ylab="", title=""){
  set.seed(12345)
  require(mgcv)
  require(dplyr)
  
  newdata <- newdata %>% mutate(dummy=0)

  #get 10th and 90th percentile of distribution to truncate plot

  
  
  
    trunc_points <- quantile(newdata$X, probs=c(0.05,0.95))
  Acuts=as.numeric(summary( subset(newdata, select=X)[,1])[c(2,3,5)])
  newdata$quartiles<-findInterval(newdata$X, Acuts)
  
 #get median of each quartiles 
  midpoints=  newdata %>% group_by(quartiles) %>%
    summarise(midpoints=median(X))
  midpoints$quartiles <- factor(midpoints$quartiles)
  
  EY <- as.data.frame(res)
  colnames(EY)[1] <- "quartiles"
  EY <- left_join(EY, midpoints, by=c("quartiles"))
  
  Wvars <- colnames(newdata)[!(colnames(newdata) %in% c("Y","X" ,"id" ,"dummy"))]

  #set covariates to the median/mode
  for(i in Wvars){
    if(class(newdata[,i])=="character"|class(newdata[,i])=="factor"){
      newdata[,i] <- Mode(newdata[,i])
    }else{
      newdata[,i] <- median(newdata[,i])
    }
  }
  
  newdata <- newdata[order(newdata$X),]
  
  Vb <- vcov(m,unconditional = TRUE)
  pred <- predict(m, newdata, se.fit = TRUE)
  fit <- pred$fit
  se.fit <- pred$se.fit
  BUdiff <- MASS::mvrnorm(n=nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg <- predict(m, newdata, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  pred <- data.frame(newdata,fit=pred$fit,se.fit=pred$se.fit)
  pred <- mutate(pred,
                 uprP = fit + (2 * se.fit),
                 lwrP = fit - (2 * se.fit),
                 uprS = fit + (crit * se.fit),
                 lwrS = fit - (crit * se.fit)
  )
  
  pred <- pred %>% arrange(X)
  
  pred_trunc <- pred %>% filter(X>trunc_points[1] & X<trunc_points[2])
  
  offset = mean(pred$fit) - mean(pred$Y)
  #offset = 0
  
  p <- ggplot(pred) + geom_ribbon(aes(x=X, ymin=lwrS-offset, ymax=uprS-offset), alpha=0.2) + 
    geom_path(aes(x=X, y=fit - offset), color="black") +
    geom_point(data=EY, aes(x=midpoints, y=emmean), size=2.5) +
    geom_linerange(data=EY, aes(x=midpoints, ymin=lower.CL, ymax=upper.CL), size=1) +
    geom_vline(aes(xintercept=Acuts[1]), linetype="dashed", color="black") +
    geom_vline(aes(xintercept=Acuts[2]), linetype="dashed", color="black") +
    geom_vline(aes(xintercept=Acuts[3]), linetype="dashed", color="black") +
    coord_cartesian(xlim=trunc_points, 
                    ylim=c(min(pred_trunc$lwrS-offset),max(pred_trunc$uprS-offset))) +
    xlab(xlab) + ylab(ylab) +
    ggtitle(title)
  
  
  return(list(p=p, pred=pred))
}


