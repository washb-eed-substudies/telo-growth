#example prescreening table:

rm(list=ls())
#renv::deactivate()
source(here::here("0-config.R"))

library(boxr)
box_auth()
d <- box_read_csv(839767614700)

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
X <- c("TS_t2")            
Y <- c("laz_t2") 

W<-c("sex","birthord", "momage", "momheight","momedu", 
             "hfiacat", "Nlt18", "Ncomp", "watmin", "floor", 
             "walls", "elec", "asset_wardrobe", "asset_table", 
             "asset_chair", "asset_clock", "asset_khat", 
             "asset_chouki", "asset_radio", "asset_tv", "asset_refrig",
             "asset_bike", "asset_moto", "asset_sewmach",
             "asset_mobile", "n_cattle", "n_goat", "n_chicken",
             "monsoon_ht2", "ageday_ht2", "tr", "cesd_sum_t2", "diar7d_t2", 
             "life_viol_any_t3")


forcedW = NULL
V = NULL
id = "clusterid"
family = "gaussian"
pval = 0.2

  if (!is.null(W)) {
    W <- subset(d, select = W)
  }
  Y <- subset(d, select = Y)
  colnames(Y) <- "Y"
  X <- subset(d, select = X)
  colnames(X) <- "X"
  id <- subset(d, select = id)
  colnames(id) <- "id"
  if (!is.null(V)) {
    Vvar <- subset(d, select = V)
    colnames(Vvar) <- "V"
  }else{
    Vvar <- data.frame(V = rep(1, nrow(d)))
  }
  if(!is.null(W)) {
    gamdat <- data.frame(Y, X, id, Vvar, W)
  }else{
    gamdat <- data.frame(Y, X, id, Vvar)
  }
  if(!is.null(W)) {
    if(sum(is.na(forcedW)) != 0) {
      colnamesW <- names(W)
    }
    else {
      if (is.null(forcedW)) {
        Wnames <- names(W)
        forcedW <- c(Wnames[Wnames == "tr" | grepl("age_", 
                                                   Wnames) | grepl("agedays_", Wnames) | 
                              grepl("ageday_", Wnames)])
      }
      cat("\nNon-prescreened covariates: ", paste(forcedW, 
                                                  sep = "", collapse = ", "), "\n")
      colnamesW <- names(W)[!(names(W) %in% forcedW)]
    }
    screenW <- subset(gamdat, select = colnamesW)
  }else {
    screenW <- NULL
  }

  
  
suppressWarnings(Wscreen <- washb_prescreen(Y = gamdat$Y, Ws = screenW, family = family, pval = pval, print = T))
    
Y = gamdat$Y
Ws = screenW
family = "gaussian"
pval = 0.2
print = TRUE
    
      require(lmtest)
      if (family[[1]] == "neg.binom") {
        require(MASS)
      }
      if (pval > 1 | pval < 0) {
        stop("P-value threshold not set between 0 and 1.")
      }
      Ws <- as.data.frame(Ws)
      dat <- data.frame(Ws, Y)
      dat <- dat[complete.cases(dat), ]
      nW <- ncol(Ws)
      LRp <- matrix(rep(NA, nW), nrow = nW, ncol = 1)
      rownames(LRp) <- names(Ws)
      colnames(LRp) <- "P-value"
      if (family[[1]] != "neg.binom") {
        for (i in 1:nW) {
          dat$W <- dat[, i]
          if (class(dat$W) == "factor" & dim(table(dat$W)) == 
              1) {
            fit1 <- fit0 <- glm(Y ~ 1, data = dat, family = family)
          }
          else {
            fit1 <- glm(Y ~ W, data = dat, family = family)
            fit0 <- glm(Y ~ 1, data = dat, family = family)
          }
          LRp[i] <- lrtest(fit1, fit0)[2, 5]
        }
      }
    
      
      
      p20 <- ifelse(LRp < pval, 1, 0)
        cat("\nLikelihood Ratio Test P-values:\n")
        print(round(LRp, 5))
          LRps <- matrix(LRp[p20 == 1, ], ncol = 1)
          rownames(LRps) <- names(Ws)[p20 == 1]
          colnames(LRps) <- "P-value"

          
  tab <- data.frame(Covariate=rownames(LRp), pval=as.vector(LRp), Selected=as.vector(ifelse(LRp<0.2,"Yes","No")))     
  tab <- tab %>% arrange(pval) %>% mutate(pval=as.character(round(pval,3)))
  tab$pval <- gsub("0.000", "<0.001", tab$pval)
  colnames(tab)[2] <-"P-value"       
  tab        
  
  tab$Covariate
  tab$Covariate <- c("Maternal height (cm)", "Household assets - working black/white or color television",
                     "Household assets - refrigerator", "Maternal education level (no education, primary secondary)",
                     "Food insecurity (4-level HFIAS categories)", "Household assets - wardrobe", "Household assets - khat",
                     "Housing materials - floor", "Household assets - watch or clock", "Household assets - electricity", 
                     "Household assets - chair or bench", "Household assets - table", "Number of individuals living in the compound",
                     "Household assets  - number of cows", "Household assets - mobile phone", "Household assets - number of chickens",
                     "Household assets - working radio", "Household assets - sewing machine", "Household assets - motorcycle",
                     "Housing materials - walls", "Maternal CESD-R Scale Score at Year 1", "Household assets - bicycle",
                     "Distance in minutes to primary drinking water source", "Household assets - number of goats", "Caregiver-reported diarrhea at Year 1",
                     "Child birth order", "Maternal age (years)", "Season of measurement at Year 1", "Household assets - chouki", 
                     "Child sex", "Number of children <18 in the household", "Maternal lifetime cumulative exposure to intimate partner violence")
  tab
  
  flextable(tab) %>% autofit() %>% save_as_image(path = here("tables/example-prescreening.png"))
          
          
   