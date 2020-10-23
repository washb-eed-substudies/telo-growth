rm(list=ls())

source(here::here("0-config.R"))
load(here("/results/telo_growth_results.Rdata"))

unadj <- rbind(h1unadj.res, h2unadj.res, h3unadj.res, h4unadj.res, 
               h5unadj.res, h6unadj.res, h7unadj.res, h8unadj.res)


adj <- rbind(h1adj.res, h2adj.res, h3adj.res, h4adj.res, 
             h5adj.res, h6adj.res, h7adj.res, h8adj.res)


# separate BH for each outcome in either unadjusted or adjusted analysis
# group by outcome and adjust BH
BH_unadj <- unadj %>% group_by(Y) %>% summarize(X=A, level=level, ATE=ATE, var=var, CI1=CI1, CI2=CI2, 
                                                Pval=Pval, compN=compN, refN=refN, meanLevel=meanLevel,
                                                meanN=meanN, meanY=meanY, mean.sd=mean.sd, mean.se=mean.se, 
                                                mean.CI1=mean.CI1, mean.CI2=mean.CI2, cutpoints=cutpoints,
                                                BH.Pval=p.adjust(Pval, method="BH"))
BH_adj <- adj %>% group_by(Y) %>% summarize(X=A, level=level, ATE=ATE, var=var, CI1=CI1, CI2=CI2, 
                                            Pval=Pval, compN=compN, refN=refN, meanLevel=meanLevel,
                                            meanN=meanN, meanY=meanY, mean.sd=mean.sd, mean.se=mean.se, 
                                            mean.CI1=mean.CI1, mean.CI2=mean.CI2, cutpoints=cutpoints,
                                            BH.Pval=p.adjust(Pval, method="BH"))


delta_unadj <- filter(BH_unadj, X=="delta_TS")
h1unadj <- delta_unadj[grep("delta", delta_unadj$Y),]
h2unadj <- delta_unadj[grep("velocity", delta_unadj$Y),]
h3unadj <- delta_unadj[grep("z_t3", delta_unadj$Y),]

t3_unadj <- filter(BH_unadj, X=="TS_t3")
h5unadj <- t3_unadj[grep("z_t3", t3_unadj$Y),]

t2_unadj <- filter(BH_unadj, X=="TS_t2")
h4unadj <- t2_unadj[grep("z_t2", t2_unadj$Y),]
h6unadj <- t2_unadj[grep("z_t3", t2_unadj$Y),]
h7unadj <- t2_unadj[grep("velocity", t2_unadj$Y),]
h8unadj <- t2_unadj[grep("delta", t2_unadj$Y),]


delta_adj <- filter(BH_adj, X=="delta_TS")
h1adj <- delta_adj[grep("delta", delta_adj$Y),]
h2adj <- delta_adj[grep("velocity", delta_adj$Y),]
h3adj <- delta_adj[grep("z_t3", delta_adj$Y),]

t3_adj <- filter(BH_adj, X=="TS_t3")
h5adj <- t3_adj[grep("z_t3", t3_adj$Y),]

t2_adj <- filter(BH_adj, X=="TS_t2")
h4adj <- t2_adj[grep("z_t2", t2_adj$Y),]
h6adj <- t2_adj[grep("z_t3", t2_adj$Y),]
h7adj <- t2_adj[grep("velocity", t2_adj$Y),]
h8adj <- t2_adj[grep("delta", t2_adj$Y),]


save(h1unadj, h2unadj, h3unadj, h4unadj, h5unadj, h6unadj, h7unadj, h8unadj,
     h1adj, h2adj, h3adj, h4adj, h5adj, h6adj, h7adj, h8adj, 
     file=here("results/telo_growth_results_BH.Rdata"))
