#### Adjust all pvalues with BH procedure ####
rm(list=ls())

source(here::here("0-config.R"))

# load all results
H1_res <- readRDS(here('results/gam_results/unadjusted/H1_res.RDS'))
H2_res <- readRDS(here('results/gam_results/unadjusted/H2_res.RDS'))
H3_res <- readRDS(here('results/gam_results/unadjusted/H3_res.RDS'))

H1_adj_res <- readRDS(here('results/gam_results/adjusted/H1_adj_res.RDS'))
H2_adj_res <- readRDS(here('results/gam_results/adjusted/H2_adj_res.RDS'))
H3_adj_res <- readRDS(here('results/gam_results/adjusted/H3_adj_res_tmleR.RDS'))

H1_res$H = 1
H2_res$H = 2
H3_res$H = 3

H1_adj_res$H = 1
H2_adj_res$H = 2
H3_adj_res$H = 3

# make splits between adjusted and unadjusted results and Z scores v. original telo length
full_res_Z <- rbind(filter(H1_res, X=="TS_t2_Z"), filter(H2_res, X=="TS_t3_Z"), 
                    filter(H3_res, X=="delta_TS_Z"))
full_res_noZ <- rbind(filter(H1_res, X=="TS_t2"), filter(H2_res, X=="TS_t3"), 
                             filter(H3_res, X=="delta_TS"))
full_adj_res_Z <- rbind(filter(H1_adj_res, X=="TS_t2_Z"), filter(H2_adj_res, X=="TS_t3_Z"), 
                    filter(H3_adj_res, X=="delta_TS_Z"))
full_adj_res_noZ <- rbind(filter(H1_adj_res, X=="TS_t2"), filter(H2_adj_res, X=="TS_t3"), 
                      filter(H3_adj_res, X=="delta_TS"))

# add BH calculations for outcomes
full_res_Z <- full_res_Z %>% group_by(Y) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

full_res_noZ <- full_res_noZ %>% group_by(Y) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

full_adj_res_Z <- full_adj_res_Z %>% group_by(Y) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

full_adj_res_noZ <- full_adj_res_noZ %>% group_by(Y) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

# combine results for saving
full_res <- rbind(full_res_Z, full_res_noZ)
full_adj_res <- rbind(full_adj_res_Z, full_adj_res_noZ)


saveRDS(full_res %>% filter(H==1) %>% select(-H), here("results/gam_results/unadjusted/H1_res_BH.RDS"))
saveRDS(full_res %>% filter(H==2) %>% select(-H), here("results/gam_results/unadjusted/H2_res_BH.RDS"))
saveRDS(full_res %>% filter(H==3) %>% select(-H), here("results/gam_results/unadjusted/H3_res_BH.RDS"))

saveRDS(full_adj_res %>% filter(H==1) %>% select(-H), here("results/gam_results/adjusted/H1_adj_res_BH.RDS"))
saveRDS(full_adj_res %>% filter(H==2) %>% select(-H), here("results/gam_results/adjusted/H2_adj_res_BH.RDS"))
saveRDS(full_adj_res %>% filter(H==3) %>% select(-H), here("results/gam_results/adjusted/H3_adj_res_BH.RDS"))
