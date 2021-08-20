




rm(list=ls())
#renv::deactivate()
source(here::here("0-config.R"))

dtelo_unadj_res <- readRDS(here("results/gam_results/unadjusted/dtelo_res.RDS")) %>% mutate(H=3)
dtelo_adj_res <- readRDS(here("results/gam_results/adjusted/dtelo_adj_res.RDS")) %>% mutate(H=3)
raw_dtelo_unadj_res <- readRDS(here("results/gam_results/unadjusted/raw_dtelo_res.RDS")) %>% mutate(H=3)
raw_dtelo_adj_res <- readRDS(here("results/gam_results/adjusted/raw_dtelo_adj_res.RDS")) %>% mutate(H=3)

d <- bind_rows(dtelo_unadj_res, raw_dtelo_unadj_res) %>% rename(A=X)
head(d)

unique(d$Y)
unique(d$A)

d$Y <- as.character(d$Y)
d$Y <- gsub("whz", "wlz", d$Y)
d$Y <- factor(d$Y, levels=c("laz_t3", "wlz_t3", "waz_t3", "hcz_t3",
                            "delta_laz_t2_t3", "delta_wlz_t2_t3", "delta_waz_t2_t3", "delta_hcz_t2_t3",    
                            "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"))
d <- d %>% filter(A != "delta_TS_Z")

d <- d %>% arrange(Y) %>%
  mutate(
    Ylab = case_when(
      Y=="delta_laz_t2_t3" ~ "Change in LAZ\nbetween years 1 and 2",
      Y=="delta_wlz_t2_t3" ~ "Change in WLZ\nbetween years 1 and 2",
      Y=="delta_waz_t2_t3" ~ "Change in WAZ\nbetween years 1 and 2",
      Y=="delta_hcz_t2_t3" ~ "Change in HCZ\nbetween years 1 and 2",
      Y=="len_velocity_t2_t3" ~ "Length velocity (cm/month)\nbetween years 1 and 2",
      Y=="wei_velocity_t2_t3" ~ "Weight velocity (kg/month)\nbetween years 1 and 2",
      Y=="hc_velocity_t2_t3" ~ "Head circumference velocity (cm/month)\nbetween years 1 and 2",
      Y=="laz_t3" ~ "LAZ at year 2",
      Y=="wlz_t3" ~ "WLZ at year 2",
      Y=="waz_t3" ~ "WAZ at year 2",
      Y=="hcz_t3" ~ "HCZ at year 2"),
    Ylab = factor(Ylab, levels=unique(Ylab)),
    Alab = case_when(
      A=="delta_TS" ~ "RTM corrected change",
      A=="raw_delta_TS" ~ "Uncorrected change"),
    anthro = case_when(
      grepl("delta_laz", Y) ~ "??? LAZ",
      grepl("delta_waz", Y) ~ "??? WAZ",
      grepl("delta_wlz", Y) ~ "??? WLZ",
      grepl("delta_hcz", Y) ~ "??? HCZ",
      grepl("laz_t3",Y) ~ "LAZ",
      grepl("waz_t3",Y) ~ "WAZ",
      grepl("wlz_t3",Y) ~ "WLZ",
      grepl("hcz_t3",Y) ~ "HCZ",
      grepl("len_",Y) ~ "LENGTH\n",
      grepl("wei_",Y) ~ "WEIGHT\n",
      grepl("hc_",Y) ~ "HEAD\nCIRCUMFERENCE"
    ),
    anthro=factor(anthro,
                  levels = c("LAZ","WAZ","WLZ","HCZ",
                             "??? LAZ", "??? WAZ", "??? WLZ", "??? HCZ",
                             "LENGTH\n","WEIGHT\n","HEAD\nCIRCUMFERENCE"))) %>%
  arrange(anthro) %>%
  mutate(Ylab = factor(Ylab, levels=unique(Ylab)))

head(d)
tail(d)

yrange <- c(-0.4,0.5)
ylabel="Unadjusted difference in mean anthropometry Z score\nbetween 10th and 90th percentile of telomere measure\n"

d <- d %>% mutate(group = case_when(
                    grepl("delta_", Y) ~ "Change in growth \nbetween Years 1 and 2\n",
                    grepl("velocity", Y) ~ "Growth velocity \n between Years 1 and 2\n",
                    grepl("z_t3",Y) ~ "Growth at Year 2\n"),
                  group = factor(group, levels = c("Growth at Year 2\n", "Change in growth \nbetween Years 1 and 2\n", 
                           "Growth velocity \n between Years 1 and 2\n")))
p <- ggplot(d, aes(x=anthro, y=point.diff)) + 
  geom_pointrange(aes(ymin=lb.diff , ymax=ub.diff, color=anthro, group=Alab, shape=Alab),
                  position = position_dodge(width = 0.4),
                  size = 1) +
  facet_grid(~group, scales = "free_x") +
  labs(y = ylabel, x =  "\nAnthropometry measure") +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=yrange) +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"),
        legend.position = "right")+
  guides(color = "none", shape=guide_legend(title="Change in Telomere\nLength"))
p




ggsave(p, file = here("figures/telo-growth-gam-differences-dtelo-unadj.tiff"), height=6, width=14)



d <- bind_rows(dtelo_adj_res, raw_dtelo_adj_res) %>% rename(A=X)
head(d)

unique(d$Y)
unique(d$A)

d$Y <- as.character(d$Y)
d$Y <- gsub("whz", "wlz", d$Y)
d$Y <- factor(d$Y, levels=c("laz_t3", "wlz_t3", "waz_t3", "hcz_t3",
                            "delta_laz_t2_t3", "delta_wlz_t2_t3", "delta_waz_t2_t3", "delta_hcz_t2_t3",    
                            "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"))
d <- d %>% filter(A != "delta_TS_Z")

d <- d %>% arrange(Y) %>%
  mutate(
    Ylab = case_when(
      Y=="delta_laz_t2_t3" ~ "Change in LAZ\nbetween years 1 and 2",
      Y=="delta_wlz_t2_t3" ~ "Change in WLZ\nbetween years 1 and 2",
      Y=="delta_waz_t2_t3" ~ "Change in WAZ\nbetween years 1 and 2",
      Y=="delta_hcz_t2_t3" ~ "Change in HCZ\nbetween years 1 and 2",
      Y=="len_velocity_t2_t3" ~ "Length velocity (cm/month)\nbetween years 1 and 2",
      Y=="wei_velocity_t2_t3" ~ "Weight velocity (kg/month)\nbetween years 1 and 2",
      Y=="hc_velocity_t2_t3" ~ "Head circumference velocity (cm/month)\nbetween years 1 and 2",
      Y=="laz_t3" ~ "LAZ at year 2",
      Y=="wlz_t3" ~ "WLZ at year 2",
      Y=="waz_t3" ~ "WAZ at year 2",
      Y=="hcz_t3" ~ "HCZ at year 2"),
    Ylab = factor(Ylab, levels=unique(Ylab)),
    Alab = case_when(
      A=="delta_TS" ~ "RTM corrected change",
      A=="raw_delta_TS" ~ "Uncorrected change"),
    anthro = case_when(
      grepl("delta_laz", Y) ~ "??? LAZ",
      grepl("delta_waz", Y) ~ "??? WAZ",
      grepl("delta_wlz", Y) ~ "??? WLZ",
      grepl("delta_hcz", Y) ~ "??? HCZ",
      grepl("laz_t3",Y) ~ "LAZ",
      grepl("waz_t3",Y) ~ "WAZ",
      grepl("wlz_t3",Y) ~ "WLZ",
      grepl("hcz_t3",Y) ~ "HCZ",
      grepl("len_",Y) ~ "LENGTH\n",
      grepl("wei_",Y) ~ "WEIGHT\n",
      grepl("hc_",Y) ~ "HEAD\nCIRCUMFERENCE"
    ),
    anthro=factor(anthro,
                  levels = c("LAZ","WAZ","WLZ","HCZ",
                             "??? LAZ", "??? WAZ", "??? WLZ", "??? HCZ",
                             "LENGTH\n","WEIGHT\n","HEAD\nCIRCUMFERENCE"))) %>%
  arrange(anthro) %>%
  mutate(Ylab = factor(Ylab, levels=unique(Ylab)))

head(d)
tail(d)

yrange <- c(-0.4,0.5)
ylabel="Adjusted difference in mean anthropometry Z score\nbetween 10th and 90th percentile of telomere measure\n"

d <- d %>% mutate(group = case_when(
  grepl("delta_", Y) ~ "Change in growth \nbetween Years 1 and 2\n",
  grepl("velocity", Y) ~ "Growth velocity \n between Years 1 and 2\n",
  grepl("z_t3",Y) ~ "Growth at Year 2\n"),
  group = factor(group, levels = c("Growth at Year 2\n", "Change in growth \nbetween Years 1 and 2\n", 
                                   "Growth velocity \n between Years 1 and 2\n")))
p <- ggplot(d, aes(x=anthro, y=point.diff)) + 
  geom_pointrange(aes(ymin=lb.diff , ymax=ub.diff, color=anthro, group=Alab, shape=Alab),
                  position = position_dodge(width = 0.4),
                  size = 1) +
  facet_grid(~group, scales = "free_x") +
  labs(y = ylabel, x =  "\nAnthropometry measure") +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=yrange) +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"),
        legend.position = "right")+
  guides(color = "none", shape=guide_legend(title="Change in Telomere\nLength"))
p




ggsave(p, file = here("figures/telo-growth-gam-differences-dtelo-adj.tiff"), height=6, width=14)

