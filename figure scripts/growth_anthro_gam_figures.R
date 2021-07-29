




rm(list=ls())
#renv::deactivate()
source(here::here("0-config.R"))

H1_adj_res <- readRDS(here("results/gam_results/adjusted/telot2_adj_res.RDS")) %>% mutate(H=1)
H2_adj_res <- readRDS(here("results/gam_results/adjusted/telot3_adj_res.RDS")) %>% mutate(H=2)
H3_adj_res <- readRDS(here("results/gam_results/adjusted/dtelo_adj_res.RDS")) %>% mutate(H=3)

d <- bind_rows(H1_adj_res, H2_adj_res, H3_adj_res) %>% rename(A=X)
head(d)

unique(d$Y)
unique(d$A)

d$Y <- as.character(d$Y)
d$Y <- gsub("whz", "wlz", d$Y)
d$Y <- factor(d$Y, levels=c("laz_t2", "wlz_t2", "waz_t2", "hcz_t2",
                            "laz_t3", "wlz_t3", "waz_t3", "hcz_t3",
                            "delta_laz_t2_t3", "delta_wlz_t2_t3", "delta_waz_t2_t3", "delta_hcz_t2_t3",    
                            "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"))

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
      Y=="laz_t2" ~ "LAZ at year 1",
      Y=="wlz_t2" ~ "WLZ at year 1",
      Y=="waz_t2" ~ "WAZ at year 1",
      Y=="hcz_t2" ~ "HCZ at year 1",
      Y=="laz_t3" ~ "LAZ at year 2",
      Y=="wlz_t3" ~ "WLZ at year 2",
      Y=="waz_t3" ~ "WAZ at year 2",
      Y=="hcz_t3" ~ "HCZ at year 2"),
    Ylab = factor(Ylab, levels=unique(Ylab)),
    Alab = case_when(
      A=="delta_TS" ~ "Change in telomere length between Years 1 and 2 (T/S ratio)",
      A=="TS_t2" ~ "Telomere length at Year 1 (T/S ratio)",
      A=="TS_t3" ~ "Telomere length at Year 2 (T/S ratio)"),
    anthro = case_when(
      grepl("laz",Y) ~ "LAZ",
      grepl("waz",Y) ~ "WAZ",
      grepl("wlz",Y) ~ "WLZ",
      grepl("hcz",Y) ~ "HCZ",
      grepl("len",Y) ~ "LEN",
      grepl("wei",Y) ~ "WEIGHT",
      grepl("hc_",Y) ~ "HCIR"
    ),
    anthro=factor(anthro,
                  levels = c("LAZ","WAZ","WLZ","HCZ",
                             "LEN","WEIGHT","HCIR"))) %>%
  arrange(anthro) %>%
  mutate(Ylab = factor(Ylab, levels=unique(Ylab)))

head(d)


#' some notes from our meeting:
#'   -Figures 2 and 3 collapsed into 1 figure
#' -Supplementary Figures updated (no density plots)
#' -1 new supp figure: density plot comparing before and after regression to the mean.
#' Thanks!
#' 


yrange <- c(-0.4,0.5)
ylabel="Adjusted difference in mean anthropometry Z score\nbetween 25th and 75th percentile of telomere measure"
plotdf <- d %>% filter(H!=3, grepl("laz",Y)| grepl("waz",Y)| grepl("wlz",Y)| grepl("hcz",Y), !grepl("_Z",A))
plotdf <- plotdf %>% mutate(
  group = case_when(
    grepl("delta_",Y) & grepl("_t2",A) ~ "Change in telomere length\nbetween Years 1 and 2 and growth\n",
    grepl("_t2",A) ~ "Telomere length at\nYear 1 and growth\n",
    grepl("_t3",Y) & grepl("_t3",A) ~ "Telomere length at\nYear 2 and growth\n"
  ),
  group=factor(group, levels=c("Telomere length at\nYear 1 and growth\n",
                               "Telomere length at\nYear 2 and growth\n",
                               "Change in telomere length\nbetween Years 1 and 2 and growth\n")),
  Y_time=ifelse(grepl("_t3",Y),"Year 2","Year 1")
)


p <- ggplot(plotdf, aes(x=anthro, y=point.diff)) + 
  geom_pointrange(aes(ymin=lb.diff , ymax=ub.diff, color=anthro, group=Y_time, shape=Y_time),
                  position = position_dodge(width = 0.4),
                  size = 1) +
  #geom_text(aes(label=ref), position = position_nudge(y = (abs(yrange[1])+abs(yrange[2]))/10)) +
  facet_grid(~group) +
  labs(y = ylabel, x =  "Anthropometry measure") +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=yrange) +
  scale_shape_manual(values=c(21,16)) +
  scale_colour_manual(values=tableau10) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        panel.spacing = unit(0, "lines"),
        legend.position = "right")+
  guides(color = "none", shape=guide_legend(title="Anthropometry\nmeasurement\ntime"))
p




ggsave(p, file = here("figures/telo-growth-gam-differences-fig1.tiff"), height=6, width=14)
