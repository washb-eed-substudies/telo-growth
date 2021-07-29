
rm(list=ls())

source(here::here("0-config.R"))

d<-haven::read_dta(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.dta"))
colnames(d)

plotdf <- bind_rows(
  data.frame(TS=d$raw_delta_TS, group="Uncorrected TS velocity"),
  data.frame(TS=d$delta_TS, group="RTM-corrected TS velocity")
  ) 
table(plotdf$group, is.na(plotdf$TS))

plotdf <- plotdf %>% filter(!is.na(TS)) %>% group_by(group) %>% mutate(Xmedian=median(TS)) %>% ungroup()


p <- ggplot(data=plotdf, aes(x=TS,group=group, fill=group)) +
  geom_density(aes(y=..density.. , alpha=0.2),color=NA) +
  #geom_vline(aes(xintercept = Xmedian)) +
  #geom_text(aes(x = Xmedian, y=lab_pos, label=Xmedian2), hjust=-0.5) +
  scale_colour_manual(values=tableau10, drop=TRUE) + 
  scale_fill_manual(values=tableau10, drop=TRUE) + 
  labs(x="Change in T/S", y="Density") +
  theme_minimal(base_size=16) +
  theme(legend.position = "bottom") +
  guides(color = "none",alpha = "none", fill=guide_legend(title=""))

p