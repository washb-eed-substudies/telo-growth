
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

ggplot(d, aes(x=raw_delta_TS, y=delta_TS), color="grey30") + geom_point(alpha=0.8) + geom_smooth(method="lm", se=F, color="grey30") +
  xlab("Uncorrected T/S") + ylab("Corrected change in T/S")

#scatter plots
test_res <- cor.test(d$TS_t2, d$delta_TS)

test_res$p.value
cor_label <- paste0("r = ",round(test_res$estimate,2), ", P-value < 0.001")
test_res$cor_label <- NA
test_res$cor_label[1] <- cor_label
p2 <- ggplot(d, aes(x=TS_t2, y=delta_TS), color="grey30") + geom_point(alpha=0.8) + geom_smooth(method="lm", se=F, color="grey30") +
  xlab("Baseline T/S") + ylab("Change in T/S") + geom_text(label=cor_label, x=2.35, y=1, size=5)

ggsave(p2, file = here("figures/telo-growth-RTM-scatter.tiff"), height=6, width=6)
