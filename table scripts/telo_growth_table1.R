rm(list=ls())
library("xtable")
source(here::here("0-config.R"))

d <- read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.csv"))

# calculating quantiles/n for each variable 
female <- round(length(d$sex[d$sex == "female"])/sum(!is.na(d$sex))* 100) #percentage
aget2 <- round(quantile(d$agemth_ht2,na.rm = T), 1)
aget3 <- round(quantile(d$agemth_ht3,na.rm = T), 1)
aged23 <- round(quantile(d$agemth_ht3-d$agemth_ht2, na.rm=T), 1)
telo1med <- round(quantile(d$TS_t2, na.rm=TRUE),2)
telo1bpmed <- round(quantile(d$ts_t2_bp, na.rm=TRUE))
telo2med <- round(quantile(d$TS_t3, na.rm=TRUE),2)
telo2bpmed <- round(quantile(d$ts_t3_bp, na.rm=TRUE))
deltatsmed <- round(quantile(d$delta_TS, na.rm=TRUE),2)
deltatsbpmed <- round(quantile(d$delta_ts_bp, na.rm=TRUE))
length_m3 <- format(round(quantile(d$laz_t1, na.rm=TRUE),2), nsmall=2)
weight_age_m3 <- format(round(quantile(d$waz_t1, na.rm=TRUE),2), nsmall=2)
weight_length_m3 <- format(round(quantile(d$whz_t1, na.rm=TRUE),2), nsmall=2)
headc_age_m3 <- format(round(quantile(d$hcz_t1, na.rm=TRUE),2), nsmall=2)
length_y1 <- format(round(quantile(d$laz_t2, na.rm=TRUE),2), nsmall=2)
weight_age_y1 <- format(round(quantile(d$waz_t2, na.rm=TRUE),2), nsmall=2)
weight_length_y1 <- format(round(quantile(d$whz_t2, na.rm=TRUE),2), nsmall=2)
headc_age_y1 <- format(round(quantile(d$hcz_t2, na.rm=TRUE),2), nsmall=2)
length_y2 <- format(round(quantile(d$laz_t3, na.rm=TRUE),2), nsmall=2)
weight_age_y2 <- format(round(quantile(d$waz_t3, na.rm=TRUE),2), nsmall=2)
weight_length_y2 <- format(round(quantile(d$whz_t3, na.rm=TRUE),2), nsmall=2)
headc_age_y2 <- format(round(quantile(d$hcz_t3, na.rm=TRUE),2), nsmall=2)
d_y1 <- round(mean(d$diar7d_t2, na.rm=TRUE) * 100) #percentage
d_y2 <- round(mean(d$diar7d_t3, na.rm=TRUE)* 100) #percentage
agem <- round(quantile(d$momage, na.rm=TRUE))
heightm <- round(quantile(d$momheight, na.rm=TRUE),1)
edumom <- round(quantile(d$momeduy, na.rm=TRUE)) 
CES_D1 <- round(quantile(d$cesd_sum_t2, na.rm=TRUE))
CES_D2 <- round(quantile(d$cesd_sum_ee_t3, na.rm=TRUE))
PPS <- round(quantile(d$pss_sum_mom_t3, na.rm=TRUE))
viol <- round(mean(d$life_viol_any_t3, na.rm=TRUE) * 100) #percentage

#create table
tbl <- data.table("1"=character(), "2"=character(), "3"=character(), "4"=numeric())
tbl <- rbind(tbl, list("Child", " ", "Female (%)", paste(length(d$sex[d$sex == "female"]), " (", female, "%)", sep="")))
tbl <- rbind(tbl, list(" ", " ", "Age (months) at Year 1", paste(aget2[3]," (", aget2[2], ", ", aget2[4], ")", sep="")))
tbl <- rbind(tbl, list(" ", " ", "Age (months) at Year 2", paste(aget3[3]," (", aget3[2], ", ", aget3[4], ")", sep="")))
tbl <- rbind(tbl, list("", " ", "Months between telomere length measurements at Year 1 and Year 2", paste(aged23[3]," (", aged23[2], ", ", aged23[4], ")", sep="")))
tbl <- rbind(tbl, list(" ", "Telomere length at Year 1", "T/S Ratio", paste(telo1med[3]," (", telo1med[2], ", ", telo1med[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Telomere length at Year 2", "3" = "T/S Ratio", 
                       "4" = paste(telo2med[3]," (", telo2med[2], ", ", telo2med[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Change in telomere length between Year 1 and Year 2", "3" = "T/S Ratio", 
                       "4" = paste(deltatsmed[3]," (", deltatsmed[2], ", ", deltatsmed[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Anthropometry (age 3 months, Month 3)", "3" = "Length-for-age Z score", 
                       "4" = paste(length_m3[3], " (", length_m3[2], ", ", length_m3[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Weight-for-age Z score", 
                       "4" = paste(weight_age_m3[3], " (", weight_age_m3[2], ", ", weight_age_m3[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Weight-for-length Z score", 
                       "4" = paste(weight_length_m3[3], " (", weight_length_m3[2], ", ", weight_length_m3[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Head circumference-for-age Z score", 
                       "4" = paste(headc_age_m3[3], " (", headc_age_m3[2], ", ", headc_age_m3[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Anthropometry (age 14 months, Year 1)", "3" = "Length-for-age Z score", 
                       "4" = paste(length_y1[3], " (", length_y1[2], ", ", length_y1[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Weight-for-age Z score", 
                       "4" = paste(weight_age_y1[3], " (", weight_age_y1[2], ", ", weight_age_y1[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Weight-for-length Z score", 
                       "4" = paste(weight_length_y1[3], " (", weight_length_y1[2], ", ", weight_length_y1[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Head circumference-for-age Z score", 
                       "4" = paste(headc_age_y1[3], " (", headc_age_y1[2], ", ", headc_age_y1[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Anthropometry (age 28 months, Year 2)", "3" = "Length-for-age Z score", 
                       "4" = paste(length_y2[3], " (", length_y2[2], ", ", length_y2[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Weight-for-age Z score", 
                       "4" = paste(weight_age_y2[3], " (", weight_age_y2[2], ", ", weight_age_y2[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Weight-for-length Z score", 
                       "4" = paste(weight_length_y2[3], " (", weight_length_y2[2], ", ", weight_length_y2[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = " ", "3" = "Head circumference-for-age Z score", 
                       "4" = paste(headc_age_y2[3], " (", headc_age_y2[2], ", ", headc_age_y2[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Diarrhea (age 14 months, Year 1)", "3" = "Caregiver-reported 7-day recall (%)", 
                       "4" = paste(sum(d$diar7d_t2, na.rm=TRUE), " (", d_y1, "%)", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Diarrhea (age 28 months, Year 2)", "3" = "Caregiver-reported 7-day recall (%)", 
                       "4" = paste(sum(d$diar7d_t3, na.rm=TRUE), " (", d_y2, "%)", sep="")))
tbl <- rbind(tbl, list("1" = "Mother", "2" = " ", "3" = "Age (years)", 
                       "4" = paste(agem[3], " (", agem[2], ", ", agem[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Anthropometry at enrollment", "3" = "Height (cm)", 
                       "4" = paste(heightm[3], " (", heightm[2], ", ", heightm[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Education", "3" = "Schooling completed (years)", 
                       "4" = paste(edumom[3], " (", edumom[2], ", ", edumom[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Depression at Year 1", "3" = "CESD-R score**", 
                       "4" = paste(CES_D1[3], " (", CES_D1[2], ", ", CES_D1[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Depression at Year 2", "3" = "CESD-R score**", 
                       "4" = paste(CES_D2[3], " (", CES_D2[2], ", ", CES_D2[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Perceived stress at Year 2", "3" = "Perceived Stress Scale score", 
                       "4" = paste(PPS[3], " (", PPS[2], ", ", PPS[4], ")", sep="")))
tbl <- rbind(tbl, list("1" = " ", "2" = "Physical, sexual, or emotional intimate partner violence", "3" = "Any lifetime exposure: number of women (%)", 
                       "4" = paste(sum(d$life_viol_any_t3, na.rm=TRUE), " (", viol, "%)", sep="")))

library(flextable)
library(officer)

tblflex <- flextable(tbl, col_keys=names(tbl))
tblflex <- set_header_labels(tblflex,
                              values = list("1" = "", "2" = "", "3" = "", "4" = "n (%) or median (IQR)"))
tblflex <- hline_top(tblflex, part="header", border=fp_border(color="black", width = 1))
tblflex <- hline_bottom(tblflex, part="all", border=fp_border(color="black", width = 1))
tblflex <- autofit(tblflex, part = "all")
tblflex <- align(tblflex, j = c(1, 2, 3), align = "left", part="all")
tblflex <- align(tblflex, j = 4, align = "center", part="all")
tblflex <- fit_to_width(tblflex, max_width=8)
tblflex <- add_footer_row(tblflex, top=F, 
                          values = "*The unit for relative telomere length is the T/S ratio. Telomere length was measured by quantitative PCR (qPCR), a method that determines relative telomere length by measuring the factor by which each DNA sample differs from a reference DNA sample in its ratio of telomere repeat copy number (T) to single-copy gene copy number (S)", colwidths = 4)
tblflex <- add_footer_row(tblflex, top=F, 
                          values = "**CESD-R = Center for Epidemiologic Studies Depression Scale Revised", colwidths = 4)


# export table as csv
write.csv(tbl, file = here("tables/main/telo_growth_table1.csv"))
write("*The unit for relative telomere length is the T/S ratio. Telomere length was measured by quantitative PCR (qPCR), a method that determines relative telomere length by measuring the factor by which each DNA sample differs from a reference DNA sample in its ratio of telomere repeat copy number (T) to single-copy gene copy number (S)",file=here("tables/main/telo_growth_table1.csv"),append=TRUE)
write("**CESD-R = Center for Epidemiologic Studies Depression Scale Revised",file=here("tables/main/telo_growth_table1.csv"),append=TRUE)

print(xtable(tbl), type="html", file=here("tables/main/telo_growth_table1.html"))
write("*The unit for relative telomere length is the T/S ratio. Telomere length was measured by quantitative PCR (qPCR), a method that determines relative telomere length by measuring the factor by which each DNA sample differs from a reference DNA sample in its ratio of telomere repeat copy number (T) to single-copy gene copy number (S)",file=here("tables/main/telo_growth_table1.html"),append=TRUE)
write("**CESD-R = Center for Epidemiologic Studies Depression Scale Revised",file=here("tables/main/telo_growth_table1.html"),append=TRUE)

sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)
save_as_docx("Table 1: Characteristics of Participants" = tblflex, path='C:/Users/Sophia/Documents/WASH/WASH Telomeres and Growth/telo-growth-enrollment.docx', 
             pr_section = sect_properties) 
