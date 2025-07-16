#!/usr/bin/env Rscript

#User: Noemi Nicole Piga

#Rationale: Compare the UK Biobank participants identified in my case definition against
#Lancet (Shrine et al. 2019) stage 1 and stage2 phenotype.

#Collate data:
#Stage 1 UKBB Lancet paper moderate-severe asthma participants:
#manually copy and paste of severe_asthma_cases+controls.csv from /rfs/TobinGroup/GWAtraits/FEV/AIRPROM/GWAS/severe_asthma/
#to /home/n/nnp5/PhD/PhD_project/Thesis_analysis/input/

#Stage 2  UKBB Lancet paper moderate-severe asthma participants:
##data file:
#awk '$NF == 1' \
#    /data/gen1/AIRPROM/assoc/severe_asthma/replication_noatopy/replication.sample | awk '{print $1}' \
#    > input/UKB_stage2_Lancet_cases_app648.txt

library(tidyverse)
library(dplyr)
library(data.table)
library(hablar)
library(VennDiagram)
#avoid log file for venn diagrams
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

bridge_app648_8389 <- fread("/data/gen1/UKBiobank/application_648/mapping_to_app8389.txt",header=T)
bridge_app648_8389$app8389 <- as.character(bridge_app648_8389$app8389)
bridge_app648_8389$app648 <- as.character(bridge_app648_8389$app648)
bridge_app648_56607 <- fread("/data/gen1/UKBiobank_500K/severe_asthma/data/bridge_app648_56607",sep=" ",header=T)
colnames(bridge_app648_56607) <- c("app648","app56607")
bridge_app648_56607$app56607 <- as.character(bridge_app648_56607$app56607)
bridge_app648_56607$app648 <- as.character(bridge_app648_56607$app648)
bridge_648_8389_56607 <- inner_join(bridge_app648_8389,bridge_app648_56607,by="app648")

#load input
#Lancet stage 1:
lancet_stage1 <- fread("input/severe_asthma_cases+controls.csv",header=T)
lancet_stage1_UKBB_cases <- lancet_stage1 %>% filter(ID_1 == "UKB_SEVEREasthma_cases") %>% rename(app8389=ID_2)
lancet_stage1_UKBB_cases_app56607 <- left_join(lancet_stage1_UKBB_cases,bridge_648_8389_56607,by=c("app8389")) %>%
                                     select(app8389,app648,app56607)

#lancet stage 2:
lancet_stage2 <- fread("input/UKB_stage2_Lancet_cases_app648.txt",header=F)
colnames(lancet_stage2)[1] <- "app648"
lancet_stage2$app648 <- as.character(lancet_stage2$app648)
lancet_stage2_UKBB_cases_app56607 <- left_join(lancet_stage2,bridge_648_8389_56607,by="app648")


#find lancet with no eid for app56607:
lancet_no_app55607_eid <- rbind(lancet_stage1_UKBB_cases_app56607,lancet_stage2_UKBB_cases_app56607) %>%
                          filter(is.na(app56607))

#Combine:
lancet_stage1_UKBB_cases_app56607 <- lancet_stage1_UKBB_cases_app56607 %>% select(app56607)
lancet_stage1_UKBB_cases_app56607$lancet_stage <- "lancet_stage1"
lancet_stage2_UKBB_cases_app56607 <-  lancet_stage2_UKBB_cases_app56607 %>% select(app56607)
lancet_stage2_UKBB_cases_app56607$lancet_stage <- "lancet_stage2"
lancet <- rbind(lancet_stage1_UKBB_cases_app56607,lancet_stage2_UKBB_cases_app56607) %>% rename(eid=app56607)
lancet$eid <- as.character(lancet$eid)
#Know how many participants have compared for stage1 and stage2:
#> dim(lancet %>% filter(lancet_stage == "lancet_stage1"))
#[1] 2992    2
#> dim(lancet %>% filter(lancet_stage == "lancet_stage2"))
#[1] 5411    2
#remove NAs:
lancet <- lancet %>% filter(!is.na(eid))

#Individuals for my severe asthma phenotype:
demo <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt",header=T,sep=" ")
demo$eid <- as.factor(demo$eid)


#lancet not present in app56607 master sample file (from which I created demographics table):
lancet_not_in_app56607_sample <- anti_join(lancet,demo,by="eid") %>% select(eid,lancet_stage)
#Know how many participants are not in the app56607 master sample file for stage1 and stage2:
#> dim(lancet_not_in_app56607_sample %>% filter(lancet_stage == "lancet_stage1"))
#[1] 8 2
#> dim(lancet_not_in_app56607_sample %>% filter(lancet_stage == "lancet_stage2"))
#[1] 105   2

#Final number of compared individuals for stage1 and stage2 (the numbers are correct):
#> dim(demo %>% filter(lancet_stage == "lancet_stage1"))
#[1] 2984   72
#> dim(demo %>% filter(lancet_stage == "lancet_stage2"))
#[1] 5306   72

#Comparison in text:
demo$lancet_or_our <- paste0(demo$cases_broad_EUR, demo$lancet_stage)
demo$lancet_or_our <- as.factor(demo$lancet_or_our)
print("Comparison between lancet stage 1/stage 2 moderate-severe participants and our severe asthma definition:")
print(summary(demo$lancet_or_our))

#Comparison in Venn Diagram:
demo_lancet1 <- demo %>% filter(lancet_stage == 'lancet_stage1') %>% rename(category=lancet_stage)
demo_lancet2 <- demo %>% filter(lancet_stage == 'lancet_stage2') %>% rename(category=lancet_stage)
demo_study <- demo %>% filter(broad_pheno_1_5_ratio == 1) %>% rename(category=lancet_stage)
venn_df <- rbind(demo_lancet1,demo_lancet2,demo_study)

#Venn diagram:
cases_venn <- venn.diagram(
   x = list(
     demo_lancet1 %>% select(eid) %>% distinct() %>% unlist(),
     demo_lancet2 %>% select(eid) %>% distinct() %>% unlist(),
     demo_study %>% select(eid) %>% distinct() %>% unlist()
    ),
   category.names = c("Lancet Stage1","Lancet stage 2","Study definition"),
   #main = "Lancet vs our cases comparison",
   #main.cex = 1,
   #main.just = c(0.5, 0.6),
   filename = "output/Lancet_vs_study_cases_comparison.png",
   output = FALSE ,
           imagetype="png" ,
           height = 1000 ,
           width = 1200 ,
           resolution = 400,
           compression = "lzw",
           lwd = 1,
           col=c("#440154ff", '#21908dff', '#fde725ff'),
           fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
           cex = 0.5,
           # Set names
           cat.pos = c(-20, 0, 20),
           cat.dist = c(0.1, 0.1, 0.1),
           fontfamily = "sans",
           cat.cex = 0.5)