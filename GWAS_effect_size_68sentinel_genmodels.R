#!/usr/bin/env Rscript
#Rationale: compare effect size of the 68 sentinel variants between the discovery GWAS and 
# the recessive and domiantn models. Update 17/07/2025 wit the sugg sentinel on chr1 that was missed.

#On bash:
#ungzipped the summary stats for noallergy GWAS:
#tar -xvzf REGENIE_assoc_31082023.tar.gz \
#home/n/nnp5/PhD/PhD_project//REGENIE_assoc/output/maf001_noallergy_pheno_betase_input_mungestat
#cp in Rdrive:
#cp home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_noallergy_pheno_betase_input_mungestat \
#/rfs/TobinGroup/nnp5/data/
#R:\TobinGroup\nnp5\PhD_DocMiscell_Bridge\sentinel_vars_severeasthma.xlsx copy and paste vars into 68sentvars_gwas
#filter for sentinel variants:
#dos2unix /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt
#remove the duplicate row for rs762279794:
#grep -w -F -f /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt /rfs/TobinGroup/nnp5/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat  | \
#    grep -v "CAA" \
#    > /scratch/gen1/nnp5/tmp_manuscript/68sentvars_gwas

#grep -w -F -f /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt /rfs/TobinGroup/nnp5/data/maf001_broad_pheno_1_5_ratio_recessive_betase_input_mungestat  | \
#    grep -v "CAA" \
#    > /scratch/gen1/nnp5/tmp_manuscript/68sentvars_recessive

#grep -w -F -f /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt /rfs/TobinGroup/nnp5/data/maf001_broad_pheno_1_5_ratio_dominant_betase_input_mungestat  | \
#    grep -v "CAA" \
#    > /scratch/gen1/nnp5/tmp_manuscript/68sentvars_dominant



#In R:
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(readxl)
library(cowplot) # it allows you to save figures in .png file
library(smplot2)
library(ggpubr)
library("ggrepel")

#input file:
gwas <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_gwas") %>% select(V1, V6)
colnames(gwas) <- c("snpid", "logodds_gwas")

##recessive:
recessive <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_recessive") %>% select(V1, V6)
colnames(recessive) <- c("snpid", "logodds_recessive")
df_plot <- gwas %>% left_join(recessive, by = "snpid")
df_plot$logOR_diff <- abs(df_plot$logodds_gwas - df_plot$logodds_recessive)


png("output/68sentinels_DiscVSRecessive_effectsize_comparison.png",units="in", width=15, height=10, res=1000)
ggplot(data = df_plot, aes(x = logodds_recessive, y = logodds_gwas)) +
  sm_statCorr(color = "black", corr_method = "pearson", text_size = 3, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 0.25, linetype = "dashed") +
  geom_smooth(method="lm") +
  geom_point(data = df_plot, aes(x = logodds_recessive,color="#D55E00",size=0.1),shape = 20 ) +
  theme_minimal() + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size = 0.25) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.5) +
  ylim(-1.5, +2.0) + xlim(-1.5, +2.5) +
  xlab("Recessive") + ylab("Additive (Discovery)") +
  geom_text_repel(
    data = subset(df_plot, logOR_diff >= 0.2 ),
    aes(label = snpid),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = 30)
dev.off()


##dominant:
dominant <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_dominant") %>% select(V1, V6)
colnames(dominant) <- c("snpid", "logodds_dominant")
df_plot <- gwas %>% left_join(dominant, by = "snpid")
df_plot$logOR_diff <- abs(df_plot$logodds_gwas - df_plot$logodds_dominant)


png("output/68sentinels_DiscVSDominant_effectsize_comparison.png",units="in", width=10, height=10, res=800)
ggplot(data = df_plot, aes(x = logodds_dominant, y = logodds_gwas)) +
  sm_statCorr(color = "black", corr_method = "pearson", text_size = 3, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 0.25, linetype = "dashed") +
  geom_smooth(method="lm") +
  geom_point(data = df_plot, aes(x = logodds_dominant,color="#D55E00",size=0.15)) +
  theme_minimal() + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size = 0.25) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.5) +
  ylim(-0.6, +0.4) + xlim(-0.6, + 0.4) +
  xlab("Dominant") + ylab("Additive (Discovery)") +
  geom_text_repel(
    data = subset(df_plot, logOR_diff >= 0.1 ),
    aes(label = snpid),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
dev.off()

##ODDS Ratio and CI plot:
gwas <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_gwas") %>% select(V1, V6, V7)
dominant <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_dominant") %>% select(V1, V6, V7)
recessive <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_recessive") %>% select(V1, V6, V7)

gwas$test <- as.factor("additive")
dominant$test <- as.factor("dominant")
recessive$test <- as.factor("recessive")

additive <- gwas %>% rename(snpid = V1, LOG_ODDS = V6, se = V7)
dominant <- dominant %>% rename(snpid = V1, LOG_ODDS = V6, se = V7)
recessive <- recessive %>% rename(snpid = V1, LOG_ODDS = V6, se = V7)

add_dom <- rbind(additive,dominant)
add_dom_rec <- rbind(add_dom,recessive)

#create odds_ratio ci_lower ci_upper
add_dom_rec$odds_ratio <- exp(add_dom_rec$LOG_ODDS)
add_dom_rec$ci_lower <- exp(add_dom_rec$LOG_ODDS - qnorm(0.975)*add_dom_rec$se)
add_dom_rec$ci_upper <- exp(add_dom_rec$LOG_ODDS + qnorm(0.975)*add_dom_rec$se)
add_dom_rec$test <- as.factor(add_dom_rec$test)

#plot OR for additive and dominant:
add_dom <- add_dom_rec %>% filter(test != "recessive")
add_dom_or_plot <- add_dom %>% ggplot(aes(x = snpid, y = odds_ratio, color = test)) +
  geom_boxplot(position=position_jitterdodge()) + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#and save plot:
ggsave("output/OR_sentinel_genetic_test_boxplot_adddom.png",add_dom_or_plot,width = 20, height = 8)


#plot all three to show how recessive has large confidence intervals:
#remove rs143191487 because it has a confidence interval too big for recessive:
add_dom_rec <- add_dom_rec %>% filter(snpid != "rs143191487")
or_plot <- add_dom_rec %>% ggplot(aes(x = snpid, y = odds_ratio, color = test)) +
  geom_boxplot(position=position_jitterdodge()) + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#and save plot:
ggsave("output/OR_sentinel_genetic_test_boxplot_adddomrec.png",or_plot,width = 20, height = 8)
