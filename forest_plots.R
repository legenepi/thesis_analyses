#!/usr/bin/env Rscript

#User: Noemi Nicole Piga

#Rational: forest plots for the 21 replicated variants in the severe asthma GWAS
#workign directory: /home/n/nnp5/PhD/PhD_project/Thesis_analysis
#in bash
#cp /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/metanalysis_input_4cohorts.txt Thesis_analysis/input/
#input/replicated_sentinels : to find replicated vars in the meta-analysis input file

library(metaviz)
library(tidyverse)
library(data.table)
library(dplyr)

#in the metanalysis_input_4cohorts.txt, this is the order of the cols:
columns_names <- c("ID", "LOG_ODDS_GERA", "LOG_ODDS_GERA_SE", "LOG_ODDS_GUU", "LOG_ODDS_SE_GUU", "LOG_ODDS_chang", "LOG_ODDS_SE_chang","LOG_ODDS_CREW","LOG_ODDS_CREW_SE")
meta_input <- fread("input/metanalysis_input_4cohorts.txt",  sep =" ", header=F)
colnames(meta_input) <- columns_names
#extract rsid from ukb to be able to fitler out replicated vars:
meta_input <- separate(data = meta_input, col = ID, into = c("ID", NA, NA, NA, NA), sep = ";")
sentinel <- fread("input/replicated_sentinels",header=F)
meta_input_replicated <- meta_input %>% filter(meta_input$ID %in% sentinel$V1)
df_fr <- meta_input_replicated

study <- c("Chang","CREW","GERA","GUU")

for (i in 1:dim(df_fr)[1]) {snp <- as.character(df_fr[i,"ID"])
logOR <- c(as.numeric(df_fr[i,"LOG_ODDS_chang"]),as.numeric(df_fr[i,"LOG_ODDS_CREW"]),as.numeric(df_fr[i,"LOG_ODDS_GERA"]),as.numeric(df_fr[i,"LOG_ODDS_GUU"]))
selogOR <- c(as.numeric(df_fr[i,"LOG_ODDS_SE_chang"]),as.numeric(df_fr[i,"LOG_ODDS_CREW_SE"]),as.numeric(df_fr[i,"LOG_ODDS_GERA_SE"]),as.numeric(df_fr[i,"LOG_ODDS_SE_GUU"]))
df <- data.frame(logOR,selogOR,study)
#df$OR <- exp(df$logOR)
#df$OR.metaL_CI95.meta <- exp(df$logOR - qnorm(0.975)*df$selogOR)
#df$OR.metaU_CI95.meta <- exp(df$logOR + qnorm(0.975)*df$selogOR)
#min_val <- min(1,as.numeric(df_fr[i,"OR.metaL_CI95.meta"]),as.numeric(df_fr[i,"L_CI95.awigen"]),as.numeric(df_fr[i,"L_CI95.ugr"]),as.numeric(df_fr[i,"L_CI95.ukbafr"]))
#min_val <- ifelse(min_val == 1, 1, 0)
#max_val_tmp <- max(1,as.numeric(df_fr[i,"OR.metaU_CI95.meta"]),as.numeric(df_fr[i,"U_CI95.awigen"]),as.numeric(df_fr[i,"U_CI95.ugr"]),as.numeric(df_fr[i,"U_CI95.ukbafr"]))
#max_val <- ifelse(max_val_tmp == 1, 1, max_val_tmp)
#print(min_val)
#print(max_val)
viz_forest(x=df[,c("logOR","selogOR")],variant = "classic",xlab = "OR", annotate_CI = T,study_labels=df[,c("study")],
x_trans_function = exp, method = "DL", text_size = 10,summary_label = "Summary effect",
col="Greys",x_limit=c(-0.5,1.05))
ggsave(paste0("output/forest_plot_",snp,".pdf"), width = 20, height = 20)
}
