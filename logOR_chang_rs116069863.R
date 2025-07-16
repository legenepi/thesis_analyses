#!/usr/bin/env Rscript

#User: Noemi Nicole Piga

#Rationale: LogOR for rs116069863 for Chang et al for meta-analysis

library(tidyverse)
library(dplyr)
library(data.table)

chang <- fread("input/chang_rs116069863")
chang$LOG_ODDS_chang <- log(chang$OR)
colnames(chang)[1] <- "SNP_chang"
colnames(chang)[3] <- "LOG_ODDS_SE_chang"
chang$SNP_UKB <- "rs116069863"
chang <- chang %>% select(SNP_UKB, SNP_chang, LOG_ODDS_chang, LOG_ODDS_SE_chang)
fwrite(chang, "input/chang_rs116069863_metasoft_input", quote = F, sep = " ")