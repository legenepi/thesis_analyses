#!/usr/bin/env Rscript

#User: Noemi Nicole Piga

#Rationale: Create metasoft input for rs11606986 from GUU and Chang

library(tidyverse)
library(dplyr)
library(data.table)


guu <- fread("input/guu_rs116069863_metasoft_input")
chang <- fread("input/chang_rs116069863_metasoft_input")
guu_chang <- guu %>% left_join(chang, by ="SNP_UKB")
guu_chang$ID <- paste0(guu_chang$SNP_UKB,";",guu_chang$SNP_GUU,";",guu_chang$SNP_chang)
guu_chang <- guu_chang %>% select(ID, LOG_ODDS_GUU, LOG_ODDS_SE_GUU, LOG_ODDS_chang, LOG_ODDS_SE_chang)
fwrite(guu_chang,"input/metasoft_changGuu_rs11606986_input",sep=" ",quote=FALSE,row.names=FALSE,col.names=F,na="NA")

#Meta-analysis for rs116069863 in Chang and GUU:
java -jar /home/n/nnp5/software/METASOFT/software-notes/docs/files/METASOFT/Metasoft.jar \
    -input input/metasoft_changGuu_rs11606986_input\
    -pvalue_table /home/n/nnp5/software/METASOFT/software-notes/docs/files/METASOFT/HanEskinPvalueTable.txt \
    -output output/sevasthma_metasoft_output_changguu_rs11606986.txt \
    -log output/sevasthma_metasoft_changguu_rs11606986_log.txt \
    -mvalue -mvalue_p_thres 5E-3
#not significant

cp output/sevasthma_metasoft_output_changguu_rs11606986.txt \
    /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/