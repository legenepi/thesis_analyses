#!/usr/bin/env Scrip

#Rationale: look up in +/- 500Kb from each suggestive sentinel region in GWAS Catalog,  (i in seq(1,length(sugg$snpid))) { on the 16/07/2025

#Download GWASCAtalog file as by 03/08/2023:
#Data in the GWAS Catalog is currently mapped to genome assembly GRCh38.p14 and dbSNP Build 154
#https://www.ebi.ac.uk/gwas/docs/file-downloads --> All associations v1.0:
#wget -d input/ https://www.ebi.ac.uk/gwas/api/search/downloads/full
#mv full input/gwas_catalog_v1_e114_r2025-07-10


#in R
library(tidyverse)
library(dplyr)
library(data.table)
library(readxl)

sugg <- read_excel("/rfs/TobinGroup/nnp5/PhD_DocMiscell_Bridge/sentinel_vars_severeasthma.xlsx", sheet = "Sheet1")
catalog <- fread("input/gwas_catalog_v1_e114_r2025-07-10",sep="\t",quote="")
catalog <- catalog %>% select("CHR_ID","CHR_POS","SNPS","DISEASE/TRAIT","MAPPED_GENE","RISK ALLELE FREQUENCY","P-VALUE","OR or BETA","95% CI (TEXT)","PUBMEDID")
catalog <- catalog %>% rename(chr=CHR_ID,pos_b38=CHR_POS)
sugg$chr <- as.numeric(sugg$chr)
sugg$pos_b38 <- as.numeric(sugg$pos_b38)
catalog$chr <- as.numeric(catalog$chr)
catalog$pos_b38 <- as.numeric(catalog$pos_b38)

sugg_catalog <- inner_join(sugg,catalog, by=c("chr","pos_b38"))
##Variants present in GWASCatalog:
##25 snpid
length(unique(sugg_catalog$snpid))
##24 genes/gene families:
length(unique(sugg_catalog$snpid))
##324 distinct traits
length(unique(sugg_catalog$'DISEASE/TRAIT'))

#sentinel in GWAS Catalog with asthma:
sugg_catalog_asthma <- sugg_catalog %>% filter(str_detect(sugg_catalog$'DISEASE/TRAIT', regex("Asthma|asthma")))
print("16 sentinel in GWAs Catalog with asthma:")
print(unique(sugg_catalog_asthma$snpid))

#save value for those identified in GWAS Catalog:
fwrite(sugg_catalog, "output/25_sentinels_inGWASCatalog", sep ="\t", quote=F)

#Associated loci +/- 500000 base pair for all the 68 suggestive sentinel variants:
columns= colnames(catalog)
vars <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(vars) = columns
width = 500000
for (i in seq(1,length(sugg$snpid))) {
sentinel_var <- as.character(sugg[i,1])
chr_sentinel <- as.numeric(sugg[i,2])
OR_sentinel_var <- sugg[i,8]
pos_sentinel <- as.numeric(sugg[i,3])
min_pos <- pos_sentinel - width
max_pos <- pos_sentinel + width
chr_sentinel <- as.numeric(sugg[i,2])
vars_tmp <- catalog %>% filter(chr == chr_sentinel) %>% filter(pos_b38 >= min_pos) %>% filter(pos_b38 <= max_pos)
vars_tmp$UKBB_sent_vars <- sentinel_var
vars_tmp$UKBB_sent_vars_OR <- OR_sentinel_var
vars = rbind(vars, vars_tmp)
}

fwrite(vars, "output/68_sentinelvars_GWASCatalog",  sep ="\t", quote=F)

#Assoc with related-asthma trait:
related_traits <- vars %>% filter(str_detect(vars$'DISEASE/TRAIT', regex("allergy|allergic|Allergy|Allergic|FEV1|fev1|FVC|lung function|Lung Function|lung function|Neutro|neutro|Eosin|eosin")))
print(unique(related_traits$UKBB_sent_vars))
#nearest proxy for a related-asthma trait for all vars:
ukbb_related <- related_traits %>% select(UKBB_sent_vars) %>% unique()
columns= colnames(catalog)
one_proxy_related_all <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(one_proxy_related_all) = columns
for (i in seq(1,length(ukbb_related$UKBB_sent_vars))) {
sentinel_var <- as.character(ukbb_related[i,1])
sentinel_proxies <- related_traits %>% filter(UKBB_sent_vars == sentinel_var)
sentinel_pos <- as.numeric(sugg %>% filter(snpid == sentinel_var) %>% select(pos_b38))
sentinel_proxies$bp_diff <- as.numeric(abs(sentinel_proxies$pos_b38 - sentinel_pos))
sentinel_proxy <- sentinel_proxies %>% filter(bp_diff == as.numeric(min(sentinel_proxies$bp_diff)))
print(sentinel_var)
one_proxy_related_all = rbind(one_proxy_related_all, sentinel_proxy)
}
fwrite(one_proxy_related_all, "output/sentinels_nearestproxy_relatedtrait_GWASCat", sep ="\t", quote=F)


#All association with asthma:
all_assoc_asthma <- vars %>% filter(str_detect(vars$'DISEASE/TRAIT', regex("Asthma|asthma")))
print(unique(all_assoc_asthma$UKBB_sent_vars))
fwrite(all_assoc_asthma, "output/68_sentinelvars_GWASCatalog_asthma",  sep ="\t", quote=F)

#nearest proxy with asthma for all 68 vars:
ukbb_asthma <- all_assoc_asthma %>% select(UKBB_sent_vars) %>% unique()
columns= colnames(catalog)
one_proxy_asthma_all <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(one_proxy_asthma_all) = columns
for (i in seq(1,length(ukbb_asthma$UKBB_sent_vars))) {
sentinel_var <- as.character(ukbb_asthma[i,1])
sentinel_proxies <- all_assoc_asthma %>% filter(UKBB_sent_vars == sentinel_var)
sentinel_pos <- as.numeric(sugg %>% filter(snpid == sentinel_var) %>% select(pos_b38))
sentinel_proxies$bp_diff <- as.numeric(abs(sentinel_proxies$pos_b38 - sentinel_pos))
sentinel_proxy <- sentinel_proxies %>% filter(bp_diff == as.numeric(min(sentinel_proxies$bp_diff)))
print(sentinel_var)
one_proxy_asthma_all = rbind(one_proxy_asthma_all, sentinel_proxy)
}
fwrite(one_proxy_asthma_all, "output/sentinels_inGWASCatalog_asthma", sep ="\t", quote=F)


#Found for physical proxy in GWAS Catalog for the other 43 sentinel variants. 
sentinel_no_catalog <- sugg %>%	filter(!sugg$snpid %in% unique(sugg_catalog$snpid))

##Associated loci +/- 500000 base pair for variants with proxy:
columns= colnames(catalog)
vars <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(vars) = columns
width = 500000
for (i in seq(1,length(sentinel_no_catalog$snpid))) {
sentinel_var <- as.character(sentinel_no_catalog[i,1])
OR_sentinel_var <- sentinel_no_catalog[i,8]
pos_sentinel <- as.numeric(sentinel_no_catalog[i,3])
min_pos <- pos_sentinel - width
max_pos <- pos_sentinel + width
chr_sentinel <- as.numeric(sentinel_no_catalog[i,2])
vars_tmp <- catalog %>% filter(chr == chr_sentinel) %>% filter(pos_b38 >= min_pos) %>% filter(pos_b38 <= max_pos)
vars_tmp$UKBB_sent_vars <- sentinel_var
vars_tmp$UKBB_sent_vars_OR <- OR_sentinel_var
vars = rbind(vars, vars_tmp)
}

columns= colnames(catalog)
one_proxy <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(one_proxy) = columns
for (i in seq(1,length(sentinel_no_catalog$snpid))) {
sentinel_var <- as.character(sentinel_no_catalog[i,1])
sentinel_proxies <- vars %>% filter(UKBB_sent_vars == sentinel_var)
sentinel_pos <- as.numeric(sentinel_no_catalog %>% filter(snpid == sentinel_var) %>% select(pos_b38))
sentinel_proxies$bp_diff <- as.numeric(abs(sentinel_proxies$pos_b38 - sentinel_pos))
sentinel_proxy <- sentinel_proxies %>% filter(bp_diff == as.numeric(min(sentinel_proxies$bp_diff)))
print(sentinel_var)
one_proxy = rbind(one_proxy, sentinel_proxy)
}

fwrite(one_proxy, "output/43_sentinels_proxyinGWASCatalog", sep ="\t", quote=F)

#proxy associated with asthma:
proxy_asthma <-	vars %>% filter(str_detect(vars$'DISEASE/TRAIT', regex("Asthma|asthma")))
print("proxy associated with asthma summarise 26 sentinel vars:")
print(unique(proxy_asthma$UKBB_sent_vars))
fwrite(proxy_asthma, "output/sentinels_proxyinGWASCatalog_withasthma", sep ="\t", quote=F)

#nearest proxy with asthma:
ukbb_wihtproxyforasthma <- proxy_asthma %>% select(UKBB_sent_vars) %>% unique()
columns= colnames(catalog)
one_proxy_asthma <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(one_proxy_asthma) = columns
for (i in seq(1,length(ukbb_wihtproxyforasthma$UKBB_sent_vars))) {
sentinel_var <- as.character(ukbb_wihtproxyforasthma[i,1])
sentinel_proxies <- proxy_asthma %>% filter(UKBB_sent_vars == sentinel_var)
sentinel_pos <- as.numeric(sentinel_no_catalog %>% filter(snpid == sentinel_var) %>% select(pos_b38))
sentinel_proxies$bp_diff <- as.numeric(abs(sentinel_proxies$pos_b38 - sentinel_pos))
sentinel_proxy <- sentinel_proxies %>% filter(bp_diff == as.numeric(min(sentinel_proxies$bp_diff)))
print(sentinel_var)
one_proxy_asthma = rbind(one_proxy_asthma, sentinel_proxy)
}
fwrite(one_proxy_asthma, "output/sentinels_proxyinGWASCatalog_withasthma_nearest", sep ="\t", quote=F)



print("All variants are found in GWAS Catalog, 43 by proxy withihn +/- 500Kb. A total of 42 variants found associated with asthma, 26 by proxy variants.")

