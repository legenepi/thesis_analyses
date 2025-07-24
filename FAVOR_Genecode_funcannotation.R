#!/usr/bin/env Rscript
#Rationale: Variant Annotation FAVOR, Genecode functional annotation.
#To obtain FAVOR: first liftOver for b38, then online FAVOR:
#liftOver online:
#input/Additional_credset_snps_March2025/Credible_sets.xlsx
#input/Additional_credset_snps_March2025/LiftOver_output_b38_additional_credset_snps.bed
###use FAVOR webtool:
#input/Additional_credset_snps_March2025/
#https://favor.genohub.org/FAVOR_input_additional_credset_SNPs.txt

sink(stderr())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


favor_file <- "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/20250310_FAVOR_output_additional_credset_SNPs_processed.csv.gz"
favor <- fread(favor_file,na.strings = c("",NA))

favor.digest <- favor %>% select("VariantVcf", "Chromosome", "Position", "Rsid", "GenecodeComprehensiveCategory",
                                 "GenecodeComprehensiveInfo", "GenecodeComprehensiveExonicCategory", "GenecodeComprehensiveExonicInfo",
                                 "CagePromoter", "CageEnhancer", "Genehancer", "SuperEnhancer", "CaddPhred", "FathmmXf",
                                 "ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability", "Clnsig",
                                 "Clnsigincl", "Clndn", "Clndnincl", "Clnrevstat", "Origin",
                                 "Clndisdb", "Clndisdbincl", "Geneinfo")

colnames(favor.digest) <- make.names(colnames(favor.digest), unique=TRUE)
functional_cat <- favor.digest %>% select("GenecodeComprehensiveCategory","Chromosome", "Position", "Rsid")
fwrite(functional_cat," "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/output/funtionalannottion_genecode_credsetvars",sep="\t",quote=F)
