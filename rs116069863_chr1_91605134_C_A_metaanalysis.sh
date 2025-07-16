#!/usr/bin/env bash

#User: Noemi Nicole Piga

#Rationale: Meta-analysis in GUU and Chang for the rs116069863 which I forgot to give to collabs.

cd /home/n/nnp5/PhD/PhD_project/Thesis_analysis/


#GUU:
cp /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/home/n/nnp5/PhD/PhD_project/Post_GWAS/input/ubio_gasp_gwas/maf001_uug_betase_input_mungestat input/
grep "chr1:91605134" input/maf001_uug_betase_input_mungestat |
    awk '{print "rs116069863", $1, $6, $7}' \
    > input/guu_rs116069863
#header:
echo "SNP_UKB SNP_GUU LOG_ODDS_GUU LOG_ODDS_SE_GUU" > input/header_guu
cat  input/header_guu input/guu_rs116069863 > input/guu_rs116069863_metasoft_input

#Chang:
cp /rfs/TobinGroup/kaf19/Public_data/Chang_et_al/Risk/shared/modsev_asthma_wgs_nodiff.tsv \
    input/
grep "rs116069863\|SNP" input/modsev_asthma_wgs_nodiff.tsv | \
    awk '{print $1, $6, $7}' - > input/chang_rs116069863

dos2unix src/logOR_chang_rs116069863.R
chmod o+x src/logOR_chang_rs116069863.R
Rscript src/logOR_chang_rs116069863.R

##merge the two sumstats:
dos2unix src/create_metasoft_input_ChangGUU_rs116069863.R
chmod o+x src/create_metasoft_input_ChangGUU_rs116069863.R
Rscript src/create_metasoft_input_ChangGUU_rs116069863.R

