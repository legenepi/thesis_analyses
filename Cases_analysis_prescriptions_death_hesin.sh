#!/usr/bin/env bash

#User: Noemi Nicole Piga

#Rationale: Compare cases if they were identified by prescription, death, or hospitalisation.

cd /home/n/nnp5/PhD/PhD_project/Thesis_analysis

#cases (take the cases id):
awk '$71 == 1 {print $1}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt \
    > input/eid_cases

#cases by prescription:
#BTS stage 4 and 5 phenotype participants:
awk '{print $1, "prescription"}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/eid_bts2019_4plus | \
    grep -F -f input/eid_cases - \
    > input/eid_asthma_bts45_eur

#cases by hesin:
awk '$4 == 1 {print $1}' \
    /scratch/gen1/nnp5/tmp_thesis/UKBiobank_datafields/output/QC_hesin_diag_asthma.txt | \
    sort -u | awk '{print $1, "hesin"}' | \
    grep -F -f input/eid_cases - \
    > input/eid_asthma_level1_hes_eur

#cases by death:
awk '{print $1}' /scratch/gen1/nnp5/tmp_thesis/UKBiobank_datafields/output/QC_asthma_PrimaryCauseOfDeath_40001.txt | \
    sort -u | awk '{print $1, "death"}' | \
    grep -F -f input/eid_cases - \
    > input/eid_asthma_primarycause_death_eur

Rscript src/venndiagram_cases.R \
    input/eid_asthma_bts45_eur \
    input/eid_asthma_level1_hes_eur \
    input/eid_asthma_primarycause_death_eur