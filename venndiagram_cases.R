#!/usr/bin/env Rscript

#User: Noemi Nicole Piga

#Rationale: Compare cases if they were identified by prescription, death, or hospitalisation.

library(tidyverse)
library(dplyr)
library(data.table)
library(hablar)
library(VennDiagram)
#avoid log file for venn diagrams
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

args = commandArgs(trailingOnly=TRUE)

prescriptions = fread(args[1])
hesin = fread(args[2])
death = fread(args[3])

#prescriptions = fread("input/eid_asthma_bts45_eur")
#hesin = fread("input/eid_asthma_level1_hes_eur")
#death = fread("input/eid_asthma_primarycause_death_eur")

#Venn diagram:
venn.diagram(
   x = list(
     prescriptions %>% select(V1) %>% distinct() %>% unlist(),
     hesin %>% select(V1) %>% distinct() %>% unlist(),
     death %>% select(V1) %>% distinct() %>% unlist()
    ),
   category.names = c("Prescriptions","Hospitalisation","Death"),
   filename = "output/Cases_definition_source_comparison.png",
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
           #cat.pos = c(-20, 0, 20),
           #cat.dist = c(0.1, 0.1, 0.1),
           #fontfamily = "sans",
           cat.cex = 0.5)