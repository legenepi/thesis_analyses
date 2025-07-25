library(RColorBrewer)
library(tidyverse)
library(pheatmap)

func <- read.table("input/functionalscore.txt",header=T)
row.names(func) <- func$SNP
func$SNP <- NULL
func <- func %>% mutate_at(c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10'), as.numeric)
colnames(func) <- c("CADD","Protein_Function","Conservation","Epigenetics_Active","Epigenetics_Repressed","Epigenetics_Transcription","Local_Nucleotide_Diversity","Mutation_Density","Transcription_Factor","Mappability")
func_matrix <- as.matrix(func)
head(func_matrix)

df <- pheatmap(func_matrix, main="My Heatmap", cluster_cols=F, border_color=NA,  cluster_rows = FALSE)
ggsave("output/funcscore_heatmap.png", df, width = 10,  height = 10, dpi = "retina")



#Person correlation coefficient to measure corr between different functional integrative scores:
library(ggplot2)
library("GGally")


func <- read.table("input/functionalscore.txt",header=T)
func$SNP <- NULL
func <- func %>% mutate_at(c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10'), as.numeric)
colnames(func) <- c("CADD","Protein_Function","Conservation","Epigenetics_Active","Epigenetics_Repressed","Epigenetics_Transcription","Local_Nucleotide_Diversity","Mutation_Density","Transcription_Factor","Mappability")

#nice one, but not ideal for my data as I have only three correlation higher thatn by chance (0.5).
ggcorr(func, geom = "blank", label = TRUE, hjust = 0.75) +
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient) > 0.5)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = FALSE, alpha = FALSE)

#function for the plot:
correlation_plot <- ggcorr(func, hjust = 0.85, size = 5, color = "grey50", nbreaks = 4, palette = "RdGy", label = TRUE, label_size = 3, label_color = "white")
ggsave("output/funcscore_corrplot.png", correlation_plot, width = 10,  height = 10, dpi = "retina")


#with all the 344 vars:
suppressMessages(library(data.table))


favor_file <- "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/20250310_FAVOR_output_additional_credset_SNPs_processed.csv.gz"
favor <- fread(favor_file,na.strings = c("",NA))

func_all <- favor %>% select("CaddPhred", "ApcProteinFunctionV3", "ApcConservationV2", "ApcEpigeneticsActive", "ApcEpigeneticsRepressed",
                                 "ApcEpigeneticsTranscription", "ApcLocalNucleotideDiversityV3", "ApcMutationDensity",
                                 "ApcTranscriptionFactor", "ApcMappability")
colnames(func_all) <- c("CADD","Protein_Function","Conservation","Epigenetics_Active","Epigenetics_Repressed","Epigenetics_Transcription","Local_Nucleotide_Diversity","Mutation_Density","Transcription_Factor","Mappability")
correlation_plot_all <-ggcorr(func_all, hjust = 0.85, size = 5, color = "grey50", nbreaks = 4, palette = "RdGy", label = TRUE, label_size = 3, label_color = "white")
ggsave("output/funcscore_corrplot_all.png", correlation_plot_all, width = 10,  height = 10, dpi = "retina")
