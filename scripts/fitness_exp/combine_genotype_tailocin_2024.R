library(DESeq2)
library(tidyverse)
library(pheatmap)

# This script reads in the rds experiments for the four experiments Done by Effie in 12/2022 and 3/2023 with Rb-Tnseq
setwd("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/tailocin/")

# Talia's tailocin tnseq results
tb_p25 <- read.csv("./all_genes_table_tailocin_6_11_2024.csv", sep=",")

#read in the metacsv with all the experiments
sample_order_filter  <- readRDS("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/full_experiments/all_sample_order_filter_6_12_2024.rds")
count_table_filter <- readRDS("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/full_experiments/all_count_table_filter_6_12_2024.rds")

# We need to just consider the p25.c2 data
indeces_p25c2 <- which(sample_order_filter$treatment == "p25c2")
sample_p25 <- sample_order_filter[indeces_p25c2,]
count_p25 <- count_table_filter[,indeces_p25c2]


#########
# This blog was helpful for getting me to understand the interactions: https://www.biostars.org/p/353618/#356111
#now with the filtered dataset, make deseq object
dds <- DESeqDataSetFromMatrix(countData = count_p25, 
                              colData = sample_p25, design = ~ plant + 
                                experiment +  plant)
dds$plant <- factor(dds$plant, levels = c("ctrl","col_0", "ey15_2"))

#perform the median of ratios method of normalization
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="plant_ey15_2_vs_ctrl" )
res_eyach <- res[,c("baseMean", "log2FoldChange", "pvalue", "padj" )]
colnames(res_eyach) <- c("baseMean_Ey", "log2FoldChange_Ey", "pvalue_Ey", "padj_Ey" )

#Now we want to convert Talia B's tnseq result file to the general ortholog mappings
ortho <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthologues_dc3000_p25_c2_7_2022.txt", header = T, sep= "\t")
