#This script was written on 9/3/2025 to take the p25.c2 plant reads from the two experiments and calculate FC

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(readr)
# This script reads in the rds experiments for the four experiments Done by Effie in 12/2022 and 3/2023 with Rb-Tnseq and does basic deseq comparison trials

# May 2025 This script takes the full counts matrix from the many trials and will do subsetting and limma voom analysis. The goal is to ask what percentage of sig genes are genetic background specific in their behavior. In August 2025 Effie determined that the original counts may have been problematic and redid the counts table.

#modified August 2025
library(limma)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(randomForest)
library(caret)
library(matrixStats)
library(grid)  # For theme text customization
library(edgeR)

# setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")
setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output")

##########################################################################
### Read in full counts table ###
# all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_5_2025.rds")
all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_8_2025.rds")
all_exp$Sample <- with(all_exp, paste(treatment, plant, time_point, position, experiment, sep = "_"))
# Move Sample column to the front
all_exp <- all_exp[, c("Sample", setdiff(names(all_exp), "Sample"))]
all_exp$plant <- tolower(all_exp$plant)

### Reassign control: when we originally did the experiment it turned out that pure culture at T0 was a more robust metric than the small amount of microbe on the plant at t0. So we need to split the pure culture controls in half and randomly assign as t0. THis should be done after removinga any actual plantxT0

# Filter: remove rows where plant ≠ Ctrl/ctrl AND time_point == "t0"
all_exp_filtered <- subset(all_exp, !(tolower(plant) != "ctrl" & time_point == "t0"))

# now for the plant controls we want to randomly assign to either ey15_2 or col_0
all_exp <- all_exp_filtered
set.seed(823)  # for reproducibility
ctrl_rows <- which(all_exp$plant == "ctrl")
all_exp$plant[ctrl_rows] <- sample(c("col_0", "ey15_2"), length(ctrl_rows), replace = TRUE)


### Subset metadata and count matrix
metadata_cols <- colnames(all_exp)[c(1:11)]
gene_cols <- colnames(all_exp)[13:dim(all_exp)[2]]

# Create metadata and counts matrix
metadata <- all_exp[,metadata_cols]
#metadata$plant <- tools::toTitleCase(metadata$plant)
rownames(metadata) <- metadata$Sample
counts <- all_exp[,gene_cols]
rownames(counts) <- all_exp$Sample

#remove samples with few reads
keep <- which(rowSums(counts, na.rm = TRUE)>=25000)
counts <- counts[keep,]
metadata_keep = metadata[keep,]
metadata <- metadata_keep
table(metadata[,c("experiment", "treatment", "time_point")])

# remove dc3000 samples
keep_p25 <- which(metadata$treatment != "dc3000")
meta_p25 <- metadata[keep_p25,]
counts_p25 <- counts[keep_p25,]
counts<- counts_p25
metadata <- meta_p25
##########################################################################
##########################################################################
### Step 2: Create interaction factors###
metadata$strain <- factor(metadata$treatment)
metadata$plant <- factor(metadata$plant)
group <- interaction(metadata$strain, metadata$plant, metadata$time_point)

# Transform counts with voom. The following takes the dge (subsetted or not, the design and a boolean for whether or not the experiment random effect should be taken into account)

run_voom <- function(dge, design, experiment=FALSE){
  v <- voom(dge, design, plot = TRUE)
  # Use duplicateCorrelation to estimate variance from random effect
  
  if(experiment==TRUE){
    corfit <- duplicateCorrelation(v, design, block = metadata$experiment)
    # Re-voom using the estimated correlation
    v <- voom(dge, design, block = metadata$experiment, correlation = corfit$consensus)
    
    # Fit model with block structure
    fit <- lmFit(v, design, block = metadata$experiment, correlation = corfit$consensus)
    
  }
  if(experiment==FALSE){
    v <- voom(dge, design)
    fit <- lmFit(v, design)
  } 
  # Empirical Bayes moderation for both
  fit <- eBayes(fit)
  return(list(v=v, fit=fit))
}



filter_meta_counts <- function(column_remove, param_remove, counts, metadata){
  # Relevel factors
  metadata$treatment <- factor(metadata$treatment, levels = c("dc3000", "p25c2"))
  metadata$time_point <- factor(metadata$time_point, levels = c("t0", "t3"))
  metadata$plant <- factor(metadata$plant)
  subset_meta <- rownames(metadata)[which(metadata[[column_remove]]==param_remove)]
  new_meta <- metadata[!rownames(metadata) %in% subset_meta, ]
  new_count <- as.data.frame(counts) %>% dplyr::select(-all_of(subset_meta))
  return(list(new_meta=new_meta, new_counts=new_count))
}

summarize_overlap_between_terms <- function(fit, coef1, coef2, pval_cutoff = 0.05) {
  # Get topTable results for both coefficients
  tt1 <- topTable(fit, coef = coef1, number = Inf, p.value = pval_cutoff, adjust.method = "BH")
  tt2 <- topTable(fit, coef = coef2, number = Inf, p.value = pval_cutoff,  adjust.method = "BH")
  
  # Gene names
  genes1 <- rownames(tt1)
  genes2 <- rownames(tt2)
  
  # Overlap
  overlap_genes <- intersect(genes1, genes2)
  
  # Summary statistics
  summary <- list(
    coef1 = coef1,
    coef2 = coef2,
    num_significant_coef1 = length(genes1),
    num_significant_coef2 = length(genes2),
    num_overlap = length(overlap_genes),
    percent_coef1_in_overlap = round(length(overlap_genes) / length(genes1) * 100, 2),
    percent_coef2_in_overlap = round(length(overlap_genes) / length(genes2) * 100, 2),
    overlap_gene_ids = overlap_genes
  )
  
}

results_table <- function(fit){
  # Build summary table for all model terms
  results_summary <- lapply(colnames(fit$coefficients), function(term) {
    tt <- topTable(fit, coef = term, number = Inf, sort.by = "none")
    
    sig_genes <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
    min_lfc   <- min(tt$logFC, na.rm = TRUE)
    max_lfc   <- max(tt$logFC, na.rm = TRUE)
    mean_lfc  <- mean(tt$logFC, na.rm = TRUE)
    
    data.frame(
      term = term,
      num_significant = sig_genes,
      mean_logFC = mean_lfc,
      min_logFC = min_lfc,
      max_logFC = max_lfc
    )
    return(results_summary)
  })
  
}

##########################################################################
##########################################################################
######### First do the analysis with all
# design matrix
#**If whole genes or samples are `NA`** → drop those rows/columns:

counts <- t(counts[,which(colSums(counts)>0) ])
design <- model.matrix(~ time_point * plant + experiment, data = metadata)
# experiment is included as a random effect only
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
res <- run_voom(dge, design, experiment = TRUE)
v <- res$v
fit <- res$fit
result_full <- summarize_overlap_between_terms(
  fit,
  coef1 = "time_pointt3",
  coef2 = "time_pointt3:plantey15_2"
)

# remove dc3000 orthologs only
lfc_time3 <- res$fit$coefficients[, "time_pointt3"]
lfc_time3_plant <- res$fit$coefficients[, "time_pointt3:plantey15_2"]

# now convert into table
## --- pick the exact coefficient names you want ---
coef_time      <- "time_pointt3"                  # example main-effect term
coef_interact  <- "time_pointt3:plantey15_2"      # example interaction term

# helper: resolve to the exact column name (strict first, then fuzzy)
resolve_coef <- function(target, cn) {
  if (target %in% cn) return(target)
  match <- grep(paste0("^", gsub("([\\[\\]\\(\\)\\^\\$\\.|?*+{}])", "\\\\\\1", target), "$"),
                cn, value = TRUE)
  if (length(match) == 1) return(match)
  # fallback: loose grep
  match <- grep(target, cn, value = TRUE)
  if (length(match) == 1) return(match)
  stop("Could not uniquely match coefficient: ", target,
       "\nAvailable:\n", paste(cn, collapse = ", "))
}

cn <- colnames(res$fit$coefficients)
coef_time     <- resolve_coef(coef_time, cn)
coef_interact <- resolve_coef(coef_interact, cn)

# extract LFCs and p-values
lfc_time  <- res$fit$coefficients[, coef_time]
p_time    <- res$fit$p.value[,      coef_time]

lfc_int   <- res$fit$coefficients[, coef_interact]
p_int     <- res$fit$p.value[,      coef_interact]

# build the output table
res_table <- data.frame(
  gene = rownames(res$fit$coefficients),
  lfc_timepoint = lfc_time,
  p_timepoint   = p_time,
  padj_timepoint = p.adjust(p_time, method = "BH"),
  lfc_plant_timepoint = lfc_int,
  p_plant_timepoint   = p_int,
  padj_plant_timepoint = p.adjust(p_int, method = "BH"),
  row.names = NULL,
  check.names = FALSE
)

# quick peek
head(res_table)

# recode table for the p25.c2 ID
orth <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", header = TRUE, row.names = 1, sep = ",")
orth[orth == ""] <- NA
lkup <- setNames(orth$p25_BJE, orth$DC3000)
res_table$p25_BJE <- ifelse(
  grepl("^WP", res_table$gene),           # only if gene starts with WP
  lkup[res_table$gene],                   # look up value (may be NA if not present)
  NA )                                     # else leave NA



#lookup <- orthologs[, c("DC3000", "p25_BJE")]
#colnames(lookup) <- c("gene", "p25_BJE")
#results table
#res_table2 <- merge(res_table, lookup, by = "gene", all.x = TRUE)
write.csv(
  res_table,
  file = "p25c2_only_tnseq_inplanta_9_3_2025_res_table.csv",
  row.names = FALSE
)

###### Let's compare with the old data
tail_fit <- read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/tailocin/tailocin_plant_fitness_612_2024.csv", sep = "\t")

both <- merge(res_table, tail_fit, by.x="p25_BJE", by.y="ID", all.x=TRUE)
plot(both$lfc_timepoint, both$log2FoldChange_Ey)

# this is terrible
cor=-0.02029946

#what went wrong?
