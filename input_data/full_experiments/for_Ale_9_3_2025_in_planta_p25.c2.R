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

setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")
#setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output")

##########################################################################
### Read in full counts table ### The 8_2025 is just wrong. Don't know why.
 all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_5_2025.rds")
#all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_8_2025.rds")
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
# remove experiment 1
metadata1 <- metadata[which(metadata$experiment == "exp_0002"),]
counts1 <- counts[which(metadata$experiment == "exp_0002"),] 
table(metadata[,c("experiment", "treatment", "time_point")])

# remove dc3000 samples
keep_p25 <- which(metadata1$treatment != "dc3000")
meta_p25 <- metadata1[keep_p25,]
counts_p25 <- counts1[keep_p25,]
counts<- counts_p25
metadata <- meta_p25

##########################################################################
##########################################################################
######### First do the analysis with all
# design matrix


counts_t <- t(counts_p25)
counts_t[is.na(counts_t)] <- 0
design <- model.matrix(~ time_point * plant, data = meta_p25)

dds <- DESeqDataSetFromMatrix(
  countData = counts_t,
  colData   = meta_p25,
  design    = ~ plant * time_point
)


dds$plant <- factor(dds$plant, levels = c("col_0", "ey15_2"))
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, test = "Wald") 

vst_mat  <- assay(vst(dds, blind = TRUE))  
p <- plotPCA(vst(dds, blind = TRUE), intgroup = c("plant","time_point"))
p
#this experiment does the LRT test to compare whether the time interaction with bacterial treatment have different effects. *** Note again that this is comparing the two models and asking which genes does it matter whether there is a different interaction term depending on which strain is used in the comparison.


# Helper to find coefficient names robustly
find_coef <- function(dds, pattern) {
  nm <- resultsNames(dds)
  hit <- grep(pattern, nm, value = TRUE)
  if (length(hit) == 0) stop("No coefficient matches: ", pattern)
  if (length(hit) > 1)  warning("Multiple matches for ", pattern, "; using first: ", hit[1])
  hit[1]
}

# 1) Identify the main time effect (T3 vs T0) and all plant×T3 interaction coefs
main_T3 <- find_coef(dds, "^time_point[_]?t3[_]?vs[_]?t0$")

# Collect interaction coefficients for each non-reference plant level
plants <- levels(dds$plant)
ref_plant <- plants[1]  # should be "ctrl" if you set it as reference
nonref_plants <- setdiff(plants, ref_plant)

int_names <- setNames(
  vapply(nonref_plants,
         function(p) find_coef(dds, paste0("^plant.*", p, ".*time_pointt3$")),
         character(1)),
  nonref_plants
)

# 2) Extract:
# (a) main T3 vs T0 (baseline plant)
res_main <- results(dds, name = main_T3) %>%
  as.data.frame() %>%
  mutate(gene = rownames(.),
         effect = "time_main",
         contrast = paste0("T3 vs T0 (baseline=", ref_plant, ")")) %>%
  select(gene, effect, contrast, log2FoldChange, lfcSE, stat, pvalue, padj)

# (b) interaction terms for each plant (how plant modifies the T3 effect vs baseline)
res_int_list <- lapply(names(int_names), function(p) {
  results(dds, name = int_names[[p]]) %>%
    as.data.frame() %>%
    mutate(gene = rownames(.),
           effect = "time_interaction",
           contrast = paste0("Interaction: ", p, " × (T3 vs T0)")) %>%
    select(gene, effect, contrast, log2FoldChange, lfcSE, stat, pvalue, padj)
})
res_int <- bind_rows(res_int_list)

# (c) plant-specific T3 vs T0 contrasts = main + interaction (or just main for baseline)
res_by_plant <- bind_rows(
  # baseline plant = main effect as-is
  results(dds, name = main_T3) %>%
    as.data.frame() %>%
    mutate(gene = rownames(.),
           effect = "time_within_plant",
           contrast = paste0(ref_plant, ": T3 vs T0")) %>%
    select(gene, effect, contrast, log2FoldChange, lfcSE, stat, pvalue, padj),
  # non-reference plants = main + interaction (use list/contrast)
  bind_rows(lapply(names(int_names), function(p) {
    results(dds, contrast = list(c(main_T3, int_names[[p]]))) %>%
      as.data.frame() %>%
      mutate(gene = rownames(.),
             effect = "time_within_plant",
             contrast = paste0(p, ": T3 vs T0")) %>%
      select(gene, effect, contrast, log2FoldChange, lfcSE, stat, pvalue, padj)
  }))
)

# 3) Combine everything for easy export/plotting
all_time_effects <- bind_rows(res_main, res_int, res_by_plant)

# Optional: log2FC shrinkage for more stable estimates (recommended for TnSeq)
# Choose which coefficients to shrink, e.g., main_T3 and the interaction coefs
# Requires apeglm: BiocManager::install("apeglm")
# shr_T3 <- lfcShrink(dds, coef = main_T3, type = "apeglm") %>% as.data.frame()

# Peek
all_time_effects %>% glimpse()
head(all_time_effects)




# Graph the top ten genes
# Get results for the main time effect
res_time <- results(dds, name = "time_point_t3_vs_t0")

# Order by adjusted p-value (or abs(log2FC, depending on what you want)
top10_genes <- rownames(head(res_time[order(res_time$padj), ], 10))
top10_genes




# Normalized counts (library size–corrected)
norm_counts <- counts(dds, normalized = TRUE)

plot_data <- as.data.frame(t(norm_counts[top10_genes, ]))
plot_data$sample     <- rownames(plot_data)
plot_data$time_point <- colData(dds)$time_point
plot_data$plant      <- colData(dds)$plant



plot_long <- plot_data %>%
  pivot_longer(cols = all_of(top10_genes), names_to = "gene", values_to = "norm_count")

# jittered scatter + boxplot per time point

plot10 <- ggplot(plot_long, aes(x = time_point, y = norm_count + 1, color = time_point)) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  labs(y = "Normalized counts (log10 scale)", title = "Top 10 time_main genes (T3 vs T0)")

pdf()
plot10
dev.off()

# compare to the original tailocin dataset
new <- as.data.frame(res_time[, c("log2FoldChange","pvalue" )])
orig <- read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/tailocin/tailocin_plant_fitness_612_2024.csv", header=TRUE, sep="\t")
keep_row <- rownames(orig[which(rownames(orig)%in% rownames(new)),])
new_keep = new[keep_row,]
orig_keep = orig[keep_row,]

plot(orig_keep$log2FoldChange_Ey, new_keep$log2FoldChange)
cor.test(orig_keep$log2FoldChange_Ey, new_keep$log2FoldChange)
# this relationship is still terrible. wtf is going on.

# let's plot the limma voom results against the deseq results.
lims <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/p25c2_only_tnseq_inplanta_9_3_2025_res_table2.csv", header = TRUE)


# merge
orig_keep$gene <- rownames(orig_keep)
new_keep$gene <- rownames(new_keep)
hm = merge(orig_keep, new_keep, by = "gene")
hm2 = merge(hm, lims, by = "gene")



plot(hm2$log2FoldChange, hm2$log2FoldChange_Ey)
plot(hm2$log2FoldChange, hm2$lfc_timepoint)
cor.test(hm2$log2FoldChange, hm2$lfc_timepoint)
#log2FoldChange and lfc_timepoint are really well correlated with one another. 

# OK the next thing we need to check is if the data tables are wildly different somehow. 
qqplot(hm2$log2FoldChange, hm2$log2FoldChange_Ey)
tail = hm2[hm2$logFC_tailocin>2,]
both = tail[tail$log2FoldChange_Ey<(-2),]
# let's compare some of our favorite genes that we have validated.
# "WP_003317533.1" 
# "WP_005771402.1"** 
# "WP_011102979.1" 
# "WP_011103413.1" ****rfbA
## "WP_011103695.1" "WP_011105261.1"

#the following are the ones in both
#"WP_005762973.1"   
#"WP_007244130.1"   "WP_007244130.1"   "WP_011103412.1.1" "WP_011103696.1"   "WP_011105261.1"   "WP_011105263.1" 

my_genes <- c("WP_007244130.1",   "WP_007244130.1",   "WP_011103412.1.1", "WP_011103696.1" ,  "WP_011105261.1",   "WP_011105263.1" )
# 2) Normalized counts subset → tidy long with a SAFE join key
norm_counts <- counts(dds, normalized = TRUE)
sub_counts  <- norm_counts[my_genes, , drop = FALSE]

plot_df <- sub_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%                    # won't collide
  pivot_longer(-gene, names_to = "sample_id", values_to = "norm_count")

# 3) Build a metadata table with a non-colliding key
meta_df <- as.data.frame(colData(dds)) %>%
  mutate(sample_id = rownames(.)) %>%                     # keep existing columns intact
  # if you ALSO have a 'sample' column already, keep both; join uses sample_id
  select(sample_id, everything())

# 4) Join safely (no duplicate column names)
plot_df <- plot_df %>%
  left_join(meta_df, by = "sample_id")

# 5) Plot
ggplot(plot_df, aes(x = time_point, y = norm_count + 1, color = time_point)) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(y = "Normalized counts (log10)",
       x = "Time point",
       title = "Selected genes: T0 vs T3")





# okay let's graph this for the old counts table. I am pulling from combine_genotype_tailocin_2024.R
# This script reads in the rds experiments for the four experiments Done by Effie in 12/2022 and 3/2023 with Rb-Tnseq
setwd("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/tailocin/")

# Talia's tailocin tnseq results
tb_p25 <- read.csv("./all_genes_table_tailocin_6_11_2024.csv", sep=",")
tb_p25 <- tb_p25[, c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
colnames(tb_p25) <- c("ID", "logFC_tailocin", "AveExpr_tailocin", "t_tailocin", "P.Value_tailocin", "adj.P.Val_tailocin")

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
mat <- fpm(dds, robust = FALSE)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="plant_ey15_2_vs_ctrl" )
res_eyach <- res[,c("baseMean", "log2FoldChange", "pvalue", "padj" )]
colnames(res_eyach) <- c("baseMean_Ey", "log2FoldChange_Ey", "pvalue_Ey", "padj_Ey" )

# need to add a column that is the initial frequency of the gene
p0_eyach_t0_p25c2<- dds[,which(dds$plant=="ctrl" & dds$treatment=="p25c2")]
res_eyach$p0_axenic <- rowMeans(fpm(p0_eyach_t0_p25c2, robust = FALSE))/1000000























# # experiment is included as a random effect only
# dge <- DGEList(counts = counts)
# dge <- calcNormFactors(dge)
# res <- run_voom(dge, design)
# v <- res$v
# fit <- res$fit
# result_full <- summarize_overlap_between_terms(
#   fit,
#   coef1 = "time_pointt3",
#   coef2 = "time_pointt3:plantey15_2"
# )
# 
# # remove dc3000 orthologs only
# lfc_time3 <- res$fit$coefficients[, "time_pointt3"]
# lfc_time3_plant <- res$fit$coefficients[, "time_pointt3:plantey15_2"]
# 
# # now convert into table
# ## --- pick the exact coefficient names you want ---
# coef_time      <- "time_pointt3"                  # example main-effect term
# coef_interact  <- "time_pointt3:plantey15_2"      # example interaction term
# 
# # helper: resolve to the exact column name (strict first, then fuzzy)
# resolve_coef <- function(target, cn) {
#   if (target %in% cn) return(target)
#   match <- grep(paste0("^", gsub("([\\[\\]\\(\\)\\^\\$\\.|?*+{}])", "\\\\\\1", target), "$"),
#                 cn, value = TRUE)
#   if (length(match) == 1) return(match)
#   # fallback: loose grep
#   match <- grep(target, cn, value = TRUE)
#   if (length(match) == 1) return(match)
#   stop("Could not uniquely match coefficient: ", target,
#        "\nAvailable:\n", paste(cn, collapse = ", "))
# }
# 
# cn <- colnames(res$fit$coefficients)
# coef_time     <- resolve_coef(coef_time, cn)
# coef_interact <- resolve_coef(coef_interact, cn)
# 
# # extract LFCs and p-values
# lfc_time  <- res$fit$coefficients[, coef_time]
# p_time    <- res$fit$p.value[,      coef_time]
# 
# lfc_int   <- res$fit$coefficients[, coef_interact]
# p_int     <- res$fit$p.value[,      coef_interact]
# 
# # build the output table
# res_table <- data.frame(
#   gene = rownames(res$fit$coefficients),
#   lfc_timepoint = lfc_time,
#   p_timepoint   = p_time,
#   padj_timepoint = p.adjust(p_time, method = "BH"),
#   lfc_plant_timepoint = lfc_int,
#   p_plant_timepoint   = p_int,
#   padj_plant_timepoint = p.adjust(p_int, method = "BH"),
#   row.names = NULL,
#   check.names = FALSE
# )
# 
# # quick peek
# head(res_table)

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
