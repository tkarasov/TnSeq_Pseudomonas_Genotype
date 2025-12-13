# May 2025 This script takes the full counts matrix from the many trials and will do subsetting and limma voom analysis. The goal is to ask what percentage of sig genes are genetic background specific in their behavior. In August 2025 Effie determined that the original counts may have been problematic and redid the counts table.
# oh no, that table must messed up. The final counts table is in 10/2025
#modified December 2025
library(limma)
library(tidyr)
library(ggplot2)
library(lme4)
library(randomForest)
library(caret)
library(matrixStats)
library(grid)  # For theme text customization
library(edgeR)
library(cowplot)
library(dplyr)

#functions
extract_timepoint_logfc_per_experiment <- function(voom_obj, metadata, experiment_col, design_formula) {
  library(limma)
  
  experiment_levels <- unique(metadata[[experiment_col]])
  logfc_list <- list()
  
  for (exp_id in experiment_levels) {
    # Subset metadata and expression
    keep <- metadata[[experiment_col]] == exp_id
    metadata_sub <- metadata[keep, ]
    expr_sub <- voom_obj$E[, keep]
    weights_sub <- voom_obj$weights[, keep]
    
    # Create clean design matrix
    design_sub <- model.matrix(design_formula, data = metadata_sub)
    colnames(design_sub) <- make.names(colnames(design_sub))  # sanitize colnames
    
    # Fit model
    fit_sub <- lmFit(expr_sub, design_sub, weights = weights_sub)
    fit_sub <- eBayes(fit_sub)
    
    # Pull the logFC directly from the sanitized time_pointt3 column
    term <- "time_pointt3"
    term_clean <- make.names(term)  # should become "time_pointt3"
    
    if (!(term_clean %in% colnames(fit_sub$coefficients))) {
      warning(paste("Coefficient", term_clean, "not found in experiment", exp_id))
      next
    }
    
    logfc_list[[exp_id]] <- fit_sub$coefficients[, term_clean]
  }
  
  logfc_df <- do.call(cbind, logfc_list)
  colnames(logfc_df) <- paste0("logFC_", names(logfc_list))
  return(as.data.frame(logfc_df))
}
run_voom <- function(dge, design, block = NULL){
  v <- voom(dge, design, plot = TRUE)
  if (!is.null(block)) {
    corfit <- duplicateCorrelation(v, design, block = block)
    v <- voom(dge, design, block = block, correlation = corfit$consensus)
    fit <- lmFit(v, design, block = block, correlation = corfit$consensus)
  } else {
    fit <- lmFit(v, design)
  }
  fit <- eBayes(fit)
  list(v = v, fit = fit)
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

# setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")
setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")

##########################################################################
### Read in full counts table ###
##########################################################################
# all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_5_2025.rds")
#all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_8_2025.rds")
all_exp <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/full_experiments/all_p25_dc_axenic_10_07_2025.csv")
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
rownames(metadata) <- metadata$Sample
counts <- all_exp[,gene_cols]
rownames(counts) <- metadata$Sample

#remove samples with few reads
keep <- which(rowSums(counts, na.rm = TRUE)>=25000)
counts <- counts[keep,]
metadata_keep = metadata[keep,]
metadata <- metadata_keep
table(metadata[,c("experiment", "treatment", "time_point")])

#counts matrix with just orthologs. Get rid of lines that dont have shared orthologs. replace NAs with 0s
ortholog_tab <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/orthology_host_run/blast_pairs_vs_orthogroups_matches_safe.csv", header = TRUE, sep = ",")
#read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", header = TRUE, row.names = 1, sep = ",")
maintain <- c(unique(ortholog_tab$DC3000))
maintain <- maintain[maintain!=""]
count_ortholog<- counts[,maintain]
count_ortholog[is.na(count_ortholog)] <- 0
count_ortholog <- count_ortholog[, colSums(count_ortholog) > 0]
count_ortholog <- t(count_ortholog)
counts <- count_ortholog

##########################################################################
##########################################################################
### Step 2: Create interaction factors###
metadata$strain <- factor(metadata$treatment)
metadata$plant <- factor(metadata$plant)
group <- interaction(metadata$strain, metadata$plant, metadata$time_point)


######### Second do the analysis with Exp 2
filtered2 <- filter_meta_counts("experiment", "exp_0001", counts, metadata)

#I want to further filter filter2 to remove genes that are not represented in exp2 in one or both strains
# 1. Identify which samples belong to each strain
is_dc   <- filtered2$new_meta$strain == "dc3000"
is_p25  <- filtered2$new_meta$strain == "p25c2"

# 2. For each gene, check if it is zero across *all* DC3000 samples
dc_all_zero <- rowSums(filtered2$new_counts[, is_dc, drop = FALSE] > 0) == 0

# 3. For each gene, check if it is zero across *all* P25.c2 samples
p25_all_zero <- rowSums(filtered2$new_counts[, is_p25, drop = FALSE] > 0) == 0

# 4. Keep genes that have **some** counts in BOTH strains
keep_genes <- !(dc_all_zero | p25_all_zero)

table(keep_genes)  # quick sanity check

# 5. Filter counts (and metadata if you want to keep them in sync)
counts_filtered <- filtered2$new_counts[keep_genes, ]

# If you want a filtered version of the full object:
filtered2_comp <- filtered2
filtered2_comp$new_counts <- counts_filtered

# OK now let's do the regression modeling
design2 <- model.matrix(~ treatment * time_point * plant, data = filtered2$new_meta)
dge2_filtered <-  DGEList(counts = filtered2_comp$new_counts)
dge2 <- calcNormFactors(dge2_filtered, method = "TMM")

v2 <- voomWithQualityWeights(dge2, design2, plot=TRUE)
fit2 <- lmFit(v2, design2)
fit2 <- eBayes(fit2)

#################
# OK I was going crazy before with an issue. 
# 1 Get the fold change values for specific coefficients in the model. 11/2025
tt_dc3000_time <- topTable(fit2,coef = "time_pointt3",number = Inf,sort.by = "none")
tt_strain_time_col0 <- topTable(fit2, coef="treatmentp25c2:time_pointt3",number = Inf,sort.by = "none")
tt_plant_affects_time <- topTable(fit2, coef="time_pointt3:plantey15_2",number = Inf,sort.by = "none")
tt_plant_strain_specific <- topTable(fit2, coef="treatmentp25c2:time_pointt3:plantey15_2",number = Inf,sort.by = "none")

# # 2. Build one combined table
result_df <- data.frame(
  gene = rownames(tt_strain_time_col0),
  
  # logFCs
  logFC_DC3000_time           = tt_dc3000_time$logFC,
  logFC_strain_time_col0      = tt_strain_time_col0$logFC,
  logFC_plant_affects_time    = tt_plant_affects_time$logFC,
  logFC_plant_strain_specific = tt_plant_strain_specific$logFC,
  
  # adjusted p-values
  padj_DC3000_time           = tt_dc3000_time$adj.P.Val,
  padj_strain_time_col0      = tt_strain_time_col0$adj.P.Val,
  padj_plant_affects_time    = tt_plant_affects_time$adj.P.Val,
  padj_plant_strain_specific = tt_plant_strain_specific$adj.P.Val,
  
  row.names = NULL
)

# ----- Plot top genes: time-specific and strain-specific -------
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)

# Safety checks
if(!exists("fit2") || !exists("v2") || !exists("result_df")){
  stop("fit2, v2 or result_df not found in the workspace. Run the voom/fit and result table first.")
}

# Use fit2 coefficients/topTables computed already:
# tt_dc3000_time and tt_strain_time_col0 were created earlier in your script.
# If not present, create them:
if(!exists("tt_dc3000_time")){
  tt_dc3000_time <- topTable(fit2, coef = "time_pointt3", number = Inf, sort.by = "none")
}
if(!exists("tt_strain_time_col0")){
  tt_strain_time_col0 <- topTable(fit2, coef = "treatmentp25c2:time_pointt3", number = Inf, sort.by = "none")
}

# Choose top 10 by adjusted p-value (smallest adj.P.Val)
top_time_genes   <- rownames(tt_dc3000_time[order(tt_dc3000_time$adj.P.Val), ])[1:10]
top_strain_genes <- rownames(tt_strain_time_col0[order(tt_strain_time_col0$adj.P.Val), ])[1:10]

# Remove any NAs if fewer than 10 available
top_time_genes   <- top_time_genes[!is.na(top_time_genes)]
top_strain_genes <- top_strain_genes[!is.na(top_strain_genes)]

# Print selected gene lists
message("Top time-specific genes (T3 vs T0):\n", paste(top_time_genes, collapse = ", "))
message("Top strain-specific-time interaction genes:\n", paste(top_strain_genes, collapse = ", "))

# Compute fraction of DC3000 time-specific genes that are strain-specific
# Use your result_df with padj columns
if(!"padj_DC3000_time" %in% colnames(result_df) || !"padj_strain_time_col0" %in% colnames(result_df)){
  stop("result_df missing expected padj columns. Make sure result_df contains padj_DC3000_time and padj_strain_time_col0.")
}
dc_time_all <- result_df %>% filter(padj_DC3000_time <= 0.05)
n_dc_time <- nrow(dc_time_all)
n_dc_time_and_strain <- sum(dc_time_all$padj_strain_time_col0 <= 0.05, na.rm = TRUE)
pct_bg_specific <- if(n_dc_time>0) 100 * n_dc_time_and_strain / n_dc_time else NA

message("Number DC3000 time-specific genes (padj <= 0.05): ", n_dc_time)
message("Of those, number with significant strain×time (padj <= 0.05): ", n_dc_time_and_strain)
message(sprintf("Percent that are genetic-background-specific: %.1f%%", pct_bg_specific))

# ---- Prepare expression table (log2-CPM) -------
# v2$E is voom log2-CPM matrix with genes as rows and samples as columns.
expr_mat <- v2$E
# Check gene names are rows
if(is.null(rownames(expr_mat))){
  stop("v2$E has no rownames (genes).")
}

# Make sample metadata: ensure it matches columns of expr_mat
meta_plot <- filtered2_comp$new_meta  # use filtered2_comp metadata used for fit2
# If your metadata has rownames as sample names, ensure they match colnames of expr_mat
if(!all(colnames(expr_mat) %in% rownames(meta_plot))){
  # try metadata assigned earlier named 'metadata' maybe:
  if(all(colnames(expr_mat) %in% rownames(metadata))){
    meta_plot <- metadata[colnames(expr_mat), ]
  } else {
    stop("Sample names in v2$E do not match rownames in metadata. Fix sample naming before plotting.")
  }
} else {
  meta_plot <- meta_plot[colnames(expr_mat), ]
}

# Create plotting function that builds a tidy data.frame for selected genes
make_tidy_for_genes <- function(gene_vector){
  genes_present <- gene_vector[gene_vector %in% rownames(expr_mat)]
  if(length(genes_present)==0) stop("No selected genes found in voom expression matrix.")
  sub_expr <- expr_mat[genes_present, , drop=FALSE]
  # transpose and melt
  sub_df <- as.data.frame(t(sub_expr))
  sub_df$Sample <- rownames(sub_df)
  sub_df <- cbind(sub_df, meta_plot[rownames(sub_df), c("treatment","plant","time_point")])
  long <- reshape2::melt(sub_df, id.vars = c("Sample","treatment","plant","time_point"),
                         variable.name = "gene", value.name = "log2CPM")
  # ensure factors
  long$time_point <- factor(long$time_point, levels = c("t0","t3"))
  long$treatment  <- factor(long$treatment)
  long$plant      <- factor(long$plant)
  return(long)
}

# Prepare tidy data
tidy_time   <- make_tidy_for_genes(top_time_genes)
tidy_strain <- make_tidy_for_genes(top_strain_genes)

# Summarise per group: mean and sd (treatment x time_point); keep plant as facet or color if you prefer
summarise_group <- function(tidy_df){
  tidy_df %>%
    group_by(gene, treatment, plant, time_point) %>%
    summarise(mean = mean(log2CPM, na.rm=TRUE),
              sd   = sd(log2CPM, na.rm=TRUE),
              n    = n(), .groups = "drop") %>%
    mutate(se = sd / sqrt(pmax(n,1)))
}

sum_time   <- summarise_group(tidy_time)
sum_strain <- summarise_group(tidy_strain)

# Plotting function: facet by gene, color by treatment (strain), shape by plant
plot_gene_set <- function(summary_df, raw_df, title_str){
  p <- ggplot(summary_df, aes(x = time_point, y = mean, group = treatment, color = treatment)) +
    geom_line(aes(linetype = treatment), position = position_dodge(width=0.1)) +
    geom_point(data = raw_df, aes(x = time_point, y = log2CPM, color = treatment, shape = plant),
               alpha = 0.5, size = 1, position = position_jitter(width = 0.05, height = 0)) +
    geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.12, position = position_dodge(width=0.1)) +
    facet_wrap(~ gene, scales = "free_y", ncol = 5) +
    theme_bw(base_size = 11) +
    labs(title = title_str,
         x = "Time point",
         y = "log2-CPM (voom)",
         color = "Strain",
         shape = "Plant genotype") +
    theme(strip.text = element_text(size=8),
          axis.text.x = element_text(size=9))
  return(p)
}

p_top_time   <- plot_gene_set(sum_time, tidy_time, "Top 10 genes: time-specific (T3 vs T0)")
p_top_strain <- plot_gene_set(sum_strain, tidy_strain, "Top 10 genes: strain-specific time effect")

# Save plots
outdir <- "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/top10_time_specific_genes_voom_log2CPM.pdf", width = 12, height = 8)
print(p_top_time)
dev.off()

pdf( "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/top10_strain_time_interaction_genes_voom_log2CPM.pdf", width = 12, height = 8)
print(p_top_strain)
dev.off()

# Optional: combined figure (stacked)
combined <- plot_grid(p_top_time + theme(legend.position = "bottom"),
                      p_top_strain + theme(legend.position = "bottom"),
                      ncol = 1, labels = c("A","B"), rel_heights = c(1,1))
pdf( "top10_time_vs_strain_time_combined.pdf", width = 12, height = 16)
print(combined)
dev.off()

message("Plots written to: ", outdir)


# Plot logFC p25.c2 time vs DC3000 time
df_time <- data.frame(
  gene = rownames(fit2$coefficients),
  
  # Time effect (T3–T0) for DC3000 in col-0
  logFC_DC3000 = fit2$coefficients[, "time_pointt3"],
  
  # Time effect (T3–T0) for P25.c2 in col-0
  logFC_P25c2 = fit2$coefficients[, "time_pointt3"] +
    fit2$coefficients[, "treatmentp25c2:time_pointt3"]
)

# Plot logFC  DC3000 in Eyach vs Col-0
df_plant <- data.frame(
  gene = rownames(fit2$coefficients),
  
  # Time effect (T3–T0) for DC3000 in col-0
  logFC_DC3000_col0 = fit2$coefficients[, "time_pointt3"],
  
  # Time effect (T3–T0) for P25.c2 in col-0
  logFC_DC3000_Eyach = fit2$coefficients[, "time_pointt3"] +
    fit2$coefficients[, "time_pointt3:plantey15_2"]
)


# Plot lfc p25.c2 and ldf DC3000
# Compute correlation
cor_val <- round(cor(df_time$logFC_DC3000, df_time$logFC_P25c2, use = "complete.obs"),2)
cor_val2 <- round(cor(df_plant$logFC_DC3000_col0, df_plant$logFC_DC3000_Eyach, use = "complete.obs"),2)
# Make the plot
strain_plot <- ggplot(df_time, aes(x = logFC_DC3000, y = logFC_P25c2)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("R = ", round(cor_val, 3)),
    hjust = 1.1, vjust = -1.1,
    size = 5
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Time response (T3–T0): DC3000 vs P25.c2",
    x = "DC3000 logFC (T3–T0)",
    y = "P25.c2 logFC (T3–T0)"
  )

plant_genot_plot <- ggplot(df_plant, aes(x = logFC_DC3000_col0, y = logFC_DC3000_Eyach)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("R = ", round(cor_val2, 3)),
    hjust = 1.1, vjust = -1.1,
    size = 5
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Time response (T3–T0): DC3000 vs P25.c2",
    x = "DC3000 in Col-0 logFC (T3–T0)",
    y = "DC3000 in Eyach logFC (T3–T0)"
  )

pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/correspondence_Strain_Plant.pdf")
plot_grid(strain_plot, plant_genot_plot)
dev.off()

# So now we have a data table result_df that we can use to look at the breakdown in significance
write.table(result_df, file.path(outdir, "sig_results_11_2025.csv"), sep=",", quote= FALSE )

# I want to make a stacked barplot for the genes associated with strain specific background vs other variables
# I want a stacked barpolot using the result_df dataframe, with three bars -- one for Strain, one for Plant and one Gene-by-Gene. In each column I want to ask whether all of the time-specific genes in DC3000 which fraction show a strain effect, which fraction show a Plant genotype effect and which show a strainxplant effect



# 1. Restrict to time-specific DC3000 genes
rownames(result_df)<-result_df$gene
dc_time <- result_df %>% 
  filter(padj_DC3000_time <= 0.05)


# 2. Summarise fractions for each effect type
# 2. Build a small summary table of fractions
#    (among *time-specific* DC3000 genes)
frac_df <- tibble(
  Bar  = c("Strain", "Plant", "Gene-by-Gene"),
  With = c(
    mean(dc_time$padj_strain_time_col0      <= 0.05, na.rm = TRUE),  # strain effect
    mean(dc_time$padj_plant_affects_time    <= 0.05, na.rm = TRUE),  # plant effect
    mean(dc_time$padj_plant_strain_specific <= 0.05, na.rm = TRUE)   # strain×plant
  )
) %>%
  mutate(Without = 1 - With) %>%
  pivot_longer(
    cols = c(With, Without),
    names_to = "EffectPresence",
    values_to = "Fraction"
  ) %>%
  mutate(
    EffectPresence = dplyr::recode(
      EffectPresence,
      With    = "Has effect",
      Without = "No effect"
    ),
    Bar = factor(Bar, levels = c("Strain", "Plant", "Gene-by-Gene"))
  )

# 3. Stacked barplot
stacked <- ggplot(frac_df, aes(x = Bar, y = Fraction, fill = EffectPresence)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Has effect" = "#3B4F8C",   # deep desaturated blue
      "No effect"  = "#C7C9CC"    # soft neutral gray
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "% of time-specific DC3000 genes",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/stacked_barplot.pdf")
stacked
dev.off()


# I want to compare (overlay) histograms for the LFC values for DC3000 vs p25.c2

library(ggplot2)
library(dplyr)
library(tidyr)


# Recreate df_time if missing (uses fit2)
if(!exists("df_time")){
  if(!exists("fit2")) stop("fit2 not found. Run voom/lmFit first.")
  df_time <- data.frame(
    gene = rownames(fit2$coefficients),
    logFC_DC3000 = fit2$coefficients[, "time_pointt3"],
    logFC_P25c2  = fit2$coefficients[, "time_pointt3"] + fit2$coefficients[, "treatmentp25c2:time_pointt3"]
  )
}

# Long format so ggplot makes a legend
df_long <- df_time %>%
  select(gene, logFC_DC3000, logFC_P25c2) %>%
  pivot_longer(cols = starts_with("logFC_"),
               names_to = "group",
               values_to = "log2FC") %>%
  mutate(group = dplyr::recode(group,
                               logFC_DC3000 = "DC3000",
                               logFC_P25c2  = "P25.c2")) %>%
  filter(!is.na(log2FC))

# Colors
pal <- c("DC3000" = "#3B4F8C", "P25.c2" = "#D95F02")

# Overlaid histograms + density lines, with legend
hist_density_plot <- ggplot(df_long, aes(x = log2FC, fill = group, color = group)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.25, bins = 60, color = NA) +
  geom_density(size = 0.9) +
  scale_fill_manual(values = pal, name = "Strain") +
  scale_color_manual(values = pal, guide = FALSE) + # density lines colored; hide duplicate legend
  labs(title = "Distribution of log2FC (T3 vs T0): DC3000 vs P25.c2",
       x = "log2 fold change (T3 - T0)",
       y = "Density") +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = "white", color = NA))

# show and optionally save
print(hist_density_plot)

# save (adjust path if needed)
outdir <- "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/plots"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
ggsave(filename = file.path(outdir, "hist_density_overlay_with_legend.pdf"),
       plot = hist_density_plot, width = 8, height = 5)
message("Saved to: ", file.path(outdir, "hist_density_overlay_with_legend.pdf"))


##############
# In this next section we will focus specifically on genes that shot time specific effects in DC3000

#I want to Graph ten genes that show time specific effects for DC3000
#I want to Graph ten genes that show strain-specific effects for DC3000

