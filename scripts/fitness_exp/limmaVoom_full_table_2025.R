# May 2025 This script takes the full counts matrix from the many trials and will do subsetting and limma voom analysis. The goal is to ask what percentage of sig genes are genetic background specific in their behavior. In August 2025 Effie determined that the original counts may have been problematic and redid the counts table.
# oh no, that table must messed up. The final counts table is in 10/2025
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
# Extract top genes for the interaction term
top_interaction <- topTable(fit, coef = "treatmentp25c2:time_pointt3", number = Inf, adjust.method = "BH")

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
})





# setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")
setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output")

##########################################################################
### Read in full counts table ###
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
rownames(counts) <- all_exp$Sample

#remove samples with few reads
keep <- which(rowSums(counts, na.rm = TRUE)>=25000)
counts <- counts[keep,]
metadata_keep = metadata[keep,]
metadata <- metadata_keep
table(metadata[,c("experiment", "treatment", "time_point")])

#counts matrix with just orthologs. Get rid of lines that dont have shared orthologs. replace NAs with 0s
ortholog_tab <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", header = TRUE, row.names = 1, sep = ",")
maintain <- c(unique(ortholog_tab[,3]))
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

# OK now let's do 
design2 <- model.matrix(~ treatment * time_point * plant, data = filtered2$new_meta)
#dge2 <- DGEList(counts = filtered2$new_counts)
#keep <- filterByExpr(dge2, design2)
dge2_filtered <-  DGEList(counts = filtered2_comp$new_counts)
dge2 <- calcNormFactors(dge2_filtered, method = "TMM")

v2 <- voomWithQualityWeights(dge2, design2, plot=TRUE)
fit2 <- lmFit(v2, design2)
fit2 <- eBayes(fit2)

#################
# OK I am going crazy with issue with the 
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


# 3. View it
head(result_df)

# Plot logFC p25.c2 time vs DC3000 time
df_time <- data.frame(
  gene = rownames(fit2$coefficients),
  
  # Time effect (T3–T0) for DC3000 in col-0
  logFC_DC3000 = fit2$coefficients[, "time_pointt3"],
  
  # Time effect (T3–T0) for P25.c2 in col-0
  logFC_P25c2 = fit2$coefficients[, "time_pointt3"] +
    fit2$coefficients[, "treatmentp25c2:time_pointt3"]
)

# Plot lfc p25.c2 and ldf DC3000
ggplot(df_time, aes(x = logFC_DC3000, y = logFC_P25c2)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_bw(base_size = 13) +
  labs(
    title = "Time response (T3–T0): DC3000 vs P25.c2",
    x = "DC3000 logFC (T3–T0)",
    y = "P25.c2 logFC (T3–T0)"
  )

# So now we have a data table result_df that we can use to look at the breakdown in significance


# I want to make a stacked barplot for the genes associated with strain specific background vs other variables


# I want to compare the LFC values for DC3000 vs p25.c2












