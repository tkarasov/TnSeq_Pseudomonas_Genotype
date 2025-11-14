# May 2025 This script takes the full counts matrix from the many trials and will do subsetting and Deseq analysis.I thought I would do deseq but now changing to edgeR

library(edgeR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(randomForest)
library(caret)
library(matrixStats)
library(DESeq2)

setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")


### Read in full counts table ###
all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_5_2025.rds")
all_exp$Sample <- with(all_exp, paste(treatment, plant, time_point, position, experiment, sep = "_"))
# Move Sample column to the front
all_exp <- all_exp[, c("Sample", setdiff(names(all_exp), "Sample"))]
all_exp$plant <- tolower(all_exp$plant)
### Reassign control: when we originally did the experiment it turned out that pure culture at T0 was a more robust metric than the small amount of microbe on the plant at t0. So we need to split the pure culture controls in half and randomly assign as t0. THis should be done after removinga any actual plantxT0

# Filter: remove rows where plant ≠ Ctrl/ctrl AND time_point == "t0"
all_exp_filtered <- subset(all_exp, !(tolower(plant) != "ctrl" & time_point == "t0"))

# now for the plant controls we want to randomly assign to either ey15_2 or col_0
all_exp <- all_exp_filtered
set.seed(444)  # for reproducibility
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
keep <- which(rowSums(counts, na.rm = TRUE)>=50000)
counts <- counts[keep,]
metadata_keep = metadata[keep,]
metadata <- metadata_keep

#counts matrix with just orthologs. Get rid of lines that dont have shared. replace NAs with 0s
ortholog_tab <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", header = TRUE, row.names = 1, sep = ",")
maintain <- c(unique(ortholog_tab[,3]))
maintain <- maintain[maintain!=""]
#count_ortholog <- counts[,grepl("^WP", colnames(counts))]
count_ortholog<- counts[,maintain]
count_ortholog[is.na(count_ortholog)] <- 0
count_ortholog <- count_ortholog[, colSums(count_ortholog) > 0]
count_ortholog <- t(count_ortholog)
counts <- count_ortholog


### Step 2: Create interaction factors###
metadata$strain <- factor(metadata$treatment)
metadata$plant <- factor(metadata$plant)
group <- interaction(metadata$strain, metadata$plant, metadata$time_point)


design <- model.matrix(~ 0 + group)  
colnames(design) <- levels(group)



######## Do edgeR with ortholog comparison
# Step 1: Match sample order between counts and metadata
#metadata <- metadata[match(colnames(counts), metadata$sample_id), ]

# Step 2: Create DGEList object
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# Estimate dispersions
dge <- estimateDisp(dge, design)




# Step 3: Filter low-expression genes
#keep <- filterByExpr(dge, group = dge$samples$group)
#dge <- dge[keep, , keep.lib.sizes = FALSE]

# Step 4: Normalize
dge <- calcNormFactors(dge)

# Check: offsets should now be finite
dge$offset <- log(dge$samples$lib.size * dge$samples$norm.factors)
summary(dge$offset)






### Create design matrix ###

###edgeR###
get_fitted_values <- function(fit, gene_name, design) {
  # Extract beta coefficients for this gene
  beta <- fit$coefficients[gene_name, ]
  
  # Get all unique combinations of design covariates
  design_matrix <- design
  
  # Create a prediction matrix for each group
  # Assumes design has no intercept (i.e., ~ 0 + group style)
  predicted_logCPM <- as.vector(design_matrix %*% beta)
  
  # Output
  data.frame(
    sample = rownames(design_matrix),
    predicted_logCPM = predicted_logCPM
  )
}

# edgeR does not allow NAs in the counts matrix so I will need to subset to orthologs shared between strains.
dge <- DGEList(counts = count_ortholog)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)

### Define contrast ###
group <- metadata$group
design <- design <- model.matrix(~ treatment * time_point * plant + experiment, data = metadata)
dge <- estimateDisp(dge, design)
### fit the model ###
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = "treatmentp25c2:time_pointt3")
topTags(lrt)


# this lrt_all is the primary model
lrt_all <- lrt

### Build a table for the coefficients of the model
# Get the coefficient names
coef_names <- colnames(design)

# Create a data frame to store results
results_summary <- data.frame(
  term = coef_names,
  num_significant = NA_integer_
)

# Loop through each coefficient and run a glmLRT
for (i in seq_along(coef_names)) {
  lrt <- glmLRT(fit, coef = i)
  res <- topTags(lrt, n = Inf)$table
  results_summary$num_significant[i] <- sum(res$FDR < 0.01)
}

# Show the summary table
print(results_summary)


#### Seems to be a strong experimental effect

# Transform counts with voom
v <- voom(dge, design, plot = TRUE)

# Use duplicateCorrelation to estimate variance from random effect
corfit <- duplicateCorrelation(v, design, block = metadata$experiment)

# Re-voom using the estimated correlation
v <- voom(dge, design, block = metadata$experiment, correlation = corfit$consensus)

# Fit model with block structure
fit <- lmFit(v, design, block = metadata$experiment, correlation = corfit$consensus)

# Empirical Bayes moderation
fit <- eBayes(fit)

# Extract top genes for the interaction term
top_interaction <- topTable(fit, coef = "treatmentp25c2:time_pointt3", number = Inf)

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

# Combine into a dataframe
results_summary <- do.call(rbind, results_summary)

# Print and export to CSV
print(results_summary)
write.csv(results_summary, "voom_limma_results_summary.csv", row.names = FALSE)


# Perform PCA on voom-normalized expression values
pca <- prcomp(t(v$E), scale. = TRUE)

# Build PCA dataframe with sample metadata
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  experiment = metadata$experiment,
  treatment = metadata$treatment,
  time_point = metadata$time_point,
  plant = metadata$plant,
  sample = metadata$Sample  # optional
)

# PCA plot: color by experiment, shape by time_point
library(ggplot2)


ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment, shape = time_point)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Voom-normalized expression",
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100), "%)")
  ) +
  theme(panel.grid = element_blank())





# Get voom-normalized expression matrix
expr <- v$E  # genes x samples

# Add experiment info
sample_exp <- metadata$experiment
colnames(expr) <- sample_exp

# Average expression per gene in each experiment
gene_means_by_exp <- t(apply(expr, 1, function(row) tapply(row, sample_exp, mean, na.rm = TRUE)))
gene_means_by_exp <- as.data.frame(gene_means_by_exp)

# Rename columns if needed (assumes two experiments)
colnames(gene_means_by_exp) <- paste0("mean_expr_", colnames(gene_means_by_exp))

# Plot correlation
cor_val <- cor(gene_means_by_exp$mean_expr_exp_0001, gene_means_by_exp$mean_expr_exp_0002, method = "pearson")
library(ggplot2)
ggplot(gene_means_by_exp, aes(x = mean_expr_exp_0001, y = mean_expr_exp_0002)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Gene-wise mean expression between experiments",
    x = "Mean expression (Experiment 1)",
    y = "Mean expression (Experiment 2)"
  ) +
  coord_equal() +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.1,
           label = paste("Pearson r =", round(cor_val, 3)), size = 5)

####### LFC correlations
# Function to run DiD analysis per experiment
run_did <- function(counts, metadata, experiment_name) {
  # Subset
  meta_sub <- metadata[metadata$experiment == experiment_name, ]
  counts_sub <- counts[, meta_sub$Sample]
  
  # Relevel factors
  meta_sub$treatment <- factor(meta_sub$treatment, levels = c("dc3000", "p25c2"))
  meta_sub$time_point <- factor(meta_sub$time_point, levels = c("t0", "t3"))
  meta_sub$plant <- factor(meta_sub$plant)
  
  # Design matrix for DiD
  design <- model.matrix(~ treatment * time_point + plant, data = meta_sub)
  
  # Voom + limma fit
  dge <- DGEList(counts = counts_sub)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Extract DiD effect
  did_results <- topTable(fit, coef = "treatmentp25c2:time_pointt3", number = Inf, sort.by = "none")
  return(did_results)
}

# Run for both experiments
did_exp1 <- run_did(counts, metadata, "exp_0001")
did_exp2 <- run_did(counts, metadata, "exp_0002")

# Ensure gene identifiers are rownames in both
common_genes <- intersect(rownames(did_exp1), rownames(did_exp2))

# Build dataframe for plotting
did_logFC_compare <- data.frame(
  gene = common_genes,
  logFC_exp1 = did_exp1[common_genes, "logFC"],
  logFC_exp2 = did_exp2[common_genes, "logFC"]
)

# Calculate Pearson correlation
cor_val <- cor(did_logFC_compare$logFC_exp1, did_logFC_compare$logFC_exp2, method = "pearson")

# Plot with ggplot2
ggplot(did_logFC_compare, aes(x = logFC_exp1, y = logFC_exp2)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Difference-in-Differences LogFC: Experiment 1 vs 2",
    x = "DiD logFC (Experiment 1)",
    y = "DiD logFC (Experiment 2)",
    subtitle = paste("Pearson r =", round(cor_val, 3))
  )

# Extract top tables for the same contrast (e.g., time_pointt3)
tt1 <- topTable(fit_exp1, coef = "time_pointt3", number = Inf, sort.by = "none")
tt2 <- topTable(fit_exp2, coef = "time_pointt3", number = Inf, sort.by = "none")

# Match on gene names
common_genes <- intersect(rownames(tt1), rownames(tt2))

# Create dataframe of logFCs
logfc_compare <- data.frame(
  gene = common_genes,
  logFC_exp1 = tt1[common_genes, "logFC"],
  logFC_exp2 = tt2[common_genes, "logFC"]
)



ggplot(logfc_compare, aes(x = logFC_exp1, y = logFC_exp2)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  labs(
    title = "LogFC comparison for time_pointt3 across experiments",
    x = "LogFC (Experiment 1)",
    y = "LogFC (Experiment 2)"
  ) +
  coord_equal()


############ Make subsets and do statistics independently


rows_dc3000_ey15_2 <- which(metadata$treatment == "dc3000" & metadata$plant == "ey15_2")
counts_dc3000_eyach <- counts[,rows_dc3000_ey15_2]
metadata_dc3000_eyach <- metadata[rows_dc3000_ey15_2,]


rows_dc3000_col <- which(metadata$treatment == "dc3000" & metadata$plant == "col_0")
counts_dc3000_col <- counts[,rows_dc3000_col]
metadata_dc3000_col <- metadata[rows_dc3000_col,]


rows_p25_ey15_2 <- which(metadata$treatment == "p25c2" & metadata$plant == "ey15_2")
counts_p25_eyach <- counts[,rows_p25_ey15_2]
metadata_p25_eyach <- metadata[rows_p25_ey15_2,]


rows_p25_col <- which(metadata$treatment == "p25c2" & metadata$plant == "col_0")
counts_p25_col <- counts[,rows_p25_col]
metadata_p25_col <- metadata[rows_p25_col,]

run_edgeR <- function(metadata, group_column, counts, output_prefix = "edgeR_result") {
  library(edgeR)
  # 3. Create group factor
  group <- factor(metadata$time_point)
  
  # 4. Create DGEList
  dge <- DGEList(counts = (counts))  # genes as rows, samples as columns
  dge$samples$group <- group
  
  # 5. Filter lowly expressed genes
  #keep <- filterByExpr(dge, group = group)
  #dge <- dge[keep,, keep.lib.sizes = FALSE]
  
  # 6. Normalize and estimate dispersion
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group)
  dge <- estimateDisp(dge, design)
  
  # 7. Fit model and run test
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit)
  
  # 8. Save results
  res <- topTags(lrt, n = Inf)$table
  out_file <- paste0(output_prefix, ".csv")
  write.csv(res, out_file, row.names = TRUE)
  
  return(lrt)
}


results_dc3000_col <- run_edgeR(
  metadata = metadata_dc3000_col,
  counts = counts_dc3000_col,
  group_column = "time_point",
  output_prefix = "dc3000_col"
)


results_dc3000_eyach <- run_edgeR(
  metadata = metadata_dc3000_eyach,
  counts = counts_dc3000_eyach,
  group_column = "time_point",
  output_prefix = "dc3000_eyach"
)

results_p25_col <- run_edgeR(
  metadata = metadata_p25_col,
  counts = counts_p25_col,
  group_column = "time_point",
  output_prefix = "p25_col"
)

results_p25_eyach <- run_edgeR(
  metadata = metadata_p25_eyach,
  counts = counts_p25_eyach,
  group_column = "time_point",
  output_prefix = "p25_eyach"
)


plot(results_dc3000_col$table$logFC, results_p25_col$table$logFC)
plot(results_dc3000_col$table$logFC, results_dc3000_eyach$table$logFC)
plot(results_p25_col$table$logFC, results_dc3000_eyach$table$logFC)
plot(results_p25_col$table$logFC, results_p25_eyach$table$logFC)

# combine results into one dataframe
combined_logFC <- data.frame(
  gene = rownames(results_dc3000_col$table),
  dc3000_col0   = results_dc3000_col$table$logFC,
  dc3000_ey15_2 = results_dc3000_eyach$table$logFC,
  p25c2_col0    = results_p25_col$table$logFC,
  p25c2_ey15_2  = results_p25_eyach$table$logFC,
  pval_fullmodel = lrt_all$table$PValue,
  fdr_fullmodel = p.adjust(lrt_all$table$PValue, method = "BH")
)

#graph

library(GGally)
library(dplyr)
df = combined_logFC
df <- df %>%
  mutate(significant = fdr_fullmodel < 0.01)
# Select only the result columns
result_cols <- df[,c("dc3000_col0", "dc3000_ey15_2", "p25c2_col0",  "p25c2_ey15_2")]

# Bind the significance column for coloring
plot_data <- bind_cols(result_cols, significant = df$significant)

# Pairwise plot colored by FDR significance
p <- ggpairs(plot_data, 
        aes(color = significant, alpha = 0.6),
        upper = list(continuous = wrap("points", size = 1.5))) + theme_classic(base_size = 12)
# Add color scale
p <- p + scale_color_manual(
  values = c("TRUE" = "navy", "FALSE" = "orange"),
  labels = c("FALSE" = "Not significant", "TRUE" = "FDR < 0.01"),
  name = "Significance"
)

#make a pdf of this plot
pdf("pairwise_comparison_logFC_effect.pdf", width = 10, height = 10, font = "ArialMT")
p
dev.off()

# Compute sign for each logFC
combined_logFC$sign_dc3000_col0   <- sign(combined_logFC$dc3000_col0)
combined_logFC$sign_p25c2_col0    <- sign(combined_logFC$p25c2_col0)
combined_logFC$sign_dc3000_ey15_2 <- sign(combined_logFC$dc3000_ey15_2)
combined_logFC$sign_p25c2_ey15_2  <- sign(combined_logFC$p25c2_ey15_2)

# Find genes with a sign flip between strains in either plant
interaction_genes <- combined_logFC[
  ((combined_logFC$sign_dc3000_col0 != 0 & combined_logFC$sign_p25c2_col0 != 0 &
      combined_logFC$sign_dc3000_col0 != combined_logFC$sign_p25c2_col0) |
     (combined_logFC$sign_dc3000_ey15_2 != 0 & combined_logFC$sign_p25c2_ey15_2 != 0 &
        combined_logFC$sign_dc3000_ey15_2 != combined_logFC$sign_p25c2_ey15_2)),
]
interaction_genes_sig <- interaction_genes[interaction_genes$fdr_fullmodel < 0.05, ]
# 85 genes show sign epistasis

# Okay time to make a table. 



