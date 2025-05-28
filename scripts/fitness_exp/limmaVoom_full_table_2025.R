# May 2025 This script takes the full counts matrix from the many trials and will do subsetting and limma voom analysis. The goal is to ask what percentage of sig genes are genetic background specific in their behavior.

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
  # Empirical Bayes moderation
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
  new_count <- as.data.frame(counts) %>% select(-all_of(subset_meta))
  return(list(new_meta=new_meta, new_counts=new_count))
}

summarize_overlap_between_terms <- function(fit, coef1, coef2, pval_cutoff = 0.05) {
  # Get topTable results for both coefficients
  tt1 <- topTable(fit, coef = coef1, number = Inf, p.value = pval_cutoff)
  tt2 <- topTable(fit, coef = coef2, number = Inf, p.value = pval_cutoff)
  
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
######### First do the analysis with all
# design matrix
design <- model.matrix(~ treatment * time_point * plant + experiment, data = metadata)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
res <- run_voom(dge, design, experiment = TRUE)
v <- res$v
fit <- res$fit
result_full <- summarize_overlap_between_terms(
  fit,
  coef1 = "time_pointt3",
  coef2 = "treatmentp25c2:time_pointt3"
)

######### Second do the analysis with Exp 2
filtered2 <- filter_meta_counts("experiment", "exp_0001", counts, metadata)
design2 <- model.matrix(~ treatment * time_point * plant, data = filtered2$new_meta)
dge2 <- DGEList(counts = filtered2$new_counts)
dge2 <- calcNormFactors(dge2)
res <- run_voom(dge2, design2, experiment= FALSE)
v2 <- res$v
fit2 <- res$fit

######### Now compare the results between the full and fit2
result_full <- summarize_overlap_between_terms(
  fit2,
  coef1 = "time_pointt3",
  coef2 = "treatmentp25c2:time_pointt3"
)


# # let's make fit1. No this doesn't work because don't have good controls in experiment 1
# filtered1 <- filter_meta_counts("experiment", "exp_0002", counts, metadata)
# design1 <- model.matrix(~ treatment * time_point * plant, data = filtered1$new_meta)
# dge1 <- DGEList(counts = filtered2$new_counts)
# dge1 <- calcNormFactors(dge1)
# res1 <- run_voom(dge, design1, experiment= FALSE)
# v1 <- res1$v
# fit1 <- res1$fit

######### The first experiment was kind of a mess, so I don't want to do it specifically 
# But let's look at the correspondence in lfc values between the two.
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

# Assuming you ran voom with a design like this:
# ~ treatment * time_point * plant + experiment
design1 <- model.matrix(~ treatment * time_point * plant + experiment, data = metadata)
res <- run_voom(dge, design1, experiment = FALSE)

# Extract LFC for interaction term "time_pointt3" for each experiment
lfc_timepoint <- extract_timepoint_logfc_per_experiment(
  voom_obj = res$v,
  metadata = metadata,
  experiment_col = "experiment",
  design_formula = ~ treatment * time_point * plant
)

# lfc_timepoint should look like: data.frame(gene_id × experiments)


# === Parameters ===
coef_name <- "time_pointt3"
fdr_cutoff <- 0.01

# === 1. Identify significant genes ===
sig_full <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
sig_full_ids <- rownames(sig_full)[sig_full$adj.P.Val < fdr_cutoff]

sig_exp2 <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")
sig_exp2_ids <- rownames(sig_exp2)[sig_exp2$adj.P.Val < fdr_cutoff]

# === 2. Merge logFC and label categories ===
sig_union <- union(sig_full_ids, sig_exp2_ids)
lfc_sig <- lfc_timepoint[rownames(lfc_timepoint) %in% sig_union, ]
lfc_sig <- na.omit(lfc_sig)

lfc_sig$status <- "Not significant"
lfc_sig$status[rownames(lfc_sig) %in% sig_full_ids] <- "Full model only"
lfc_sig$status[rownames(lfc_sig) %in% sig_exp2_ids] <- "Exp2 only"
lfc_sig$status[rownames(lfc_sig) %in% intersect(sig_full_ids, sig_exp2_ids)] <- "Both"
lfc_sig$status <- factor(lfc_sig$status, levels = c("Both", "Full model only", "Exp2 only"))

# === 3. Subset for "Both" genes ===
lfc_both <- lfc_sig[lfc_sig$status == "Both", ]

# Compute Pearson correlation for "Both"
cor_both <- cor.test(lfc_both$logFC_exp_0001, lfc_both$logFC_exp_0002)
R_label <- paste0("R = ", round(cor_both$estimate, 2))

# === 4. Define colorblind-friendly palette ===
color_map <- c(
  "Both" = "#0072B2",           # Blue
  "Full model only" = "#E69F00",  # Orange
  "Exp2 only" = "#CC79A7" )        # Magenta
  
# Load required packages
  library(limma)
  library(ggplot2)
  
  # === Parameters ===
  coef_name <- "time_pointt3"
  fdr_cutoff <- 0.01
  
  # === 1. Identify significant genes ===
  sig_full <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
  sig_full_ids <- rownames(sig_full)[sig_full$adj.P.Val < fdr_cutoff]
  
  sig_exp2 <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")
  sig_exp2_ids <- rownames(sig_exp2)[sig_exp2$adj.P.Val < fdr_cutoff]
  
  # === 2. Merge logFC and label significance categories ===
  sig_union <- union(sig_full_ids, sig_exp2_ids)
  lfc_sig <- lfc_timepoint[rownames(lfc_timepoint) %in% sig_union, ]
  lfc_sig <- na.omit(lfc_sig)
  
  lfc_sig$status <- "Not significant"
  lfc_sig$status[rownames(lfc_sig) %in% sig_full_ids] <- "Full model only"
  lfc_sig$status[rownames(lfc_sig) %in% sig_exp2_ids] <- "Exp2 only"
  lfc_sig$status[rownames(lfc_sig) %in% intersect(sig_full_ids, sig_exp2_ids)] <- "Both"
  lfc_sig$status <- factor(lfc_sig$status, levels = c("Both", "Full model only", "Exp2 only"))
  
  # === 3. Subset for correlation + regression ===
  lfc_both <- lfc_sig[lfc_sig$status == "Both", ]
  cor_both <- cor.test(lfc_both$logFC_exp_0001, lfc_both$logFC_exp_0002)
  R_label <- paste0("R = ", round(cor_both$estimate, 2))
  
  # === 4. Define colorblind-friendly colors ===
  color_map <- c(
    "Both" = "#0072B2",            # Blue
    "Full model only" = "#E69F00", # Orange
    "Exp2 only" = "#CC79A7"        # Magenta
  )
  
  # === 5. Save PDF ===
  
  library(extrafont)
  
  # Register system fonts (only needs to be done once)
  extrafont::font_import(prompt = FALSE)
  loadfonts(device = "pdf")

pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/logFC_comparison_DC3000.pdf", width = 3.5, height = 3.5, family = "Arial")
  # 89 mm × 89 mm
  
  ggplot(lfc_sig, aes(x = logFC_exp_0001, y = logFC_exp_0002, color = status)) +
    geom_point(alpha = 0.85, size = 1.6) +
    geom_smooth(
      data = lfc_both,
      method = "lm", se = FALSE, color = "gray40", linewidth = 0.8
    ) +
    annotate("text", x = Inf, y = -Inf, label = R_label, hjust = 1.1, vjust = -0.5, size = 3.5, family = "Arial") +
    scale_color_manual(values = color_map) +
    labs(
      x = expression(log[2]*"FC (DC3000, Experiment 1)"),
      y = expression(log[2]*"FC (DC3000, Experiment 2)"),
      color = "Significant in"
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(size = 9),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.position = "right"
    )
  
dev.off()

#### OK now we have a subset of genes we can follow. Those that show a significant fitness effect in DC3000 over time. What I want to do with them is ask questions about how many of them are sensitive to the strain genetic background and how many of them are sensitive to the plant genetic background. And I want good graphics to relay this information.

summarize_background_sensitivity_with_3way <- function(
    fit,
    coef_time = "time_pointt3",
    coef_strain_time = "treatmentp25c2:time_pointt3",
    coef_plant_time = "time_pointt3:plantey15_2",
    coef_3way = "treatmentp25c2:time_pointt3:plantey15_2",
    fdr_cutoff = 0.01,
    plot_type = c("barplot", "upset"),
    plot = TRUE
  ) {
    suppressPackageStartupMessages({
      library(limma)
      library(dplyr)
      library(tidyr)
      library(ggplot2)
    })
    
    plot_type <- match.arg(plot_type)
    upset_available <- requireNamespace("ComplexUpset", quietly = TRUE)
    if (plot_type == "upset" && !upset_available) {
      warning("ComplexUpset not installed, defaulting to barplot.")
      plot_type <- "barplot"
    }
    
    # Extract topTables
    tt_time <- topTable(fit, coef = coef_time, number = Inf, adjust.method = "BH")
    tt_st <- topTable(fit, coef = coef_strain_time, number = Inf, adjust.method = "BH")
    tt_pt <- topTable(fit, coef = coef_plant_time, number = Inf, adjust.method = "BH")
    tt_3w <- topTable(fit, coef = coef_3way, number = Inf, adjust.method = "BH")
    
    sig_genes <- rownames(tt_time)[tt_time$adj.P.Val < fdr_cutoff]
    
    df <- data.frame(
      Gene = sig_genes,
      Time = TRUE,
      Strain_Time = sig_genes %in% rownames(tt_st)[tt_st[sig_genes, "adj.P.Val"] < fdr_cutoff],
      Plant_Time = sig_genes %in% rownames(tt_pt)[tt_pt[sig_genes, "adj.P.Val"] < fdr_cutoff],
      Three_Way = sig_genes %in% rownames(tt_3w)[tt_3w[sig_genes, "adj.P.Val"] < fdr_cutoff]
    )
    
    if (plot) {
      if (plot_type == "barplot") {
        df_counted <- df %>%
          pivot_longer(
            cols = c(Strain_Time, Plant_Time, Three_Way),
            names_to = "Effect",
            values_to = "Significant"
          ) %>%
          filter(Significant == TRUE) %>%
          group_by(Effect) %>%
          summarise(n = n(), .groups = "drop")
        
        total_genes <- nrow(df)  # total number of genes with a significant time effect
        
        p <- ggplot(df_counted, aes(x = Effect, y = n, fill = Effect)) +
          geom_col() +
          geom_hline(yintercept = total_genes, linetype = "dashed", color = "gray40") +
          annotate("text", x = 1.5, y = total_genes + 5, label = paste("Total:", total_genes),
                   hjust = 0, size = 3, family = "Arial", color = "gray20") +
          scale_fill_manual(values = c(
            "Strain_Time" = "#E69F00",
            "Plant_Time" = "#56B4E9",
            "Three_Way" = "#009E73"
          )) +
          labs(
            x = NULL,
            y = "Genes with background-sensitive fitness effect (subset of time-effect genes)",
            title = "Genetic background sensitivity among DC3000 time-responsive genes"
          ) +
          theme_minimal(base_family = "Arial") +
          theme(
            legend.position = "none",
            axis.text = element_text(size = 9),
            axis.title.y = element_text(size = 9)
          )
        print(p)
        
      } else {
        library(ComplexUpset)
        
        # Create the base theme separately
        custom_theme <- theme_minimal(base_family = "Arial") +
          theme(
            axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
            axis.title.x = element_text(size = 9)
          )
        
        # Generate the UpSet plot
        # Generate the UpSet plot
        ComplexUpset::upset(
          df,
          intersect = c("Strain_Time", "Plant_Time", "Three_Way"),
          name = "Background Effect",
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size(
              text_mapping = aes(label = ..count..)  # <- corrected here
            )
          ),
          themes = ComplexUpset::upset_modify_themes(
            theme_minimal(base_family = "Arial") +
              theme(
                axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
                axis.title.x = element_text(size = 9)
              )
          )
        )
       # print(p)
      }
    }
    
    return(df)
  }

gene_summary <- summarize_background_sensitivity_with_3way(
  fit = fit,
  coef_time = "time_pointt3",
  coef_strain_time = "treatmentp25c2:time_pointt3",
  coef_plant_time = "time_pointt3:plantey15_2",
  coef_3way = "treatmentp25c2:time_pointt3:plantey15_2",
  fdr_cutoff = 0.01,
  plot_type = "upset"  # or "upset" if you have ComplexUpset
  )

# the plotting of the upset plot isn't working well. But only the strain genetic background had an effect.
# So I really just want to make a bar plot
tot_sig <- length(which(gene_summary$Time==TRUE))
strain_sig <- length(which(gene_summary$Strain_Time==TRUE))
plant_sig <- length(which(gene_summary$Plant_Time==TRUE))
threeway_sig <- length(which(gene_summary$Three_Way==TRUE))

# Now make a stacked barplot showing the proportion that are significant
# Compute non-significant counts
strain_nonsig <- tot_sig - strain_sig
plant_nonsig <- tot_sig - plant_sig
threeway_nonsig <- tot_sig - threeway_sig

# Build dataframe for plotting
bar_data <- data.frame(
  Background = rep(c("Strain", "Plant", "Gene-by-Gene"), each = 2),
  Significance = rep(c("Significant", "Not significant"), 3),
  Count = c(strain_sig, strain_nonsig, plant_sig, plant_nonsig, threeway_sig, threeway_nonsig)
)

# Factor ordering
bar_data$Background <- factor(bar_data$Background, levels = c("Strain", "Plant", "Gene-by-Gene"))
bar_data$Significance <- factor(bar_data$Significance, levels = c("Not significant", "Significant"))

# Plot
pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/genetic_background_sensitivity_barplot.pdf", width = 3.5, height = 3.5, family = "Arial")

ggplot(bar_data, aes(x = Background, y = Count, fill = Significance)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Not significant" = "gray85", "Significant" = "#0072B2")) +
  labs(
    y = "Number of Time-Sensitive Genes",
    x = NULL,
    title = "Genetic background sensitivity of time-responsive genes"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

dev.off()
  
# OK let's also graph the model with the logFC p25.c2 and logFC DC3000
# Step 1: Prepare metadata and DGEList
metadata$treatment <- factor(metadata$treatment, levels = c("dc3000", "p25c2"))
metadata$time_point <- factor(metadata$time_point, levels = c("t0", "t3"))
metadata$plant <- factor(metadata$plant)

dge <- DGEList(counts = counts)

dge <- calcNormFactors(dge)

# Step 2: Create design matrix with interaction
design <- model.matrix(~ treatment * time_point + plant + experiment, data = metadata)

# Step 3: Run voom with duplicateCorrelation to adjust for experiment
v <- voom(dge, design, plot = TRUE)
corfit <- duplicateCorrelation(v, design, block = metadata$experiment)
v <- voom(dge, design, block = metadata$experiment, correlation = corfit$consensus)
fit <- lmFit(v, design, block = metadata$experiment, correlation = corfit$consensus)
fit <- eBayes(fit)

# Step 4: Extract logFCs
lfc_mat <- fit$coefficients
dc3000_lfc <- lfc_mat[, "time_pointt3"]
p25c2_lfc <- dc3000_lfc + lfc_mat[, "treatmentp25c2:time_pointt3"]

# Step 5: Identify significant genes in DC3000 (main effect)
tt <- topTable(fit, coef = "time_pointt3", number = Inf, adjust.method = "BH")
sig_genes <- rownames(tt)[tt$adj.P.Val < 0.01]

# Step 6: Assemble dataframe
lfc_df <- data.frame(
  Gene = rownames(fit),
  logFC_DC3000 = dc3000_lfc,
  logFC_P25C2 = p25c2_lfc,
  Significant = ifelse(rownames(fit) %in% sig_genes, "Yes", "No")
)

# Step 7: Plot
pdf("logFC_P25C2_vs_DC3000.pdf", width = 3.5, height = 3.5, family = "Helvetica")  # Nature Eco Evo format
ggplot(lfc_df, aes(x = logFC_DC3000, y = logFC_P25C2, color = Significant)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_smooth(data = subset(lfc_df, Significant == "Yes"),
              method = "lm", se = FALSE, color = "gray40", linewidth = 0.7) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("R = ", round(cor.test(lfc_df$logFC_DC3000[lfc_df$Significant == "Yes"],
                                                 lfc_df$logFC_P25C2[lfc_df$Significant == "Yes"])$estimate, 2)),
           hjust = 1.2, vjust = -0.5, size = 3.5, family = "Helvetica") +
  scale_color_manual(values = c("Yes" = "#D55E00", "No" = "gray70")) +
  labs(
    x = expression(log[2]*"FC (DC3000, t3 vs t0)"),
    y = expression(log[2]*"FC (P25.C2, t3 vs t0)"),
    color = "Significant\nin DC3000"
  ) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.grid = element_blank()  # optional: removes background grid
  )
dev.off()
  
  
  
  
  
  
  
  
  


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
  return(v, fit)
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



