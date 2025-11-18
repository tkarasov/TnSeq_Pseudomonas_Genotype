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

# Transform counts with voom. The following takes the dge (subsetted or not, the design and a boolean for whether or not the experiment random effect should be taken into account)

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

design <- model.matrix(~ treatment * time_point * plant, data = metadata)
# experiment is included as a random effect only
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
res  <- run_voom(dge, design, block = metadata$experiment)
v <- res$v
fit <- res$fit
result_full <- summarize_overlap_between_terms(
  fit,
  coef1 = "time_pointt3",
  coef2 = "treatmentp25c2:time_pointt3"
)
# 172 significant coef1, 130 significant coef2, num_overlap 64
#238 genes total are significant

# If I include experiment as a fixed effect, get coef1=177, coef2=122, num_overlap=65
######### Second do the analysis with Exp 2
filtered2 <- filter_meta_counts("experiment", "exp_0001", counts, metadata)
design2 <- model.matrix(~ treatment * time_point * plant, data = filtered2$new_meta)
dge2 <- DGEList(counts = filtered2$new_counts)
dge2 <- calcNormFactors(dge2)
res2 <- run_voom(dge2, design2, block = NULL)
v2 <- res2$v
fit2 <- res2$fit

######### Now compare the results between the full and fit2
result_full2 <- summarize_overlap_between_terms(
  fit2,
  coef1 = "time_pointt3",
  coef2 = "treatmentp25c2:time_pointt3"
)
# 147 coef1, 84 coef2, 47 overlap



######### The first experiment was kind of a mess, so I don't want to do it specifically because doesn't have good controls
# But let's look at the correspondence in lfc values between the full model and just exp2.
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
design_full <- design
# design <- model.matrix(~ treatment * time_point * plant, data = metadata)
res <- run_voom(dge, design_full, block = metadata$experiment)
v <- res$v
fit <- res$fit

# Extract LFC for interaction term "time_pointt3" for each experiment
lfc_timepoint <- extract_timepoint_logfc_per_experiment(
  voom_obj = res$v,
  metadata = metadata,
  experiment_col = "experiment",
  design_formula = ~ treatment * time_point * plant
)


# === Parameters ===
coef_name <- "time_pointt3"
fdr_cutoff <- 0.05

# === 1. Identify significant genes ===
sig_full <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
sig_full_ids <- rownames(sig_full)[sig_full$adj.P.Val < fdr_cutoff]
length(sig_full_ids)
#this is 273. I ran again and now it's 139. OK it's 139 again. OK now it's 172. Ive changed the FDR to 0.05. But 172 is consistent with above 

sig_exp2 <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")
sig_exp2_ids <- rownames(sig_exp2)[sig_exp2$adj.P.Val < fdr_cutoff]
length(sig_exp2_ids)
#this is 126. Still 126. Nope now 147 wtf

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
  fdr_cutoff <- 0.05
  
  # === 1. Identify significant genes ===
  sig_full <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
  sig_full_ids <- rownames(sig_full)[sig_full$adj.P.Val < fdr_cutoff]
  length(sig_full_ids)
  # 139 genes are in this sig_full_ids. Now it's 172 after correcting for FDR=0.05 BH
  
  
  
  sig_exp2 <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")
  sig_exp2_ids <- rownames(sig_exp2)[sig_exp2$adj.P.Val < fdr_cutoff]
  length(sig_exp2_ids)
  #126 genes are in this group. Now it's 147 in this group after correcting for FDR=0.05.
  
  
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

pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/logFC_comparison_DC3000.pdf", width = 3.5, height = 3.5, family = "Arial")
  # 89 mm × 89 mm
  
replication<-  ggplot(lfc_sig, aes(x = logFC_exp_0001, y = logFC_exp_0002, color = status)) +
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
replication 
dev.off()

#### OK now we have a subset of genes we can follow. Those that show a significant fitness effect in DC3000 over time. We will continue to look at the full model genes. What I want to do with them is ask questions about how many of them are sensitive to the strain genetic background and how many of them are sensitive to the plant genetic background. And I want good graphics to relay this information.

summarize_background_sensitivity_with_3way <- function(
    fit,
    coef_time = "time_pointt3",
    coef_strain_time = "treatmentp25c2:time_pointt3",
    coef_plant_time = "time_pointt3:plantey15_2",
    coef_3way = "treatmentp25c2:time_pointt3:plantey15_2",
    fdr_cutoff = 0.05,
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
  fit = fit2,
  coef_time = "time_pointt3",
  coef_strain_time = "treatmentp25c2:time_pointt3",
  coef_plant_time = "time_pointt3:plantey15_2",
  coef_3way = "treatmentp25c2:time_pointt3:plantey15_2",
  fdr_cutoff = 0.05,
  plot_type = "upset"  # or "upset" if you have ComplexUpset
  )

# the plotting of the upset plot isn't working well. But only the strain genetic background had an effect.
# So I really just want to make a bar plot
tot_sig <- length(which(gene_summary$Time==TRUE))
#139 genes. Now 172 after FDR=0.05
strain_sig <- length(which(gene_summary$Strain_Time==TRUE))
# 45 genes. Now 68 genes after FDR=0.05
plant_sig <- length(which(gene_summary$Plant_Time==TRUE))
# 0 genes. Still 0 after FDR=0.05

threeway_sig <- length(which(gene_summary$Three_Way==TRUE))
#0 genes. Still 0 after FDR=0.05


# write the gene_summary table to file to be used in the gene_ontology_contrasts_script
write.csv(gene_summary, file = "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/gene_sig_summary_11_2025.csv", row.names = FALSE,  quote = FALSE)

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
pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/genetic_background_sensitivity_barplot.pdf", width = 3.5, height = 3.5, family = "Arial")

sig_bar <- ggplot(bar_data, aes(x = Background, y = Count, fill = Significance)) +
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
sig_bar
dev.off()
  
# OK let's also graph the model with the logFC p25.c2 and logFC DC3000
# # Step 1: Prepare metadata and DGEList
# metadata$treatment <- factor(metadata$treatment, levels = c("dc3000", "p25c2"))
# metadata$time_point <- factor(metadata$time_point, levels = c("t0", "t3"))
# metadata$plant <- factor(metadata$plant)
# 
# dge <- DGEList(counts = counts)
# 
# dge <- calcNormFactors(dge)
# 
# # Step 2: Create design matrix with interaction
# design <- model.matrix(~ treatment * time_point * plant , data = metadata)
# 
# # Step 3: Run voom with duplicateCorrelation to adjust for experiment. We've done this already many times. Just making sure we are using the right model
# v <- voom(dge, design, plot = TRUE)
# corfit <- duplicateCorrelation(v, design, block = metadata$experiment)
# v <- voom(dge, design, block = metadata$experiment, correlation = corfit$consensus)
# fit <- lmFit(v, design, block = metadata$experiment, correlation = corfit$consensus)
# fit <- eBayes(fit)
# 
# ############# Newon 8/4/2025
# # Clean up design matrix
# # Ensure correct factor levels
# metadata$treatment <- factor(metadata$treatment, levels = c("dc3000", "p25c2"))
# metadata$time_point <- factor(metadata$time_point, levels = c("t0", "t3"))
# metadata$plant <- factor(metadata$plant, levels = c("col_0", "ey15_2"))
# 
# # Create design matrix and sanitize column names
# design <- model.matrix(~ treatment * time_point * plant, data = metadata)
# colnames(design) <- make.names(colnames(design))
# 
# # Normalize counts and run voom with duplicateCorrelation
# dge <- DGEList(counts = counts)
# dge <- calcNormFactors(dge)
# 
# v <- voom(dge, design, plot = TRUE)
# corfit <- duplicateCorrelation(v, design, block = metadata$experiment)
# v <- voom(dge, design, block = metadata$experiment, correlation = corfit$consensus)
# fit <- lmFit(v, design, block = metadata$experiment, correlation = corfit$consensus)
# fit <- eBayes(fit)

# Define contrasts for t3 vs t0 within each strain × plant condition
colnames(design2) <- make.names(colnames(design2))
# 1. DC3000 in col_0 (baseline)
dc_col0_contrast <- makeContrasts(
  logFC_P25c2_col0 = 
    `time_pointt3` + `treatmentp25c2.time_pointt3`,
  levels = design2
)

# 2. P25c2 in col_0
p25_col0_contrast <- makeContrasts(
  logFC_P25c2_col0 = time_pointt3 + treatmentp25c2.time_pointt3,
  levels = design2
)

# 3. DC3000 in ey15_2
dc_ey15_contrast <- makeContrasts(
  logFC_DC3000_ey15 = time_pointt3 + time_pointt3.plantey15_2,
  levels = design2
)

# 4. P25c2 in ey15_2
p25_ey15_contrast <- makeContrasts(
  logFC_P25c2_ey15 = time_pointt3 + time_pointt3.plantey15_2 +
    treatmentp25c2.time_pointt3 + treatmentp25c2.time_pointt3.plantey15_2,
  levels = design2
)

# Apply contrasts and extract logFCs

lfc_df <- data.frame(gene = rownames(fit2$coefficients))

fit_dc_col0 <- eBayes(contrasts.fit(fit2, dc_col0_contrast))
lfc_df$logFC_DC3000_col0 <- fit_dc_col0$coefficients[, 1]

fit_p25_col0 <- eBayes(contrasts.fit(fit2, p25_col0_contrast))
lfc_df$logFC_P25c2_col0 <- fit_p25_col0$coefficients[, 1]

fit_dc_ey <- eBayes(contrasts.fit(fit2, dc_ey15_contrast))
lfc_df$logFC_DC3000_ey15 <- fit_dc_ey$coefficients[, 1]

fit_p25_ey <- eBayes(contrasts.fit(fit2, p25_ey15_contrast))
lfc_df$logFC_P25c2_ey15 <- fit_p25_ey$coefficients[, 1]

# Preview result
head(lfc_df)







# Step 5: Identify significant genes in DC3000 (main effect)
tt <- topTable(fit, coef = "time_pointt3", number = Inf, adjust.method = "BH")
sig_genes <- rownames(tt)[tt$adj.P.Val < 0.05]
#sig_genes now has 172

# Step 6: Assemble dataframe
lfc_df$Significant = ifelse(rownames(fit) %in% sig_genes, "Yes", "No")


#Write out dataframe with the lfc for DC3000 and for p25.c2
write.table(lfc_df, "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/logFCDC3000logFCp25c2_10_2025.csv", quote=FALSE, col.names =TRUE, sep=",", row.names=FALSE)


# Step 7: Plot
pdf("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/logFC_P25C2_vs_DC3000_col_eyach.pdf", width = 3.5, height = 7, family = "Helvetica")  # Nature Eco Evo format
dc3000_p25_col0 <- ggplot(lfc_df, aes(x = logFC_DC3000_col0, y = logFC_P25c2_col0, color = Significant)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_smooth(data = subset(lfc_df, Significant == "Yes"),
              method = "lm", se = FALSE, color = "gray40", linewidth = 0.7) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("R = ", round(cor.test(lfc_df$logFC_DC3000_col0[lfc_df$Significant == "Yes"],
                                                 lfc_df$logFC_P25c2_col0[lfc_df$Significant == "Yes"])$estimate, 2)),
           hjust = 1.2, vjust = -0.5, size = 3.5, family = "Helvetica") +
  scale_color_manual(values = c("Yes" = "#D55E00", "No" = "gray70")) +
  labs(
    x = expression(log[2]*"FC (DC3000)"),
    y = expression(log[2]*"FC (P25.C2)"),
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

  

col_eyach <- ggplot(lfc_df, aes(x = logFC_DC3000_col0 , y = logFC_DC3000_ey15, color = Significant)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_smooth(data = subset(lfc_df, Significant == "Yes"),
              method = "lm", se = FALSE, color = "gray40", linewidth = 0.7) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("R = ", round(cor.test(lfc_df$logFC_DC3000_col0[lfc_df$Significant == "Yes"],
                                                 lfc_df$logFC_DC3000_ey15[lfc_df$Significant == "Yes"])$estimate, 2)),
           hjust = 1.2, vjust = -0.5, size = 3.5, family = "Helvetica") +
  scale_color_manual(values = c("Yes" = "#D55E00", "No" = "gray70")) +
  labs(
    x = expression(log[2]*"FC (DC3000 in Col-0)"),
    y = expression(log[2]*"FC (DC3000 in Ey1.5)"),
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
  
  
  plot_grid(col_eyach,dc3000_p25_col0,  nrow=2)
  
  
  dev.off()
pdf("overlaid_hist_logFC_P25C2_vs_DC3000_col.pdf", width = 3.5, height =3.5, family = "Helvetica")  # Nature Eco Evo format
  # Plotting overlaid histograms of logFC (t3 - t0) for DC3000 and p25.c2 in col_0
  
  # Subset relevant columns
  df_plot <- lfc_df[, c("logFC_DC3000_col0", "logFC_P25c2_col0")]
  colnames(df_plot) <- c("DC3000", "p25.c2")
  
  # Convert to long format
  df_long <- df_plot %>%
    pivot_longer(cols = everything(), names_to = "Strain", values_to = "logFC")
  
  # Create smoothed density plot
  ggplot(df_long, aes(x = logFC, fill = Strain, color = Strain)) +
    geom_density(alpha = 0.4, linewidth = 1.2) +
    theme_minimal(base_family = "Arial") +
    scale_fill_manual(values = c("DC3000" = "#0072B2", "p25.c2" = "#E69F00")) +
    scale_color_manual(values = c("DC3000" = "#0072B2", "p25.c2" = "#E69F00")) +
    labs(
      title = "Smoothed logFC (t3 - t0) in Col-0",
      x = "log2(Fitness fold change)",
      y = "Density",
      fill = "Strain",
      color = "Strain"
    ) +
    theme(
      text = element_text(size = 12, family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )
  
  
 dev.off() 
  
    
  
  


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

# Combine into a dataframe
results_summary <- do.call(rbind, results_summary)
# I think the calculations in results summary are not correct.

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


# 

######### starting over 11/14/2025
## --- Restrict analysis to experiment 2 only ---

metadata_exp2 <- metadata[metadata$experiment == "exp_0002", , drop = FALSE]
metadata_exp2 <- droplevels(metadata_exp2)

# counts is genes x samples (after your transpose)
# subset columns to samples in experiment 2
counts_exp2 <- counts[, rownames(metadata_exp2), drop = FALSE]

# Make sure factors are correctly set
metadata_exp2$treatment  <- factor(metadata_exp2$treatment, levels = c("dc3000", "p25c2"))
metadata_exp2$time_point <- factor(metadata_exp2$time_point, levels = c("t0", "t3"))
metadata_exp2$plant      <- factor(metadata_exp2$plant,     levels = c("col_0", "ey15_2"))

# Design matrix for experiment 2 ONLY
design_exp2 <- model.matrix(~ treatment * time_point * plant, data = metadata_exp2)
colnames(design_exp2) <- make.names(colnames(design_exp2))

dge_exp2 <- DGEList(counts = counts_exp2)
dge_exp2 <- calcNormFactors(dge_exp2)

v_exp2   <- voom(dge_exp2, design_exp2, plot = TRUE)
fit_exp2 <- lmFit(v_exp2, design_exp2)
fit_exp2 <- eBayes(fit_exp2)

# 1. DC3000 in col_0 (baseline)
dc_col0_contrast <- makeContrasts(
  logFC_DC3000_col0 = time_pointt3,
  levels = design_exp2
)

# 2. P25c2 in col_0
p25_col0_contrast <- makeContrasts(
  logFC_P25c2_col0 = time_pointt3 + treatmentp25c2.time_pointt3,
  levels = design_exp2
)

# 3. DC3000 in ey15_2
dc_ey15_contrast <- makeContrasts(
  logFC_DC3000_ey15 = time_pointt3 + time_pointt3.plantey15_2,
  levels = design_exp2
)

# 4. P25c2 in ey15_2
p25_ey15_contrast <- makeContrasts(
  logFC_P25c2_ey15 = time_pointt3 + time_pointt3.plantey15_2 +
    treatmentp25c2.time_pointt3 + treatmentp25c2.time_pointt3.plantey15_2,
  levels = design_exp2
)

lfc_df_exp2 <- data.frame(gene = rownames(fit_exp2$coefficients))

fit_dc_col0_exp2 <- eBayes(contrasts.fit(fit_exp2, dc_col0_contrast))
lfc_df_exp2$logFC_DC3000_col0 <- fit_dc_col0_exp2$coefficients[, 1]

fit_p25_col0_exp2 <- eBayes(contrasts.fit(fit_exp2, p25_col0_contrast))
lfc_df_exp2$logFC_P25c2_col0 <- fit_p25_col0_exp2$coefficients[, 1]

fit_dc_ey_exp2 <- eBayes(contrasts.fit(fit_exp2, dc_ey15_contrast))
lfc_df_exp2$logFC_DC3000_ey15 <- fit_dc_ey_exp2$coefficients[, 1]

fit_p25_ey_exp2 <- eBayes(contrasts.fit(fit_exp2, p25_ey15_contrast))
lfc_df_exp2$logFC_P25c2_ey15 <- fit_p25_ey_exp2$coefficients[, 1]


# Significance
tt_exp2 <- topTable(fit_exp2, coef = "time_pointt3", number = Inf, adjust.method = "BH")
sig_genes_exp2 <- rownames(tt_exp2)[tt_exp2$adj.P.Val < 0.05]

lfc_df_exp2$Significant <- ifelse(rownames(fit_exp2) %in% sig_genes_exp2, "Yes", "No")





