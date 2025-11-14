# This script takes the gene conservation information from the panX output and combines it with the in planta growth information to ask which variables are the best predictors of the function of a gene. 

library(car)
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
library(iml)
library(patchwork)

setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")

my_tufte_plot <- function(p) {
  p + 
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      plot.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    )
}
########## READ IN COUNTS TABLE
# May 2025 This script takes the full counts matrix from the many trials and will do subsetting and limma voom analysis. The goal is to ask what percentage of sig genes are genetic background specific in their behavior.
# setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/full_experiments/")

setwd("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")

##########################################################################
### Read in full counts table ###
#all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_8_2025.rds")
#all_exp <- readRDS("../full_experiments/all_p25_dc_axenic_5_2025.rds")
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

#remove experiments with dc3000
remove_dc3000 <- which(metadata$treatment=="dc3000")
p25_counts <- counts[-c(remove_dc3000),]
p25_metadata <- metadata[-c(remove_dc3000),]

#remove empty counts columns
p25_counts<-p25_counts[, !(colSums(p25_counts, na.rm = TRUE) == 0 | is.na(colSums(p25_counts, na.rm = TRUE)))]

# Now run the limmavoom model to get an estimate for the logFC for each of the genes over time 
ey <- which(p25_metadata$plant=="ey15_2")
p25_ey_metadata <- p25_metadata[ey,]
p25_ey_counts <- t(p25_counts[ey,])
lib_sizes <- colSums(p25_ey_counts)
p25_ey_counts_clean <- p25_ey_counts[, lib_sizes > 0]
design <- model.matrix(~ time_point, data = p25_ey_metadata)



# 
# v <- voom(dge, design, block = p25_ey_metadata$experiment)
# 
# fit <- lmFit(v, design, block = p25_ey_metadata$experiment, correlation = corfit$consensus)
# fit <- eBayes(fit)
# top_df <- topTable(fit, coef = "time_pointt3", number = Inf)
# 

###### THIS IS WHERE I GENERATE THE LOGFC DATA######### Actually I pulled in the lfc datatable from limmaVoom
lfc_df <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/logFCDC3000logFCp25c2_10_2025.csv")

#this is the output from running the python evolutionary genetics scripts on the panX output
conservation = read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd.csv", row.names = 1)

#this is the output from Effie's blasting (from which the orthologs were originally called). This is where we pull the sequence divergence.
blast_div <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit_2.txt")
blast_div <- blast_div[,c("V1", "V2", "V11")]
colnames(blast_div) <- c("DAK","WP", "perc_div_aa")
rownames(blast_div) <- blast_div$WP

#this is data file for the output from the tnseq experiment for p25.c2 in the tailocin plus also some of the plant genotype info
#eyach_tailocin <- read.table("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/tailocin/tailocin_plant_fitness_selection_717_2024.csv", header = TRUE)

#dc3000_eyach <- read.table("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/in_planta_rbtnseq_p25c2_dc3000/in_planta_may23/dc3000_may23/dc3000_sample_info_may23.txt", header = TRUE)

#I need to rename the gene names in conservation file
orthologs <- read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", 
   header = TRUE, sep = ",", row.names = 2)

#And assign NA to empty DC3000 orthology
orthologs[orthologs$DC3000=="",]$DC3000<-NA

#this was the direct output file of Effie.
other_ortholog <- read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit_2.txt",  sep="\t")

# change conservation DAKCFMEO to BJE identifier
conservation$BJE <- orthologs[rownames(conservation),]$p25_BJE
conservation$WP <- orthologs[rownames(conservation),]$DC3000

# remove those genes that are duplicated in the BJE genome
# get numbers of duplicated rows. Keeping original number
conservation2<-conservation %>%
  arrange(BJE) %>%
  filter(duplicated(BJE) == FALSE)
conservation2$DAC <- rownames(conservation2)
# remove genes with NA for BJE
conservation2 <- conservation2 %>% filter(is.na(BJE)==FALSE)
rownames(conservation2)<- conservation2$BJE

#Now I need to rename in terms of WP (the DC3000 mapping)
i=1
dupl_WP<-which(duplicated(conservation2$WP))
conservation2 <- conservation2[-dupl_WP,]
for(rowname in rownames(conservation2)){
   if(!is.na(conservation2[rowname,]$WP)){
    rownames(conservation2)[i]=conservation2[rowname,]$WP
   }
  i=i+1
}

# Now combine the diversity and tnseq files
rownames(lfc_df) <- lfc_df$gene
fit_con <- merge(conservation2, lfc_df, by=0, all=TRUE)
rownames(fit_con) <- fit_con$Row.names
fit_con <- merge(fit_con, blast_div, by=0, all=TRUE)
# There are ~6 genes in blast_div not in fit_con before merge

# now we can ask how the time in tree relates to the fitness effect
model1 <- lm(data=fit_con, logFC_DC3000_col0 ~ Time.in.tree + Genetic_Diversity + Num_gene_events)
names(fit_con) <- make.unique(names(fit_con), sep = "__dup")
#############

# 1. Compute correlation in each bin
df <- fit_con %>%
  mutate(div_bin = cut(perc_div_aa, breaks = seq(50, 100, by = 5), include.lowest = TRUE))
df <- df %>%
  mutate(div_bin = as.character(div_bin))
cor_by_bin <- df %>%
  filter(!is.na(div_bin)) %>%          # <--- drop NA bin
  group_by(div_bin) %>%
  summarise(
    cor_val = cor(logFC_DC3000_col0,
                  logFC_P25c2_col0,
                  use = "complete.obs"),
    .groups = "drop"
  )
# 2: Convert bin labels (e.g., "[60,65)") into numeric midpoints
cor_by_bin2 <- df %>%
  filter(!is.na(div_bin)) %>%
  group_by(div_bin) %>%
  summarise(
    n_pairs = sum(!is.na(logFC_DC3000_col0) &
                    !is.na(logFC_P25c2_col0)),
    cor_val = ifelse(
      n_pairs >= 2,
      cor(logFC_DC3000_col0, logFC_P25c2_col0, use = "complete.obs"),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  # --- Fisher CI only if n >= 10 (recommended threshold) ---
  mutate(
    z      = ifelse(n_pairs >= 10, atanh(cor_val), NA_real_),
    se_z   = ifelse(n_pairs >= 10, 1 / sqrt(n_pairs - 3), NA_real_),
    z_low  = ifelse(is.finite(se_z), z - 1.96 * se_z, NA_real_),
    z_high = ifelse(is.finite(se_z), z + 1.96 * se_z, NA_real_),
    cor_low  = tanh(z_low),
    cor_high = tanh(z_high)
  ) %>%
  mutate(
    div_range = gsub("\\[|\\(|\\)|\\]|\\s", "", div_bin)
  ) %>%
  separate(
    div_range, into = c("lower", "upper"),
    sep = ",", convert = TRUE
  ) %>%
  mutate(div_mid = (lower + upper) / 2) %>%
  select(div_bin, div_mid, n_pairs, cor_val, cor_low, cor_high)

# 3. Plot
div_pear <- ggplot(cor_by_bin2, aes(x = div_mid, y = cor_val)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_errorbar(
    aes(ymin = cor_low, ymax = cor_high),
    width = 0,
    linewidth = 0.6
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_pairs), show.legend = FALSE) +  # size ~ n_pairs (optional)
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(
    x = "% Amino Acid Divergence",
    y = "Correlation of fitness effect (DC3000 vs P25.c2 logFC)"
  )

pdf("divergence_correlatedFC.pdf", width = 3.5, height = 3, family = "Helvetica")  # Nature Eco
my_tufte_plot(div_pear)
dev.off()

################
#Now Let's move into the section where we look at the explanatory impact of different variables. read in the Tsuda expression data
################
exp <- read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/Tsuda_PNAS_expression_data/tk_modified_files/expression_gene_mapping.txt", sep="\t", header = TRUE)

# subset_expression to a few of interest. WP, RS, ID, name, annotation, Col.0_Pto
exp_sub <- exp %>% dplyr::select(c("WP", "RS", "ID", "name", "annotation", "Col.0_Pto", "KB_Pto",
                            "MM_Pto", "SA_Pto"))

var_exp <-exp %>% dplyr::select(c("Col.0_Pto", "KB_Pto",
                           "MM_Pto", "SA_Pto"))%>% as.matrix() %>% rowSds()
  #exp %>% select(c(4:35))%>% as.matrix() %>% rowSds()
exp_sub$variance <- var_exp
ratio<- exp_sub$Col.0_Pto/exp_sub$KB_Pto
exp_sub$ratio <- ratio

# merge expression info with fit_con
fit_exp <- merge(fit_con, exp_sub, by.x="WP.x", by.y="WP", all.x = TRUE)

# I would like to add the log2FC data from the limmaVoom_full_table script
#lfc_both_df <- read.csv("logFCDC3000logFCp25c2_7_2025.csv", header=TRUE)
#hm <- merge(fit_exp, lfc_both_df, by.x= "WP", by.y = "Gene")
#fit_exp <- hm

model2 <- lm(data=fit_exp,  logFC_DC3000_col0 ~ logFC_P25c2_col0 +Time.in.tree + KB_Pto +Genetic_Diversity + Col.0_Pto + Num_gene_events )
summary(model2)
anova_table <-Anova(model2, type=3)
anova_table$Term <- rownames(anova_table)
write.table(anova_table, "ANOVA_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write an output file that will then be fed into a random forest model

# The expression level in planta seems to be the best predictor of the selection coefficient in planta. 
# What about if I look at genes that are differentially expressed
cor.test(exp$Col.0_Pto, exp$KB_Pto)
cor.test(exp$Col.0_Pto, exp$MM_Pto)
# The minimal media did not seem to be a better predictor of gene expression in the plant than did the  KB
# it is not clear to me why the ratio of Col0/KB was not a better predictor than the overall variance. 


# Let's do modeling to understand which predictors matter most.At first I did random forest but I am not liking how I am handling the interaction betwen divergence and expression so will try glmnet




library(glmnet)

# Step 0: Define relevant variables
vars <- c("logFC_P25C2", "logFC_DC3000", "perc_div_aa", 
          "Num_gene_events", "Genetic_Diversity", "Col.0_Pto", "KB_Pto")

# Step 1: Clean data and add interaction term
# fit_exp_clean <- fit_exp %>%
#   dplyr::select(all_of(vars)) %>%
#   mutate(div_x_dc3000 = perc_div_aa * logFC_DC3000) %>%
#   na.omit()
fit_exp$scaled_div <- scale(fit_exp$perc_div_aa)
fit_exp$scaled_logFC <- scale(fit_exp$logFC_DC3000)
fit_exp$interaction <- fit_exp$scaled_div * fit_exp$scaled_logFC
fit_exp_clean <- fit_exp %>% 
  na.omit()

#first let's graph against all the variables of interest
# Define variables
predictor_vars <- c("perc_div_aa", "Num_gene_events", "Genetic_Diversity", "Col.0_Pto", "KB_Pto",
                    "interaction", "logFC_P25C2", "Time.in.tree")

# Reshape to long format
plot_data <- fit_exp_clean %>%
  dplyr::select(logFC_DC3000, all_of(predictor_vars)) %>%
  pivot_longer(cols = -logFC_DC3000, names_to = "Variable", values_to = "Value")

ggplot(plot_data, aes(x = Value, y = logFC_DC3000)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, linetype = "dashed", color = "black") +
  facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
  labs(
    x = NULL,
    y = expression(log[2] * "FC in DC3000"),
    title = expression("Predictors vs. log"[2] * "FC in DC3000")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3)
  )






# Step 2: Define response and predictor matrix
Y <- fit_exp_clean$logFC_P25C2
X <- model.matrix(~ logFC_DC3000 + 
                    perc_div_aa +
                    Genetic_Diversity +
                    Col.0_Pto + 
                    KB_Pto +
                    interaction +
                    Time.in.tree, data = fit_exp_clean[, -1])[, -1]  # Remove intercept

# Step 3: Fit elastic net model with cross-validation
set.seed(123)
fit <- cv.glmnet(X, Y, alpha = 0.5)

# Step 4: View coefficients at optimal lambda
coef(fit, s = "lambda.min")

# Step 5 (optional): Calculate R²
y_pred <- predict(fit, newx = X, s = "lambda.min")
rsq <- 1 - sum((Y - y_pred)^2) / sum((Y - mean(Y))^2)
cat("R² =", round(rsq, 3), "\n")


# Predict values at lambda.min

# Build data frame for plotting
  library(grid)  # for unit()
  
  # Predict values at lambda.min
  y_pred <- as.numeric(predict(fit, newx = X, s = "lambda.min"))
  
  # Build data frame for plotting
  plot_df <- data.frame(
    Observed = Y,
    Predicted = y_pred
  )
  
  # Compute R²
  r2 <- round(cor(plot_df$Observed, plot_df$Predicted)^2, 2)
  
  # Plot
  ggplot(plot_df, aes(x = Observed, y = Predicted)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "black", linetype = "dashed") +
    annotate("text", x = min(plot_df$Observed), y = max(plot_df$Predicted),
             label = paste0("R² = ", r2), hjust = 0, vjust = 1.5, size = 3.5) +
    theme_minimal(base_size = 10) +
    labs(
      x = expression("Observed " * log[2] * "FC"),
      y = expression("Predicted " * log[2] * "FC")
    ) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.3),
      plot.margin = unit(c(5, 5, 5, 5), "mm"),
      aspect.ratio = 1
    )


# Those analyses are fine. 
  # Linear model with interaction term
  model_interact <- lm(logFC_P25c2_col0 ~ logFC_DC3000_col0 * perc_div_aa, data = fit_con)
  model_no_interact <- lm(logFC_P25c2_col0 ~ logFC_DC3000_col0 , data = fit_con)
  # Get adjusted R²
  no_adj_r2 <- summary(model_no_interact)$adj.r.squared
  adj_r2 <- summary(model_interact)$adj.r.squared
  
  # Report
  cat(sprintf("Adjusted R² with interaction = %.1f%%\n", adj_r2 * 100))
  
  
  
  
  

# https://www.guru99.com/r-random-forest-tutorial.html
# Define the control
# library(caret)
# 
# # Assume fit_exp is your dataframe
# # Columns: log2FC, KB_Pto, Col.o_Pto, variance, Time.in.tree, Genetic_Diversity, Num_gene_events
# 
# # Remove rows with missing values (RF can handle NAs with some packages, but safest to clean)
# fit_exp_clean <- na.omit(fit_exp)
# fit_exp_clean$interaction <- interaction(fit_exp_clean$logFC_P25C2, fit_exp_clean$perc_div_aa)
# model_formula <- logFC_DC3000 ~logFC_P25C2 + perc_div_aa + interaction + KB_Pto + Col.0_Pto +  Time.in.tree + Genetic_Diversity + Num_gene_events
# 
# set.seed(123)  # for reproducibility
# 
# ctrl <- trainControl(method = "cv", number = 5)  # 5-fold CV
# 
# tuneGrid <- expand.grid(.mtry = c(2, 3, 4))  # Tune mtry between 2 and 4
# 
# rf_model <- train(
#   model_formula,
#   data = fit_exp_clean,
#   method = "rf",
#   trControl = ctrl,
#   tuneGrid = tuneGrid,
#   importance = TRUE
# )
# 
# 
# print(rf_model)
# varImp(rf_model)
# plot(varImp(rf_model), top = 10)
# 
# 
# 
# predicted <- predict(rf_model, newdata = fit_exp_clean)
# plot_df <- data.frame(
#     Observed = fit_exp_clean$logFC_DC3000,
#   Predicted = predicted
# )
# r2 <- cor(plot_df$Observed, plot_df$Predicted)^2
# observed_predict_plot <- ggplot(plot_df, aes(x = Observed, y = Predicted)) +
#   geom_point(alpha = 0.6) +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   annotate("text", x = min(plot_df$Observed), y = max(plot_df$Predicted),
#            label = paste0("R² = ", round(r2, 3)), hjust = 0, vjust = 1, size = 5) +
#   labs(
#     title = "Observed vs Predicted log2FC",
#     x = "Observed log2FC",
#     y = "Predicted log2FC"
#   ) +
#   theme_minimal(base_size = 14)
# 
# 
# # Holy shit. Why is it so good at predicting the positive fold change? Why is there nothing in the wrong place?
# #Let's look at some specific datapoints. Why does the positive get predicted positive? And the really negative really neative. why does the KB Pto bump around so much? This doesn't make sense to me. OK apparently I was looking at apparent performance rather than CV performance
# 
# #this is the CV performance
# # Use all cv_preds (since mtry is constant). Actually no. I want to filter the cv_preds to those with a |log2FC|>1
# cv_preds_best <- cv_preds  %>% filter(abs(obs)>1)
# 
# # Compute CV R²
# r2_cv <- cor(cv_preds_best$obs, cv_preds_best$pred)^2
# 
# # Plot
# cv_plot <- ggplot(cv_preds_best, aes(x = obs, y = pred)) +
#   geom_point(alpha = 0.6, color = "grey50") +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   annotate("text", 
#            x = min(cv_preds_best$obs, na.rm = TRUE), 
#            y = max(cv_preds_best$pred, na.rm = TRUE),
#            label = paste0("CV R² = ", round(r2_cv, 3)), 
#            hjust = 0, vjust = 1, size = 5) +
#   labs(
#     title = "Observed vs CV Predicted log2FC",
#     x = "Observed log2FC",
#     y = "Predicted log2FC"
#   ) +
#   theme_minimal(base_size = 14)
# 
# print(cv_plot)
# ggsave("observed_predict_plot.pdf", cv_plot, width = 10, height = 8)
# 
# ## can I compare predictive capacity for different ranges?
# # Example: Define bins
# cv_preds_best <- cv_preds_best %>%
#   mutate(log2FC_bin = cut(obs, 
#                           breaks = c(-Inf, -1.5, 0, 1.5, Inf), 
#                           labels = c("Low", "Medium-Neg", "Medium-Pos", "High")))
# 
# # Compute R² in each bin
# bin_r2 <- cv_preds_best %>%
#   group_by(log2FC_bin) %>%
#   summarise(
#     n = n(),
#     R2 = cor(obs, pred)^2
#   )
# 
# print(bin_r2)
# #now print observed vs predicted:
# ggplot(cv_preds_best, aes(x = obs, y = pred)) +
#   geom_point(alpha = 0.5) +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   facet_wrap(~ log2FC_bin) +
#   labs(title = "Observed vs CV Predicted by log2FC bin",
#        x = "Observed log2FC", y = "Predicted log2FC") +
#   theme_minimal()
# 
# 
# #######
# # Define predictors
# predictors <- c("KB_Pto", "Col.0_Pto",  "Time.in.tree", "Genetic_Diversity", "Num_gene_events")
# 
# # Create plots
# plot_list <- lapply(predictors, function(var) {
#   
#   # Compute R2
#   model <- lm(as.formula(paste("log2FC ~", var)), data = fit_exp_clean)
#   r2_val <- summary(model)$r.squared
#   
#   ggplot(fit_exp_clean, aes_string(x = var, y = "log2FC")) +
#     geom_point(alpha = 0.3, color = "grey50") +  # dull points
#     geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.2) +  # trend line with SE
#     annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, 
#              label = paste0("R² = ", round(r2_val, 3)), size = 4) +
#     labs(
#       title = var,
#       x = var,
#       y = "log2FC"
#     ) +
#     theme_minimal(base_size = 12)
# })
# 
# # Combine into grid
# combined_plot <- wrap_plots(plot_list, ncol = 2)
# 
# # Display
# print(combined_plot)
# ggsave("log2FC_predictor_scatterplots_withR2_SE.pdf", combined_plot, width = 10, height = 8)
