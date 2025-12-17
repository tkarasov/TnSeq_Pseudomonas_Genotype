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

source("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/scripts/isme_fig_settings.R")


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

######### Actually I pulled in the lfc datatable from limmaVoom
#lfc_df <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/logFCDC3000logFCp25c2_10_2025.csv")

# On 12/15/2025 I replace the lfc_df file with sig_results_12_2025.csv
sig_results <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/sig_results_12_2025.csv",header=TRUE)


#this is the output from running the python evolutionary genetics scripts on the panX output
conservation = read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd.csv", row.names = 1)

#this is the output from Effie's blasting (from which the orthologs were originally called). This is where we pull the sequence divergence.
blast_div <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit_2.txt")
blast_div <- blast_div[,c("V1", "V2", "V11")]
colnames(blast_div) <- c("DAK","WP", "perc_div_aa")
rownames(blast_div) <- blast_div$WP

#I need to rename the gene names in conservation file. They currently have the DAK naming but need the BJE naming.
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

# What is the output for the estimated LFC of P25.C2 and DC3000.
# Now combine the diversity and tnseq files
lfc_df <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000logFC_time_DC3000_vs_P25c2_12_2025.csv", header=TRUE)


rownames(lfc_df) <- lfc_df$gene
fit_con <- merge(conservation2, lfc_df, by=0, all=TRUE)
rownames(fit_con) <- fit_con$Row.names
fit_con$gene <- rownames(fit_con)
blast_div$gene <- rownames(blast_div)
# I want to give lfc_df a new column if the gene  is significant in DC3000 in the padj_DC3000_time of sig_results (<0.01)
lfc_df$Significant <- ifelse(lfc_df$gene %in% sig_results[sig_results$padj_DC3000_time<0.01,]$gene, "Yes", "No")
# remove if present
if ("Row.names" %in% names(fit_con)) fit_con$Row.names <- NULL

# now safe to merge
merged <- merge(fit_con, blast_div, by = "gene", all = TRUE)

# sanity check
any(duplicated(names(merged)))   # should be FALSE
fit_con <- merge(fit_con, blast_div, by = "gene")

fit_con <- merge(fit_con, blast_div, by=0, all=TRUE)

# now we can ask how the time in tree relates to the fitness effect
model1 <- lm(data=fit_con, logFC_DC3000 ~ Time.in.tree + Genetic_Diversity + Num_gene_events)
names(fit_con) <- make.unique(names(fit_con), sep = "__dup")
#############

# 1. Compute correlation in each bin
df <- fit_con %>%
  mutate(div_bin = cut(perc_div_aa.x, breaks = seq(50, 100, by = 5), include.lowest = TRUE))
df <- df %>%
  mutate(div_bin = as.character(div_bin))
cor_by_bin <- df %>%
  filter(!is.na(div_bin)) %>%          # <--- drop NA bin
  group_by(div_bin) %>%
  summarise(
    cor_val = cor(logFC_DC3000,
                  logFC_P25c2,
                  use = "complete.obs"),
    .groups = "drop"
  )
# 2: Convert bin labels (e.g., "[60,65)") into numeric midpoints
library(dplyr)
library(tidyr)

cor_by_bin2 <- df %>%
  filter(!is.na(div_bin)) %>%
  group_by(div_bin) %>%
  summarise(
    n_pairs = as.numeric(sum(!is.na(logFC_DC3000) & !is.na(logFC_P25c2))),
    cor_val = ifelse(n_pairs >= 2,
                     cor(logFC_DC3000, logFC_P25c2, use = "complete.obs"),
                     NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    # Fisher z only where we have enough pairs
    z = ifelse(n_pairs >= 10, atanh(cor_val), NA_real_)
  ) %>%
  # --- safe vectorized computation of se_z: no sqrt() on bad values ---
  { tmp <- .
  tmp$se_z <- NA_real_                    # initialize
  good_idx <- which(!is.na(tmp$n_pairs) & tmp$n_pairs >= 10 & (tmp$n_pairs - 3) > 0)
  if (length(good_idx) > 0) {
    tmp$se_z[good_idx] <- 1 / sqrt(tmp$n_pairs[good_idx] - 3)
  }
  tmp
  } %>%
  mutate(
    z_low  = ifelse(!is.na(se_z), z - 1.96 * se_z, NA_real_),
    z_high = ifelse(!is.na(se_z), z + 1.96 * se_z, NA_real_),
    cor_low  = ifelse(!is.na(z_low), tanh(z_low), NA_real_),
    cor_high = ifelse(!is.na(z_high), tanh(z_high), NA_real_)
  ) %>%
  mutate(
    div_range = gsub("\\[|\\(|\\)|\\]|\\s", "", div_bin)
  ) %>%
  separate(div_range, into = c("lower", "upper"), sep = ",", convert = TRUE) %>%
  mutate(div_mid = (lower + upper) / 2) %>%
  select(div_bin, div_mid, n_pairs, cor_val, cor_low, cor_high)

# Inspect
cor_by_bin2
str(cor_by_bin2)



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

#pdf("divergence_correlatedFC.pdf", width = 3.5, height = 3, family = "Helvetica")  # Nature Eco
p1 <- my_tufte_plot(div_pear)
#dev.off()
save_figure("divergence_correlatedFC2.pdf",p1)

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

model2 <- lm(data=fit_exp,  logFC_DC3000 ~ logFC_P25c2 +Time.in.tree + KB_Pto +Genetic_Diversity + Col.0_Pto + Num_gene_events )
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

# # Step 0: Define relevant variables
# vars <- c("logFC_P25C2", "logFC_DC3000", "perc_div_aa", 
#           "Num_gene_events", "Genetic_Diversity", "Col.0_Pto", "KB_Pto")
# 
# # Step 1: Clean data and add interaction term
# # fit_exp_clean <- fit_exp %>%
# #   dplyr::select(all_of(vars)) %>%
# #   mutate(div_x_dc3000 = perc_div_aa * logFC_DC3000) %>%
# #   na.omit()
# fit_exp$scaled_div <- scale(fit_exp$perc_div_aa.x)
# fit_exp$scaled_logFC <- scale(fit_exp$logFC_DC3000)
# fit_exp$interaction <- fit_exp$scaled_div * fit_exp$scaled_logFC
# fit_exp_clean <- fit_exp %>% 
#   na.omit()
# 
# #first let's graph against all the variables of interest
# # Full code: faceted predictor plots with R / p labels in upper-right corner
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(purrr)
# 
# 
# # --- Ensure required columns exist; create numeric interaction if missing ---
# if (!"interaction" %in% names(fit_exp_clean)) {
#   if (all(c("logFC_P25c2_col0", "perc_div_aa") %in% names(fit_exp_clean))) {
#     fit_exp_clean <- fit_exp_clean %>%
#       mutate(interaction = logFC_P25c2 * perc_div_aa)
#   } else {
#     stop("Neither 'interaction' column exists nor the necessary columns to create it (logFC_P25c2_col0, perc_div_aa) are present.")
#   }
# }



library(dplyr); library(tidyr); library(ggplot2); library(broom)

# ---------- 0) Choose predictors you want to plot / model --------------
predictor_vars <- c("perc_div_aa.x", "Num_gene_events", "Genetic_Diversity",
                    "Col.0_Pto", "logFC_P25c2", "Time.in.tree")   # edit as needed

# ---------- 1) Quick availability report (useful summary) -------------
avail <- sapply(c("logFC_DC3000", predictor_vars), function(col) {
  if (col %in% names(fit_exp)) {
    c(non_na = sum(!is.na(fit_exp[[col]])), unique_examples = paste(head(unique(na.omit(fit_exp[[col]])),3), collapse = ", "))
  } else {
    c(non_na = 0, unique_examples = NA_character_)
  }
})
avail <- as.data.frame(t(avail), stringsAsFactors = FALSE)
avail$non_na <- as.numeric(avail$non_na)
cat("Availability summary (non-NA counts):\n")
print(avail)

# ---------- 2) Normalize names if there are alternate names -------------
# (you referenced both logFC_P25c2 and logFC_P25c2_col0 in earlier code)
if (!("logFC_P25c2" %in% names(fit_exp)) & ("logFC_P25c2" %in% names(fit_exp))) {
  fit_exp <- fit_exp %>% mutate(logFC_P25c2 = logFC_P25c2)
}
if (!("perc_div_aa.x" %in% names(fit_exp)) & ("perc_div_aa" %in% names(fit_exp))) {
  fit_exp <- fit_exp %>% mutate(perc_div_aa.x = perc_div_aa)
}

# refresh list of actually usable predictors after normalization
use_predictors <- intersect(predictor_vars, names(fit_exp))
cat("Using predictors:", paste(use_predictors, collapse = ", "), "\n")

# ---------- 3) Build a cleaned table for flexible plotting -------------
# Keep rows with response present and at least one predictor present.
fit_exp_clean2 <- fit_exp %>%
  # ensure response exists
  filter(!is.na(logFC_DC3000)) %>%
  # require at least one predictor to be non-NA so that a facet will have data
  filter(if_any(all_of(use_predictors), ~ !is.na(.x))) %>%
  # keep original identifiers for troubleshooting
  mutate(rowid = row_number())

cat("Rows kept for plotting (response present + >=1 predictor):", nrow(fit_exp_clean2), "\n")

# ---------- 4) Optionally create model-ready dataset (all predictors non-NA) ----------
model_ready <- fit_exp_clean2 %>%
  filter(complete.cases(select(., all_of(use_predictors))))

cat("Rows with ALL chosen predictors present (model-ready):", nrow(model_ready), "\n")
if (nrow(model_ready) == 0) {
  cat("No rows have all predictors present. Consider dropping predictors or imputing.\n")
}

# ---------- 5) Create interaction column safely if needed ---------------
if (!"interaction" %in% names(fit_exp_clean2)) {
  if (all(c("logFC_P25c2", "perc_div_aa.x") %in% names(fit_exp_clean2))) {
    fit_exp_clean2 <- fit_exp_clean2 %>%
      mutate(interaction = logFC_P25c2 * perc_div_aa.x)
    cat("Created 'interaction' = logFC_P25c2 * perc_div_aa.x\n")
  } else {
    cat("Did NOT create 'interaction' (missing logFC_P25c2 or perc_div_aa.x)\n")
  }
}

# ---------- 6) Build plot_data (long format) for facetted plotting -----
plot_data <- fit_exp_clean2 %>%
  select(any_of(c("gene.x", "gene", "logFC_DC3000", use_predictors))) %>%
  # keep an ID column if present
  rename(gene_id = any_of(c("gene.x","gene"))) %>%
  pivot_longer(cols = all_of(use_predictors), names_to = "Variable", values_to = "Value") %>%
  mutate(Value = as.numeric(Value), logFC_DC3000 = as.numeric(logFC_DC3000)) %>%
  filter(!is.na(Value) & !is.na(logFC_DC3000))

cat("Rows in plot_data (after pivot and filtering):", nrow(plot_data),
    "unique Variables:", length(unique(plot_data$Variable)), "\n")

if (nrow(plot_data) == 0) stop("No rows to plot after filtering. Check predictor names or data availability.")

# ---------- 7) Compute labels_df with n, R, p and positions ----------
labels_df <- plot_data %>%
  group_by(Variable) %>%
  summarise(
    n_obs = sum(!is.na(Value) & !is.na(logFC_DC3000)),
    R_val = ifelse(n_obs >= 2, cor(Value, logFC_DC3000, use = "complete.obs"), NA_real_),
    p_val = if (n_obs >= 3) broom::tidy(cor.test(Value, logFC_DC3000))$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    label_main = paste0("n=", n_obs, ifelse(!is.na(p_val), paste0("\n", "p=", signif(p_val,2)), "")),
    label_R = ifelse(!is.na(R_val), paste0("R=", signif(R_val,2)), NA_character_)
  )

# sensible label positions (top-right inside each facet)
label_pos <- plot_data %>%
  group_by(Variable) %>%
  summarise(x_pos = max(Value, na.rm = TRUE),
            y_max = max(logFC_DC3000, na.rm = TRUE),
            y_min = min(logFC_DC3000, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(y_span = y_max - y_min,
         y_pos_main = ifelse(y_span>0, y_max - 0.02*y_span, y_max),
         y_pos_R = ifelse(y_span>0, y_max - 0.06*y_span, y_max - 0.03))

labels_df <- labels_df %>% left_join(label_pos, by = "Variable") %>% filter(!is.na(x_pos))

# ---------- 8) Plot (R in red) ----------
p1 <- ggplot(plot_data, aes(x = Value, y = logFC_DC3000)) +
  geom_point(size = 1.0, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, linetype = "dashed", color = "black") +
  facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
  labs(x = NULL, y = expression(log[2] * "FC in DC3000"),
       title = "Predictors vs log2 FC in DC3000") +
  theme_minimal(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3))

# add labels: main (black) and R (red)
if (nrow(labels_df) > 0) {
  p1 <- p1 +
    geom_text(data = labels_df, aes(x = x_pos, y = y_pos_main, label = label_main),
              hjust = 1, vjust = 1, size = 3.0, color = "black", inherit.aes = FALSE) +
    geom_text(data = labels_df %>% filter(!is.na(label_R)), aes(x = x_pos, y = y_pos_R, label = label_R),
              hjust = 1, vjust = 1, size = 3.0, color = "red", inherit.aes = FALSE)
}

print(p1)











########## STOPPED HERE ON 12/15

outdir <- "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
ggsave(file.path(outdir, "predictors_vs_logFC_DC3000_facet_Rlabels.pdf"), plot = p1, width = 9, height = 8)


##################
# This next part takes the genes in DC3000 that are found to be important over time (are signficant), and asks about predictors for those ones.


lfc_dcsig <- lfc_df %>% dplyr::filter(Significant=="Yes")





# 2) Subset fit_exp to genes that are DC3000-significant.
#    In your merges, gene identifier in fit_exp is in column "WP.x" (from earlier merge).
sig_fit <- fit_exp %>%
  filter(!is.na(WP.x) & WP.x %in% lfc_dcsig$gene) %>% 
  # keep only rows with at least one response present
  filter(!is.na(logFC_DC3000) | !is.na(logFC_P25c2))

cat("Number of DC3000-significant genes found in fit_exp:", nrow(sig_fit), "\n")

# Quick table of missingness
sapply(sig_fit[,c("logFC_DC3000","logFC_P25c2","perc_div_aa.x",
                  "Num_gene_events","Genetic_Diversity","Col.0_Pto","Time.in.tree")],
       function(x) sum(is.na(x)))

# -------------------------
# A. Linear model (within DC3000-significant genes)
#    Predicting p25.c2 effect from your candidate predictors
# -------------------------
lm_p25 <- lm(data = sig_fit, logFC_P25c2 ~
               perc_div_aa.x + Num_gene_events + Genetic_Diversity +
               Col.0_Pto + KB_Pto + Time.in.tree  + logFC_DC3000)
summary(lm_p25)
anova_p25 <- car::Anova(lm_p25, type = 3)
print(anova_p25)


# #############################################
# # FACET PLOT OF ALL PREDICTORS vs RESPONSE for the SNPs that are significant in DC3000
# # (fixed: no ambiguous unname() call; robust to all-NA groups)
# #############################################
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(purrr)
# 
# ### Choose dataset:
# # df_to_plot <- fit_exp_clean        # full dataset
# df_to_plot <- sig_fit               # DC3000-significant genes
# 
# ### Choose response variable ("logFC_DC3000_col0" or "logFC_P25c2_col0")
# response <- "logFC_DC3000"
# 
# predictors <- c(
#   "perc_div_aa",
#   "Num_gene_events",
#   "Genetic_Diversity",
#   "Col.0_Pto",
#   "KB_Pto",
#   "Time.in.tree",
#   "interaction",
#   "logFC_P25c2_col0"
# )
# 
# # Keep only available predictors
# predictors <- predictors[predictors %in% names(df_to_plot)]
# 
# # Build long-format data for plotting
# plot_df <- df_to_plot %>%
#   dplyr::select(all_of(c(response, predictors))) %>%
#   filter(!is.na(.data[[response]])) %>%
#   pivot_longer(cols = all_of(predictors),
#                names_to = "Variable",
#                values_to = "Value")
# 
# # Compute per-variable Pearson correlations (robust)
# corr_labels <- plot_df %>%
#   group_by(Variable) %>%
#   summarize(
#     n = sum(!is.na(Value) & !is.na(.data[[response]])),
#     ct = list(if (n >= 3) {
#       tryCatch(cor.test(Value, .data[[response]], use = "complete.obs"),
#                error = function(e) NULL)
#     } else NULL),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     R = map_dbl(ct, ~ if (is.null(.x)) NA_real_ else as.numeric(.x$estimate)),
#     p = map_dbl(ct, ~ if (is.null(.x)) NA_real_ else as.numeric(.x$p.value)),
#     label = map2_chr(R, p, ~ if (is.na(.x)) "n<3" else paste0("R = ", formatC(.x, 2, format = "f"), "\np = ", signif(.y, 2)))
#   )
# 
# # Position labels in upper-right of each facet (robust to all-NA)
# positions <- plot_df %>%
#   group_by(Variable) %>%
#   summarize(
#     x = if (all(is.na(Value))) NA_real_ else max(Value, na.rm = TRUE),
#     y = if (all(is.na(.data[[response]]))) NA_real_ else max(.data[[response]], na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# corr_labels <- left_join(corr_labels, positions, by = "Variable")
# 
# # If position is NA (all NA), set a fallback (so geom_text doesn't error)
# corr_labels <- corr_labels %>%
#   mutate(
#     x = ifelse(is.na(x), 0, x),
#     y = ifelse(is.na(y), 0, y)
#   )
# 
# # ---- Final Faceted Plot ----
# p_all <- ggplot(plot_df, aes(Value, .data[[response]])) +
#   geom_point(alpha = 0.75, size = 1) +
#   geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
#   facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
#   geom_text(data = corr_labels,
#             aes(x = x, y = y, label = label),
#             hjust = 1.05, vjust = 1.1, size = 3.1, na.rm = TRUE) +
#   labs(
#     x = NULL,
#     y = response,
#     title = paste("Predictors vs", response)
#   ) +
#   theme_minimal(base_size = 11) +
#   theme(
#     panel.grid = element_blank(),
#     strip.text = element_text(face = "bold"),
#     axis.line = element_line(color = "black", linewidth = 0.3)
#   )
# 
# print(p_all)
# 
# # Save to file
# ggsave("ALL_predictors_vs_response_facet_plot_fixed.pdf",
#        plot = p_all, width = 9, height = 8)


#Function that returns a sentence describing variance explained

variance_sentence <- function(model, response_name, predictor_name) {
  r2 <- summary(model)$r.squared
  pct <- round(100 * r2, 2)
  
  sentence <- paste0(
    "Variation in ", predictor_name, 
    " explains ", pct, "% of the variance in ", response_name, "."
  )
  
  return(sentence)
}


############################################################
# 1. Variance explained in P25.c2 by DC3000
############################################################

model_p25_by_dc <- lm(logFC_P25c2 ~ logFC_DC3000, data = fit_exp_clean2)

sentence_p25 <- variance_sentence(
  model_p25_by_dc,
  response_name = "logFC_P25c2",
  predictor_name = "logFC_DC3000"
)

cat(sentence_p25, "\n\n")

# ############################################################
# # 2. Variance explained between Col-0 and Ey15 in DC3000
# ############################################################
# # Update column names here if different:
# col0_col <- "logFC_DC3000"
# ey15_col <- "logFC_DC3000"
# 
# model_ey15_by_col0 <- lm(fit_exp_clean2[[ey15_col]] ~ fit_exp_clean2[[col0_col]])
# model_col0_by_ey15 <- lm(fit_exp_clean2[[col0_col]] ~ fit_exp_clean2[[ey15_col]])
# 
# sentence_ey15 <- variance_sentence(
#   model_ey15_by_col0,
#   response_name = "logFC_DC3000_Ey15",
#   predictor_name = "logFC_DC3000"
# )
# 
# sentence_col0 <- variance_sentence(
#   model_col0_by_ey15,
#   response_name = "logFC_DC3000_Col0",
#   predictor_name = "logFC_DC3000_Ey15"
# )
# 
# cat(sentence_ey15, "\n")
# cat(sentence_col0, "\n")
# 
# 

############################################################
# 3. Make a graph of how the local Pearson changes with the FC of DC3000
############################################################

# Moderately-smoothed local curves WITH 95% CI ribbons (LOESS se-based)
# Requires: dplyr, ggplot2
# Paste into R and run (uses sig_fit if present else fit_exp_clean else fit_exp)

library(dplyr)
library(ggplot2)

# ---- USER TUNABLES ----
window_width <- 1.2   # sliding window half-width (moderate smoothing)
loess_span   <- 0.4   # LOESS smoothing span (moderate)
min_points    <- 10   # minimum points in a window to compute local stats
n_grid        <- 250  # resolution for evaluation along DC3000 axis
z_alpha       <- 1.96 # 95% CI ~ 1.96 * se
# ------------------------

# --- select dataset (prefers sig_fit) ---
if (exists("sig_fit")) {
  df <- sig_fit
} else if (exists("fit_exp_clean2")) {
  df <- fit_exp_clean2
} else if (exists("fit_exp")) {
  df <- fit_exp
} else stop("No dataset found: define 'sig_fit' or 'fit_exp_clean' or 'fit_exp'")

# required columns
req <- c("logFC_DC3000", "logFC_P25c2")
miss <- setdiff(req, names(df))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

# prepare data (drop NA)
df2 <- df %>% dplyr::select(all_of(req)) %>% na.omit()

if (nrow(df2) < 30) warning("Fewer than 30 complete points — interpret results with caution.")

# define DC3000 grid
dc_vals <- seq(min(df2$logFC_DC3000, na.rm = TRUE),
               max(df2$logFC_DC3000, na.rm = TRUE),
               length.out = n_grid)

# storage
local_R2    <- rep(NA_real_, n_grid)
local_slope <- rep(NA_real_, n_grid)
local_r     <- rep(NA_real_, n_grid)
local_n     <- rep(0L, n_grid)

# compute local stats using sliding window
for (i in seq_along(dc_vals)) {
  x0 <- dc_vals[i]
  window <- df2 %>% filter(abs(logFC_DC3000 - x0) <= window_width)
  local_n[i] <- nrow(window)
  if (nrow(window) >= min_points) {
    m <- lm(logFC_P25c2 ~ logFC_DC3000, data = window)
    local_R2[i]    <- summary(m)$r.squared
    local_slope[i] <- coef(m)[2]
    local_r[i]     <- suppressWarnings(cor(window$logFC_DC3000, window$logFC_P25c2))
  } else {
    local_R2[i] <- NA_real_
    local_slope[i] <- NA_real_
    local_r[i] <- NA_real_
  }
}

local_df <- data.frame(
  DC3000 = dc_vals,
  n_points = local_n,
  local_R2 = local_R2,
  local_slope = local_slope,
  local_r = local_r
)

# Identify valid rows per metric
ok_R2 <- which(!is.na(local_df$local_R2))
ok_slope <- which(!is.na(local_df$local_slope))
ok_r <- which(!is.na(local_df$local_r))

# ---- Fit LOESS with se=TRUE and compute 95% CI for each metric ----
safe_loess_predict <- function(x, y, newx, span) {
  # returns a list with fit, se.fit (vector of length newx)
  if (length(x) < 10 || all(is.na(y))) {
    return(list(fit = rep(NA_real_, length(newx)), se = rep(NA_real_, length(newx))))
  }
  lo <- tryCatch(loess(y ~ x, span = span, control = loess.control(surface = "interpolate")),
                 error = function(e) NULL)
  if (is.null(lo)) {
    return(list(fit = rep(NA_real_, length(newx)), se = rep(NA_real_, length(newx))))
  }
  pr <- tryCatch(predict(lo, newdata = data.frame(x = newx), se = TRUE),
                 error = function(e) NULL)
  if (is.null(pr)) {
    return(list(fit = rep(NA_real_, length(newx)), se = rep(NA_real_, length(newx))))
  }
  # predict.loess with se=TRUE returns a list with fit and se.fit; handle both possibilities
  fit <- as.numeric(pr$fit)
  se  <- as.numeric(pr$se.fit %||% pr$se) # fallbacks
  list(fit = fit, se = se)
}

# LOESS for R2
if (length(ok_R2) >= 5) {
  res_R2 <- safe_loess_predict(local_df$DC3000[ok_R2], local_df$local_R2[ok_R2], local_df$DC3000, loess_span)
  local_df$R2_fit <- res_R2$fit
  local_df$R2_se <- res_R2$se
  local_df$R2_lo <- local_df$R2_fit - z_alpha * local_df$R2_se
  local_df$R2_hi <- local_df$R2_fit + z_alpha * local_df$R2_se
} else {
  local_df$R2_fit <- NA; local_df$R2_se <- NA; local_df$R2_lo <- NA; local_df$R2_hi <- NA
}

# LOESS for slope
if (length(ok_slope) >= 5) {
  res_slope <- safe_loess_predict(local_df$DC3000[ok_slope], local_df$local_slope[ok_slope], local_df$DC3000, loess_span)
  local_df$slope_fit <- res_slope$fit
  local_df$slope_se <- res_slope$se
  local_df$slope_lo <- local_df$slope_fit - z_alpha * local_df$slope_se
  local_df$slope_hi <- local_df$slope_fit + z_alpha * local_df$slope_se
} else {
  local_df$slope_fit <- NA; local_df$slope_se <- NA; local_df$slope_lo <- NA; local_df$slope_hi <- NA
}

# LOESS for r
if (length(ok_r) >= 5) {
  res_r <- safe_loess_predict(local_df$DC3000[ok_r], local_df$local_r[ok_r], local_df$DC3000, loess_span)
  local_df$r_fit <- res_r$fit
  local_df$r_se <- res_r$se
  local_df$r_lo <- local_df$r_fit - z_alpha * local_df$r_se
  local_df$r_hi <- local_df$r_fit + z_alpha * local_df$r_se
} else {
  local_df$r_fit <- NA; local_df$r_se <- NA; local_df$r_lo <- NA; local_df$r_hi <- NA
}

# clamp R2 CI to sensible range [0,1]
local_df$R2_lo <- pmax(0, local_df$R2_lo)
local_df$R2_hi <- pmin(1, local_df$R2_hi)

# ---- PLOT 1: Local R^2 with CI ribbon ----
p_r2 <- ggplot(local_df, aes(x = DC3000)) +
  geom_line(aes(y = local_R2), color = "grey80", size = 0.5, alpha = 0.6) +
  geom_ribbon(aes(ymin = R2_lo, ymax = R2_hi), fill = "steelblue", alpha = 0.18, inherit.aes = TRUE) +
  geom_line(aes(y = R2_fit), color = "steelblue", size = 1.05) +
  geom_point(data = local_df %>% filter(n_points < min_points), aes(x = DC3000, y = 0), shape = 4, color = "grey50", size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Local R²: DC3000 → P25.c2 (moderately smoothed, 95% CI)",
       subtitle = paste0("window_width = ", window_width, ", LOESS span = ", loess_span, ", min_points = ", min_points),
       x = expression(log[2]*"FC in DC3000"),
       y = expression("Local R"^2)) +
  theme_minimal(base_size = 13)

ggsave("local_R2_moderately_smoothed_CI.pdf", p_r2, width = 7, height = 4.5)
print(p_r2)

# ---- PLOT 2: Local slope with CI ribbon ----
p_slope <- ggplot(local_df, aes(x = DC3000)) +
  geom_line(aes(y = local_slope), color = "grey80", size = 0.5, alpha = 0.6) +
  geom_ribbon(aes(ymin = slope_lo, ymax = slope_hi), fill = "firebrick3", alpha = 0.18) +
  geom_line(aes(y = slope_fit), color = "firebrick3", size = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Local slope of P25.c2 ~ DC3000 (moderately smoothed, 95% CI)",
       subtitle = paste0("window_width = ", window_width, ", LOESS span = ", loess_span),
       x = expression(log[2]*"FC in DC3000"),
       y = "Local slope") +
  theme_minimal(base_size = 13)

ggsave("local_slope_moderately_smoothed_CI.pdf", p_slope, width = 7, height = 4.5)
print(p_slope)

# ---- PLOT 3: Local Pearson r with CI ribbon ----
p_r <- ggplot(local_df, aes(x = DC3000)) +
  geom_line(aes(y = local_r), color = "grey80", size = 0.5, alpha = 0.6) +
  geom_ribbon(aes(ymin = r_lo, ymax = r_hi), fill = "darkgreen", alpha = 0.18) +
  geom_line(aes(y = r_fit), color = "darkgreen", size = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Local Pearson r: DC3000 vs P25.c2 (moderately smoothed, 95% CI)",
       subtitle = paste0("window_width = ", window_width, ", LOESS span = ", loess_span),
       x = expression(log[2]*"FC in DC3000"),
       y = "Local Pearson r") +
  theme_minimal(base_size = 13)

ggsave("local_correlation_moderately_smoothed_CI.pdf", p_r, width = 7, height = 4.5)
print(p_r)

cat("Saved:\n - local_R2_moderately_smoothed_CI.pdf\n - local_slope_moderately_smoothed_CI.pdf\n - local_correlation_moderately_smoothed_CI.pdf\n")
############################################################
# I think that previous section may be useless. This next section uses a GAM-predicted effect modeling.
# Full GAM + derivative pipeline (robust to gratia output column names)
# Uses data frame `fit_exp`
# Predictor: logFC_DC3000 (effect_A)
# Response:  logFC_P25c2 (response_B)
# Covariates: perc_div_aa.x, Num_gene_events, Genetic_Diversity, Col.0_Pto, Time.in.tree

# Combined: GAM predicted curve + parallel bootstrap local R^2 (Option B)
# Assumes `fit_exp` exists and has required columns.
# Packages
# install.packages(c("mgcv","ggplot2","dplyr","parallel"))
library(mgcv)
library(ggplot2)
library(dplyr)
library(parallel)

# ---------- 0) Sanity checks & data prep ----------
stopifnot(exists("fit_exp"))
req_cols <- c("logFC_DC3000","logFC_P25c2","perc_div_aa.x",
              "Num_gene_events","Genetic_Diversity","Col.0_Pto","Time.in.tree")
stopifnot(all(req_cols %in% names(fit_exp)))

# Subset + drop NA and rename
df <- fit_exp %>%
  dplyr::select(all_of(req_cols)) %>%
  na.omit() %>%
  rename(effect_A = logFC_DC3000,
         response_B = logFC_P25c2) %>%
  mutate(effect_A_pos = pmax(effect_A, 0),
         effect_A_neg = -pmin(effect_A, 0))

# Convert Col.0_Pto to factor if it appears categorical
if(!is.factor(df$Col.0_Pto)) {
  if(length(unique(df$Col.0_Pto)) <= 10) df$Col.0_Pto <- factor(df$Col.0_Pto)
  else df$Col.0_Pto <- as.numeric(df$Col.0_Pto)
}

# ---------- 1) Fit the GAM (separate smooths for + and - sides) ----------
fit_gam <- gam(
  response_B ~ s(effect_A_pos, k = 10) + s(effect_A_neg, k = 10) +
    perc_div_aa.x + Num_gene_events + Genetic_Diversity + Col.0_Pto + Time.in.tree,
  data = df,
  family = gaussian(),
  method = "REML"
)

# quick diagnostics (print)
print(summary(fit_gam))
gam.check(fit_gam)

# ---------- 2) Predicted curve plot ----------
grid_n <- 400
grid <- tibble(effect_A = seq(min(df$effect_A, na.rm = TRUE),
                              max(df$effect_A, na.rm = TRUE),
                              length.out = grid_n)) %>%
  mutate(effect_A_pos = pmax(effect_A, 0),
         effect_A_neg = -pmin(effect_A, 0))

# set covariates to typical values
grid$perc_div_aa.x   <- mean(df$perc_div_aa.x, na.rm = TRUE)
grid$Num_gene_events <- mean(df$Num_gene_events, na.rm = TRUE)
grid$Genetic_Diversity <- mean(df$Genetic_Diversity, na.rm = TRUE)
grid$Time.in.tree <- mean(df$Time.in.tree, na.rm = TRUE)
if(is.factor(df$Col.0_Pto)) {
  grid$Col.0_Pto <- factor(names(which.max(table(df$Col.0_Pto))), levels = levels(df$Col.0_Pto))
} else {
  grid$Col.0_Pto <- mean(df$Col.0_Pto, na.rm = TRUE)
}

pred <- predict(fit_gam, newdata = grid, se.fit = TRUE, type = "response")
grid <- grid %>% mutate(fit = pred$fit, se_fit = pred$se.fit,
                        fit_upper = fit + 1.96 * se_fit, fit_lower = fit - 1.96 * se_fit)

p_resp <- ggplot(grid, aes(x = effect_A, y = fit)) +
  geom_ribbon(aes(ymin = fit_lower, ymax = fit_upper), alpha = 0.25) +
  geom_line(size = 1.05) +
  labs(x = "log2 Fold Change in DC3000 (logFC_DC3000)",
       y = "Predicted log2 Fold Change in P25c2 (logFC_P25c2)",
       title = "GAM prediction: response_B vs effect_A") +
  theme_minimal()
print(p_resp)

# ---------- 3) Prepare for faster bootstraped local R^2 (Option B) ----------
# We'll use rolling windows in terms of observation counts (nearest neighbors).
df2 <- df %>% arrange(effect_A) %>% mutate(idx = row_number())
# store model predictions per observation (used in R^2 calc)
df2$pred <- predict(fit_gam, newdata = df2, type = "response")

n_total <- nrow(df2)
message("Total observations: ", n_total)

# TUNABLE PARAMETERS (adjust for speed/precision)
window_n <- min(120, max(10, floor(n_total/8)))  # number of observations per window (example default)
step <- max(1, floor(window_n / 5))              # step in rows between reported centers (saves time)
B <- 200                                         # bootstrap reps per window (reduce to 100 to speed up)
cores <- max(1, detectCores() - 1)               # parallel workers

message("window_n = ", window_n, "; step = ", step, "; B = ", B, "; cores = ", cores)

# center indices for windows (avoid edges too-short)
start_idx <- ceiling(window_n / 2)
end_idx   <- n_total - floor(window_n / 2)
if(end_idx < start_idx) {
  stop("Dataset too small for the chosen window size. Reduce window_n.")
}
centers_idx <- seq(start_idx, end_idx, by = step)

# ---------- helper: compute bootstrap local R2 for a given center index ----------
bootstrap_center <- function(center_idx, df_local, window_n_local, B_local) {
  lo <- center_idx - floor(window_n_local / 2)
  hi <- center_idx + floor(window_n_local / 2)
  sub <- df_local[lo:hi, , drop = FALSE]
  n <- nrow(sub)
  if(n < 8) return(data.frame(center = NA_real_, r2 = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n = n))
  y <- sub$response_B
  p <- sub$pred
  ss_res <- sum((y - p)^2); ss_tot <- sum((y - mean(y))^2)
  r2_obs <- ifelse(ss_tot == 0, NA_real_, 1 - ss_res/ss_tot)
  r2_boot <- numeric(B_local)
  for(b in seq_len(B_local)) {
    idx <- sample.int(n, n, replace = TRUE)
    yb <- sub$response_B[idx]
    pb <- sub$pred[idx]
    ss_res_b <- sum((yb - pb)^2); ss_tot_b <- sum((yb - mean(yb))^2)
    r2_boot[b] <- ifelse(ss_tot_b == 0, NA_real_, 1 - ss_res_b/ss_tot_b)
  }
  r2_boot <- r2_boot[!is.na(r2_boot)]
  if(length(r2_boot) < 10) {
    ci <- c(NA_real_, NA_real_)
  } else {
    ci <- quantile(r2_boot, probs = c(0.025, 0.975), names = FALSE)
  }
  center_effect <- sub$effect_A[floor(n/2) + 1 - ifelse(n %% 2 == 0, 0, 1)] # approximate center effect
  data.frame(center = center_effect, r2 = r2_obs, ci_low = ci[1], ci_high = ci[2], n = n)
}

# ---------- 4) Run bootstrap in parallel ----------
# Export minimal objects to cluster to avoid heavy copying
cl <- makeCluster(cores)
clusterExport(cl, varlist = c("df2","window_n","B","bootstrap_center"), envir = environment())
clusterEvalQ(cl, library(dplyr))

# run
res_list <- parLapply(cl, centers_idx, function(ci) {
  bootstrap_center(ci, df2, window_n, B)
})
stopCluster(cl)

res_df <- bind_rows(res_list) %>% filter(!is.na(center))

# Remove extreme NA rows if any
res_df <- res_df %>% arrange(center)

# ---------- 5) Plot local R^2 with bootstrap CI ----------
p_localR2 <- ggplot(res_df, aes(x = center, y = r2)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey80", alpha = 0.8) +
  geom_line(color = "black") +
  geom_point(size = 0.9) +
  labs(x = "log2 Fold Change in DC3000 (effect_A)",
       y = expression(Local~R^2),
       title = paste0("Local R^2 (window_n=", window_n, ", step=", step, ", B=", B, ")")) +
  theme_minimal()
print(p_localR2)

# optional: plot sample size per window
p_n <- ggplot(res_df, aes(x = center, y = n)) + geom_line() +
  labs(x = "effect_A", y = "N in window", title = "Window sample size") + theme_minimal()
print(p_n)

# ---------- 6) Save results (optional) ----------
# write.csv(res_df, "local_R2_results.csv", row.names = FALSE)

message("Done. Adjust window_n, step or B to trade accuracy vs speed.")

#############################################
# Now I want variable importance for the GAM model  
# Variable explanatory-power table + importance plot for GAM covariates
# Requires: mgcv, dplyr, ggplot2, broom
# Assumes df and fit_gam already exist (df contains response_B, effect_A, covariates)
# install.packages(c("mgcv","dplyr","ggplot2","broom"))

# Covariate importance including DC3000 logFC (effect_A)
# Assumes df (with response_B, effect_A, covariates) and fit_gam exist
# Packages required: mgcv, dplyr, ggplot2, broom
# install.packages(c("mgcv","dplyr","ggplot2","broom"))

library(mgcv)
library(dplyr)
library(ggplot2)
library(broom)

# ---------- sanity ----------
stopifnot(exists("df"))
stopifnot(exists("fit_gam"))

# Define covariates to evaluate; include "effect_A" as a single item representing both smooths
covars <- c("effect_A", "perc_div_aa.x","Num_gene_events","Genetic_Diversity","Col.0_Pto","Time.in.tree")

# Full model info
full <- fit_gam
full_deviance <- deviance(full)
full_aic <- AIC(full)
full_preds <- predict(full, newdata = df, type = "response")
full_rmse <- sqrt(mean((df$response_B - full_preds)^2))

# Helper: remove terms from formula robustly
# var can be:
#  - "effect_A": remove both s(effect_A_pos) and s(effect_A_neg)
#  - other covariate names: remove top-level terms that match the var name exactly
remove_terms_for_var <- function(model, var) {
  f <- formula(model)
  terms_full <- attr(terms(f), "term.labels")
  if(var == "effect_A") {
    # terms to remove: any smooth or term containing effect_A_pos or effect_A_neg
    to_remove <- terms_full[grepl("effect_A_pos|effect_A_neg|s\\(effect_A", terms_full)]
  } else {
    # remove exact matches of variable as a top-level term (word boundary)
    to_remove <- terms_full[grepl(paste0("\\b", var, "\\b"), terms_full)]
  }
  # if nothing to remove, return NULL
  if(length(to_remove) == 0) return(NULL)
  to_keep <- setdiff(terms_full, to_remove)
  # If smoothing terms remain but we removed everything (i.e., empty), keep intercept-only model
  if(length(to_keep) == 0) {
    # create intercept-only formula
    as.formula(paste(as.character(f[[2]]), "~ 1"))
  } else {
    reformulate(to_keep, response = as.character(f[[2]]))
  }
}

results <- list()

for(v in covars) {
  reduced_formula <- remove_terms_for_var(full, v)
  if(is.null(reduced_formula)) {
    warning("Variable '", v, "' not found as a top-level term (or special handling) in model; skipping.")
    next
  }
  # Fit reduced model with same family and method
  reduced <- gam(reduced_formula, data = df, family = family(full)$family, method = "REML")
  
  # Partial R2 based on deviance (RSS for Gaussian)
  rss_full <- deviance(full)
  rss_red  <- deviance(reduced)
  partial_R2 <- ifelse(rss_red <= 0, NA_real_, (rss_red - rss_full) / rss_red)
  
  # AIC difference
  aic_full <- AIC(full)
  aic_red  <- AIC(reduced)
  delta_aic <- aic_red - aic_full   # positive => full better
  
  # RMSE change
  pred_red <- predict(reduced, newdata = df, type = "response")
  rmse_red <- sqrt(mean((df$response_B - pred_red)^2))
  delta_rmse <- rmse_red - full_rmse   # positive => removing v made RMSE worse
  
  # ANOVA comparison (F test)
  an <- try(anova(reduced, full, test = "F"), silent = TRUE)
  pval <- NA_real_; fstat <- NA_real_
  if(!inherits(an, "try-error")) {
    an_tbl <- broom::tidy(an)
    if(nrow(an_tbl) >= 1) {
      if("p.value" %in% names(an_tbl)) pval <- an_tbl$p.value[nrow(an_tbl)]
      else if("Pr..F." %in% names(an_tbl)) pval <- an_tbl$Pr..F.[nrow(an_tbl)]
      if("statistic" %in% names(an_tbl)) fstat <- an_tbl$statistic[nrow(an_tbl)]
    }
  }
  
  results[[v]] <- data.frame(
    variable = v,
    partial_R2 = partial_R2,
    delta_AIC = delta_aic,
    rmse_full = full_rmse,
    rmse_reduced = rmse_red,
    delta_RMSE = delta_rmse,
    anova_p = pval,
    anova_F = fstat,
    stringsAsFactors = FALSE
  )
}

res_df <- bind_rows(results)

# Formatting and sorting
res_df <- res_df %>%
  arrange(desc(partial_R2)) %>%
  mutate(
    partial_R2 = round(partial_R2, 4),
    delta_AIC = round(delta_AIC, 3),
    delta_RMSE = round(delta_RMSE, 4),
    anova_p = signif(anova_p, 3)
  )

# Print and save table
print(res_df)
write.csv(res_df, "covariate_explanatory_power_table_with_effectA.csv", row.names = FALSE)

# ---------- importance barplot including effect_A ----------
# Use partial R2 as main metric; include delta_AIC as label
p_bar <- ggplot(res_df, aes(x = reorder(variable, partial_R2), y = partial_R2)) +
  geom_col(fill = "#2b8cbe") +
  coord_flip() +
  geom_text(aes(label = paste0("ΔAIC=", delta_AIC)), hjust = -0.05, size = 3) +
  labs(x = "", y = "Partial R² (variance uniquely explained)",
       title = "Covariate explanatory power (GAM) — including DC3000 logFC (effect_A)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Expand y limit a bit to fit labels
ylim_max <- max(res_df$partial_R2, na.rm = TRUE) * 1.25
if(is.finite(ylim_max) && ylim_max > 0) p_bar <- p_bar + ylim(0, ylim_max)

print(p_bar)
ggsave("covariate_importance_barplot_with_effectA.png", p_bar, width = 7, height = 4, dpi = 300)

# done
message("Saved results to covariate_explanatory_power_table_with_effectA.csv and covariate_importance_barplot_with_effectA.png")
