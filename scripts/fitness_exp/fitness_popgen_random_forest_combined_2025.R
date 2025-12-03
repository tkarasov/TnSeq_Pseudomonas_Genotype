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

######### Actually I pulled in the lfc datatable from limmaVoom
lfc_df <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/logFCDC3000logFCp25c2_10_2025.csv")

#this is the output from running the python evolutionary genetics scripts on the panX output
conservation = read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd.csv", row.names = 1)

#this is the output from Effie's blasting (from which the orthologs were originally called). This is where we pull the sequence divergence.
blast_div <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit_2.txt")
blast_div <- blast_div[,c("V1", "V2", "V11")]
colnames(blast_div) <- c("DAK","WP", "perc_div_aa")
rownames(blast_div) <- blast_div$WP

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
fit_exp$scaled_logFC <- scale(fit_exp$logFC_DC3000_col0)
fit_exp$interaction <- fit_exp$scaled_div * fit_exp$scaled_logFC
fit_exp_clean <- fit_exp %>% 
  na.omit()

#first let's graph against all the variables of interest
# Full code: faceted predictor plots with R / p labels in upper-right corner
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Clean, fixed, end-to-end plotting block
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# --- Ensure required columns exist; create numeric interaction if missing ---
if (!"interaction" %in% names(fit_exp_clean)) {
  if (all(c("logFC_P25c2_col0", "perc_div_aa") %in% names(fit_exp_clean))) {
    fit_exp_clean <- fit_exp_clean %>%
      mutate(interaction = logFC_P25c2_col0 * perc_div_aa)
  } else {
    stop("Neither 'interaction' column exists nor the necessary columns to create it (logFC_P25c2_col0, perc_div_aa) are present.")
  }
}

# --- Define predictor variables (as you provided) ---
predictor_vars <- c("perc_div_aa", "Num_gene_events", "Genetic_Diversity",
                    "Col.0_Pto", "KB_Pto", "interaction",
                    "logFC_P25c2_col0", "Time.in.tree")

# --- Prepare long-format plotting data (drop rows lacking the selected columns) ---
needed_cols <- unique(c("logFC_DC3000_col0", predictor_vars))
missing_cols <- setdiff(needed_cols, names(fit_exp_clean))
if (length(missing_cols) > 0) {
  stop("Missing required columns in fit_exp_clean: ", paste(missing_cols, collapse = ", "))
}

plot_data <- fit_exp_clean %>%
  dplyr::select(all_of(needed_cols)) %>%
  # keep rows that have at least the response present (we'll let cor.test handle cases per variable)
  filter(!is.na(logFC_DC3000_col0)) %>%
  pivot_longer(cols = -logFC_DC3000_col0,
               names_to = "Variable",
               values_to = "Value")

# --- Compute per-Variable Pearson correlation tests safely ---
labels_df <- plot_data %>%
  group_by(Variable) %>%
  summarize(
    n_obs = sum(!is.na(Value) & !is.na(logFC_DC3000_col0)),
    cor_test = list(if (n_obs >= 3) {
      tryCatch(cor.test(Value, logFC_DC3000_col0, use = "complete.obs"),
               error = function(e) NULL)
    } else NULL),
    .groups = "drop"
  ) %>%
  mutate(
    R = map_dbl(cor_test, ~ if (is.null(.x)) NA_real_ else as.numeric(.x$estimate)),
    p = map_dbl(cor_test, ~ if (is.null(.x)) NA_real_ else as.numeric(.x$p.value)),
    label = map2_chr(R, p, ~ {
      if (is.na(.x)) {
        "insufficient\npoints"
      } else {
        paste0("R = ", formatC(.x, digits = 2, format = "f"), "\n", "p = ", signif(.y, 2))
      }
    })
  )

# --- Robust upper-right label positions per facet (5% inset) ---
pos_df <- plot_data %>%
  group_by(Variable) %>%
  summarize(
    x_max = if (all(is.na(Value))) NA_real_ else max(Value, na.rm = TRUE),
    x_min = if (all(is.na(Value))) NA_real_ else min(Value, na.rm = TRUE),
    y_max = if (all(is.na(logFC_DC3000_col0))) NA_real_ else max(logFC_DC3000_col0, na.rm = TRUE),
    y_min = if (all(is.na(logFC_DC3000_col0))) NA_real_ else min(logFC_DC3000_col0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_span = x_max - x_min,
    y_span = y_max - y_min,
    # if span is zero/NA, fall back to a small absolute offset so the label can be placed
    x_span = ifelse(is.na(x_span) | x_span == 0,
                    ifelse(is.na(x_max), NA_real_, max(abs(x_max), 1) * 0.01 + 1e-6),
                    x_span),
    y_span = ifelse(is.na(y_span) | y_span == 0,
                    ifelse(is.na(y_max), NA_real_, max(abs(y_max), 1) * 0.01 + 1e-6),
                    y_span),
    inset_frac = 0.05,  # 5% inset from the top-right corner
    x_pos = ifelse(is.na(x_max), NA_real_, x_max - inset_frac * x_span),
    y_pos = ifelse(is.na(y_max), NA_real_, y_max - inset_frac * y_span)
  ) %>%
  select(Variable, x_pos, y_pos)

# --- Join label text with positions ---
labels_df <- labels_df %>% left_join(pos_df, by = "Variable")

# --- Build and display the faceted plot ---
p <- ggplot(plot_data, aes(x = Value, y = logFC_DC3000_col0)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, linetype = "dashed", color = "black") +
  facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
  geom_text(
    data = labels_df,
    aes(x = x_pos, y = y_pos, label = label),
    hjust = 1, vjust = 1, size = 3.1,
    na.rm = TRUE
  ) +
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

print(p)


# outdir <- "path/to/outdir"
# if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
# ggsave(file.path(outdir, "predictors_vs_logFC_DC3000_facet_Rlabels.pdf"), plot = p, width = 9, height = 8)


















