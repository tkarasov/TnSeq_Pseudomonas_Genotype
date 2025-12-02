# The goal of this script is to take lists of gene names and do gene ontology enrichments.
#Conclusion as of 7/2025 is that there are enrichments of real basal processes, but not any story that makes a lot of sense#Coming back to this on 11/20/2025.
# first step is to write a function that takes two lists and gets enrichments
setwd("~/Google\ Drive/My\ Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")

#sig_genes is the output from the limmavoom file for which genes were significant in which interaction.
result_df = read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/sig_results_11_2025.csv", header = TRUE)

# ortholog table for possible genes
ortholog_tab <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", header = TRUE, row.names = 1, sep = ",")

#!/usr/bin/env Rscript
# go_and_heatmap_top50.R
# Combined GO ORA + heatmap pipeline. Heatmaps limited to top 50 genes per set.

# ---------------------- User settings ----------------------
lfc_path  <- "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/sig_results_11_2025_with_uniprot.csv"
go_map_path <- "/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/DC3000_genome_mappings/gene_ontology_mapping_uniprotkb_proteome_UP000002515_2025_07_01 (2).tsv"

# Parameters
fdr_thresh      <- 0.05
minGSSize       <- 3
showCategory    <- 20
top_n_heatmap   <- 50      # <<-- LIMIT: top 50 genes per set for heatmaps
scale_rows      <- TRUE
wrap_width_labels <- 40
output_prefix   <- "ORA_results_top50"

# ---------------------- Libraries & conflict handling ----------------------
suppressPackageStartupMessages({
  if(!requireNamespace("conflicted", quietly = TRUE)) install.packages("conflicted", repos = "https://cloud.r-project.org")
  if(!requireNamespace("png", quietly = TRUE)) install.packages("png", repos = "https://cloud.r-project.org")
  library(conflicted)
  library(dplyr)
  library(stringr)
  library(clusterProfiler)
  library(enrichplot)
  library(GO.db)
  library(pheatmap)
  library(ggplot2)
  library(readr)
  library(grid)
  library(png)
# Parameters
top_n_genes     <- 25     # top 25 genes chosen by padj_strain_time_col0
scale_rows      <- TRUE   # row-wise z-score for heatmaps
max_label_length <- 60    # truncate label length for readability
output_prefix   <- "strainxtime_top25"

# padj / logFC column names used in your dataset
padj_strain_time_col <- "padj_strain_time_col0"
logfc_strain_time_col <- "logFC_strain_time_col0"

# Other logFC columns will be detected automatically (columns starting with "logFC_")
# but we'll reorder them so strain_time comes first (if found)
# ---------------------- Libraries & conflict handling ----------------------
suppressPackageStartupMessages({
  if(!requireNamespace("conflicted", quietly = TRUE)) install.packages("conflicted", repos = "https://cloud.r-project.org")
  if(!requireNamespace("png", quietly = TRUE)) install.packages("png", repos = "https://cloud.r-project.org")
  library(conflicted)
  library(dplyr)
  library(stringr)
  library(pheatmap)
  library(ggplot2)
  library(grid)
  library(png)
  library(readr)
})
# enforce dplyr verbs and base intersect
conflict_prefer("select",  "dplyr")
conflict_prefer("filter",  "dplyr")
conflict_prefer("mutate",  "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("rename",  "dplyr")
conflict_prefer("intersect", "base")

# ---------------------- Read input LFC table ----------------------
message("Reading LFC table from: ", lfc_path)
lfc_file <- tryCatch(read.csv(lfc_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
                     error = function(e) stop("Failed to read lfc_file: ", e$message))
if(!"uniprot" %in% colnames(lfc_file)) stop("lfc_file must contain a column named 'uniprot' with UniProt Entry IDs.")

# try reading go_mappings (only for row labels if present)
go_mappings <- NULL
if(file.exists(go_map_path)) {
  go_mappings <- tryCatch(read.delim(go_map_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE),
                          error = function(e) { message("Could not read go_map_path: ", e$message); return(NULL) })
}

# ---------------------- Detect logFC columns and reorder so strain_time is first ----------------------
logfc_cols_all <- grep("^logFC_", colnames(lfc_file), value = TRUE)
if(length(logfc_cols_all) == 0) stop("No columns starting with 'logFC_' found in lfc_file.")

# find the strain×time logFC column and place it first if present
if(logfc_strain_time_col %in% logfc_cols_all) {
  logfc_cols_ordered <- c(logfc_strain_time_col, setdiff(logfc_cols_all, logfc_strain_time_col))
} else {
  # If expected column name not found, attempt to find something similar
  candidate <- grep("strain.*time|strain_time|strain.time", logfc_cols_all, ignore.case = TRUE, value = TRUE)
  if(length(candidate) >= 1) {
    logfc_cols_ordered <- c(candidate[1], setdiff(logfc_cols_all, candidate[1]))
    logfc_strain_time_col <- candidate[1]
    message("Using detected strain×time logFC column: ", candidate[1])
  } else {
    # no strain_time column found — keep default ordering
    logfc_cols_ordered <- logfc_cols_all
    message("No strain×time logFC column found; columns will not be reordered.")
  }
}
message("logFC column order (first column intended as strain×time):\n", paste(logfc_cols_ordered, collapse = ", "))

# ---------------------- Select top 25 genes by padj_strain_time_col ----------------------
if(!padj_strain_time_col %in% colnames(lfc_file)) {
  # try to find a similar padj column name if exact name not found
  cand_padj <- grep("padj.*strain.*time|padj.*strain_time|padj_strain_time", colnames(lfc_file), ignore.case = TRUE, value = TRUE)
  if(length(cand_padj) > 0) {
    padj_strain_time_col <- cand_padj[1]
    message("Using detected padj column for Strain×Time: ", padj_strain_time_col)
  } else {
    stop("Cannot find padj column for Strain×Time. Expected '", padj_strain_time_col, "' or similar.")
  }
}

# filter rows with available padj and uniprot, then pick top 25 by smallest padj
sig_df <- lfc_file %>%
  filter(!is.na(.data[[padj_strain_time_col]]), !is.na(uniprot), uniprot != "") %>%
  arrange(.data[[padj_strain_time_col]]) %>%
  distinct(uniprot, .keep_all = TRUE)  # ensure unique uniprot

if(nrow(sig_df) == 0) stop("No rows with padj for Strain×Time found in your data.")

top25_df <- sig_df %>% slice_head(n = top_n_genes)
top25_uniprot <- top25_df$uniprot
cat("Selected top", length(top25_uniprot), "genes by", padj_strain_time_col, "\n")
print(top25_df[, c("uniprot", padj_strain_time_col)][1:min(10, nrow(top25_df)), ])

# ---------------------- Build matrix for the two heatmaps (same genes)
#  - Matrix columns ordered with Strain×Time logFC first
#  - Heatmap A: strain×time top25 (shows the same matrix but emphasises the strain×time column)
#  - Heatmap B: Time view for the same top25 (same matrix)
# ------------------------------------------------------------------------------

# extract matrix rows for top25, columns = logfc_cols_ordered
mat_df <- lfc_file %>%
  filter(uniprot %in% top25_uniprot) %>%
  distinct(uniprot, .keep_all = TRUE) %>%
  select(uniprot, all_of(logfc_cols_ordered))

# ensure numeric
mat_df[ , logfc_cols_ordered] <- lapply(mat_df[ , logfc_cols_ordered, drop = FALSE], function(x) as.numeric(as.character(x)))
# align rows in the same order as top25_uniprot
mat_df <- mat_df[match(top25_uniprot, mat_df$uniprot), ]
rownames_mat <- mat_df$uniprot
mat <- as.matrix(mat_df[, logfc_cols_ordered, drop = FALSE])
rownames(mat) <- rownames_mat

# ---------------------- Row labels: Protein name (UniProt) if available; else UniProt ID
row_labels <- rownames(mat)
if(!is.null(go_mappings) && "Entry" %in% colnames(go_mappings)) {
  prot_col <- NULL
  if("Protein names" %in% colnames(go_mappings)) prot_col <- "Protein names"
  if(is.null(prot_col) && "Protein.names" %in% colnames(go_mappings)) prot_col <- "Protein.names"
  if(!is.null(prot_col)) {
    prot_lookup <- go_mappings %>% select(Entry, !!sym(prot_col)) %>% distinct(Entry, .keep_all = TRUE)
    prot_vec <- prot_lookup[[prot_col]]
    names(prot_vec) <- prot_lookup$Entry
    lab <- prot_vec[rownames(mat)]
    lab[is.na(lab) | lab == ""] <- rownames(mat)[is.na(lab) | lab == ""]
    row_labels <- make.unique(paste0(lab, " (", rownames(mat), ")"))
  }
}
# truncate long labels and wrap so they're readable
row_labels <- sapply(row_labels, function(x) {
  x2 <- if(nchar(x) > max_label_length) paste0(substr(x,1,max_label_length-3), "...") else x
  stringr::str_wrap(x2, width = 60)
}, USE.NAMES = FALSE)

# ---------------------- Optionally scale rows (z-score) ----------------------
mat_to_plot <- mat
if(scale_rows) {
  mat_scaled <- t(scale(t(mat_to_plot)))
  mat_scaled[is.na(mat_scaled)] <- 0
  mat_to_plot <- mat_scaled
}

# ---------------------- Build pheatmap objects (silent) ----------------------
breaks <- seq(-3, 3, length.out = 101)
colors <- colorRampPalette(c("#67001F", "white", "#053061"))(100)

ph_time_for_top25 <- pheatmap::pheatmap(mat_to_plot,
                                        cluster_rows = TRUE,
                                        cluster_cols = TRUE,
                                        show_rownames = TRUE,
                                        labels_row = row_labels,
                                        fontsize_row = if(nrow(mat_to_plot) <= 25) 10 else 8,
                                        fontsize_col = 10,
                                        color = colors,
                                        breaks = breaks,
                                        main = "Top25 (by Strain×Time padj) — All contrasts",
                                        border_color = NA,
                                        silent = TRUE)

# We'll keep the same pheatmap for strain×time (same matrix) but change title (emphasize strain×time column)
ph_strain_time_top25 <- ph_time_for_top25  # same gtable; we will save twice with different titles if needed

# ---------------------- Save PNG + rasterized PDF helper ----------------------
save_ph_as_png_and_raster_pdf <- function(ph_object, png_file, pdf_file, width_in = 10, height_in = 7, dpi = 300) {
  if(is.null(ph_object)) { message("ph_object NULL, skipping: ", png_file); return(NULL) }
  png(png_file, width = width_in, height = height_in, units = "in", res = dpi)
  grid::grid.newpage()
  grid::grid.draw(ph_object$gtable)
  dev.off()
  img <- png::readPNG(png_file)
  pdf(pdf_file, width = width_in, height = height_in)
  grid::grid.newpage()
  grid::grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
  dev.off()
  message("Wrote PNG:", png_file, " (", round(file.info(png_file)$size/1024), "KB )")
  message("Wrote PDF:", pdf_file, " (", round(file.info(pdf_file)$size/1024), "KB )")
  invisible(list(png = png_file, pdf = pdf_file))
}

# Determine height depending on number of rows for readability
height_in <- max(5, 0.32 * max(1, nrow(mat_to_plot)) + 2)

# Save "Strain×Time" heatmap (same matrix but title emphasizes Strain×Time)
save_ph_as_png_and_raster_pdf(ph_strain_time_top25,
                              png_file = paste0(output_prefix, "_strainxtime_top25.png"),
                              pdf_file = paste0(output_prefix, "_strainxtime_top25.pdf"),
                              width_in = 10, height_in = height_in, dpi = 300)

# Save "Time" heatmap for the same genes (title adjusted)
save_ph_as_png_and_raster_pdf(ph_time_for_top25,
                              png_file = paste0(output_prefix, "_time_for_top25.png"),
                              pdf_file = paste0(output_prefix, "_time_for_top25.pdf"),
                              width_in = 10, height_in = height_in, dpi = 300)

# ---------------------- Save TSV of selected genes (with key columns) ----------------------
# Write the top25 table with selected key columns (including padj for Strain×Time and any key logFCs)
save_cols <- c("uniprot", padj_strain_time_col, logfc_strain_time_col)
# include other detected logFCs if present
other_lfc <- setdiff(logfc_cols_ordered, logfc_strain_time_col)
save_cols <- unique(c(save_cols, other_lfc))
# keep only those columns that exist in lfc_file
save_cols <- save_cols[save_cols %in% colnames(lfc_file)]
top25_out <- lfc_file %>% filter(uniprot %in% top25_uniprot) %>% distinct(uniprot, .keep_all = TRUE) %>% select(all_of(save_cols))
# arrange by padj
if(padj_strain_time_col %in% colnames(top25_out)) top25_out <- top25_out %>% arrange(.data[[padj_strain_time_col]])
write.table(top25_out, paste0(output_prefix, "_top25_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote top25 table: ", paste0(output_prefix, "_top25_table.tsv"))

message("Done. Files produced:")
message("- ", paste0(output_prefix, "_strainxtime_top25.png"))
message("- ", paste0(output_prefix, "_strainxtime_top25.pdf"))
message("- ", paste0(output_prefix, "_time_for_top25.png"))
message("- ", paste0(output_prefix, "_time_for_top25.pdf"))
message("- ", paste0(output_prefix, "_top25_table.tsv"))

