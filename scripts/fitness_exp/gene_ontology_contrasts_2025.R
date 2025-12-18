# The goal of this script is to take lists of gene names and do gene ontology enrichments.
#Conclusion as of 7/2025 is that there are enrichments of real basal processes, but not any story that makes a lot of sense#Coming back to this on 11/20/2025.
# first step is to write a function that takes two lists and gets enrichments
setwd("~/Google\ Drive/My\ Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")

#sig_genes is the output from the limmavoom file for which genes were significant in which interaction.
#result_df = read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/sig_results_11_2025.csv", header = TRUE)
sig_results <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/sig_results_12_2025.csv",header=TRUE)


# go_and_heatmap_top50.R
# Combined GO ORA + heatmap pipeline. Heatmaps limited to top 50 genes per set.

# ---------------------- User settings ----------------------
# lfc_path  <- "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/sig_results_11_2025_with_uniprot.csv"
go_map_path <- "/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/DC3000_genome_mappings/gene_ontology_mapping_uniprotkb_proteome_UP000002515_2025_07_01 (2).tsv"

# On 12/17/2025, switching to this LFC file but need to triple check this file has the right mappings.
lfc_df <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000logFC_time_DC3000_vs_P25c2_12_2025.csv", header=TRUE)

# Parameters
fdr_thresh      <- 0.05
minGSSize       <- 3
showCategory    <- 20
top_n_heatmap   <- 50      # <<-- LIMIT: top 50 genes per set for heatmaps
scale_rows      <- TRUE
wrap_width_labels <- 40
output_prefix   <- "ORA_results_top50"


  # ---------------------- User settings ----------------------
  lfc_path  <- "/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_8_2025/output/sig_results_11_2025_with_uniprot.csv"

  go_map_path <- "/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/DC3000_genome_mappings/gene_ontology_mapping_uniprotkb_proteome_UP000002515_2025_07_01 (2).tsv"
  
  # Parameters
  fdr_thresh       <- 0.05
  minGSSize        <- 3
  showCategory     <- 20
  top_n_heatmap    <- 25      # <<-- LIMIT: top 25 genes per set for heatmaps
  scale_rows       <- TRUE
  max_label_length <- 60      # truncate labels to this many chars for readability
  output_prefix    <- "ORA_results_top25_readable_setdiff_fixed"
  
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
  })
  # prefer dplyr verbs and base::intersect / base::setdiff
  conflict_prefer("select",  "dplyr")
  conflict_prefer("filter",  "dplyr")
  conflict_prefer("mutate",  "dplyr")
  conflict_prefer("arrange", "dplyr")
  conflict_prefer("rename",  "dplyr")
  conflict_prefer("intersect", "base")
  conflict_prefer("setdiff", "base")   # <--- FIX: prefer base::setdiff
  
  # ---------------------- Read inputs ----------------------
  message("Reading LFC table from: ", lfc_path)
  lfc_file <- tryCatch(read.csv(lfc_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
                       error = function(e) stop("Failed to read lfc_file: ", e$message))
  if(!"uniprot" %in% colnames(lfc_file)) stop("lfc_file must contain a column named 'uniprot' with UniProt Entry IDs.")
  
  message("Reading GO mapping from: ", go_map_path)
  go_mappings <- tryCatch(read.delim(go_map_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE),
                          error = function(e) { warning("Failed to read go_mappings: ", e$message); return(NULL) })
  
  # ---------------------- Build TERM2GENE (auto-detect GO columns) ----------------------
  message("Scanning GO mapping columns for 'GO:' occurrences (if go_mappings available)...")
  if(!is.null(go_mappings)) {
    go_col_scores <- sapply(colnames(go_mappings), function(col) {
      vals <- as.character(go_mappings[[col]])
      sum(grepl("GO:\\d{4,7}", vals, perl = TRUE), na.rm = TRUE)
    })
    go_col_scores_df <- data.frame(column = names(go_col_scores), matches = as.integer(go_col_scores), stringsAsFactors = FALSE)
    go_col_scores_df <- go_col_scores_df[order(-go_col_scores_df$matches), ]
    print(head(go_col_scores_df, 20))
    candidate_cols <- names(go_col_scores)[go_col_scores > 0]
    if(length(candidate_cols) == 0) {
      candidate_cols <- colnames(go_mappings)[grepl("Gene Ontology|Gene.Ontology|Gene Ontology IDs|Gene Ontology \\(", colnames(go_mappings), ignore.case = TRUE)]
      if(length(candidate_cols) == 0) candidate_cols <- character(0)
      if(length(candidate_cols)>0) warning("Using name-based GO candidate columns: ", paste(candidate_cols, collapse = ", "))
    }
  } else {
    candidate_cols <- character(0)
  }
  message("GO columns detected: ", ifelse(length(candidate_cols)>0, paste(candidate_cols, collapse = ", "), "<none detected>"))
  
  # Build TERM2GENE only if go_mappings present
  term2gene <- data.frame(GO=character(0), uniprot=character(0), stringsAsFactors = FALSE)
  if(!is.null(go_mappings) && length(candidate_cols) > 0) {
    entry_col_name <- if("Entry" %in% colnames(go_mappings)) "Entry" else colnames(go_mappings)[1]
    extract_go_ids <- function(x) {
      if(is.na(x) || x == "") return(character(0))
      s <- as.character(x)
      s <- gsub("\\|", ";", s)
      s <- gsub(",", ";", s)
      parts <- unlist(strsplit(s, ";\\s*"))
      parts <- trimws(parts)
      parts[grepl("^GO:\\d+", parts)]
    }
    go_list_combined <- vector("list", nrow(go_mappings))
    names(go_list_combined) <- as.character(go_mappings[[entry_col_name]])
    for(i in seq_len(nrow(go_mappings))) {
      collected <- character(0)
      for(col in candidate_cols) collected <- c(collected, extract_go_ids(go_mappings[[col]][i]))
      collected <- unique(collected)
      go_list_combined[[i]] <- collected
    }
    keep <- lengths(go_list_combined) > 0
    go_list_filtered <- go_list_combined[keep]
    if(length(go_list_filtered) > 0) {
      go_df_combined <- stack(go_list_filtered)
      colnames(go_df_combined) <- c("GO", "uniprot")
      term2gene <- as.data.frame(go_df_combined, stringsAsFactors = FALSE) %>% filter(!is.na(GO), GO != "", !is.na(uniprot), uniprot != "") %>% distinct(GO, uniprot)
      cat("TERM2GENE built:", nrow(term2gene), "rows; unique GO terms:", length(unique(term2gene$GO)), "; unique uniprot entries:", length(unique(term2gene$uniprot)), "\n")
    } else {
      message("No GO IDs found in go_mappings candidate columns.")
    }
  } else {
    message("Skipping TERM2GENE build (no go_mappings or no candidate GO columns).")
  }
  
  # ---------------------- ORA wrapper (safe) ----------------------
  run_enricher_with_desc_safe <- function(gene_vec, universe, term2gene_df, minGSSize = minGSSize, qvalueCutoff = fdr_thresh) {
    if(length(gene_vec) < 1) { message("Empty gene set: returning NULL"); return(NULL) }
    if(nrow(term2gene_df) == 0) { message("Empty TERM2GENE: cannot run enricher"); return(NULL) }
    e <- tryCatch({
      clusterProfiler::enricher(gene = gene_vec, universe = universe, TERM2GENE = term2gene_df,
                                pAdjustMethod = "BH", minGSSize = minGSSize, qvalueCutoff = qvalueCutoff, maxGSSize = 5000)
    }, error = function(err) { message("enricher() failed: ", err$message); return(NULL) })
    if(is.null(e)) return(NULL)
    resdf <- as.data.frame(e)
    if(nrow(resdf) == 0) return(NULL)
    desc_vec <- sapply(resdf$ID, function(gid) { term <- GOTERM[[gid]]; if(is.null(term)) NA_character_ else Term(term) }, USE.NAMES = TRUE)
    desc_vec[is.na(desc_vec)] <- names(desc_vec)[is.na(desc_vec)]
    resdf$description <- as.character(desc_vec)
    resdf$description_wrapped <- stringr::str_wrap(resdf$description, width = 60)
    resdf$ontology <- sapply(resdf$ID, function(gid) { t <- GOTERM[[gid]]; if(is.null(t)) NA_character_ else Ontology(t) }, USE.NAMES = FALSE)
    e_result <- e@result
    match_idx <- match(e_result$ID, resdf$ID)
    desc_full <- rep(NA_character_, nrow(e_result))
    matched <- !is.na(match_idx)
    if(any(matched)) desc_full[matched] <- resdf$description_wrapped[match_idx[matched]]
    desc_full[is.na(desc_full)] <- e_result$ID[is.na(desc_full)]
    e@result$Description <- desc_full
    return(list(enrich_result = e, df = resdf))
  }
  
  # ---------------------- Prepare universe & sig lists ----------------------
  lfc_file$uniprot <- as.character(lfc_file$uniprot)
  unique_uniprot_in_lfc <- unique(na.omit(lfc_file$uniprot))
  if(nrow(term2gene) > 0) {
    universe_uniprot <- base::intersect(unique_uniprot_in_lfc, unique(term2gene$uniprot))
  } else {
    universe_uniprot <- unique_uniprot_in_lfc
  }
  cat("Unique uniprot in lfc_file:", length(unique_uniprot_in_lfc), " - universe size:", length(universe_uniprot), "\n")
  if(length(universe_uniprot) == 0) stop("No universe uniprot IDs available. Check inputs.")
  
  padj_time_col <- "padj_DC3000_time"
  padj_strain_time_col <- "padj_strain_time_col0"
  padj_plant_time_col <- "padj_plant_affects_time"
  padj_three_way_col <- "padj_plant_strain_specific"
  
  get_sig_uniprot <- function(df, padj_col, thresh = fdr_thresh) {
    if(!padj_col %in% colnames(df)) { warning(sprintf("Column '%s' not found — returning empty vector.", padj_col)); return(character(0)) }
    df %>% filter(!is.na(.data[[padj_col]])) %>% filter(.data[[padj_col]] <= thresh) %>% pull(uniprot) %>% na.omit() %>% unique()
  }
  
  sig_time_uniprot        <- get_sig_uniprot(lfc_file, padj_time_col, fdr_thresh)
  sig_strain_time_uniprot <- get_sig_uniprot(lfc_file, padj_strain_time_col, fdr_thresh)
  sig_plant_time_uniprot  <- get_sig_uniprot(lfc_file, padj_plant_time_col, fdr_thresh)
  sig_three_way_uniprot   <- get_sig_uniprot(lfc_file, padj_three_way_col, fdr_thresh)
  cat("Significant counts (FDR <= ", fdr_thresh, "): Time:", length(sig_time_uniprot), " Strain×Time:", length(sig_strain_time_uniprot), " Plant×Time:", length(sig_plant_time_uniprot), " Three-way:", length(sig_three_way_uniprot), "\n", sep="")
  
  # ---------------------- Run ORA for contrasts (if TERM2GENE present) ----------------------
  if(nrow(term2gene) > 0) {
    message("Running ORA for contrasts...")
    ora_time        <- run_enricher_with_desc_safe(sig_time_uniprot,        universe_uniprot, term2gene, minGSSize = minGSSize)
    ora_strain_time <- run_enricher_with_desc_safe(sig_strain_time_uniprot, universe_uniprot, term2gene, minGSSize = minGSSize)
    ora_plant_time  <- run_enricher_with_desc_safe(sig_plant_time_uniprot,  universe_uniprot, term2gene, minGSSize = minGSSize)
    ora_three_way   <- run_enricher_with_desc_safe(sig_three_way_uniprot,   universe_uniprot, term2gene, minGSSize = minGSSize)
    if(!is.null(ora_time)) write.table(ora_time$df, paste0(output_prefix, "_time.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    if(!is.null(ora_strain_time)) write.table(ora_strain_time$df, paste0(output_prefix, "_strain_time.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    if(!is.null(ora_plant_time)) write.table(ora_plant_time$df, paste0(output_prefix, "_plant_time.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    if(!is.null(ora_three_way)) write.table(ora_three_way$df, paste0(output_prefix, "_three_way.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    message("Skipping ORA (TERM2GENE not built).")
    ora_time <- ora_strain_time <- ora_plant_time <- ora_three_way <- NULL
  }
  
  # ---------------------- Dotplots (if ORA results exist) ----------------------
  plot_dot_manual <- function(ora_obj, outfile, title = "", top_n = showCategory, width = 10, height = 7, dpi = 300, wrap_width = 60, fontsize = 12) {
    if(is.null(ora_obj) || is.null(ora_obj$df)) { message("No results to plot for ", title); return(NULL) }
    df <- ora_obj$df %>%
      mutate(label_raw = ifelse(!is.na(description_wrapped) & description_wrapped != "", description_wrapped, ifelse(!is.na(description) & description != "", description, as.character(ID))),
             label = stringr::str_wrap(label_raw, width = wrap_width),
             p.adjust = as.numeric(p.adjust),
             Count = as.numeric(Count)) %>%
      arrange(p.adjust) %>%
      slice_head(n = top_n) %>%
      mutate(negLog10FDR = -log10(p.adjust + 1e-300))
    if(nrow(df) == 0) { message("No rows for ", title); return(NULL) }
    cat("Plotting top", nrow(df), "terms for:", title, "\n")
    print(df$label)
    p <- ggplot(df, aes(x = Count, y = reorder(label, negLog10FDR))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_continuous(name = "FDR", trans = "reverse") +
      scale_size_continuous(name = "Count") +
      labs(x = "Gene count", y = "", title = title) +
      theme_minimal(base_size = fontsize) +
      theme(axis.text.y = element_text(size = fontsize - 1), plot.title = element_text(face = "bold", size = fontsize + 1), legend.position = "right")
    print(p)
    ggsave(outfile, p, width = width, height = height, dpi = dpi)
    invisible(p)
  }
  
  if(!is.null(ora_time)) plot_dot_manual(ora_time, paste0("dot_", output_prefix, "_time.png"), "GO ORA — Time (descriptions)", top_n = showCategory)
  if(!is.null(ora_strain_time)) plot_dot_manual(ora_strain_time, paste0("dot_", output_prefix, "_strain_time.png"), "GO ORA — Strain×Time (descriptions)", top_n = showCategory)
  
  # ---------------------- Heatmap generation for Time & Strain×Time (top 25, readable labels) ----------------------
  logfc_cols <- grep("^logFC_", colnames(lfc_file), value = TRUE)
  if(length(logfc_cols) == 0) stop("No columns starting with 'logFC_' found in lfc_file.")
  
  build_heatmat <- function(gene_uniprot_vec, padj_col_name, label_name = "set") {
    gene_uniprot_vec <- unique(na.omit(as.character(gene_uniprot_vec)))
    if(length(gene_uniprot_vec) == 0) { message("No genes for ", label_name); return(NULL) }
    subdf <- lfc_file %>% filter(uniprot %in% gene_uniprot_vec) %>% distinct(uniprot, .keep_all = TRUE)
    if(nrow(subdf) == 0) { message("No matching rows in lfc_file for ", label_name); return(NULL) }
    if(!is.null(padj_col_name) && padj_col_name %in% colnames(subdf) && top_n_heatmap > 0) {
      subdf <- subdf %>% arrange(.data[[padj_col_name]]) %>% slice_head(n = top_n_heatmap)
    } else if(top_n_heatmap > 0) {
      subdf <- subdf %>% slice_head(n = top_n_heatmap)
    }
    mat <- as.data.frame(subdf[, logfc_cols, drop = FALSE])
    mat[] <- lapply(mat, function(x) as.numeric(as.character(x)))
    rownames(mat) <- subdf$uniprot
    
    # build readable row labels: Protein name (UniProt) if available; else UniProt
    row_labels <- rownames(mat)
    if(!is.null(go_mappings) && "Entry" %in% colnames(go_mappings)) {
      prot_col <- NULL
      if("Protein names" %in% colnames(go_mappings)) prot_col <- "Protein names"
      if(is.null(prot_col) && "Protein.names" %in% colnames(go_mappings)) prot_col <- "Protein.names"
      if(!is.null(prot_col)) {
        prot_lookup <- go_mappings %>% select(Entry, !!sym(prot_col)) %>% distinct(Entry, .keep_all = TRUE)
        prot_vec <- prot_lookup[[prot_col]]; names(prot_vec) <- prot_lookup$Entry
        lab <- prot_vec[rownames(mat)]
        lab[is.na(lab) | lab == ""] <- rownames(mat)[is.na(lab) | lab == ""]
        row_labels <- make.unique(paste0(lab, " (", rownames(mat), ")"))
      }
    }
    # truncate long labels and wrap if necessary
    row_labels <- sapply(row_labels, function(x) {
      x2 <- if(nchar(x) > max_label_length) paste0(substr(x,1,max_label_length-3), "...") else x
      stringr::str_wrap(x2, width = 60)
    }, USE.NAMES = FALSE)
    
    mat_to_plot <- as.matrix(mat)
    if(scale_rows) {
      mat_scaled <- t(scale(t(mat_to_plot)))
      mat_scaled[is.na(mat_scaled)] <- 0
      mat_to_plot <- mat_scaled
    }
    return(list(mat = mat_to_plot, labels = row_labels, subdf = subdf))
  }
  
  time_mat_obj <- build_heatmat(sig_time_uniprot, padj_time_col, label_name = "Time")
  strain_time_mat_obj <- build_heatmat(sig_strain_time_uniprot, padj_strain_time_col, label_name = "Strain×Time")
  
  # Helper: create pheatmap object (silent)
  create_pheatmap_gtable <- function(mat_obj, title="heatmap", breaks=seq(-3,3,length.out=101),
                                     colors = colorRampPalette(c("#67001F","white","#053061"))(100)) {
    if(is.null(mat_obj)) return(NULL)
    ph <- pheatmap::pheatmap(mat_obj$mat,
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             show_rownames = TRUE,
                             labels_row = mat_obj$labels,
                             fontsize_row = if(nrow(mat_obj$mat) <= 25) 10 else if(nrow(mat_obj$mat) <= 40) 8 else 6,
                             fontsize_col = 10,
                             color = colors,
                             breaks = breaks,
                             main = title,
                             border_color = NA,
                             silent = TRUE)
    return(ph)
  }
  
  # Helper: save PNG + rasterized PDF
  save_heatmap_png_to_pdf <- function(ph_object, png_file, pdf_file, width_in = 10, height_in = 8, dpi = 300) {
    if(is.null(ph_object)) { message("ph_object is NULL; skipping ", png_file); return(NULL) }
    png(png_file, width = width_in, height = height_in, units = "in", res = dpi)
    grid::grid.newpage()
    grid::grid.draw(ph_object$gtable)
    dev.off()
    img <- png::readPNG(png_file)
    pdf(pdf_file, width = width_in, height = height_in)
    grid::grid.newpage()
    grid::grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
    dev.off()
    png_size <- if(file.exists(png_file)) file.info(png_file)$size else NA
    pdf_size <- if(file.exists(pdf_file)) file.info(pdf_file)$size else NA
    message("Wrote PNG: ", png_file, " (", ifelse(is.na(png_size), "missing", paste0(round(png_size/1024), " KB")), ")")
    message("Wrote PDF: ", pdf_file, " (", ifelse(is.na(pdf_size), "missing", paste0(round(pdf_size/1024), " KB")), ")")
    invisible(list(png = png_file, pdf = pdf_file))
  }
  
  # Generate and save heatmaps (top 25)
  if(!is.null(time_mat_obj)) {
    ph_time <- create_pheatmap_gtable(time_mat_obj, title = "Time-significant genes (top 25, scaled rows)")
    if(!is.null(ph_time)) {
      save_heatmap_png_to_pdf(ph_time,
                              png_file = paste0(output_prefix, "_time_heatmap_top25.png"),
                              pdf_file = paste0(output_prefix, "_time_heatmap_top25.pdf"),
                              width_in = 10,
                              height_in = max(5, 0.28 * max(1, nrow(time_mat_obj$mat)) + 2),
                              dpi = 300)
      write.table(time_mat_obj$subdf, paste0(output_prefix, "_time_gene_table_top25.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  } else message("No Time-significant genes to plot.")
  
  if(!is.null(strain_time_mat_obj)) {
    ph_strain_time <- create_pheatmap_gtable(strain_time_mat_obj, title = "Strain×Time-significant genes (top 25, scaled rows)")
    if(!is.null(ph_strain_time)) {
      save_heatmap_png_to_pdf(ph_strain_time,
                              png_file = paste0(output_prefix, "_strain_time_heatmap_top25.png"),
                              pdf_file = paste0(output_prefix, "_strain_time_heatmap_top25.pdf"),
                              width_in = 10,
                              height_in = max(5, 0.28 * max(1, nrow(strain_time_mat_obj$mat)) + 2),
                              dpi = 300)
      write.table(strain_time_mat_obj$subdf, paste0(output_prefix, "_strain_time_gene_table_top25.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  } else message("No Strain×Time-significant genes to plot.")
  
  message("Done. Check PNGs, rasterized PDFs and TSVs in the working directory.")
  

