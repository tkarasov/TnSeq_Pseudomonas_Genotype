# The goal of this script is to take lists of gene names and do gene ontology enrichments.
#Conclusion as of 7/2025 is that there are enrichments of real basal processes, but not any story that makes a lot of sense.
# first step is to write a function that takes two lists and gets enrichments
setwd("~/Google\ Drive/My\ Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")

#sig_genes is the output from the limmavoom file for which genes were significant in which interaction.
sig_genes = read.csv("gene_sig_summary.csv", header = TRUE)


# ortholog table for possible genes
ortholog_tab <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", header = TRUE, row.names = 1, sep = ",")

# Load biomart
# Load biomaRt
library(biomaRt)
library(ggplot2)
library(GO.db)
library(tidyr)
library(dplyr)
# I downloaded the uniprot mapping for the DC3000 genome
uniprot <- read.table(
  "uniprotkb_proteome_UP000002515_2025_07_01.tsv",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)
df <- uniprot[,c("Entry", "RefSeq", "Gene.Names")]

# Split RefSeq on ";"
refseq_split <- strsplit(df$RefSeq, ";")

# Extract first and second RefSeq (NP_ and WP_)
df$NP_RefSeq <- sapply(refseq_split, function(x) {
  np <- grep("^NP_", x, value = TRUE)
  if (length(np) > 0) np[1] else NA
})

df$WP_RefSeq <- sapply(refseq_split, function(x) {
  wp <- grep("^WP_", x, value = TRUE)
  if (length(wp) > 0) wp[1] else NA
})

# Extract PSPTO from Gene Names
df$PSPTO <- sub(".*(PSPTO_[0-9]+).*", "\\1", df$Gene.Names)
# remove the "_" in the PSPTO
# Remove underscore
df$PSPTO <- gsub("_", "", df$PSPTO)

# Select columns
df_tidy <- df[, c("Entry", "NP_RefSeq", "WP_RefSeq", "PSPTO")]


# Build the named vector (WP ID = name, Q5 = value)
pspto_to_q5 <- setNames(df$Entry, df$PSPTO)

# Remove any NA names if no WP found (just in case)
#pspto_to_q5 <- pspto_to_q5[!is.na(names(pspto_to_q5))]
# this narrowed the protein number from 3179 to 2335. Not clear why a bunch are missing WPs. TBD

#trying again
tk_remapping <- read.table("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/DC3000_genome_mappings/tk_remapping.txt",
                           sep = "\t", 
                           header = TRUE, 
                           fill = TRUE,   # allows blank fields by filling missing columns
                           quote = "",    # avoids issues with unexpected quotes
                           comment.char = "")  # disables comment parsing if you want to keep lines starting with #)

tk_wp_pspto <- tk_remapping[, c("non.redundant_refseq","single_old_locus")]
colnames(tk_wp_pspto) <- c("WP_RefSeq", "PSPTO")
refseq_dict <- setNames(tk_wp_pspto$PSPTO, tk_wp_pspto$WP_RefSeq)

tk_wp_pspto$uniprot <- pspto_to_q5[tk_wp_pspto$PSPTO]
print(tk_wp_pspto$uniprot)
write.table(tk_wp_pspto,"/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/DC3000_genome_mappings/tk_remapping_with_uniprot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE )


# ok now I can take my gene lists and try to do something interesting I hope. Hallelujah that worked!
sig_genes$PSPTO <-  refseq_dict[sig_genes$Gene]
sig_genes$uniprot <- pspto_to_q5[sig_genes$PSPTO]

# sig_genes are those that were significant and the orthology ones were the universe
universe <- data.frame(Gene=names(table(ortholog_tab$DC3000)))
universe$PSPTO <- refseq_dict[universe$Gene]
universe <- universe[!(is.na(universe$PSPTO) | universe$PSPTO== ""), ]
universe$uniprot <- pspto_to_q5[universe$PSPTO]
# 3852 genes in the universe

# can we get a list of all genes under consideration
go_mappings = read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/DC3000_genome_mappings/gene_ontology_mapping_uniprotkb_proteome_UP000002515_2025_07_01 (2).tsv", header = TRUE, sep="\t")
# Example with your data
go_dict <- lapply(go_mappings$Gene.Ontology.IDs, function(x) unlist(strsplit(x, "; ")))
names(go_dict) <- go_mappings$Entry


# Convert the list to a long dataframe
go_df <- stack(go_dict)
colnames(go_df) <- c("GO", "uniprot")
universe_q8 <- universe$uniprot


#an example of significant genes
sig_q8 <- universe %>% filter(Gene %in% sig_genes$Gene) %>% pull(uniprot)



# 6️⃣ (OPTIONAL) Add GO descriptions if you have them
# Example: go_descriptions <- data.frame(GO=c("GO:0001234"), description=c("Example function"))
# go_sig <- go_sig %>% left_join(go_descriptions, by="GO")


go_sig$description <- Term(GOTERM[go_sig$GO])
# Example list of GO IDs
go_ids <- go_sig$GO  # or whatever your GO ID column is called

# Get GO term names
go_terms <- Term(GOTERM[go_ids])

# Get GO ontology category (BP, MF, CC)
go_ont <- Ontology(GOTERM[go_ids])

# Combine into a dataframe
go_sig <- go_sig %>%
  mutate(
    description = go_terms,
    ontology = go_ont
  )



compute_go_enrichment <- function(sig_q8, universe_q8, GO_dict) {
  # this function takes a list of target genes and a universe of genes and detects go_enrichments
  go_df <- stack(GO_dict)
  colnames(go_df) <- c("GO", "uniprot")
  
  go_df <- go_df[go_df$uniprot %in% universe_q8, ]
  
  go_counts <- go_df %>%
    group_by(GO) %>%
    summarize(
      n_universe = n(),
      n_sig = sum(uniprot %in% sig_q8),
      .groups = "drop"
    ) %>%
    mutate(
      n_bg = length(universe_q8) - n_universe,
      pval = phyper(n_sig - 1, n_universe, n_bg, length(sig_q8), lower.tail = FALSE),
      pval_adj = p.adjust(pval, method = "fdr"),
      description = sapply(GO, function(go) {
        term <- GOTERM[[go]]
        if (is.null(term)) NA_character_ else Term(term)
      }),
      ontology = sapply(GO, function(go) {
        term <- GOTERM[[go]]
        if (is.null(term)) NA_character_ else Ontology(term)
      })
    )
  
  return(go_counts)
}


plot_go_enrichment <- function(go_results, fdr_thresh = 0.05, top_n = 20) {
  go_sig <- go_results %>%
    filter(pval_adj < fdr_thresh) %>%
    arrange(pval_adj) %>%
    slice_head(n = top_n)
  
  ggplot(go_sig, aes(x = reorder(description, -pval_adj), y = -log10(pval_adj), fill = ontology)) +
    geom_col(width = 0.7, color = "black") +
    coord_flip() +
    scale_fill_manual(
      values = c(
        BP = "#4D648D",   # soft blue
        MF = "#E27D60",   # soft coral
        CC = "#85DCBA"    # soft green
      ),
      na.value = "grey70"
    ) +
    labs(
      title = "GO enrichment of significant genes",
      x = "GO term",
      y = expression(-log[10]~"(FDR-adjusted p-value)"),
      fill = "Ontology"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# Run enrichment
go_results <- compute_go_enrichment(sig_q8, universe_q8, go_dict)

# Plot
all_plot<-plot_go_enrichment(go_results)

# let's do subplots of genes that are important strain-specifically
# sig_time_strain
time_strain<- sig_genes %>% 
  filter(Strain_Time==TRUE) %>%
  pull(uniprot)

Plant_Time<- sig_genes %>% 
  filter(Plant_Time==TRUE) %>%
  pull(uniprot)

time<- sig_genes %>% 
  filter(Time==TRUE) %>%
  pull(uniprot)


go_time_results <- compute_go_enrichment(time, universe_q8, go_dict)
go_time_strain_results <- compute_go_enrichment(time_strain, universe_q8, go_dict)
go_Plant_results <- compute_go_enrichment(Plant_Time, universe_q8, go_dict)
time <- plot_go_enrichment(go_time_results)
time_strain <- plot_go_enrichment(go_time_strain_results)
# plot_go_enrichment(go_Plant_results) None are significant



# OK time_strain has a number of significcant results, seemingly many in central metabolism.
# Combine without overlap
combined_plot <- plot_grid(
  time, time_strain,
  ncol = 2,
  align = "h",
  label_fontfamily = "Arial",
  rel_widths = c(1, 1)
)
font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")
ggsave("go_enrichment_comparison.tiff", combined_plot, width = 12, height = 6, dpi = 300, compression = "lzw")

#############
# There is some sort of weird peak in the histogram of fitness values for DC3000 around +1.0. Which genes are these?
lfc_df <- read.csv("logFCDC3000logFCp25c2_7_2025.csv", header = TRUE)
lfc_dc_1 <- lfc_df %>% filter(logFC_DC3000>.75) %>% filter(logFC_DC3000<1.75) 

lfc1_Q <- universe %>% filter(Gene %in% lfc_dc_1$Gene) %>% pull(uniprot)
lfc_universe <- universe %>% pull(uniprot)
lfc1_GO_enrich <- compute_go_enrichment(lfc1_Q, lfc_universe, go_dict)
plot_go_enrichment(lfc1_GO_enrich)
# only thing that is significant is cytoplasm. ??

#I would like to make a heatmap for the genes associated with virulence
# Search GO.db for terms containing "virulence" or "pathogenesis"
all_go <- Term(GOTERM)
virulence_go <- grep("virulence|pathogenesis|hrp|tox|effector|coronatine|avr|secretion|type III secretion system", all_go, value = TRUE, ignore.case = TRUE)

# Extract GO IDs
virulence_go_ids <- names(virulence_go)
# Add some known T3SS GO terms
t3ss_go_ids <- c("GO:0030257", "GO:0009289", "GO:0005576")

# Combine with your original virulence list
combined_go_ids <- union(virulence_go_ids, t3ss_go_ids)
# Subset genes that are classified with one of these go ids
# go_mappings Entry column is the q value and the Gene.Ontology.IDs is the GO ids separated by semicolons
library(stringr)
filtered_go_mappings <- go_mappings %>%
  filter(str_detect(Gene.Ontology.IDs, str_c(combined_go_ids, collapse = "|")))

# now I want to do a heatmap for those of those genes that are in the lfc dataframe
vir_genes <- universe %>% filter(uniprot %in% filtered_go_mappings$Entry) 
vir_lfc <- lfc_df %>% filter(Gene %in% vir_genes$Gene)
vir_lfc$uniprot <- vir_genes[vir_genes$Gene %in% vir_lfc$Gene, c("uniprot")]
#vir_lfc$Protein.names <- filtered_go_mappings %>% filter(uniprot %in% vir_lfc$uniprot ) %>% pull(Protein.names)
vir_lfc <- left_join(
  vir_lfc,
  filtered_go_mappings[, c("Entry", "Protein.names")],
  by = c("uniprot" = "Entry")
)



library(pheatmap)
library(grid)

# Prepare matrix
mat <- as.matrix(vir_lfc[, c("logFC_DC3000", "logFC_P25C2")])
colnames(mat) <- c("DC3000 in Col-0", "P25.c2 in Col-0")
rownames(mat) <- make.unique(vir_lfc$Protein.names)

# Color scale
breaks <- seq(-2, 2, length.out = 100)
colors <- colorRampPalette(c("#67001F", "white", "#053061"))(99)

# Draw the heatmap and capture the gtable object
pdf("virulence_gene_log2FC_heatmap.pdf", width = 7, height = 6)  # ~89 mm × variable height
heatmap_obj <- pheatmap(mat,
                        cluster_rows = TRUE,
                        cluster_cols = FALSE,
                        color = colors,
                        breaks = breaks,
                        fontsize = 8,
                        fontsize_row = 6,
                        fontsize_col = 8,
                        border_color = NA,
                        main = "")

# Add custom legend title (log₂FC)
grid.text("Virulence Gene log₂FC (T3–T0)", x = 0.5, y = 0.97, gp = gpar(fontsize = 10))
grid.text(expression(log[2] * "FC"), x = 0.92, y = 0.75, gp = gpar(fontsize = 10))
dev.off()
