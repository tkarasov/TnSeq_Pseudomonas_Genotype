# This script compares the output from Effie's blast best hits to the orthofinder analysis done with steps1 and 2
#!/usr/bin/env Rscript
# step3_compare_blast_orthofinder_safe.R
# Safe, verbose comparison: for BLAST rows where DC3000 is present
# check whether (p25_DAK, DC3000) occur in the same row of Orthogroups.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(rlang)
})

# ---------- Edit paths here ----------
#blast_fp <- "/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv"
#new diamond blast mapping from 12/11/2025
blast_fp="/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/ortholog_diamond_12_2025/rbh_20251211_165439/rbh.csv"

orthogroups_tsv <- "/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/orthology_host_run/orthofinder_pair_20251210_165321/orthofinder_run/Results_Dec10/Orthogroups/one2one_pairs/one2one_DC3000_protein_12_2025.faa_vs_plate25.C2.annotation.faa.tsv"
out_csv <- "/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/orthology_host_run/blast_pairs_vs_orthogroups_matches_safe.csv"
# -------------------------------------

blast=read.csv(blast_fp, stringsAsFactors = FALSE)
colnames(blast)[1:2]<-c("DC3000", "p25_DAK")
orthogroups=read.delim(orthogroups_tsv, stringsAsFactors = FALSE)

# Create a lookup table in blast where p25_DAK column is the lookup key and DC3000 is the value
blast_lookup <- blast %>%
  select(p25_DAK, DC3000) %>%
  filter(!is.na(DC3000) & DC3000 != "") %>%
  distinct()

# Now iterate through the keys in the lookup table and look to see if there are rows in orthogroups in which that string is found. If so, see if the value in the lookup table is found in that same row of the orthogroups file
results <- data.frame(p25_DAK=character(), DC3000=character(), Orthogroup_Match=logical(), stringsAsFactors=FALSE)
for (i in 1:nrow(blast_lookup)) {
  p25_dak_key <- blast_lookup$p25_DAK[i]
  dc3000_value <- blast_lookup$DC3000[i]
  
  # Find rows in orthogroups where p25_dak_key is found
  matching_rows <- orthogroups %>%
    filter(if_any(everything(), ~ str_detect(., fixed(p25_dak_key))))
  
  # Check if dc3000_value is found in any of those matching rows
  orthogroup_match <- any(apply(matching_rows, 1, function(row) {
    any(str_detect(row, fixed(dc3000_value)))
  }))
  
  # Append result
  results <- rbind(results, data.frame(p25_DAK=p25_dak_key, DC3000=dc3000_value, Orthogroup_Match=orthogroup_match, stringsAsFactors=FALSE))
}
# Write results to CSV in the same directory as the orthogroup file in a file name "blast_pairs_vs_orthogroups_matches_safe.csv"
out_csv="/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/orthology_host_run/blast_pairs_vs_orthogroups_matches_safe.csv"
write.csv(results, out_csv, row.names = FALSE)

# 4011 positive matches and 24 negative matches
ortholog_subset <- results[!is.na(results$DC3000) &
                                  trimws(results$DC3000) != "", ]

# find columns in which have duplicated DC3000 values
duplicated_dc3000 <- ortholog_subset$DC3000[duplicated(ortholog_subset$DC3000) | duplicated(ortholog_subset$DC3000, fromLast = TRUE)]
duplicated_rows <- ortholog_subset %>%
  filter(DC3000 %in% duplicated_dc3000) %>%
  arrange(DC3000)

# OK there are 85 rows of duplicated DC3000 values. No longer. After I used the diamond pipeline on 12/12/2025 I no longer had the issue with duplicates
