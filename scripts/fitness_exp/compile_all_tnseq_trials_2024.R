# In December 2022 and May 2023 Effie completed two experiments with DC3000 and p25.c2 Rb-Tnseq grown in two ecotypes. This file is building a list for each experiment which is then written to an rds
library(gsubfn)
library(dplyr)
setwd("~/Google\ Drive/My\ Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")
orthologs <- read.table("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/orthologues_dc3000_p25_c2_7_2022.txt", 
                        header = TRUE,
                        sep = "\t")

##################################################
# In this first section I am compiling Effie's experiments together. This corresponds to exp_0001 and exp_0002
##################################################

compile_experiment <- function(path_files, gff, sample_info){
  sample_order_p = read.table(paste(path_files, sample_info, sep=""),header = TRUE)
  rownames(sample_order_p) <- sample_order_p$sample
  counts = read.delim(paste(path_files, gff, sep=''), header = T)
  count_list <- names(which(colSums(counts[,c(10:dim(counts)[2])])>1))
  count_table = counts[,count_list]
  #reorder sample_order
  sample_order_p = sample_order_p[count_list,]
  sample_order_p<- na.omit(sample_order_p)
  rownames(count_table) = counts$gene_id
  count_table_p <- t(count_table[, colnames(count_table) %in% sample_order_p$sample])
  #count_table_p is the genes in rows and samples in columns
  #sample_order_p is the samples in rows and columns are all of the variables
  experiment0 <- list(  sample_order_p, count_table_p)
  return(experiment0)
}




assign_gene_name <- function(gene_record){
  #The goal of this function is to take a geneID from a Tnseq result file and output its ortholog name in the DC3000 genome and in the p25.c2 genome
  Gene_name1 <- strsplit(gene_record, "Name=")[2]
  Gene_name2 <- strsplit(Gene_name1, ";")[1]
}


#Paths to all experiments
dc_22_path <- ("./in_planta_dec22/dc3000_dec22/")
p25_22 <- ("./in_planta_dec22/p25c2_dec22/")
dc_23 <- ("./in_planta_may23/dc3000_may23/")
p25_23 <- ("./in_planta_may23/p25c2_may23/")

#Pull meta data and make count table
#DC3000 2022
path_files <- dc_22_path
gff <- "in_planta_dec22_dc3000_header.gff"
sample_info <-"sample_info_dc3000_0323.txt"
dc_22_counts <- compile_experiment(path_files, gff, sample_info)
genelist_dc3000 <- colnames(dc_22_counts[2][[1]]) 
genelist_dc30002 <- genelist_dc3000[grep("^ID=cds", genelist_dc3000)]
genelist_dc30003 <- gsub('.*;(.*);locus_tag=','',genelist_dc30002)
genelist_dc22 <-substr(genelist_dc30003 , 1, regexpr("\\;", genelist_dc30003 )-1)

  
saveRDS(dc_22_counts, "../full_experiments/dc22.rds")

#p25c2 2022
path_files <- p25_22
gff <- "in_planta_dec22_p25c2_counts_header.gff"
sample_info <-"sample_info_p25c2_0323.txt"
p25_22_counts <- compile_experiment(path_files, gff, sample_info)
saveRDS(p25_22_counts, "../full_experiments/p25c2_22.rds")

#DC3000 2023
path_files <- dc_23
gff <- "count_table_dc3000_in_planta_may23_header.gff"
sample_info <-"dc3000_sample_info_may23.txt"
dc_23_counts <- compile_experiment(path_files, gff, sample_info)
saveRDS(dc_23_counts, "../full_experiments/dc23.rds")

#p25c2 2023
path_files <- p25_23
gff <- "in_planta_may23_p25c2_header.gff"
sample_info <-"in_planta_may23_info_p25c2.txt"
p25_23_counts <- compile_experiment(path_files, gff, sample_info)
saveRDS(p25_23_counts, "../full_experiments/p25c2_23.rds")

# Create a file to merge on gene name
ortho <- read.table("../orthologues_dc3000_p25_c2_7_2022.txt", sep = "\t", header = TRUE)
# the first column is the WP ID in DC3000, the second column is the name in p25_c2
#add the ortholog data to eash rds

gene_dc22 <- colnames(dc_23_counts[[2]])

pull_gene_name <- function(counts_table, start_gene){
  # relabel dc3000 counts table 
  # start_gene is the variable that corresponds to the column where the gene listing starts
  gene_list <- colnames(counts_table[[2]])
  poo <- strsplit(gene_list, "ID=")
  poo2 <- sapply( poo, "[", 2 )
  poo3 <- strsplit(poo2, ";")
  poo4 <- sapply(poo3, "[[", 1)
  print(length(poo4))
  df_ct <- data.frame(counts_table)
  print(dim(df_ct))
  colnames(df_ct)[start_gene:dim(df_ct)[2]] <- poo4
  return(df_ct)
}

# counts_table = p25_22_counts
# ortho_dict <- ortho[,3]
# names(ortho_dict) <- ortho[,1]
# reverse_dict <- names(ortho_dict)
# names(reverse_dict) <- ortho_dict

pull_gene_name_ortho <- function(counts_table, ortho, start_gene){
  # relabel p25_c2 counts table 
  gene_list <- colnames(counts_table[[2]])
  poo <- strsplit(gene_list, "ID=") #split the gene annotation string
  poo2 <- sapply( poo, "[", 2 ) #take only the second element
  poo3 <- strsplit(poo2, ";")
  poo4 <- sapply(poo3, "[[", 1) # take only the first element (the gene or cds name)
  df_ct <- data.frame(counts_table)
  print(length(poo4))
  colnames(df_ct)[start_gene:dim(df_ct)[2]] <- poo4
  
  #now we need to rename the orthologs of p25.c2 to DC3000. Iterate through the colnames of df_ct, and ask if the value is in reverse_dict 
  new_poo4 <- poo4
  for(i in 1:length(poo4)){
    new_poo4[i] <- reverse_dict[poo4[i]]
    if(is.na(new_poo4[i])) new_poo4[i] = poo4[i]
    
    }
  colnames(df_ct)[start_gene:dim(df_ct)[2]] <- new_poo4
  return(df_ct)
}

dc22_relabel <- pull_gene_name(dc_22_counts, start_gene = 10)
dc23_relabel <- pull_gene_name(dc_23_counts, start_gene = 9)
p25_22_relabel <- pull_gene_name_ortho(p25_22_counts, ortho, start_gene = 10)
p25_23_relabel <- pull_gene_name_ortho(p25_23_counts, ortho, start_gene = 8)

# for dc22_relabel and dc23_relable we need to remove columns that dont have cds in them
remove_gene<-function(relabel_file, start_gene){
  keep_gene <- relabel_file[,stri_detect_fixed(colnames(relabel_file), "cds")]
  new_name = sapply( strsplit(colnames(keep_gene),"cds-"), "[", 2 )
  keep_file <- keep_gene
  colnames(keep_file) <- new_name
  keep_file<-cbind(relabel_file[,1:(start_gene -1)], keep_file)
  return(keep_file)
}

dc22_relabel <- remove_gene(dc22_relabel, start_gene = 11)
dc23_relabel <- remove_gene(dc23_relabel, start_gene = 9)
p22_limit <- p25_22_relabel[,c(10:dim(p25_22_relabel)[2])]
p23_limit <- p25_23_relabel[,c(8:dim(p25_23_relabel)[2])]
p25_22_new<- t(rowsum(t(p22_limit), group = colnames(p22_limit), na.rm = T))
p25_23_new<- t(rowsum(t(p23_limit), group = colnames(p23_limit), na.rm = T))



p25_22_new2<- cbind(p25_22_relabel[,c(1:9)],p25_22_new )
p25_23_new2<- cbind(p25_23_relabel[,c(1:7)],p25_23_new )

# experiment column as first column to the four dataframes
p25_22_new2 <- cbind("exp_0001", p25_22_new2)
p25_23_new2 <- cbind("exp_0002", p25_23_new2)
dc22_relabel <- cbind("exp_0001", dc22_relabel)
dc23_relabel <- cbind("exp_0002", dc23_relabel)
colnames(p25_22_new2)[1]<-"experiment"
colnames(p25_23_new2)[1]<-"experiment"
colnames(dc22_relabel)[1]<-"experiment"
colnames(dc23_relabel)[1]<-"experiment"

#p25_22_new2<- sapply(split.default(p22_limit, names(p22_limit)), rowSums, na.rm = TRUE)
#end of day on 6/4/2024 there is a mismatch in the number of columns between the 22 and 23 dates. This mismatch occurs because the 22 experiment had two extra columns, "cycles" and "reads". We need to do merge with all.

dc_exp <- full_join(x=dc22_relabel, y=dc23_relabel)

#78 of the p25 columns are not unique (it seems like two of the p25 genes map to the same DC3000 ortholog. What can we do? how about merge colulmns with the same column name)
p25_exp <- full_join(x=p25_22_new2, y=p25_23_new2)

# all experiment join
all_exp <- full_join(x=dc_exp, y=p25_exp)
saveRDS(all_exp, "../full_experiments/all_p25_dc_axenic_6_2024.rds")
write.csv(all_exp, "../full_experiments/all_p25_dc_axenic_6_2024.csv", quote = FALSE, row.names =FALSE)

