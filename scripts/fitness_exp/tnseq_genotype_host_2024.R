library(DESeq2)
library(tidyverse)
library(pheatmap)
library(readr)
# This script reads in the rds experiments for the four experiments Done by Effie in 12/2022 and 3/2023 with Rb-Tnseq and does basic deseq comparison trials

setwd("~/Google\ Drive/My\ Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000")

# #read in the rds files generated by compile_all_tnseq_trials_2024.R
# dc22 <- readRDS( "../full_experiments/dc22.rds")
# dc23 <- readRDS( "../full_experiments/dc23.rds")
# p25_22 <- readRDS("../full_experiments/p25c2_22.rds")
# p25_23 <- readRDS("../full_experiments/p25c2_23.rds")


#read in the metacsv with all the experiments. I'm confused. I think we have a better one. I don't know where this came from. 
#all_exp <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data//full_experiments/all_p25_dc_axenic_6_2024.csv", header = TRUE)
all_exp <- read.csv("/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/compiled_trials_3_2024/data/in_planta_rbtnseq_p25c2_dc3000/", header = TRUE)

# I think it's best if we just compare internally within the experiments and then look at how replicable the results are

# 2023 and 2023 experiments. Let's analyze the orthologs shared by all.
###########
#there is inconsistency in naming. everything needs to be made lowercase

all_orth <- all_exp[, colSums(is.na(all_exp[1:7321])) == 0]
count_table_p <- t(all_orth[,9:dim(all_orth)[2]])
sample_order_p <- (all_orth[,1:8])
sample_order_p <- sample_order_p %>% mutate(plant = tolower(plant))

# Let's filter out all samples that have fewer than 25,000 reads
high_enough <- which(colSums(count_table_p)>=25000)
count_table_filter <-count_table_p[,high_enough]
sample_order_filter <- sample_order_p[high_enough,]

# Let's write the count_table and sample_order to file
write_rds(sample_order_filter, "~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/full_experiments/all_sample_order_filter_6_12_2024.rds")
write_rds(count_table_filter, "~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/full_experiments/all_count_table_filter_6_12_2024.rds")


#########
# This blog was helpful for getting me to understand the interactions: https://www.biostars.org/p/353618/#356111


# first I will randomly assign the controls to either eyach or col0
plants <- c(rep("ey15_2",24), rep("col_0",24))
samps <- sample( plants, size=48, replace = F)
sample_order_filter[sample_order_filter$plant == "ctrl",]$plant <- samps

#now with the filtered dataset, make deseq object
dds <- DESeqDataSetFromMatrix(countData = count_table_filter, 
                              colData = sample_order_filter, design = ~ plant + 
                                experiment + treatment + plant:treatment +plant:time_point+treatment:time_point)
#dds$plant <- factor(dds$plant, levels = c("ctrl","col_0", "ey15_2"))
dds <- estimateSizeFactors(dds)


#this experiment does the LRT test to compare whether the time interaction with bacterial treatment have different effects. *** Note again that this is comparing the two models and asking which genes does it matter whether there is a different interaction term depending on which strain is used in the comparison.
dds_red <- DESeq(dds, test="LRT", reduced = ~ plant + 
                   experiment + treatment + plant:treatment +plant:time_point)
res_red <- results(dds_red)

#Let's do some basic visualization. a variance stabilizing transformation is useful for when doing clustering methods
ntd <- normTransform(dds_red)
select <- order(rowMeans(counts(dds_red,normalized=TRUE)),
                decreasing=TRUE)[1:134]
df <- as.data.frame(colData(dds_red)[,c( "treatment")])
pdf("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/figs/heatmap_strain_background_effect.pdf")
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()


# Let's watch the counts change
fiss <- plotCounts(dds_red, which.min(res_red$padj), 
             intgroup = c("time_point", "treatment"), returnData = TRUE)
#most significant gene is WP_005771391.1 which is a glycosyl transferase family 2
# WP_005771391.1	DAKCFMEO_03457	BJEIHDPM_00932	PSPTO_RS05640
#BJEIHDPM_00932;product=putative glycosyltransferase
#in 2022 experiment is absent across all samples. In 2023, has more appreciable numbers but goes to zero at the t3 timepoint.  

pdf("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/figs/WP_005771391.1_genetic_background.pdf")
ggplot(fiss, aes(x=time_point, y=count,group=treatment, color = treatment)) + 
  geom_point() + stat_summary(fun.y=median, geom="line") + scale_y_log10()
dev.off()

# this vignette (#9 https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments) was useful for figuring out which genes are important for the genetic background.

#perform the median of ratios method of normalization
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="treatmentp25c2.time_pointt3")

# or to shrink log fold changes association with condition:
res_bac <- lfcShrink(dds, coef="treatmentp25c2.time_pointt3", type="apeglm")
plotMA(res_bac, ylim = c(-5,5))
table(res$padj<0.01)



############# From here we are going to narrow our dds to the p25.c2 samples and look at the dynamics of p25.c2 because this is the one for which we have the best data.
keep<-which(sample_order_filter$treatment!="dc3000")
count_table_p25<-count_table_filter[,keep]
sample_order_p25 <- sample_order_filter[keep,]
dds_p25 <- DESeqDataSetFromMatrix(countData = count_table_p25, 
                              colData = sample_order_p25, design = ~ plant + 
                                experiment  + plant:time_point)

#perform the median of ratios method of normalization
dds_p25 <- estimateSizeFactors(dds_p25)
dds_p25 <- DESeq(dds_p25)
resultsNames(dds_p25) # lists the coefficients


### What percentage of genes are differentially important depending on the host genotype
# This is a contrast between Ey1.5 and Col-0 asking for genes that changed in importance based on host genetic background.
contrast <- c("plantcol_0.time_pointt3", "plantey15_2.time_pointt3")
res_plant <- results(dds_p25, contrast=list(contrast))
table(res_plant$padj<0.01)

# FALSE  TRUE 
# 3317   243
#6.8% of genes show a genotypic-specific effect

# ### What percentage of genes are important dependent on the strain background
# res_plant <- results(dds, name = "treatment_p25c2_vs_dc3000"). This analysis was completely wrong based off of not doing the correct change from T0 to T3. Differences in starting concentrations in the two libraries undoubtably had a strong effect. 
# table(res_plant$padj<0.01)
# 
# # FALSE  TRUE 
# # 1685  2030 
# # 54.6% of the genome

##################
# Now let's do a comparison with orthologs that are unique to DC3000
uniq_dc <- ortho[sapply(ortho[,5], function(x) all( x == '' )),1]
shared_dc <- ortho[sapply(ortho[,5], function(x) all( x != '' )),1]
shared_dc <- shared_dc[which(shared_dc %in% rownames(res_col0))]
uniq_dc <- uniq_dc[which(uniq_dc %in% rownames(res_col0))]
#For some reason 75 of the DC3000 genes are not in the res_col0 at all. 

# Let's compare the genes in their importance in colonizing Col-0
res_col0<- results(dds, name = "plant_col_0_vs_ctrl")
# 417 genes important 
shared<-res_col0[shared_dc,]
shared$shared <- c(rep("Shared", dim(shared)[1]))
unique<-res_col0[uniq_dc,]
unique$shared <- "Unique"
All <- rbind(shared,unique)

#Now let's plot shared vs unique to see if they differ in their distributions of fold changes.
# The unique genes are definitely less likely to be deleterious when knocked out. 
hist_shared_unique <- ggplot(All, aes(log2FoldChange, fill = shared)) + geom_density(alpha = 0.4) +
  theme_bw() + scale_fill_viridis_d()


# What the hell is going on with that second peak in the unique?
unique_pos <- data.frame(All) %>% filter(shared=="Unique") 
unique_pos <- unique_pos %>% filter(log2FoldChange>0.75)
# 0.775443563482869  0.93896798461601  1.02126669109762  1.20684935637317  2.13313372879183 
#2                 2                 2               198                 2 
table(ortho[which(ortho[,1]=="WP_003409202.1"),][,2])

# let's create an expanded FC file
resultsNames(dds)
# [1] "Intercept"                       "plant_col_0_vs_ctrl"            
# [3] "plant_ey15_2_vs_ctrl"            "experiment_exp_0002_vs_exp_0001"
# [5] "treatment_p25c2_vs_dc3000"       "plantcol_0.treatmentp25c2"      
# [7] "plantey15_2.treatmentp25c2"   

# the contrasts I want for my paper are the FC for Eyach vs ctrl, Col vs ctrl, p25c3 vs dc3000
res_col0<- results(dds, name = "plant_col_0_vs_ctrl")
res_col0$contrast <- "plant_col_0_vs_ctrl"

res_eyach<- results(dds, name = "plant_ey15_2_vs_ctrl")
res_eyach$contrast <- "plant_ey15_2_vs_ctrl"

res_dc_p25<- results(dds, name = "treatment_p25c2_vs_dc3000")
res_dc_p25$contrast<- "treatment_p25c2_vs_dc3000"

# add p0 to each of these for the control
# need to add a column that is the initial frequency of the gene
p0_eyach_t0_p25c2<- dds[,which(dds$plant=="ctrl" & dds$treatment=="p25c2")]

res_col0$p0_axenic <- rowMeans(fpm(p0_eyach_t0_p25c2, robust = FALSE))/1000000
res_eyach$p0_axenic <- rowMeans(fpm(p0_eyach_t0_p25c2, robust = FALSE))/1000000
res_dc_p25$p0_axenic<- rowMeans(fpm(p0_eyach_t0_p25c2, robust = FALSE))/1000000

comb_res <- rbind(res_col0, res_eyach)
comb_res <- rbind(comb_res, res_dc_p25)

#write table to file
write.table(comb_res, "~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/tailocin/all_plant_contrasts_deseq_7_3_2024.csv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep="\t")
