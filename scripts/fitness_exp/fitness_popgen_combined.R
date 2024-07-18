# This script takes the gene conservation information from the panX output and combines it with the in planta growth information.
library(tidyr)
library(dplyr)
library(gpplot2)
library(lme4)



#this is the output from running the python evolutionary genetics scripts on the panX output
conservation = read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd.csv", row.names = 1)

#this is data file for the output from the tnseq experiment. 
eyach_tailocin <- read.table("~/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/tailocin/tailocin_plant_fitness_selection_717_2024.csv", header = TRUE)

#I need to rename the gene names in conservation file
orthologs <- read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/p25c2_dc3000_ortholog_7_2_2024/p25c2_to_dc3000_noReps.csv", 
                        header = TRUE,
                        sep = ",", row.names = 2)
orthologs[orthologs$DC3000=="",]$DC3000<-NA



# change conservation DAKCFMEO to BJE identifier
conservation$BJE <- orthologs[rownames(conservation),]$p25_BJE
conservation$WP <- orthologs[rownames(conservation),]$DC3000
# remove those genes that are duplicated in the BJE genome
# get numbers of duplicated rows. Keeping original number
conservation2<-conservation %>%
  arrange(BJE) %>%
  filter(duplicated(BJE) == FALSE)

conservation2$DAC <- rownames(conservation2)
# remove genes with 
conservation2 <- conservation2 %>% filter(is.na(BJE)==FALSE)
rownames(conservation2)<- conservation2$BJE

#Now I need to rename in terms of WP (the DC3000 mapping)
i=1
dupl_WP<-which(duplicated(conservation2$WP))
conservation2 <- conservation2[-dupl_WP,]
#rownames(conservation2)<- conservation2$WP
for(rowname in rownames(conservation2)){
   if(!is.na(conservation2[rowname,]$WP)){
    rownames(conservation2)[i]=conservation2[rowname,]$WP
   }
  i=i+1
}

# Now combine the diversity and tnseq files
fit_con <- merge(conservation2, eyach_tailocin, by=0, all=TRUE)

# now we can ask how the time in tree relates to the fitness effect
model1 <- lm(data=fit_con, log2FoldChange_Ey ~ Time.in.tree)






