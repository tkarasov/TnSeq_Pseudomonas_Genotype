# This script takes the gene conservation information from the panX output and combines it with the in planta growth information.
library(tidyr)
library(dplyr)
library(gpplot2)
library(lme4)
library(randomForest)
library(caret)


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


################
# Let's read in the Tsuda expression data
################
exp <- read.csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/Tsuda_PNAS_expression_data/tk_modified_files/expression_gene_mapping.txt", sep="\t", header = TRUE)

# subset_expression to a few of interest. WP, RS, ID, name, annotation, Col.0_Pto
exp_sub <- exp %>% select(c("WP", "RS", "ID", "name", "annotation", "Col.0_Pto", "KB_Pto",
                            "MM_Pto", "SA_Pto"))
var_exp <-exp %>% select(c("Col.0_Pto", "KB_Pto",
                           "MM_Pto", "SA_Pto"))%>% as.matrix() %>% rowSds()
  #exp %>% select(c(4:35))%>% as.matrix() %>% rowSds()
exp_sub$variance <- var_exp
ratio<- exp_sub$Col.0_Pto/exp_sub$KB_Pto
exp_sub$ratio <- ratio

# merge expression info with fit_con
fit_exp <- merge(fit_con, exp_sub, by.x="WP", by.y="WP", all.x = TRUE)

model2 <- lm(data=fit_exp, selection_Eyach ~ Time.in.tree + KB_Pto +Genetic_Diversity+ Col.0_Pto*variance )
summary(model2)
anova(model2)


# The expression level in planta seems to be the best predictor of the selection coefficient in planta. 
# What about if I look at genes that are differentially express
cor.test(exp$Col.0_Pto, exp$KB_Pto)
cor.test(exp$Col.0_Pto, exp$MM_Pto)
# The minimal media did not seem to be a better predictor of gene expression in the plant than did the 
# it is not clear to me why the ratio of Col0/KB was not a better predictor than the overall variance. 


# Let's do random forest modeling to look at the best predictors
# https://www.guru99.com/r-random-forest-tutorial.html
# Define the control
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")


set.seed(1234)
# Run the model
df <- fit_exp %>% select(selection_Eyach, KB_Pto, MM_Pto, SA_Pto, variance, Time.in.tree, Genetic_Diversity, Num_gene_events)
rownames(df) <- fit_exp$Row.names


rf_default <- train(selection_Eyach~.,
                    data = data_train,
                    method = "rf",
                    metric = "Accuracy",
                    trControl = trControl)
# Print the results
print(rf_default)

# search best mtry
set.seed(1234)
tuneGrid <- expand.grid(.mtry = c(1: 10))
rf_mtry <- train(selection_Eyach~.,
                 data = data_train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 14,
                 ntree = 300)
print(rf_mtry)

# search the best maxnodes
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(5: 15)) {
  set.seed(1234)
  rf_maxnode <- train(survived~.,
                      data = data_train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)

# search the best ntrees
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(5678)
  rf_maxtrees <- train(survived~.,
                       data = data_train,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = 24,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)

# evaluate the model
prediction <-predict(fit_rf, data_test)
confusionMatrix(prediction, data_test$survived)

#visualize result
varImpPlot(fit_rf)



