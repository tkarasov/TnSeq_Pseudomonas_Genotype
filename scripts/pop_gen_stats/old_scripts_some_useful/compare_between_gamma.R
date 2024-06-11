#the goal of this script is to compare genes across gammproteobacteria

find_shared_conditions <-{}





















distance_matrix = read.table("/ebio/abt6_projects9/tnseq/data/across_genera/distance_matrix_across_strains.txt", sep = ",")
orthology = read.table("/ebio/abt6_projects9/tnseq/data/pan_genome_datasets/gammaproteobacteria_mappings.txt", sep = ',', fill = TRUE)
conditions = read.table("/ebio/abt6_projects9/tnseq/data/fitness_datasets/all_conditions.txt", header = TRUE, sep = "\t")


file.names <- list.files(path = "/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/")
# read in each file in the directory naming it with the interesting bit of the filename

genome_list = c(length(file.names))
for(i in 1:length(file.names)){
  start.stripped.i <- unlist(strsplit(x = file.names[i], split = 'organism_'))[2]
  obj.name.i <- unlist(strsplit(x = start.stripped.i, split = '\\.'))[1] # escape character before . so it's not treated as a wildcard 
  X.i <- read.csv(file.names[i], sep="\t")
  assign(x = obj.name.i, value = X.i)
  genome_list[i]=obj.name.i
  rm(list = c('start.stripped.i', 'obj.name.i', 'X.i'))
  gc()
}


