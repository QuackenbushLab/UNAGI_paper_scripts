# Load the network zoo.
library("netZooR")

# Read in the files.
filteredCoexpressionDir <- NULL
randomSetDir <- NULL
skeletal_muscle <- readRDS(paste0(filteredCoexpressionDir, "Skeletal_Muscle_sm_S.RDS"))
subcutaneous_adipose <- readRDS(paste0(filteredCoexpressionDir, "Adipose_Subcutaneous_sm_S.RDS"))
skin <- readRDS(paste0(filteredCoexpressionDir, "Skin_sm_S.RDS"))
lung <- readRDS(paste0(filteredCoexpressionDir, "Lung_sm_S.RDS"))
aorta <- readRDS(paste0(filteredCoexpressionDir, "Artery_Aorta_sm_S.RDS"))
print("read all files")

# Generate random sets.
num_rand_sets <- 100
num_genes <- 33
sharedGenes <- unique(c(skeletal_muscle$Var1, skeletal_muscle$Var2))
for(tissue in list(skeletal_muscle, subcutaneous_adipose, skin, lung, aorta)){
  for(i in 1:length(tissue)){
    sharedGenes <- intersect(sharedGenes, unique(c(tissue$Var1, tissue$Var2)))
  }
}
str(sharedGenes)
randsets <- lapply(1:num_rand_sets, function(i){
  return(sample(sharedGenes, num_genes))
})
str(randsets)
saveRDS(randsets, paste0(randomSetDir, "randSets.RDS"))
