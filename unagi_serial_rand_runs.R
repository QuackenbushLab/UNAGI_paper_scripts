# This install the devel branch of netZooR.
if(!require(netZooR)){
  devtools::install_github("netZoo/netZooR@devel")
}
library(netZooR)

# Read in the files.
coexpressionDir <- NULL
randomSetDir <- NULL
randResultDir <- NULL
dir.create(randResultDir)

cellTypes <- c("Adipose_Subcutaneous", "Artery_Aorta", "Lung", "Skeletal_Muscle", "Skin")
randsetCount <- 100

for(cellType in cellTypes){
  for(i in 1:randsetCount){
    tissue <- readRDS(paste0(coexpressionDir, cellType, "_sm_S.RDS"))
    randsets <- readRDS(paste0(randomSetDir, "randSets.RDS"))
    
    # Run UNAGI on each randomization, for each tissue.
    tissue[,1] <- as.character(tissue[,1])
    tissue[,2] <- as.character(tissue[,2])
    colnames(tissue) <- c("source", "target", "score")
    unagi <- netZooR::RunUNAGI(nodeSet = randsets[[i]],
                               network = tissue, hopConstraint = 2,
                               verbose = TRUE)
    saveRDS(unagi, paste0(randResultDir, "/", cellType, "_", i, ".RDS"))
  }
}
