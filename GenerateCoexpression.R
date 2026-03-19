# Install and load packages.
library(Matrix)
library(parallel)
if(!require("reshape2")){
  install.packages("reshape2")
}
library("propagate")
if(!require("propagate")){
  install.packages("propagate")
}
library("propagate")
library(qlcMatrix)

# Change to directory where expression data are stored.
expressionDir <- NULL
coexpressionDir <- NULL
dir.create(coexpressionDir)

# Function to compute correlation
computeCor <- function(data){
  corMat <- propagate::bigcor(data, fun = c("cor"), size = 2000)
  meltedCorMat <- reshape2::melt(corMat)
  absCorMat <- abs(meltedCorMat[,3])
  decile <- quantile(absCorMat, probs = 0.9)
  corMatFilt <- meltedCorMat[which(absCorMat > decile),]
  return(corMatFilt)
}

# ===== Load data =====
ocDataSkin <- read.csv(paste0(expressionDir, "/Skin.csv"), row.names = 1)
ocDataAdipose <- read.csv(paste0(expressionDir, "/Adipose_Subcutaneous.csv"), row.names = 1)
ocDataAorta <- read.csv(paste0(expressionDir, "/Artery_Aorta.csv"), row.names = 1)
ocDataLung <- read.csv(paste0(expressionDir, "/Lung.csv"), row.names = 1)
ocDataMuscle <- read.csv(paste0(expressionDir, "/Skeletal_Muscle.csv"), row.names = 1)

doCorAndSave <- function(ocData, outFileName){
  ocData <- as.matrix(ocData)
  mode(ocData) <- "numeric"
  ocResult <- computeCor(ocData)
  
  # Use colnames since we're transposing the data
  S <- threshold_sparse_cov_fast_mc(t(ocData), gene_names = rownames(ocData), parallel = TRUE)
  
  # Save result
  saveRDS(S, paste0(coexpressionDir, "/", outFileName))
}
doCorAndSave(ocDataSkin, "skin_S.RDS")
doCorAndSave(ocDataAdipose, "adipose_S.RDS")
doCorAndSave(ocDataAorta, "artery_S.RDS")
doCorAndSave(ocDataLung, "lung_S.RDS")
doCorAndSave(ocDataMuscle, "skeletal_S.RDS")

