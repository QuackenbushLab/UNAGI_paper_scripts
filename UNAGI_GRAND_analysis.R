# Load the package to flatten the matrix.
library(reshape2)

# Load the network zoo (BLOBFISH package).
library(netZooR)

# Load igraph to perform graph-based analysis.
library(igraph)

# Convert Symbols to Ensembl IDs.
library("org.Hs.eg.db")

# Directories
unagiDir <- NULL
randResultDir <- NULL
tableDir <- NULL
plotDir <- NULL
dir.create(tableDir)
dir.create(plotDir)

# Read the networks. These are exclusive.
skeletal_muscle <- read.csv(paste0(unagiDir, "skeletal_muscle_only_updated.csv"), row.names = 1)
subcutaneous_adipose <- read.csv(paste0(unagiDir, "subcutaneous_adipose_only_updated.csv"), row.names = 1)
skin <- read.csv(paste0(unagiDir, "skin_only_updated.csv"), row.names = 1)
lung <- read.csv(paste0(unagiDir, "lung_only_updated.csv"), row.names = 1)
aorta <- read.csv(paste0(unagiDir, "aorta_only_updated.csv"), row.names = 1)

# Label outlying genes.
muscle_dev_genes <- c("MB", "MYH2", "MYL2", "DES", "TNNC1", "TNNC2", "ENO3", "MYL3",
                      "TTN", "TPM1", "TCAP", "MYL1", "TPM3", "TPM4", "COX5B",
                      "COX5A", "CKMT2", "TUBA1A", "TUBA1B", "TUBA4A", "TUBA1C", 
                      "TUBA3C", "TUBA8", "TUBA3D", "TUBA3E", "TUBA4B")
muscle_dev_ensembl <- mapIds(org.Hs.eg.db, keys = muscle_dev_genes, keytype = "SYMBOL", column="ENSEMBL")

# Adipogenic markers from 36647068
adipogenic_genes <- c("PPARG", "FASN", "SREBF1", "SCD", "CEBPA", "ADIPOQ", "FABP4")
adipogenic_ensembl <- mapIds(org.Hs.eg.db, keys = adipogenic_genes, keytype = "SYMBOL", column="ENSEMBL")
geneColorMapping <- data.frame(gene = c(muscle_dev_genes, adipogenic_genes), 
                               color = c(rep("hotpink", length(muscle_dev_genes)), 
                                         rep("goldenrod1", length(adipogenic_genes)))) 
colnames(skeletal_muscle) <- c("Var1", "Var2", "value")
colnames(subcutaneous_adipose) <- c("Var1", "Var2", "value")
colnames(skin) <- c("Var1", "Var2", "value")
colnames(lung) <- c("Var1", "Var2", "value")
colnames(aorta) <- c("Var1", "Var2", "value")

# Get all genes from each network (Var1 and Var2 combined)
get_all_genes <- function(network) {
  unique(as.character(c(network$Var1, network$Var2)))
}

muscle_genes_in_network <- get_all_genes(skeletal_muscle)
adipose_genes_in_network <- get_all_genes(subcutaneous_adipose)
skin_genes_in_network <- get_all_genes(skin)
lung_genes_in_network <- get_all_genes(lung)
aorta_genes_in_network <- get_all_genes(aorta)

# Check the average number of genes connecting each pair.
geneCountAllPairs <- function(genes, network){
  network[,1] <- as.character(network[,1])
  network[,2] <- as.character(network[,2])
  counts <- unlist(lapply(genes, function(i){
    print(i)
    return(unlist(lapply(setdiff(genes, i), function(j){
      tfsConnectedToi <- c(unique(network[which(network[,1] == i), 2]),
                           unique(network[which(network[,2] == i), 1]))
      tfsConnectedToj <- c(unique(network[which(network[,1] == j), 2]),
                           unique(network[which(network[,2] == i), 1]))
      return(length(intersect(tfsConnectedToi, tfsConnectedToj)))
    })))
  }))
  return(list(mean = mean(counts), sd = sd(counts)))
}
tfCountAllPairsMuscleMuscle <- geneCountAllPairs(muscle_dev_ensembl, skeletal_muscle)
tfCountAllPairsMuscleAdipose <- geneCountAllPairs(adipogenic_ensembl, skeletal_muscle)
tfCountAllPairsAdiposeMuscle <- geneCountAllPairs(muscle_dev_ensembl, subcutaneous_adipose)
tfCountAllPairsAdiposeAdipose <- geneCountAllPairs(adipogenic_ensembl, subcutaneous_adipose)
tfCountAllPairsSkinMuscle <- geneCountAllPairs(muscle_dev_ensembl, skin)
tfCountAllPairsSkinAdipose <- geneCountAllPairs(adipogenic_ensembl, skin)
tfCountAllPairsLungMuscle <- geneCountAllPairs(muscle_dev_ensembl, lung)
tfCountAllPairsLungAdipose <- geneCountAllPairs(adipogenic_ensembl, lung)
tfCountAllPairsAortaMuscle <- geneCountAllPairs(muscle_dev_ensembl, aorta)
tfCountAllPairsAortaAdipose <- geneCountAllPairs(adipogenic_ensembl, aorta)

# Do randomized analyses.
# Get the TF counts.
tissues <- c("adipose", "aorta", "lung", "muscle", "skin")
tissueSpecificMeans <- list()
tissueSpecificMeans[["adipose"]] <- list()
tissueSpecificMeans[["aorta"]] <- list()
tissueSpecificMeans[["lung"]] <- list()
tissueSpecificMeans[["muscle"]] <- list()
tissueSpecificMeans[["skin"]] <- list()
tissueSpecificSds <- list()
tissueSpecificSds[["adipose"]] <- list()
tissueSpecificSds[["aorta"]] <- list()
tissueSpecificSds[["lung"]] <- list()
tissueSpecificSds[["muscle"]] <- list()
tissueSpecificSds[["skin"]] <- list()

for(randsetIndex in 1:100){
  unagis <- list()
  unagis[["adipose"]] <- readRDS(paste0(randResultDir, "/Adipose_Subcutaneous_", randsetIndex, ".RDS"))
  unagis[["aorta"]] <- readRDS(paste0(randResultDir, "/Artery_Aorta_", randsetIndex, ".RDS"))
  unagis[["lung"]] <- readRDS(paste0(randResultDir, "/Lung_", randsetIndex, ".RDS"))
  unagis[["muscle"]] <- readRDS(paste0(randResultDir, "/Skeletal_Muscle_", randsetIndex, ".RDS"))
  unagis[["skin"]] <- readRDS(paste0(randResultDir, "/Skin_", randsetIndex, ".RDS"))
  
  for(tissue in tissues){
    allOtherunagis <- unlist(lapply(setdiff(tissues, tissue), function(tis){return(rownames(unagis[[tis]]))}))
    unagiRes <- unagis[[tissue]][setdiff(rownames(unagis[[tissue]]), allOtherunagis),]
    tfCounts <- geneCountAllPairs(genes = randsets[[randsetIndex]], network = unagiRes)
    tissueSpecificMeans[[tissue]][[randsetIndex]] <- tfCounts$mean
    tissueSpecificSds[[tissue]][[randsetIndex]] <- tfCounts$sd
  }
}

# Compute the mean values for mean and sd per tissue.
meanTissueSpecificMeans <- lapply(tissues, function(tissue){
  return(mean(unlist(tissueSpecificMeans[[tissue]])))
})
names(meanTissueSpecificMeans) <- tissues

meanTissueSpecificSds <- lapply(tissues, function(tissue){
  return(mean(unlist(tissueSpecificSds[[tissue]])))
})
names(meanTissueSpecificSds) <- tissues

# Make the table.
meanTable <- data.frame(muscle = c(tfCountAllPairsMuscleMuscle$mean, tfCountAllPairsMuscleAdipose$mean, meanTissueSpecificMeans[["muscle"]]),
                        adipose = c(tfCountAllPairsAdiposeMuscle$mean, tfCountAllPairsAdiposeAdipose$mean, meanTissueSpecificMeans[["adipose"]]),
                        aorta = c(tfCountAllPairsAortaMuscle$mean, tfCountAllPairsAortaAdipose$mean, meanTissueSpecificMeans[["aorta"]]),
                        lung = c(tfCountAllPairsLungMuscle$mean, tfCountAllPairsLungAdipose$mean, meanTissueSpecificMeans[["lung"]]),
                        skin = c(tfCountAllPairsSkinMuscle$mean, tfCountAllPairsSkinAdipose$mean, meanTissueSpecificMeans[["skin"]]))
rownames(meanTable) <- c("muscleGenes", "adiposeGenes", "randomGenes")
sdTable <- data.frame(muscle = c(tfCountAllPairsMuscleMuscle$sd, tfCountAllPairsMuscleAdipose$sd, meanTissueSpecificSds[["muscle"]]),
                        adipose = c(tfCountAllPairsAdiposeMuscle$sd, tfCountAllPairsAdiposeAdipose$sd, meanTissueSpecificSds[["adipose"]]),
                        aorta = c(tfCountAllPairsAortaMuscle$sd, tfCountAllPairsAortaAdipose$sd, meanTissueSpecificSds[["aorta"]]),
                        lung = c(tfCountAllPairsLungMuscle$sd, tfCountAllPairsLungAdipose$sd, meanTissueSpecificSds[["lung"]]),
                        skin = c(tfCountAllPairsSkinMuscle$sd, tfCountAllPairsSkinAdipose$sd, meanTissueSpecificSds[["skin"]]))
rownames(sdTable) <- c("muscleGenes", "adiposeGenes", "randomGenes")
write.csv(meanTable, paste0(tableDir, "/meanConnectorCount.csv"))
write.csv(sdTable, paste0(tableDir, "/sdConnectorCount.csv"))

zscoreTable <- data.frame(muscle = c((meanTable$muscle[1] - sdTable$muscle[3]) / meanTable$muscle[3],
                                     (meanTable$muscle[2] - sdTable$muscle[3]) / meanTable$muscle[3]),
                          adipose = c((meanTable$adipose[1] - sdTable$adipose[3]) / meanTable$adipose[3],
                                     (meanTable$adipose[2] - sdTable$adipose[3]) / meanTable$adipose[3]),
                          aorta = c((meanTable$aorta[1] - sdTable$aorta[3]) / meanTable$aorta[3],
                                      (meanTable$aorta[2] - sdTable$aorta[3]) / meanTable$aorta[3]),
                          lung = c((meanTable$lung[1] - sdTable$lung[3]) / meanTable$lung[3],
                                    (meanTable$lung[2] - sdTable$lung[3]) / meanTable$lung[3]),
                          skin = c((meanTable$skin[1] - sdTable$skin[3]) / meanTable$skin[3],
                                    (meanTable$skin[2] - sdTable$skin[3]) / meanTable$skin[3]))
write.csv(zscoreTable, paste0(tableDir, "/zConnectorCount.csv"))

# Check for direct connections.
percentConnected <- function(genes, network){
  network[,1] <- as.character(network[,1])
  network[,2] <- as.character(network[,2])
  edgeCount <- sum(unlist(lapply(genes, function(i){
    geneSrcConnections <- network[intersect(which(network[,1] == i),
                                            which(network[,2] %in% setdiff(genes, i))), 2]
    geneTgtConnections <- network[intersect(which(network[,2] == i),
                                            which(network[,1] %in% setdiff(genes, i))), 1]
    return(length(unique(c(geneSrcConnections, geneTgtConnections))))
  })))

  # We allow for double-counting.
  numEdgesInClique <- (length(genes) * (length(genes) - 1))
  print(edgeCount / numEdgesInClique)
  return(edgeCount / numEdgesInClique)
}
percentConnectedMuscleMuscle <- percentConnected(muscle_dev_ensembl, skeletal_muscle)
percentConnectedMuscleAdipose <- percentConnected(adipogenic_ensembl, skeletal_muscle)
percentConnectedAdiposeMuscle <- percentConnected(muscle_dev_ensembl, subcutaneous_adipose)
percentConnectedAdiposeAdipose <- percentConnected(adipogenic_ensembl, subcutaneous_adipose)
percentConnectedSkinMuscle <- percentConnected(muscle_dev_ensembl, skin)
percentConnectedSkinAdipose <- percentConnected(adipogenic_ensembl, skin)
percentConnectedLungMuscle <- percentConnected(muscle_dev_ensembl, lung)
percentConnectedLungAdipose <- percentConnected(adipogenic_ensembl, lung)
percentConnectedAortaMuscle <- percentConnected(muscle_dev_ensembl, aorta)
percentConnectedAortaAdipose <- percentConnected(adipogenic_ensembl, aorta)

meanPercentConnected <- list()
sdPercentConnected <- list()
unagis <- list()
unagis[["adipose"]] <- lapply(1:100, function(randsetIndex){
  return(readRDS(paste0(randResultDir, "/Adipose_Subcutaneous_", randsetIndex, ".RDS")))
})
unagis[["aorta"]] <- lapply(1:100, function(randsetIndex){
  return(readRDS(paste0(randResultDir, "/Artery_Aorta_", randsetIndex, ".RDS")))
})
unagis[["lung"]] <- lapply(1:100, function(randsetIndex){
  return(readRDS(paste0(randResultDir, "/Lung_", randsetIndex, ".RDS")))
})
unagis[["muscle"]] <- lapply(1:100, function(randsetIndex){
  return(readRDS(paste0(randResultDir, "/Skeletal_Muscle_", randsetIndex, ".RDS")))
})
unagis[["skin"]] <- lapply(1:100, function(randsetIndex){
  return(readRDS(paste0(randResultDir, "/Skin_", randsetIndex, ".RDS")))
})

for(tissue in tissues){
  perc <- unlist(lapply(1:100, function(randsetIndex){
    return(percentConnected(genes = randsets[[randsetIndex]], network = unagis[[tissue]][[randsetIndex]]))
  }))
  meanPercentConnected[[tissue]] <- mean(perc)
  sdPercentConnected[[tissue]] <- sd(perc)
}

resTable <- data.frame(muscle = c(percentConnectedMuscleMuscle, percentConnectedMuscleAdipose,
                                  (percentConnectedMuscleMuscle - meanPercentConnected[["muscle"]]) / sdPercentConnected[["muscle"]],
                                  (percentConnectedMuscleAdipose - meanPercentConnected[["muscle"]]) / sdPercentConnected[["muscle"]],
                                   meanPercentConnected[["muscle"]], sdPercentConnected[["muscle"]]),
                        adipose = c(percentConnectedAdiposeMuscle, percentConnectedAdiposeAdipose, 
                                    (percentConnectedAdiposeMuscle - meanPercentConnected[["adipose"]]) / sdPercentConnected[["adipose"]],
                                    (percentConnectedAdiposeAdipose - meanPercentConnected[["adipose"]]) / sdPercentConnected[["adipose"]],
                                    meanPercentConnected[["adipose"]], sdPercentConnected[["adipose"]]),
                        aorta = c(percentConnectedAortaMuscle, percentConnectedAortaAdipose,
                                  (percentConnectedAortaMuscle - meanPercentConnected[["aorta"]]) / sdPercentConnected[["aorta"]],
                                  (percentConnectedAortaAdipose - meanPercentConnected[["aorta"]]) / sdPercentConnected[["aorta"]],
                                  meanPercentConnected[["aorta"]], sdPercentConnected[["aorta"]]),
                        lung = c(percentConnectedLungMuscle, percentConnectedLungAdipose, 
                                 (percentConnectedLungMuscle - meanPercentConnected[["lung"]]) / sdPercentConnected[["lung"]],
                                 (percentConnectedLungAdipose - meanPercentConnected[["lung"]]) / sdPercentConnected[["lung"]],
                                 meanPercentConnected[["lung"]], sdPercentConnected[["lung"]]),
                        skin = c(percentConnectedSkinMuscle, percentConnectedSkinAdipose, 
                                 (percentConnectedSkinMuscle - meanPercentConnected[["skin"]]) / sdPercentConnected[["skin"]],
                                 (percentConnectedSkinAdipose - meanPercentConnected[["skin"]]) / sdPercentConnected[["skin"]],
                                 meanPercentConnected[["skin"]], sdPercentConnected[["skin"]]))
rownames(resTable) <- c("muscleGenes", "adiposeGenes", "muscleZ", "adiposeZ", "randomGenesMean", "randomGenesSd")
write.csv(resTable, paste0(tableDir, "/percentConnected.csv"))

zscoreTable <- data.frame(muscle = c((meanTable$muscle[1] - sdTable$muscle[3]) / meanTable$muscle[3],
                                     (meanTable$muscle[2] - sdTable$muscle[3]) / meanTable$muscle[3]),
                          adipose = c((meanTable$adipose[1] - sdTable$adipose[3]) / meanTable$adipose[3],
                                      (meanTable$adipose[2] - sdTable$adipose[3]) / meanTable$adipose[3]),
                          aorta = c((meanTable$aorta[1] - sdTable$aorta[3]) / meanTable$aorta[3],
                                    (meanTable$aorta[2] - sdTable$aorta[3]) / meanTable$aorta[3]),
                          lung = c((meanTable$lung[1] - sdTable$lung[3]) / meanTable$lung[3],
                                   (meanTable$lung[2] - sdTable$lung[3]) / meanTable$lung[3]),
                          skin = c((meanTable$skin[1] - sdTable$skin[3]) / meanTable$skin[3],
                                   (meanTable$skin[2] - sdTable$skin[3]) / meanTable$skin[3]))
write.csv(zscoreTable, paste0(tableDir, "/zConnectorCount.csv"))

# Plot networks.
filterNetAndConvertSymbols <- function(net, cutoff){
  filtNetDirect <- net[intersect(which(net[,1] %in% c(muscle_dev_ensembl, adipogenic_ensembl)),
                             which(net[,2] %in% c(muscle_dev_ensembl, adipogenic_ensembl))),]
  filtNetDirect[,1] <- mapIds(org.Hs.eg.db, keys = as.character(filtNetDirect[,1]), 
                              keytype = "ENSEMBL", column="SYMBOL")
  filtNetDirect[,2] <- mapIds(org.Hs.eg.db, keys = as.character(filtNetDirect[,2]), 
                              keytype = "ENSEMBL", column="SYMBOL")
  print(dim(filtNetDirect))
  filtNetComplete <- filtNetDirect[intersect(which(!is.na(filtNetDirect[,1])),
                                       which(!is.na(filtNetDirect[,2]))),]
  print(dim(filtNetComplete))
  toRemove <- unlist(lapply(1:nrow(filtNetComplete), function(i){
    whichSrc <- which(filtNetComplete[,2] == filtNetComplete[i,1])
    whichTgt <- which(filtNetComplete[,1] == filtNetComplete[i,2])
    whichSorted <- which(filtNetComplete[,1] < filtNetComplete[,2])
    return(Reduce(intersect, list(whichSrc, whichTgt, whichSorted)))
  }))
  filtNetDedup <- filtNetComplete[setdiff(1:nrow(filtNetComplete), toRemove),]
  rownames <- paste(filtNetDedup[,1], filtNetDedup[,2], sep = "_")
  toRemove2 <- unlist(lapply(unique(rownames), function(pair){
    retval <- NA
    whichPair <- which(rownames == pair)
    if(length(whichPair) > 1){
      retval <- whichPair[2:length(whichPair)]
    }
    return(retval)
  }))
  filtNetDedup <- filtNetDedup[setdiff(1:nrow(filtNetDedup), toRemove2),]
  print(dim(filtNetDedup))
  return(filtNetDedup)
}
skeletal_muscle_filt <- filterNetAndConvertSymbols(skeletal_muscle)
subcutaneous_adipose_filt <- filterNetAndConvertSymbols(subcutaneous_adipose)
lung_filt <- filterNetAndConvertSymbols(lung)
aorta_filt <- filterNetAndConvertSymbols(aorta)
skin_filt <- filterNetAndConvertSymbols(skin)
write.csv(skeletal_muscle_filt, paste0(plotDir, "skeletalMuscle.csv"))
write.csv(subcutaneous_adipose_filt, paste0(plotDir, "adipose.csv"))
write.csv(lung_filt, paste0(plotDir, "lung.csv"))
write.csv(aorta_filt, paste0(plotDir, "aorta.csv"))
write.csv(skin_filt, paste0(plotDir, "skin.csv"))

PlotNetworkU(skeletal_muscle_filt, geneColorMapping = geneColorMapping,
             nodeSize = 3)
PlotNetworkU(subcutaneous_adipose_filt, geneColorMapping = geneColorMapping,
             nodeSize = 3)
PlotNetworkU(lung_filt, geneColorMapping = geneColorMapping,
             nodeSize = 3)
PlotNetworkU(aorta_filt, geneColorMapping = geneColorMapping,
             nodeSize = 3)
PlotNetworkU(skin_filt, geneColorMapping = geneColorMapping,
             nodeSize = 3)

sa <- subcutaneous_adipose
sa[,1] <- as.character(subcutaneous_adipose[,1])
sa[,2] <- as.character(subcutaneous_adipose[,2])