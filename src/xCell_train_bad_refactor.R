if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("celldex", quietly = TRUE))
  BiocManager::install("celldex")
if (!requireNamespace("singscore", quietly = TRUE))
  BiocManager::install("singscore")
if (!requireNamespace("GSEABase", quietly = TRUE))
  BiocManager::install("GSEABase")
if (!requireNamespace("lubridate", quietly = TRUE))
  install.packages("lubridate")
if (!requireNamespace("here", quietly = TRUE))
  install.packages("here")
if (!requireNamespace("minpack.lm", quietly = TRUE))
  install.packages("minpack.lm")

library(GSEABase)
library(celldex)
library(singscore)
library(ggplot2)
library(lubridate)
library(here)
library(readxl)
library(minpack.lm)
library(corrplot)

#/home/gidi/repos/xCell2/xCell/Dev_scripts
working.dir= paste0(here(), "/src")
filePath = paste0(working.dir, "/xCell_train_functions.R")
source(filePath)

# ----------------- Preparing the data -----------------
# Read database- should be from user.
ref <- BlueprintEncodeData()
# get samples names
samples= ref$label.main
types=unique(samples)
ntypes= length(types)
singleCellData= create.single.cell.data.mat(ref, types) # single cell data for all cell types- maybe change to average?

# ----------------- Generating Signatures -----------------
gsc <- create.signatures(ref, types, samples, NA)

# Create in sillico mixtures.

inSillicoMix= create.inSillicoMixtures(ref, gsc, types, singleCellData)

# get mix score.
scores= create.signatureScores(ref, gsc,  inSillicoMix, types)

# Choose best signatures.
bestSig<- choose.signatures(inSillicoMix, types, scores)


# ----------------- Corrlation testing and Linear Transformation -----------------
gscSignatures= gsc[c(bestSig)]

percentage= create.random.percentage(types)
mixtures= create.corrlationMixtures(percentage, ref, singleCellData)
mixScores= testing.correlation(mixtures, types, precentage, gscSignatures, bestSig)
plot.linear.regression(working.dir, mixScores, percentage)
parameterMatrix= create.parameter.matrix(working.dir, types, mixScores, percentage)

transformedMixScore= transform.mix.score(parameterMatrix, mixScores)

corr= cor(t(transformedMixScore), t(percentage))
# Plot the corralation matrix before imp
corrMatrixPath= file.path(working.dir,"corrlation matrix before spillover.jpg")
jpeg(file=corrMatrixPath)
plot.new()
corrplot(corr, order = 'hclust', main="corrlation matrix before spillover")
dev.off()
jpeg(file = NULL)

# -----------------Correlation and spillover matrix--------------

refmixExp=create.ref.mix(types, singleCellData, corr)
refmixRankedData<-rankGenes(refmixExp)
refMixScore= score.ref.mix(types, bestSig, gscSignatures)


transformedRefmix= transform.mix.score(parameterMatrix, refMixScore)

spilloverMat= 0.2*((transformedRefmix)) / diag((transformedRefmix))
corrplot(spilloverMat)

spillMatPath= file.path(working.dir,"spillover matrix.jpg")
jpeg(file=spillMatPath)
plot.new()
corrplot(spilloverMat,  main="spillover matrix")
dev.off()
jpeg(file = NULL)
