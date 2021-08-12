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
if (!requireNamespace("MatrixGenerics", quietly = TRUE))
  BiocManager::install("MatrixGenerics")

library(GSEABase)
library(celldex)
library(singscore)
library(ggplot2)
library(lubridate)
library(here)
library(readxl)
library(minpack.lm)
library(corrplot)

working.dir= paste0(here(), "/src")
filePath = paste0(working.dir, "/xCell_train.R")
source(filePath)
train.function.path= paste0(working.dir, "/xCell_train_functions.R")
source(train.function.path)

xCell.path= paste0(working.dir, "/xCell.R")
source(xCell.path)


sdy = readRDS(paste0(working.dir, '/sdy420.rds'))
genes.to.use= rownames(sdy$expr)
xCell.my.data= xCellReferenceGenerate(save.figures=FALSE, genes.to.use)




sed= SummarizedExperiment(assays=list(counts=(as.matrix(sdy$expr))))

scores= weighted.score(sed, xCell.my.data$signaturesRank, xCell.my.data$signatures, xCell.my.data$types)
transfomed.score= transform.mix.score(xCell.my.data$fitValues, scores)
after.spill= spillOver(transformedScores = transfomed.score, xCell.my.data$spill)


fcs = sdy$fcs#[rownames(after.spill),colnames(sdy$expr)]
fcs <- fcs[ , order(colnames(fcs))]
after.spill <- after.spill[ , order(colnames(after.spill))]
res = cor(t(after.spill),t(fcs))
# res= res[order(rownames(res)), ]
# res= res[, order(colnames(res))]
maxxr <- apply(res, 1, max)
result <- res[order(-maxxr),]
maxxc <- apply(result, 2, max)
result <- result[,order(-maxxc)]

orgnizedRes= matrix(0, nrow =dim(res)[1] )

for (i in 1:dim(res)[1]){
  
  
}

plot.new()
jpeg(file="TrainedOnMonacoImmuneDataFine.jpg")
corrplot(result, is.corr = F)
dev.off()
heatmap(res)
