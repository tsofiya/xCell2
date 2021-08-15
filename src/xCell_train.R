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
library(ggpubr)

xCellReferenceGenerate= function(ref= NULL, save.figures= TRUE, genes.to.use){
  # ----------------- Preparing the data -----------------
  # Read database- should be from user.
  working.dir= paste0(here(), "/src")
  filePath = paste0(working.dir, "/xCell_train_functions.R")
  source(filePath)
  
  if (is.null(ref))
    ref <- BlueprintEncodeData()
      
  genes= intersect(genes.to.use,rownames(ref))
  ref= ref[genes, ]
  
  # get samples names
  samples= ref$label.main
  types=unique(samples)
  ntypes= length(types)
  singleCellData= create.single.cell.data.mat(ref, types) # single cell data for all cell types- maybe change to average?
  
  
  # ----------------- Generating Signatures -----------------
  # Generating signatures. 
  gsc <- create.signatures(ref, types, samples, NA)
  
  # Generating percentage matrix in preperation for mixtures simulations
  percentage= create.percentage.matrix(types, (1:250)/1250)
  
  # Generating in-sillico mixtures
  mixtures= create.mixtures(singleCellData, percentage)
  
  # scoring mixtures according to the generated signatures. 
  scoreMatrix= score.signature(mixtures, gsc)
  
  # Testing the corrletion of the score and the real percentages and returns the IDs of the best signatures
  bestSignatures= choose.signatures(percentage, types, scoreMatrix)
  bestGsc= gsc[c(bestSignatures)]
  
  # ----------------- Corrlation testing and Linear Transformation -----------------
  # Generating percentage matrix for linear regression
  percentage= create.percentage.matrix(types, (1:20)/100)
  
  # Generating in-sillico mixtures
  mixtures= create.mixtures(singleCellData, percentage)
  
  # scoring mixtures according to the generated signatures. 
  weightedScoreMat= weighted.score(mixtures, bestSignatures, bestGsc, types)
  # if (save.figures){
  #   for (type in types){
  #     score = weightedScoreMat[type, ]
  #     df = data.frame(CellType=types, Score = score)
  #     ggplot(df,aes(x=CellType,y=score,fill=CellType))+geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ggtitle(paste0(type, " signature"))
  #     ggsave( paste0(paste0(type, " signature"),".png"))
  #   }
  # }
  
  if (save.figures)
    plot.linear.regression(working.dir, weightedScoreMat, percentage, types, save.figures)
  #working.dir, mixScores, percentage, types, save.graph= FALSE
  
  parameterMatrix= create.parameter.matrix(working.dir, types, weightedScoreMat, percentage, save.figures)
  
  
  transformedMixScore= transform.mix.score(parameterMatrix, weightedScoreMat)
  
  
  
  corr= cor(t(transformedMixScore), (percentage), method = "spearman")
  # Plot the corralation matrix before im
  if (save.figures){
    corrMatrixPath= file.path(working.dir,"corrlation matrix before spillover.jpg")
    jpeg(file=corrMatrixPath)
    plot.new()
    corrplot(corr, order = 'hclust', main="corrlation matrix before spillover")
    dev.off()
    #jpeg(file = NULL)
  }
  
  
  # -----------------Correlation and spillover matrix--------------
  corr= cor(t(transformedMixScore), t(transformedMixScore), method = "spearman")
  controls=apply(corr,1,which.min)
  refPercentage= create.ref.percentage(types, corr, controls)
  refmixExp=create.mixtures(singleCellData, refPercentage)
  refMixScore= weighted.score(refmixExp, bestSignatures, bestGsc, types)
  transformedRefmix= transform.mix.score(parameterMatrix, refMixScore)
  colnames(transformedRefmix)= types
  spilloverMat= 0.2*((transformedRefmix)) / diag((transformedRefmix))
  spilloverMat[cbind(1:ntypes, controls)]=0
  #spilloverMat=spilloverMat[row.names(spilloverMat) != "Erythrocytes", , drop = FALSE] 
  #spilloverMat=spilloverMat[,-c(9)]#percentageMat[cbind(1:length(types), controls)]= 0
  
  if (save.figures){
    corrplot(spilloverMat, is.corr = FALSE, main="spillover matrix")
  #with(mtcars, corrplot(spilloverMat, is.corr = FALSE, main="spillover matrix"), cl.lim = c(-0.2, 0.2))
  
  spillMatPath= file.path(working.dir,"spillover matrix.jpg")
  jpeg(file=spillMatPath)
  plot.new()
  corrplot(spilloverMat,  is.corr = FALSE, main="spillover matrix")
  dev.off()
  
  
  #jpeg(file = NULL)
  }
  
  xCell.data = list(spill=spilloverMat,fitValues= parameterMatrix,signatures=bestGsc, signaturesRank= bestSignatures, genes=ref@NAMES, types= types)
  save(xCell.data,file=file.path(here(),'xCell.data.Rdata'))
  xCell.data
  
}
