if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GSEABase", quietly = TRUE))
  BiocManager::install("GSEABase")

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
filePath = paste0(working.dir, "/xCell_train_functions_old.R")
source(filePath)

# ----------------- Preparing the data -----------------
# Read database- should be from user.
ref <- BlueprintEncodeData()
# get samples names
samples= ref$label.main
types=unique(samples)


# ----------------- Generating Signatures -----------------
gsc <- create.signatures(ref, types, samples, NA)

# Create in sillico mixtures.
inSillicoMix= create.inSillicoMixtures(ref, gsc, types)

# get mix score.
scores= create.signatureScores(ref, gsc,  inSillicoMix, types)

# Choose best signatures.
bestSig<- choose.signatures(inSillicoMix, types, scores)


# ----------------- Corrlation testing -----------------
gscSignatures= gsc[c(bestSig)]

# Create percentage for in-Sillico mix.
percentage= matrix(NA,length(types),20*ntypes)
colnames(percentage)<-1:(20*ntypes)
rownames(percentage)= types


ntypes= length(types)
i=1
for (type in types){
for (k in 1:20){
    randPrec= sample(100, ntypes, replace=T)
    sumPrec= sum(randPrec)-randPrec[i]
    prec= (randPrec/sumPrec)*((100-k)/100)
    prec[i]= k/100
    percentage[,(i-1)*20+k]= prec

}
  i=i+1
}

mixtures= create.corrlationMixtures(percentage, ref)
mixScores= testing.correlation(mixtures, types, precentage, gscSignatures, bestSig)

dir.create(file.path(working.dir,"LinearRegressionGraphs"), showWarnings = FALSE)
for (type in types){
  x= percentage[type,]
  y= mixScores[type, ]
  df = data.frame(x, y)
  line = lm(y~x)
  pred_line = predict(line)
  LinerRegressionPath= file.path(working.dir,"LinearRegressionGraphs",
                                 paste0(type,"LinearRegression.jpg"))
  jpeg(file=LinerRegressionPath)
  plot.new()

  plot(y, x,main=paste0(type," Linear Regression"),
       xlab="percentage", ylab="Score", pch = 11, col = "blue")
  #abline(line, col="red")
  lines(pred_line,x,col="red")
  dev.off()


}


parameterMatrix= matrix(NA, ntypes, 2)
rownames(parameterMatrix)=types
colnames(parameterMatrix)<-c("b", "calib")
dir.create(file.path(working.dir,"LinearTransformGraphs"), showWarnings = FALSE)
for (i in 1:ntypes){
  type= types[i]
  df = data.frame(x=mixScores[i,],y=percentage[i,])
  df$x=df$x-min(df$x)
  orderdInd= order(df$x)
  df$x= df$x[orderdInd]
  df$y= df$y[orderdInd]

  z = nlsLM(formula = y~(a*x^b),start = list(a=1,b=1),data = df)
  b = coef(z)[2]
  calib = coef(lm(df$x^b~df$y))[2]
  parameterMatrix[i, ]= c(b, calib)

  LinerTransformPath= file.path(working.dir,"LinearTransformGraphs",
                                 paste0(type,"LinerTransform.jpg"))
  jpeg(file=LinerTransformPath)
  plot.new()
  plot((df$x^b)/calib,df$y, main=paste0(type," Linear Transform"),
      xlab="percentage", ylab="Score", col="blue")

  x = (df$x^b)/calib
  y = df$y
  line = lm(x ~ y)
  pred_line = predict(line)
  lines(pred_line, y, col="red")
  dev.off()
}


mix_scores_2 = mixScores
for(i in 1:dim(mixScores)[1]){
  vec = mix_scores_2[i,]
  vec = vec + abs(min(vec))
  vec = vec^parameterMatrix[i,1]
  vec = vec/parameterMatrix[i,2]
  mix_scores_2[i,] = vec
}


corr= cor(t(mix_scores_2), t(percentage))
# Plot the corralation matrix before imp
corrMatrixPath= file.path(working.dir,"corrlation matrix before spillover.jpg")
jpeg(file=corrMatrixPath)
plot.new()
corrplot(corr, order = 'hclust', main="corrlation matrix before spillover")
dev.off()
jpeg(file = NULL)



# -----------------Correlation and spillover matrix--------------

refmixData= matrix(NA, length(rownames(ref)),ntypes)
rownames(refmixData)= rownames(ref)
colnames(refmixData)= types


refMixPerc= matrix(0, ntypes, ntypes)
rownames(refMixPerc)=types
colnames(refMixPerc)<-1:ntypes
control.type<-c()

for(i in 1:dim(corr)[2]){
  corrPercenatages = abs(corr[,i])
  minInd = which(corrPercenatages==min(corrPercenatages))
  control.type<-c(control.type,minInd)
  contorlTypeName = colnames(corr)[minInd]
  typeName = colnames(corr)[i]
  print(c("base type:",typeName," min corr type:", contorlTypeName))

  locationControl= match(contorlTypeName, ref$label.main)
  controlData= assays(ref)$logcounts[,locationControl]*0.8
  locationType= match(typeName, ref$label.main)
  typeData= assays(ref)$logcounts[,locationType]*0.2
  refmixData[,i] <- controlData+typeData
  refMixPerc[i, i]= 0.2
  refMixPerc[i, minInd]=0.8
}


refmixExp= SummarizedExperiment(assays=list(counts=refmixData))
refmixRankedData<-rankGenes(refmixExp)

refMixScore= matrix(NA, ntypes, ntypes)
colnames(refMixScore)<-types
rownames(refMixScore)=1:ntypes

for (i in 1:ntypes){
  sigs= bestSig[i,]
  nsigs= length(sigs)
  sumWeight= sum(1:nsigs)
  j= nsigs
  totScore= rep(0,ntypes)
  for (sig in sigs){
    gs= gscSignatures[[sig]]
    totScore <-totScore+ (j/sumWeight)*(simpleScore(refmixRankedData, geneIds(gs), centerScore = TRUE)$TotalScore)
    j= j-1
  }

  refMixScore[i, ]= totScore
}


transformedRfmix= refMixScore
# Use transformation
for(i in 1:dim(transformedRfmix)[1]){
  vec = transformedRfmix[i,]
  vec = vec - (min(vec))
  vec = vec^parameterMatrix[i,1]
  vec = vec/parameterMatrix[i,2]
  transformedRfmix[i,] = vec
}

transformedRefmix[cbind(1:dim(transformedRefmix)[1],controls)] = 0

spilloverMat= 0.2*((transformedRfmix)) / diag((transformedRfmix))

spill_corr= cor(t(spilloverMat), (refMixPerc))
corrMatrixPath= file.path(working.dir,"spillovermat.jpg")
jpeg(file=corrMatrixPath)
plot.new()
corrplot(spilloverMat, order = 'hclust', main="spillovermat.jpg")
dev.off()
jpeg(file = NULL)
corrplot(spill_corr)



