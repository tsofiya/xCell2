library(ggplot2)
library(ggpubr)

create.single.cell.data.mat= function(ref, types, is.main){
  ntypes= length(types)
  singleCellData= matrix(NA, ntypes, length(rownames(ref)))
  rownames(singleCellData)= types
  colnames(singleCellData)<- rownames(ref)
  if ( is.main)
    samples= ref$label.main
  else
    samples=ref$label.fine
  
  for (type in types){
    location= match(type, samples)
    singleCellData[type,] <- assays(ref)$logcounts[,location]
  }
  
  singleCellData
}

create.signatures = function(ref, types, samples, dependencies) {
  probs = c(.1,.25,.33333333,.5,.6666666,.75,.9)
  diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2,3,4,5)
  
  #TODO: add dependencies
  message('get quantiles...')
  
  
  q = lapply(types,function(y){
    A=samples==y
    if (sum(A)==1) {
      ex = cbind(assays(ref)$logcounts[,samples==y],assays(ref)$logcounts[,samples==y])
    } else {
      ex = assays(ref)$logcounts[,A]
    }
    q = apply(ex,1,function(z) quantile(z,probs,na.rm=TRUE))
    q
  })
  
  message('create all signatures...')
  ntypes= length(types)
  #signature= data.frame(ref@colData)
  listIn= list()
  rankData <- rankGenes(ref)
  counter <-0
  set_names= c()
  for (i in 1:ntypes) {
    for (diff in diff_vals) {
      for (j in 1:round(length(probs)/2+0.5)) {
        diff_genes = lapply(q,function(x) q[[i]][j,]>x[length(probs)-j,]+diff)
          output <- matrix(unlist(diff_genes), nrow = ntypes, byrow = TRUE)
          for (n in (ntypes-1):(ntypes-3)) {
            g = colSums(output)>=n
            if (sum(g)>7 & sum(g)<201){
              set_name= sprintf("%s-%g-%g-%g-%g",types[i],round(probs[j]*100),diff,ntypes-n, i)
              set_names<- c(set_names, set_name)
              gs=GeneSet(rownames(ref)[g],
                         setName=set_name)
              counter<-counter+1
              listIn<-c(listIn, c(gs))
          
          }
        }
      }
    }
  }
  gsc= GeneSetCollection(listIn)
  gsc
}



#' Title
#'
#' @param types - a vector of the unique cell type name
#' @param cellPercentage - a vector of the wanted percentage (percentage in fractions [0,1])
#'
#' @return a matrix of percentage. Each cell type will have an expirement in which he has cellPercentage fraction, and the other are chosen randomly.  
#' @export
#'
#' @examples
create.percentage.matrix= function(types, cellPercentage){
  ntypes= length(types)
  experiemts= ntypes*length(cellPercentage)
  percentageMat= matrix(NA, experiemts, length(types)) 
  rownames(percentageMat)= 1:experiemts
  colnames(percentageMat) <- types
  
  message(c("creating precentage mat for: ", cellPercentage))
  
  i=1
  k=1
  for (type in types){
    for (cellprec in cellPercentage){
      randPrec= sample(100, ntypes, replace=T)
      sumPrec= sum(randPrec)-randPrec[i]
      prec= (randPrec/sumPrec)*(1-cellprec)
      prec[i]= cellprec
      percentageMat[k,]= prec
      message(c("current type: ", type, ". current percentage: ", cellprec, ". row sum: ", sum(prec), "."))
      k=k+1
    }
    i=i+1
  }
  
  percentageMat
}

#' Title
#'
#' @param singleCellData - RNA single cell data of all cell types. one sample per type.
#' @param percentage - the wanted percentage of the cell type in the mixtures. 
#'
#' @return a matrix containing RNA mixed data in the wanted percentage
#' @export
#'
#' @examples
create.mixtures= function(singleCellData, percentage){
  experiments= length(rownames(percentage))
  countMix= matrix(NA, length(colnames(singleCellData)), experiments)
  rownames(countMix)= colnames(singleCellData)
  colnames(countMix)= 1:experiments
  
  for (i in 1:experiments){
    perc= percentage[i, ]
    #mix= rep(0, length(ref@NAMES))
    mix = colSums(perc * (singleCellData)) # My Gideon's line
    countMix[,i]= mix
  }
  
  mixtures= SummarizedExperiment(assays=list(counts=countMix))
  mixtures
}

#' Title
#'
#' @param mixtures 
#' @param signatures 
#'
#' @return a matrix containing simplescore result of the signatures over the mixtures. 
#' @export
#'
#' @examples
score.signature= function(mixtures, signatures){
  message("Scoring signatures")
  rankData= rankGenes(mixtures)
  experiments= length(colnames(mixtures))
  scoreMatrix= matrix(NA, length(signatures), experiments)
  rownames(scoreMatrix)= names(signatures)
  colnames(scoreMatrix)= 1:experiments
  
  for (gs in signatures){
    score <- simpleScore(rankData, geneIds(gs), centerScore = TRUE)$TotalScore
    scoreMatrix[gs@setName,]= score
  }
  scoreMatrix
}

choose.signatures = function(percentage, types, scores){
  message("Chosing best signatures")
  ntypes= length(types)
  bestSig= matrix(NA, ntypes, 5)
  colnames(bestSig)<- 1:5
  rownames(bestSig)= types
  
  for (type in types){
    typeScores= scores[startsWith(rownames(scores), type),]
    perc= (percentage[, type])
    corrlation= cor(t(typeScores), perc, method= "spearman")
    top5= head( order(corrlation, decreasing = TRUE), 5)
    top5Sig= rownames(corrlation)[top5]
    bestSig[type,]=top5Sig
  }
  bestSig
}

weighted.score= function(mixtures, bestSignatures, bestGsc, types){
  message("Calculating weighted score")
  mixRankedData<-rankGenes(mixtures)
  ntypes= length(types)
  experiments= length(colnames(mixtures))
  
  # Score matrix
  mixScore= matrix(NA, ntypes, experiments)
  rownames(mixScore)=types
  colnames(mixScore)= colnames(mixtures)
  
  for (type in types){
    sigs= bestSignatures[type,]
    nsigs= length(sigs)
    sumWeight= sum(1:nsigs)
    j= nsigs
    totScore= rep(0,experiments)
    for (sig in sigs){
      gs= bestGsc[[sig]]
      totScore <-totScore+ (j/sumWeight)*(simpleScore(mixRankedData, geneIds(gs), centerScore = TRUE)$TotalScore)
      j= j-1
    }
    
    mixScore[type, ]= totScore
  }
  mixScore
}

plot.linear.regression= function(working.dir, mixScores, percentage, types, save.graph= FALSE){
  message("Plotting linear regression")
  dir.create(file.path(working.dir,"LinearRegressionGraphs"), showWarnings = FALSE)
  
  panel= matrix(0, length(types)*dim(mixScores)[2], 3)
  colnames(panel)=c("Score", "Percentage", "cell.type")
  
  panel[, 1]= as.vector(t(mixScores))
  panel[,2 ]= as.vector((percentage))
  panel[,3]= rep(types, each = dim(mixScores)[2])
  
  paneldf= as.data.frame(panel, stringsAsFactors = FALSE)
  paneldf$Score= as.numeric(paneldf$Score)
  paneldf$Percentage= as.numeric(paneldf$Percentage)
  ggplot(paneldf,mapping= aes(Score,Percentage))+
    geom_smooth(method=lm, se=FALSE, color="red")+
    geom_point(size=0.1)+facet_wrap(.~cell.type)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name	="rho")
  if (save.graph)
    ggsave(filename="Linear Regression Panel.jpg")
}

create.parameter.matrix= function(working.dir, types, mixScores, percentage, save.graph= FALSE){
  message("Calculate fit values")
  ntypes= length(types)
  parameterMatrix= matrix(NA, ntypes, 2)
  rownames(parameterMatrix)=types
  colnames(parameterMatrix)<-c("b", "calib")
  dir.create(file.path(working.dir,"LinearTransformGraphs"), showWarnings = FALSE)
  for (type in types){
    df = data.frame(x=mixScores[type,],y=percentage[,type])
    df$x=df$x-min(df$x)
    orderdInd= order(df$x)
    df$x= df$x[orderdInd]
    df$y= df$y[orderdInd]
    
    z = nlsLM(formula = y~a*x^b,start = list(a=1,b=1),data = df)
    b = coef(z)[2]
    calib = coef(lm(df$x^b~df$y))[2]
    parameterMatrix[type, ]= c(b, calib)
    
    
  }
  
  if (save.graph){
    mixScores= transform.mix.score(parameterMatrix, mixScores)
    panel= matrix(0, length(types)*dim(mixScores)[2], 3)
    colnames(panel)=c("Score", "Percentage", "cell.type")
    panel[, 1]= as.vector(t(mixScores))
    panel[,2 ]= as.vector((percentage))
    panel[,3]= rep(types, each = dim(mixScores)[2])
    
    paneldf= as.data.frame(panel, stringsAsFactors = FALSE)
    paneldf$Score= as.numeric(paneldf$Score)
    paneldf$Percentage= as.numeric(paneldf$Percentage)
    ggplot(paneldf,mapping= aes(Score,Percentage))+
      geom_smooth(method=lm, se=FALSE, color="red")+
      geom_point(size=0.1)+facet_wrap(.~cell.type)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())+
      stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name	="rho")
    ggsave(filename="Transform Panel Panel.jpg")
  }
  
  parameterMatrix
}

transform.mix.score= function(parameterMatrix, mixScores){
  message("Transform score to percentage")
  transformedMix= mixScores
  for(i in 1:dim(mixScores)[1]){
    vec = transformedMix[i,]
    vec = vec + abs(min(vec))
    vec = vec^parameterMatrix[i,1]
    vec = vec/parameterMatrix[i,2]
    transformedMix[i,] = vec
  }
  
  
  transformedMix
}

create.ref.percentage= function(types, corr, controls){
  message("Creating percentage matrix for refmix (spillover matrix")
  percentageMat= diag(0.2, length(types), length(types))
  rownames(percentageMat)= types
  colnames(percentageMat)= 1:length(types)
  
  percentageMat[cbind(controls,1:length(types))]= 0.8
  percentageMat
  
}

