library(GSVA)
library(GSEABase)
library(pracma)
library(RColorBrewer)
library(pheatmap)
library(singscore)


read.types.dependencies = function(file.name) {
  con  <- file(file.path(file.name), open = "r")
  out <- list()
  i = 1
  while (length(oneLine <-
                readLines(con, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(oneLine, "\t"))
    out$types[i] = vec[1]
    n = max(which(!(vec == "")))
    if (n < 2)
      out$dep[i] = list(types = "")
    else
      out$dep[i] = list(types = vec[2:n])
    i = i + 1
  }
  close(con)
  out
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

create.single.cell.data.mat= function(ref, types){
  ntypes= length(types)
  singleCellData= matrix(NA, ntypes, length(rownames(ref)))
  rownames(singleCellData)= types
  colnames(singleCellData)<- rownames(ref)
  
  i=1
  for (type in types){
    location= match(type, ref$label.main)
    singleCellData[i,] <- assays(ref)$logcounts[,location]
    i=i+1
  }
  
  singleCellData
}

create.inSillicoMixtures= function(ref, gsc, types, singleCellData){
  ntypes= length(types)
  inSillicoMix= c()
  for (type in types){
    # Create mixtures
    countMix= matrix(NA,length(ref@NAMES), 1000)
    colnames(countMix)<-(1:1000)/50
    rownames(countMix)=rownames(ref)

    # Create in sillico mixtures
    for (i in 1:1000){
      randPerc= sample(1:20, ntypes-1, replace=T)
      sumPerc= sum(randPerc)
      perc= (randPerc/sumPerc)*((100-i/50)/100)
      perc= c(perc, i/(100*50))
      mix= rep(0, length(ref@NAMES))
      mix = colSums(perc * singleCellData) # My Gideon's line
      countMix[,i]= mix
    }

    # SimpleScore all signatures
    checkSED= SummarizedExperiment(assays=list(counts=countMix))
    inSillicoMix<- c(inSillicoMix, checkSED)
  }

  inSillicoMix
}

create.signatureScores= function(ref, gsc, inSillicoMix, types){
  scores= list()
  k=1
  for (type in types){
    checkSED= inSillicoMix[[k]]
    mixRankedData<-rankGenes(checkSED)
    gscSignatures= gsc[startsWith( names(gsc), type)] # get only the signature of the certain cell type we are interested in.
    mixScores= matrix(NA, length(names(gscSignatures)), 1000)
    colnames(mixScores)<-1:1000 # magic number, change or replace with global variable
    rownames(mixScores)= names(gscSignatures)
    j=1
    for (gs in gscSignatures){
      score <- simpleScore(mixRankedData, geneIds(gs), centerScore = TRUE)$TotalScore
      score.sing = matrix(NA,ncol(mixRankedData),length(gs))
      for (i in 1:length(gs)){
        score.sing[,i] <- simpleScore(mixRankedData, geneIds(gs), centerScore = TRUE)$TotalScore
      }
      score = t(score)
      rownames(score) = names(gs)
      colnames(score) = colnames(mixRankedData)
      score[is.na(score.sing)] = 0
      mixScores[j,] <- score
      j= j+1
      # df = data.frame(CellType=ref$label.main, Score = t(score.sing[1,]))
      # ggplot(df,aes(x=CellType,y=score,fill=CellType))+geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ggtitle(gs@setName)
      # ggsave( paste0(toString(gs@setName),".png"), path="~/graphs")
    }
    scores[[k]]= mixScores
    k= k+1
  }
  scores
}


choose.signatures = function(inSillicoMix, types, scores){
  ntypes= length(types)
  bestSig= matrix(NA, ntypes, 5)
  colnames(bestSig)<- 1:5
  rownames(bestSig)= types
  k=1

  for (type in types){
    # Choose the signatures with the best corrlation
    mixScores= scores[[k]]
    corrlation= cor(t(mixScores), (1:1000)/50, method= "spearman")
    top5= head( order(corrlation, decreasing = TRUE), 5)
    top5Sig= rownames(corrlation)[top5]
    bestSig[k,]=top5Sig
    k= k+1
  }
  bestSig
}

create.corrlationMixtures= function(percentage, ref, singleCellData){

  ntypes= length(types)

  # Create mixtures
  countMix= matrix(NA,length(ref@NAMES), length(colnames(percentage)))
  colnames(countMix)<-1:length(colnames(percentage))
  rownames(countMix)=rownames(ref)

  for (k in 1:length(colnames(percentage))){
    perc= unname(percentage[,k])
    mix= rep(0, length(ref@NAMES))
    mix = colSums(perc * singleCellData) # My Gideon's line
    countMix[,k]= mix
  }

  checkSED= SummarizedExperiment(assays=list(counts=countMix))
  checkSED
}

testing.correlation = function(mix, types,precentage, signatures, bestRank){
  mixRankedData<-rankGenes(mix)
  ntypes= length(types)

  # Score matrix
  mixScore= matrix(NA, ntypes, ntypes*20)
  colnames(mixScore)<-1:(ntypes*20)
  rownames(mixScore)=types

  for (i in 1:ntypes){
    sigs= bestRank[i,]
    nsigs= length(sigs)
    sumWeight= sum(1:nsigs)
    j= nsigs
    totScore= rep(0,ntypes*20)
    for (sig in sigs){
      gs= signatures[[sig]]
      totScore <-totScore+ (j/sumWeight)*(simpleScore(mixRankedData, geneIds(gs), centerScore = TRUE)$TotalScore)
      j= j-1
    }

    mixScore[i, ]= totScore
  }
  mixScore
}

create.random.percentage = function(types){
  ntypes= length(types)
  # Create percentage for in-Sillico mix.
  percentage= matrix(NA,length(types),20*ntypes)
  colnames(percentage)<-1:(20*ntypes)
  rownames(percentage)= types
  
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
  percentage
}

plot.linear.regression= function(working.dir, mixScores, percentage){
  
  dir.create(file.path(working.dir,"LinearRegressionGraphs"), showWarnings = FALSE)
  for (type in types){
    x= percentage[type,] 
    y=mixScores[type, ]
    df = data.frame(x, y)
    line = lm(y~x)
    pred_line = predict(line)
    LinerRegressionPath= file.path(working.dir,"LinearRegressionGraphs",
                                   paste0(type,"LinearRegression.jpg"))
    jpeg(file=LinerRegressionPath)
    plot.new()
    
    plot(y, x,main=paste0(type," Linear Regression"),
         xlab="percentage", ylab="Score", pch = 11, col = "blue")
    lines(pred_line,x,col="red")
    dev.off()
  }
}

create.parameter.matrix= function(working.dir, types, mixScores, percentage){
  ntypes= length(types)
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
    
    z = nlsLM(formula = y~a*x^b,start = list(a=1,b=1),data = df)
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
  
  parameterMatrix
}

transform.mix.score= function(parameterMatrix, mixScores){
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

create.ref.mix= function(types, singleCellData, control.type){
  ntypes= length(types)
  refmixData= matrix(NA, length(colnames(singleCellData)),ntypes)
  rownames(refmixData)= rownames(ref)
  colnames(refmixData)= types
  
  refMixPerc= matrix(0, ntypes, ntypes)
  rownames(refMixPerc)=types
  colnames(refMixPerc)<-1:ntypes
  
  
  for(i in 1:ntypes){
    contorlTypeName= types[control.type[i]]
    typeName= types[i]
    locationControl= match(contorlTypeName, rownames(singleCellData))
    controlData= assays(ref)$logcounts[,locationControl]*0.8
    locationType= match(typeName, rownames(singleCellData))
    typeData= assays(ref)$logcounts[,locationType]*0.2
    refmixData[,i] <- controlData+typeData
    refMixPerc[i, i]= 0.2
    refMixPerc[i, control.type[i]]=0.8
  }
  
  refmixExp= SummarizedExperiment(assays=list(counts=refmixData))
  refmixExp
}

score.ref.mix= function(types, bestSig, gscSignatures){
  ntypes= length(types)
  refMixScore= matrix(NA, ntypes, ntypes)
  colnames(refMixScore)<-1:ntypes
  rownames(refMixScore)=types
  
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
  refMixScore
}

create.refmix.controls= function(types, cormat){
  ntypes= length(types)
  control.type<-c()
  for(i in 1:ntypes){
    corrPercenatages = abs(cormat[,i])
    minInd = which(corrPercenatages==min(corrPercenatages))
    control.type<-c(control.type,minInd)
    contorlTypeName = colnames(cormat)[minInd]
    typeName = colnames(cormat)[i]
    print(c("base type:",typeName," min corr type:", contorlTypeName))
  }
  
  control.type
}