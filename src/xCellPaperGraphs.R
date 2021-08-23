library(here)
library(celldex)
library(corrplot)
library(xCell)
library(psych)
library(ggplot2)


working.dir= here()
spillover.dir= paste0(working.dir,'/Graphs/Spillover')
sdy.dir= paste0(working.dir,'/Graphs/sdy420graphs')
dir.create(paste0(working.dir, '/Graphs'), showWarnings = FALSE)
dir.create(spillover.dir, showWarnings = FALSE)
dir.create(sdy.dir, showWarnings = FALSE)

code.dir= paste0(here(), "/src")
train.dir = paste0(code.dir, "/xCell_train.R")
source(train.dir)
train.function.path= paste0(code.dir, "/xCell_train_functions.R")
source(train.function.path)
xCell.path= paste0(code.dir, "/xCell.R")
source(xCell.path)


sdy = readRDS(paste0(code.dir, '/sdy420.rds'))
genes.to.use= rownames(sdy$expr)

# ----------------- Blueprint main -----------------
xCell.my.data= xCellReferenceGenerate(ref= BlueprintEncodeData(), save.figures= TRUE, genes.to.use=genes.to.use, is.main=TRUE)

plot.new()
jpeg(file=paste0(spillover.dir, "/BlueprintMainSpillover.jpg"))
maxx= max(xCell.my.data$spill)
corrplot(xCell.my.data$spill, is.corr = T, col.lim=c(-maxx, maxx), tl.cex=0.75)
dev.off()


sed= SummarizedExperiment(assays=list(counts=(as.matrix(sdy$expr))))
scores= weighted.score(sed, xCell.my.data$signaturesRank, xCell.my.data$signatures, xCell.my.data$types)
transfomed.score= transform.mix.score(xCell.my.data$fitValues, scores)
after.spill= spillOver(transformedScores = transfomed.score, xCell.my.data$spill)


fcs = sdy$fcs#[rownames(after.spill),colnames(sdy$expr)]
fcs <- fcs[ , order(colnames(fcs))]
after.spill <- after.spill[ , order(colnames(after.spill))]
res1 = cor(t(after.spill),t(fcs))

plot.new()
jpeg(file=paste0(sdy.dir,"/BlueprintMainSDY420.jpg"))
corrplot(res1, is.corr = F)
dev.off()

eq= intersect(rownames(res1), colnames(res1))
test1= res1[rownames(res1) %in% eq, ]
test1= test1[,colnames(test1) %in% eq ]
test1= test1[order(rownames(test1)),]
test1= test1[,order(colnames(test1))]

plot.new()
jpeg(file=paste0(sdy.dir,"/BlueprintMainSDY420OnlySame.jpg"))
corrplot(test1, is.corr = F)
dev.off()


# ----------------- Real xCell -----------------

cell.types.use = intersect(colnames(xCell.data$spill$K),
                           rownames(sdy$fcs))
scores = xCellAnalysis(sdy$expr, rnaseq=F,
                       cell.types.use = cell.types.use)

fcs = sdy$fcs[rownames(scores),colnames(scores)]
res = corr.test(t(scores),t(fcs),adjust='none')
qplot(x=rownames(res$r),y=diag(res$r),
      fill=diag(res$p)<0.05,geom='col',
      main='SDY420 association with immunoprofiling',
      ylab='Pearson R', xlab="") + labs(fill = "p-value<0.05") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


res = cor(t(scores),t(fcs))

test= res[rownames(res) %in% eq, ] #use blueprint eq
test= test[,colnames(test) %in% eq ]
test= test[order(rownames(test)),]
test= test[,order(colnames(test))]


plot.new()
jpeg(file=paste0(sdy.dir,"/xCellSDY420OnlySame.jpg"))
corrplot(test, is.corr = F)
dev.off()


# ----------------- Blueprint fine -----------------
xCell.my.data= xCellReferenceGenerate(ref= BlueprintEncodeData(), save.figures= FALSE, genes.to.use=genes.to.use, is.main=FALSE)

plot.new()
jpeg(file=paste0(spillover.dir, "/BlueprintFineSpillover.jpg"))
maxx= max(xCell.my.data$spill)
corrplot(xCell.my.data$spill, is.corr = T, col.lim=c(-maxx, maxx), tl.cex=0.75)
dev.off()


sed= SummarizedExperiment(assays=list(counts=(as.matrix(sdy$expr))))
scores= weighted.score(sed, xCell.my.data$signaturesRank, xCell.my.data$signatures, xCell.my.data$types)
transfomed.score= transform.mix.score(xCell.my.data$fitValues, scores)
after.spill= spillOver(transformedScores = transfomed.score, xCell.my.data$spill)


fcs = sdy$fcs#[rownames(after.spill),colnames(sdy$expr)]
fcs <- fcs[ , order(colnames(fcs))]
after.spill <- after.spill[ , order(colnames(after.spill))]
res = cor(t(after.spill),t(fcs))

# result <- result[,order(-maxxc)]

# eq= intersect(rownames(res), colnames(res))
# test= result[rownames(result) %in% eq, ]
# test= test[,colnames(result) %in% eq ]

plot.new()
jpeg(file=paste0(sdy.dir,"/BlueprintFineSDY420.jpg"))
corrplot(res, is.corr = F)
dev.off()

# ----------------- MonacoImmuneData main -----------------
xCell.my.data= xCellReferenceGenerate(ref= MonacoImmuneData(), save.figures= FALSE, genes.to.use=genes.to.use, is.main=TRUE)

plot.new()
jpeg(file=paste0(spillover.dir, "/MonacoImmuneDataMainSpillover.jpg"))
maxx= max(xCell.my.data$spill)
corrplot(xCell.my.data$spill, is.corr = T, col.lim=c(-maxx, maxx), tl.cex=0.75)
dev.off()


sed= SummarizedExperiment(assays=list(counts=(as.matrix(sdy$expr))))
scores= weighted.score(sed, xCell.my.data$signaturesRank, xCell.my.data$signatures, xCell.my.data$types)
transfomed.score= transform.mix.score(xCell.my.data$fitValues, scores)
after.spill= spillOver(transformedScores = transfomed.score, xCell.my.data$spill)


fcs = sdy$fcs#[rownames(after.spill),colnames(sdy$expr)]
fcs <- fcs[ , order(colnames(fcs))]
after.spill <- after.spill[ , order(colnames(after.spill))]
res = cor(t(after.spill),t(fcs))

# result <- result[,order(-maxxc)]

# eq= intersect(rownames(res), colnames(res))
# test= result[rownames(result) %in% eq, ]
# test= test[,colnames(result) %in% eq ]

plot.new()
jpeg(file=paste0(sdy.dir,"/MonacoImmuneDataMainSDY420.jpg"))
corrplot(res, is.corr = F)
dev.off()



# ----------------- MonacoImmuneData fine -----------------
xCell.my.data= xCellReferenceGenerate(ref= MonacoImmuneData(), save.figures= FALSE, genes.to.use=genes.to.use, is.main=FALSE)

plot.new()
jpeg(file=paste0(spillover.dir, "/MonacoImmuneDataFineSpillover.jpg"))
maxx= max(xCell.my.data$spill)
spill= xCell.my.data$spill
spill= spill[!(row.names(spill) %in% c("Classical monocytes")),]
spill= spill[,!(colnames(spill) %in% c("Classical monocytes"))]
maxx= max(spill)
corrplot(spill, is.corr = T, col.lim=c(-maxx, maxx), tl.cex=0.75)
dev.off()


sed= SummarizedExperiment(assays=list(counts=(as.matrix(sdy$expr))))
scores= weighted.score(sed, xCell.my.data$signaturesRank, xCell.my.data$signatures, xCell.my.data$types)
transfomed.score= transform.mix.score(xCell.my.data$fitValues, scores)
after.spill= spillOver(transformedScores = transfomed.score, xCell.my.data$spill)


fcs = sdy$fcs#[rownames(after.spill),colnames(sdy$expr)]
fcs <- fcs[ , order(colnames(fcs))]
after.spill <- after.spill[ , order(colnames(after.spill))]
res = cor(t(after.spill),t(fcs))


# eq= intersect(rownames(res), colnames(res))
# test= res[rownames(res) %in% eq, ]
# test= test[,colnames(res) %in% eq ]

plot.new()
jpeg(file=paste0(sdy.dir,"/MonacoImmuneDataFineSDY420.jpg"))
corrplot(res, is.corr = F)
dev.off()
