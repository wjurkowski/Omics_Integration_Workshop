#########################################################################################
##########       HANDS ON: INTEGRATION EXAMPLE 6 - PLS and related methods     ##########
#########################################################################################

## By Sonia
## Dec-2016


setwd("~/Dropbox/cursos/IntegrationCourseCIPF/")


library(mixOmics)
library(STATegRa)




# Data -------------------------------------

myexpression = read.delim("HandsOn/IE6_PLS/EcoliData/expressionEcoli.txt", sep = " ", dec = ",")
head(myexpression); dim(myexpression)   # 85 31

mymetabolites = read.delim("HandsOn/IE6_PLS/EcoliData/metabolitesEcoli.txt", sep = " ", dec = ",")
head(mymetabolites); dim(mymetabolites)  ## missing values!   261 31

myproteins = read.delim("HandsOn/IE6_PLS/EcoliData/proteinsEcoli.txt", sep = " ", dec = ",")
head(myproteins); dim (myproteins)  ## missing values!  67 31

mydesign = read.delim("HandsOn/IE6_PLS/EcoliData/EcoliDesign.txt", as.is = TRUE, row.names = 1)
head(mydesign); dim(mydesign)
table(mydesign)

myexpression = myexpression[, rownames(mydesign)]
mymetabolites = mymetabolites[, rownames(mydesign)]
myproteins = myproteins[, rownames(mydesign)]



# Exploring data ----------------------------------------------------------

summary(myexpression)
summary(myproteins)

# Are there proteins with missing values in all the samples?
numNA = apply(myproteins, 1, function (x) sum(is.na(x)))
summary(numNA)

# Removing proteins with missing values in 25 samples or more (80%)
myproteins = myproteins[names(which(numNA < 25)),]  # 57 proteins

# Are there samples with missing values for all the proteins?
numNA = apply(myproteins, 2, function (x) sum(is.na(x)))
summary(numNA)

# Remove samples with missing values for all proteins from expression and protein data
myproteins = myproteins[, names(which(numNA < 57))]  # 30 samples
myexpression = myexpression[, colnames(myproteins)]
mydesign = mydesign[colnames(myproteins),]

# Boxplots
boxplot(myexpression, las = 2)
boxplot(myproteins, las = 2)

min(myexpression); min(myproteins[!is.na(myproteins)])



# PCA
X <- t(myexpression)
Y <- t(myproteins)

boxplot(X)
boxplot(scale(X, center = TRUE, scale = TRUE))

pcaX <- pca(X, ncomp = 30, center = TRUE, scale = TRUE)
plot(pcaX, ncomp = 30)

pcaX <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)

mycolors = c(rainbow(24), rep(1,6))

plotIndiv(object = pcaX, comp = 1:2, ind.names = apply(mydesign, 1, paste, collapse = "-"),
          group = apply(mydesign, 1, paste, collapse = "-"),
          col = mycolors) 


boxplot(Y)
boxplot(scale(Y, center = TRUE, scale = TRUE))
pcaY <- pca(Y, ncomp = 3, center = TRUE, scale = TRUE) 

plotIndiv(object = pcaY, comp = 1:2, ind.names = apply(mydesign, 1, paste, collapse = "-"),
          group = apply(mydesign, 1, paste, collapse = "-"),
          col = mycolors) 



# PLS ---------------------------------------------------------------------

X = scale(X, center = TRUE, scale = FALSE)
Y = scale(Y, center = TRUE, scale = FALSE)

# first choose the maximum number of components and then plot explained variance to decide the optimal ncomp
myresult <- pls(X, Y, ncomp = 30, mode = "regression", scale= FALSE)  
par(mfrow = c(1,2))
barplot(myresult$explained_variance$X, las = 2, main = "X")
barplot(myresult$explained_variance$Y, las = 2, main = "Y")

myresult <- pls(X, Y, ncomp = 2, mode = "regression", scale = FALSE)

# scores
plotIndiv(object = myresult, comp = 1:2, rep.space = "XY-variate", 
          ind.names = apply(mydesign, 1, paste, collapse = "-"),
          group = apply(mydesign, 1, paste, collapse = "-"), 
          col = mycolors) 


# loadings
plotVar(object = myresult, comp = 1:2, cex = c(4, 4))




# Sparse PLS --------------------------------------------------------------

mySPresult <- spls(X, Y, ncomp = 3, mode = 'regression', scale = FALSE,
                   keepX = c(20, 20, 20), keepY = c(20, 20, 20))

# LOO CV --> Validation is not possible with missing values


plotIndiv(object = mySPresult, comp = c(2,3),
          ind.names = apply(mydesign, 1, paste, collapse = "-"),
          group = apply(mydesign, 1, paste, collapse = "-"), 
          col = mycolors)

plotVar(object = mySPresult, comp = 1:2, cex = c(4, 4))







# N-PLS -------------------------------------------------------------------

library(Multiway)
library(plsdepot)
library(sfsmisc)

## Loading gene data
load("HandsOn/IE6_PLS/ToxicogenomicsData/Toxicogenomics.RData", verbose = TRUE)


## Preparing metabolomic data

metabolites = read.delim("HandsOn/IE6_PLS/ToxicogenomicsData/toxico_metabolomics.txt", 
                         as.is = TRUE, header = TRUE, row.names = 1, dec = ",")
dim(metabolites)
metabolites[1:5,1:5]
colnames(metabolites)

tratam = colnames(metabolites)
tratam = sapply(tratam, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1:2])
dim(tratam); head(t(tratam))
tratam = apply(tratam, 2, paste, collapse = "_")
tratam

# Averaging replicates
metaAvg = t(apply(metabolites, 1, tapply, INDEX = tratam, mean))
head(metaAvg)

dosis = c("UT", "CO", "LO", "ME", "HI")
tiempos = c(6,24,48)

# Converting matrix into a 3d array
Metab = array(0, dim = c(length(dosis), nrow(metaAvg), length(tiempos)))

for (i in 1:length(dosis)) {
  for (j in 1:length(tiempos)) {
    Metab[i,,j] = metaAvg[,paste(dosis[i], tiempos[j], sep = "_")]
  }
}

dimnames(Metab)[[1]] = c("UT", "CO", "LO", "ME", "HI")
dimnames(Metab)[[2]] = rownames(metaAvg)
dimnames(Metab)[[3]] = c("6h", "24h", "48h")


## Centering 
GeneExprCent = GeneExpr
MetabCent = Metab
for (i in 1:dim(GeneExpr)[3]) {
  GeneExprCent[,,i] = scale(GeneExpr[,,i], center = TRUE, scale = FALSE)
  MetabCent[,,i] = scale(Metab[,,i], center = TRUE, scale = FALSE)
}


## NPLS model
myNpls = npls(XN = GeneExprCent, YN = MetabCent, F = 2)

pdf("HandsOn/IE6_PLS/ToxicogenomicsData/NPLSplots_task.pdf", width = 10, height = 10)
plotSpace(myNpls, what = "X", PCs = 1:2, cutoff = 20)
plotSpace(myNpls, what = "Y", PCs = 1:2, cutoff = 20)
dev.off()

myloadingsX = myNpls$FactorsX$Mode2
myloadingsY = myNpls$FactorsY$Mode2

plot(density(myloadingsX[,2]))
plot(density(myloadingsY[,2]))

summary(myloadingsX[,2])
summary(myloadingsY[,2])






# O2PLS -------------------------------------------------------------------

data(breast.TCGA)

names(breast.TCGA$data.train)
sapply(breast.TCGA$data.train, dim)
subtypes = data.frame("classname" = breast.TCGA$data.train$subtype)
rownames(subtypes) = rownames(breast.TCGA$data.train$mrna)


# Block1 - gene expression data
B1 <- createOmicsExpressionSet(Data = t(breast.TCGA$data.train$mrna), pData = subtypes,
                               pDataDescr=c("classname"))

# Block2 - miRNA expression data
B2 <- createOmicsExpressionSet(Data = t(breast.TCGA$data.train$mirna), pData = subtypes,
                               pDataDescr=c("classname"))

# Block3 - Protein expression data
B3 <- createOmicsExpressionSet(Data = t(breast.TCGA$data.train$protein), pData = subtypes,
                               pDataDescr=c("classname"))



## Model selection 1
ms <- modelSelection(Input=list(B1,B2), Rmax=4, fac.sel="single%",
                     varthreshold=0.03)
ms
# $common
# [1] 3
# 
# $dist
# [1] 3 5        

## Model selection 2
ms <- modelSelection(Input=list(B1,B3), Rmax=4, fac.sel="single%",
                     varthreshold=0.03)
ms
# $common
# [1] 3
# 
# $dist
# [1] 3 3



## O2-PLS model1
o2plsRes1 <- omicsCompAnalysis(Input=list(B1, B2),Names=c("expr", "mirna"),
                               method="O2PLS", Rcommon=3, Rspecific=c(3, 5),
                               center=TRUE, scale=TRUE, weight=TRUE)

## O2-PLS model2
o2plsRes2 <- omicsCompAnalysis(Input=list(B1, B3),Names=c("expr", "prot"),
                               method="O2PLS", Rcommon=3, Rspecific=c(3, 3),
                               center=TRUE, scale=TRUE, weight=TRUE)


# Scatterplot of scores variables associated to common components MODEL 1
# Associated to first block
p1 <- plotRes(object=o2plsRes1, comps=c(1, 2), what="scores", type="common",
              combined=FALSE, block="expr", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to second block
p2 <- plotRes(object=o2plsRes1, comps=c(1, 2), what="scores", type="common",
              combined=FALSE, block="mirna", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Combine both plots
legend <- grid_arrange_shared_legend(p1,p2)   


# Scatterplot of scores variables associated to common components MODEL 2
# Associated to first block
p1 <- plotRes(object=o2plsRes2, comps=c(1, 2), what="scores", type="common",
              combined=FALSE, block="expr", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to second block
p2 <- plotRes(object=o2plsRes2, comps=c(1, 2), what="scores", type="common",
              combined=FALSE, block="prot", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Combine both plots
legend <- grid_arrange_shared_legend(p1,p2)  




# O2-PLS scores scatterplot associated to individual components MODEL 1
# Associated to first block
p1 <- plotRes(object=o2plsRes1, comps=c(1, 2), what="scores", type="individual",
              combined=FALSE, block="expr", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to second block
p2 <- plotRes(object=o2plsRes1, comps=c(1, 2), what="scores", type="individual",
              combined=FALSE, block="mirna", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Combine plots
legend <- grid_arrange_shared_legend(p1,p2) 



# O2-PLS scores scatterplot associated to individual components MODEL 2
# Associated to first block
p1 <- plotRes(object=o2plsRes2, comps=c(1, 2), what="scores", type="individual",
              combined=FALSE, block="expr", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Associated to second block
p2 <- plotRes(object=o2plsRes2, comps=c(1, 2), what="scores", type="individual",
              combined=FALSE, block="prot", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
# Combine plots
legend <- grid_arrange_shared_legend(p1,p2) 



# O2-PLS scores combined plot for individual components
plotRes(object=o2plsRes1, comps=c(1, 1), what="scores", type="individual",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)

plotRes(object=o2plsRes2, comps=c(1, 1), what="scores", type="individual",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)


# O2-PLS combined plot of scores for common and individual components MODEL 1
p1 <- plotRes(object=o2plsRes1, comps=c(1, 1), what="scores", type="both",
              combined=TRUE, block="expr", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
p2 <- plotRes(object=o2plsRes1, comps=c(1, 1), what="scores", type="both",
              combined=TRUE, block="mirna", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
legend <- grid_arrange_shared_legend(p1,p2) 


# O2-PLS combined plot of scores for common and individual components MODEL 2
p1 <- plotRes(object=o2plsRes2, comps=c(1, 1), what="scores", type="both",
              combined=TRUE, block="expr", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
p2 <- plotRes(object=o2plsRes2, comps=c(1, 1), what="scores", type="both",
              combined=TRUE, block="prot", color="classname", shape=NULL,
              labels=NULL, background=TRUE, palette=NULL, pointSize=4,
              labelSize=NULL, axisSize=NULL, titleSize=NULL)
legend <- grid_arrange_shared_legend(p1,p2) 




# Loadings plot for common components --> Combined plot
plotRes(object=o2plsRes1, comps=c(1, 2), what="loadings", type="common",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
plotRes(object=o2plsRes2, comps=c(1, 2), what="loadings", type="common",
        combined=TRUE, block="", color="classname", shape=NULL,
        labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)




# Biplot common part. O2PLS MODEL 1
p1 <- biplotRes(object=o2plsRes1, type="common", comps=c(1, 2),
                block="expr", title=NULL, colorCol="classname",
                sizeValues=c(2, 4), shapeValues=c(17, 0),
                background=TRUE, pointSize=4, labelSize=NULL,
                axisSize=NULL, titleSize=NULL)
p2 <- biplotRes(object=o2plsRes1, type="common", comps=c(1, 2),
                block="mirna", title=NULL, colorCol="classname",
                sizeValues=c(2, 4), shapeValues=c(17, 0),
                background=TRUE, pointSize=4, labelSize=NULL,
                axisSize=NULL, titleSize=NULL)
legend <- grid_arrange_shared_legend(p1,p2) 


# Biplot common part. O2PLS MODEL 2
p1 <- biplotRes(object=o2plsRes2, type="common", comps=c(1, 2),
                block="expr", title=NULL, colorCol="classname",
                sizeValues=c(2, 4), shapeValues=c(17, 0),
                background=TRUE, pointSize=4, labelSize=NULL,
                axisSize=NULL, titleSize=NULL)
p2 <- biplotRes(object=o2plsRes2, type="common", comps=c(1, 2),
                block="prot", title=NULL, colorCol="classname",
                sizeValues=c(2, 4), shapeValues=c(17, 0),
                background=TRUE, pointSize=4, labelSize=NULL,
                axisSize=NULL, titleSize=NULL)
legend <- grid_arrange_shared_legend(p1,p2) 






# DIABLO ------------------------------------------------------------------

data("STATegRa_S3")

# this is the X data as a list of mRNA, miRNA and proteins
data = list(expr = t(Block1.PCA), miRNA = t(Block2.PCA))
sapply(data, dim)

# scaling data
data = lapply(data, scale, center = TRUE, scale = FALSE)

# 120 samples for train data
set.seed(123)
trainSamples = sample(1:169, 120, replace = FALSE)
testSamples = setdiff(1:169, trainSamples)

data.train = lapply(data, function(x) x[trainSamples,])
data.test = lapply(data, function(x) x[testSamples,])


# cancer subtypes
mysubtypes = ed.PCA$classname  # Y

# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
diag(design) =  0
design


# set number of components 
mydiablo = block.splsda(X = data.train, Y = mysubtypes[trainSamples], scale = FALSE, 
                        ncomp = 4, design = design)
perf.numcomp <- perf(mydiablo, validation = 'Mfold', folds = 10, progressBar = TRUE, nrepeat = 50) 
plot(perf.numcomp, overlay = 'measure', sd=TRUE)

ncomp = 3


# set number of variables to select, per component and per data set (this is set arbitrarily)
list.keepX = list(expr = seq(5, 50, 5), miRNA = seq(5, 50, 5))


# Tuning DIABLO model
tuningDiablo = tune.block.splsda(X = data.train, Y = mysubtypes[trainSamples], scale = FALSE, ncomp = ncomp, 
                                 test.keepX = list.keepX, design = design)
tuningDiablo$choice.keepX
# $expr
# [1]  5 20  5
# 
# $miRNA
# [1]  5 10 25


# DIABLO model
mydiablo = block.splsda(X = data.train, Y = mysubtypes[trainSamples], scale = FALSE, ncomp = ncomp, 
                        keepX = tuningDiablo$choice.keepX, design = design)

plotIndiv(mydiablo, ind.names = FALSE, legend = TRUE)

# illustrates coefficient weights in each block
plotLoadings(mydiablo, ncomp = 1, contrib = 'max')
plotLoadings(mydiablo, ncomp = 2, contrib = 'max')
plotLoadings(mydiablo, ncomp = 3, contrib = 'max')

plotVar(mydiablo, style = 'graphics', legend = TRUE)

plotDiablo(mydiablo, ncomp = 1)
plotDiablo(mydiablo, ncomp = 3)

# Circos plots
circosPlot(mydiablo, comp = 1:2, cutoff = 0.8, size.variables = 0.5)


X11()
cimDiablo(mydiablo, size.legend = 1)



# prediction
myDiabloPrediction = predict(mydiablo, data.test)

head(myDiabloPrediction$class$mahalanobis.dist)
comparison = data.frame(myDiabloPrediction$MajorityVote$mahalanobis.dist, mysubtypes[testSamples])

table(comparison[,c(4,1)])
table(comparison[,c(4,2)])
table(comparison[,c(4,3)])






