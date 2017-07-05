#########################################################################################
##########         INTEGRATION EXAMPLE 6 - PLS and related methods             ##########
#########################################################################################

## By Sonia
## Dec-2016


setwd("~/Dropbox/cursos/IntegrationCourseCIPF/")



install.packages("mixOmics")
library(mixOmics)






# PLS example from mixOmics R package -------------------------------------

data(liver.toxicity)
head(liver.toxicity$gene[,1:10]); dim(liver.toxicity$gene)
head(liver.toxicity$clinic[,1:10]); dim(liver.toxicity$clinic)
table(liver.toxicity$treatment[,3:4])
liver.toxicity$treatment
# reordering rats according to dose and time
liver.toxicity$treatment = liver.toxicity$treatment[order(liver.toxicity$treatment$Dose.Group,
                                                          liver.toxicity$treatment$Time.Group),]

X <- liver.toxicity$gene[paste0("ID", liver.toxicity$treatment$Animal.Number),]
Y <- liver.toxicity$clinic[as.character(liver.toxicity$treatment$Animal.Number),]


# first choose the maximum number of components and then plot explained variance to decide the optimal ncomp
myresult <- pls(X, Y, ncomp = 10, mode = "regression")  
par(mfrow = c(1,2))
barplot(myresult$explained_variance$X, las = 2, main = "X")
barplot(myresult$explained_variance$Y, las = 2, main = "Y")

# final model
myresult <- pls(X, Y, ncomp = 3, mode = "regression", scale = TRUE) # scale = TRUE means autoscaling

dim(myresult$loadings$X)
dim(myresult$loadings$Y)

# scores
dim(myresult$variates$X)
dim(myresult$variates$Y)

plotIndiv(object = myresult)

plotIndiv(object = myresult, comp = 1:2, rep.space = "XY-variate", 
          ind.names = apply(liver.toxicity$treatment[,3:4], 1, paste, collapse = "-"),
          group = liver.toxicity$treatment$Treatment.Group, 
          col = rep(rainbow(16), each = 4), style = "ggplot2", ellipse = FALSE, ellipse.level = 0.95, centroid = FALSE) 

plotIndiv(object = myresult, comp = c(1,2),
          ind.names = apply(liver.toxicity$treatment[,3:4], 1, paste, collapse = "-"),
          group = liver.toxicity$treatment$Treatment.Group, 
          col = rep(rainbow(16), each = 4), style = "ggplot2", ellipse = FALSE, 
          ellipse.level = 0.95, centroid = FALSE,
          star = FALSE, legend = FALSE, abline = TRUE, alpha = 0)

# loadings
plotVar(object = myresult, comp = 1:2, cex = c(2, 5))


# VIP: Variable Importance in Projection
myVIP = vip(myresult)  
head(myVIP) 

mySelectedGenes = rownames(myVIP)[myVIP[,1] > 1]
mySelectedGenes = rownames(myVIP)[myVIP[,1] > 2]


# Validation

myperfMfold = perf(myresult, validation = "Mfold", folds = 10, progressBar = TRUE)
myperfLoo = perf(myresult, validation = "loo", progressBar = TRUE)

myperfMfold$MSEP  # Mean Square Error Prediction for each Y variable
myperfMfold$R2
myperfMfold$Q2
myperfMfold$Q2.total
myperfMfold$RSS
myperfMfold$PRESS

myperfLoo$MSEP

plot(x = myperfLoo,
     criterion = "MSEP",
     xlab = "number of components",
     ylab = NULL,
     LimQ2 = 0.0975,
     LimQ2.col = "darkgrey")





# Sparse PLS from mixOmics package ----------------------------------------

mySPresult <- spls(X, Y, ncomp = 3, mode = 'regression', keepX = c(50, 50, 50), keepY = c(10, 10, 10), scale = TRUE)

# LOO CV
liver.spls.loo <- perf(mySPresult, ncomp = 3, mode = 'regression', keepX = c(50, 50, 50), validation = 'loo')      

palette(rainbow(10))
par(mfrow = c(2, 5))
for(i in 1:10){
  spls.rmsep <- sqrt(liver.spls.loo$MSEP[i, ])
  pls.rmsep <- sqrt(myperfLoo$MSEP[i, ])
  matplot(cbind(spls.rmsep, pls.rmsep), type = 'l', col = i, ylab = 'RMSEP', lwd = 2,
          xlab = 'dim', lty = c(1, 2), axes = FALSE)
  axis(1, 1:3, labels = 1:3)
  axis(2)
  title(main = paste(rownames(liver.spls.loo$MSEP)[i]))
}
palette("default")


# Alternative plot
plot(liver.spls.loo, criterion = 'RMSEP', type = 'l', layout = c(3, 4))


# Samples plot (scores)
plotIndiv(object = mySPresult, comp = c(1,2),
          ind.names = apply(liver.toxicity$treatment[,3:4], 1, paste, collapse = "-"),
          group = liver.toxicity$treatment$Treatment.Group, 
          col = rep(rainbow(16), each = 4), style = "ggplot2", ellipse = FALSE, 
          ellipse.level = 0.95, centroid = FALSE,
          star = FALSE, legend = FALSE, abline = TRUE, alpha = 0)


# Loadings plot
plotVar(mySPresult, comp = 1:2, cex = c(3, 5))



# Network
color.edge <- colorRampPalette(c("red4", "white", "darkgreen"))
pdf("HandsOn/IE6_PLS/network_sPLS.pdf")
network(mySPresult, comp = 1:2, shape.node = c("rectangle", "rectangle"),
        color.node = c("white", "lightblue"), color.edge = color.edge(10))
dev.off()


# Heatmap 
X11()
cim(mySPresult, comp = 1:2, xlab = "", ylab = "", margins = c(5, 6))



# PLS-DA from mixOmics ------------------------------------------------------------------

data(srbct)
names(srbct)

## Data
X <- srbct$gene

dim(X)
X[1:10,1:5]

Y <- srbct$class 
summary(Y)

## Number of components for PLS-DA
set.seed(32)
srbct.plsda.perf <- plsda(X, Y, ncomp = 10, scale = TRUE)

perf.plsda <- perf(srbct.plsda.perf, validation = 'Mfold', folds = 5,
                   progressBar = TRUE, nrepeat = 10) 
# better with folds = 10 and higher nrepeat value, but more time-consuming

plot(perf.plsda, overlay = 'measure', sd=TRUE)


srbct.plsda <- plsda(X, Y, ncomp = 3, scale = TRUE)

plotIndiv(srbct.plsda , comp = c(1,2),
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'SRBCT, PLSDA comp 1 - 2')

plotVar(srbct.plsda, comp = c(1,2), cex = 2) 




### Sparse PLS-DA


## Tuning the model

set.seed(32)
# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 100, 10))

tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, validation = 'Mfold', folds = 5, 
                                 progressBar = TRUE, dist = 'mahalanobis.dist',
                                 test.keepX = list.keepX, nrepeat = 10) # nrepeat 50-100

head(tune.splsda.srbct$error.rate)
tune.splsda.srbct$choice.keepX

pdf("pls-da.pdf")
plot(tune.splsda.srbct, optimal = TRUE, sd = TRUE)
dev.off()


## Final model

# optimal number of variables to select on 3 comps:
select.keepX = c(10,50,90) #from tuning step

splsda.srbct <- splsda(X, Y, ncomp = 3, keepX = select.keepX, scale = TRUE) 

plotIndiv(splsda.srbct, comp = c(1,2),
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'SRBCT, sPLSDA comp 1 - 2')

plotVar(splsda.srbct, comp = c(1,2), 
        var.names = list(substr(srbct$gene.name[, 2], 1, 5)), cex = 3)

plotLoadings(splsda.srbct, comp = 2, method = 'mean', contrib = 'max')



## Evaluating sPLS-DA

set.seed(32)  
perf.splsda <- perf(splsda.srbct, folds = 5, validation = "Mfold", 
                    dist = "max.dist", progressBar = TRUE, nrepeat = 10)
# perf.srbct  # lists the different outputs
perf.splsda$error.rate

perf.splsda$error.rate.class

# The same variables are selected across the different CV folds
head(perf.splsda$features$stable[[1]])
plot(perf.splsda$features$stable[[1]], type = 'h', 
     xlab = 'variables selected across CV folds', ylab = 'Stability frequency', 
     main='Feature stability for comp = 1')

# Selected variables for component 1 and their loadings
selectVar(splsda.srbct, comp = 1)$value  


## Prediction
Xnew = X[1,,drop = FALSE] + rnorm(ncol(X), 0, 0.3) 
Xnew[,1:3, drop = FALSE] # EWS

myprediction = predict(splsda.srbct, Xnew, dist = "mahalanobis.dist")
myprediction$class




# DIABLO (mixOmics): PLS-DA for multiomics --------------------------------

data("breast.TCGA")
names(breast.TCGA)
names(breast.TCGA$data.train)
sapply(breast.TCGA$data.train[1:3], dim) # train = 150
names(breast.TCGA$data.test)
sapply(breast.TCGA$data.test[1:2], dim) # test = 70


# this is the X data as a list of mRNA, miRNA and proteins
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
            protein = breast.TCGA$data.train$protein)

# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
diag(design) =  0
design

# set number of component per data set
ncomp = c(2)

# set number of variables to select, per component and per data set (this is set arbitrarily)
list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2), protein = rep(10, 2))


TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype, scale = TRUE,
                                 ncomp = ncomp, keepX = list.keepX, design = design)
TCGA.block.splsda

plotIndiv(TCGA.block.splsda, ind.names = FALSE, legend = TRUE)

# illustrates coefficient weights in each block
pdf("diablo.pdf")
plotLoadings(TCGA.block.splsda, ncomp = 1, contrib = 'max')
dev.off()

plotVar(TCGA.block.splsda, style = 'graphics', legend = TRUE)

plotDiablo(TCGA.block.splsda, ncomp = 1)


# Circos plots
circosPlot(TCGA.block.splsda, comp = 1:2, cutoff = 0.7, size.variables = 0.5)


X11()
cimDiablo(TCGA.block.splsda, size.legend = 1)



# prediction
myDiabloPrediction = predict(TCGA.block.splsda, 
                             list(mrna = breast.TCGA$data.test$mrna, 
                                  mirna = breast.TCGA$data.test$mirna))
                             
                                  
head(myDiabloPrediction$class$mahalanobis.dist)
head(data.frame(myDiabloPrediction$MajorityVote$mahalanobis.dist, breast.TCGA$data.test$subtype))

# confusion matrix to compare the real subtypes with the predicted subtypes for a 2 component model
table(myDiabloPrediction$MajorityVote$mahalanobis.dist[,2],
      breast.TCGA$data.test$subtype)
