
# This code is derived from "AllSQLitePCA.R" in this folder:
#  C:\Users\mooredf\Box\gp9\NJacts\investigators\MariaGennaro\targetedMSexpand\progs


# from AllSQLite.R:

#setwd("C:\\Users\\mooredf\\Box\\gp9\\NJacts\\investigators\\MariaGennaro\\targetedMSexpand\\dataUse")
#save(healthyData, healthyNewData, kawasakiData, mildAsympData, miscData,
#     pneumoniaData, geneProt, file="AllSQLite.RData")

setwd("C:\\Users\\mooredf\\Box\\gp9\\NJacts\\investigators\\MariaGennaro\\targetedMSexpand\\dataUse")
load(file="AllSQLite.RData")

ncol(healthyData)
#  22
ncol(healthyNewData)
#  6
ncol(kawasakiData)
#  13
ncol(mildAsympData)
# 10
ncol(miscData)
# 19
ncol(pneumoniaData)
# 17


#x

#allDataX <- data.frame(healthyData, healthyNewData, kawasakiData, mildAsympData,
#                      miscData, pneumoniaData)
allDataX <- data.frame(healthyData, healthyNewData, miscData, kawasakiData,
                       pneumoniaData, mildAsympData)
rownames(allDataX) <- paste(geneProt$Genes, geneProt$Protein.Group, sep="_")
dim(allDataX)
#[1] 1391   87


# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #




# drop the 28 "healthy" patients ("healthy" and "healthy new")
allDataXb <- allDataX[,-c(1:28)]   # same thing, but with row names as genes
#allDataXb <- data.frame(kawasakiData, mildAsympData,
#                       miscData, pneumoniaData)
dim(allDataXb)
# [1] 1391   59
col.sample.b <- c( rep("purple", 19), rep("blue", 13),
                 rep("brown", 17), rep("green", 10))
pch.sample.b <- c(rep(3,19), rep(1,13), rep(4,17), rep(2,10))
let.sample.b <- c(rep("M", 19), rep("K", 13),  rep("P", 17), rep("A", 10))
length(col.sample.b)
#[1] 59
table(col.sample.b, pch.sample.b)   # not diagonal, but checks out
table(col.sample.b, let.sample.b)
let.sample.f <- factor(let.sample.b, levels=c("M", "K", "P", "A"))
# remove proteins with too many missing values;
# do this separately for each of the four subgroups
# use the "completeX" function
round(0.3*c(19, 13, 17, 10))
# [1] 6 4 5 3
# allow these numbers of missing (30% rounded):
# M: 6 / 19
# K: 4 / 13
# P: 5 / 17
# A: 3 / 10

k.ind <- completeX(datax=kawasakiData, X=4)
table(k.ind)
#FALSE  TRUE
#1087   304
kawasakiU <- kawasakiData[k.ind,]
a.ind <- completeX(datax=mildAsympData, X=3)
table(a.ind)
#FALSE  TRUE
#1034   357
m.ind <- completeX(datax=miscData, X=6)
table(m.ind)
#FALSE  TRUE
#1048   343
p.ind <- completeX(datax=pneumoniaData, X=5)
table(p.ind)
#FALSE  TRUE
#976   415

# use this if setting up the data as analyzed in the paper:
#all.ind <- {k.ind & a.ind & m.ind & p.ind}  # not excessive missing in any

# get everything:
all.ind <- rep(TRUE, nrow(miscData))


##allData.b <- allDataXb[complete.cases(allDataXb),]
##dim(allData.b)
### [1] 185  59
## We started with 1391 genes, reduced to 190 genes with complete data

allData.b <- allDataXb[all.ind,]
dim(allData.b)
#[1] 284  59   # this way (June 10, 2024) we have 284 proteins, up from 185
allGene.b <- rownames(allData.b)

# remove immunoglobulins
tempID <- grep("^IG", allGene.b)
allGene.b[tempID]
# per Marila Gennaro (June 20, 2024) the following are not immunoglobulins:
# "IGF2R_P11717"    "IGFALS_P35858;P35858-2"      "IGF2_P01344-3"
# These are numbers 12, 13, and 18
tempIDMod <- tempID[-c(12, 13, 18)]
allGene.b[tempIDMod]
# successful!
allData.c <- allData.b[-tempIDMod,]
allGene.c <- allGene.b[-tempIDMod]
dim(allData.c)
# [1] 269  59
allGene.b[tempIDMod]  # immunoglobulins that were removed

dim(allData.b)  # before removing immunoglobulins
# [1] 1391   59
dim(allData.c)  # after removing immunoglobulins (38 of them)
# [1] 1353   59

#log2out.c <- log2((allData.c)) # wrong!!
log2out.c <- allData.c
dim(log2out.c)
# [1]   269 59
sum(is.na(log2out.c))
# [1] 337
geneProteinNames <- row.names(log2out.c)
temp <- unlist(strsplit(geneProteinNames, "_"))
geneProteinNames <- matrix(temp, ncol=2, byrow=TRUE)

geneNames <- matrix(temp, ncol=2, byrow=TRUE)[,1] # drop protein identifiers
# # # # # # # # # # # # # # # # # #
#


grpInd <- let.sample.f  # factor with levels in new order
table(grpInd)
#  M  K  P  A
# 19 13 17 10

# M: MIS-C
# K: Kawasaki
# P: Pneumonia
# A: Asymptomatic/mild

n.row <- nrow(log2out.c)

names(log2out.c) <- paste(grpInd, ".", 1:ncol(log2out.c), sep="")
head(log2out.c)

#setwd("C:\\Users\\mooredf\\Box\\gp10\\Res\\MariaGennaroMISC\\manuscript\\RCode\\data")
write.csv(log2out.c, file="mis-c_dataSet2.csv", row.names=TRUE)
# data set complete

# # # # # # # # # # # # # # # # # #
# do one vs all others comparisons
# # # # # # # # # # # # # # # # # #

X="M"

#grpselect <- {{let.sample.f == X}} # TRUE for X="M"
log2out.d <- log2out.c[complete.cases(log2out.c),]
dim(log2out.d)
#[1] 176  59
#log2protData.X.Y <- log2out.c
log2protData.X.Y <- log2out.d  # complete cases only
grpIndNames <- let.sample.f
grpInd <- as.numeric(grpIndNames == X)
table(grpInd)
# grpInd
# 0  1
# 40 19
resultOut <- protOut(log2protData=log2protData.X.Y, grpInd=grpInd)
names(resultOut) <- paste(names(resultOut), X, "vs rest", sep=".")

pvalCrit <- 0.05 / nrow(resultOut)
ordProt <- order(resultOut[,3])
resultOut.ord <- resultOut[ordProt,]
geneProtNames <- rownames(resultOut)
#tempGene <- matrix(unlist(strsplit(geneProtNames, "_")), ncol=2, byrow=TRUE)
geneName <- matrix(unlist(strsplit(geneProtNames, "_")), ncol=2, byrow=TRUE)[,1]
if (TRUE) {
  pvalCrit <- mean(c(2.043919e-04 , 2.314744e-04))
  volcanoPlot(log2ratio=resultOut[,2], pvalue=resultOut[,3],
              proteinName=rownames(resultOut), pvalCrit=pvalCrit, pvalCrit2=pvalCrit,
              mainTitle=paste(X,"vs all others"))
  windows()
  par(mar=c(5,5.5,4,2))
  volcanoPlot(log2ratio=resultOut[,2], pvalue=resultOut[,3],
              proteinName="", pvalCrit=pvalCrit, pvalCrit2=pvalCrit,
              mainTitle=paste(X,"vs all others"))
  identify(resultOut[,2], -log10(resultOut[,3]), labels=geneName, cex=0.7)
  box()
}
# extract gene name
geneList <- matrix(unlist(strsplit(x=rownames(resultOut), split ="_")),
               ncol=2, byrow=TRUE)[,1]
identify(-log10(resultOut[,3]) ~ resultOut[,2] , labels=geneList, cex=0.6)
# #
# gene list for all with Holm p-values < 0.05
genesAll <- rownames(resultOut)
indHolm <- {resultOut[,4] < 0.05}
genes.X.Y <- genesAll[indHolm]
#X.Y.grp <- rep(paste(X,Y,sep="."), length(genes.X.Y) )
#genes.list <- data.frame(X.Y.grp, genes.X.Y)
#result <- data.frame(result, resultOut)
resultReturn <- result[,-1]  # drop first column
geneHolmList <- rbind(geneHolmList, genes.list)

resultAll_M_vs_others <- data.frame(dataOutAllSQLiteA, resultOut)
head(resultAll_M_vs_others)
resultAll_M_vs_others_ord <- resultAll_M_vs_others[pvalOrd,]
head(resultAll_M_vs_others_ord)
write.csv(resultAll_M_vs_others_ord, file="resultAll_M_vs_others_ord.csv",
          row.names=FALSE)

pvalOrdM <- order(resultOut[,3])
resultAll_M_vs_others_ordM <- resultAll_M_vs_others[pvalOrdM,]

head(resultAll_M_vs_others_ordM)
# remove those with missing values ("nm" means "not missing")
resultAll_M_vs_others_ordM_nm <- resultAll_M_vs_others_ordM[resultAll_M_vs_others_ordM$sumMiss == 0,]
dim(resultAll_M_vs_others_ordM_nm)


tempList <- resultOut[pvalOrdM,][1:30,1:4]
ordIntercept <- order(tempList[,1], decreasing=TRUE)
round(tempList[ordIntercept,][1:20,], digits=4)
# # # # # # # # # # # # # # # # # # # #
# SVM prediction
# # # # # # # # # # # # # # # # # # # #
#head(resultAll_M_vs_others_ord)
#
head(geneList[pvalOrd], n=20)   # proteins with smallest p-values

#genes3 <- rownames(resultAll_M_vs_others_ordM)[1:3]
geneListOrd <- geneList[pvalOrdM]
genes3 <- geneListOrd[1:3]
genes4 <- geneListOrd[1:4]
genes5 <- geneListOrd[1:5]
genes6 <- geneListOrd[1:6]
genes20 <- geneListOrd[1:20]
genes7 <- geneListOrd[1:7]

# gene list from Gennaro, with positive coefficients:
#genes3alt <- c("VWF_P04275", "FCGBP_Q9Y6R7", "SERPINA3_P01011")
#grpselect <- {rownames(log2out.c) %in% genes3alt}

# just use gene names
genes3alt <- c("VWF", "FCGBP", "SERPINA3")
genes3alt.vwf.serpina3 <- c("VWF", "FCGBP")
genes3alt.vwf.fcgbp <- c("VWF", "SERPINA3")
genes4alt <- c("FCGR3A", "LCP1", "SERPINA3", "BCHE")
# genes with high intercept and low p-values (Holm < 0.05)
# note: IGF2R has 3 missing values, so don't include it
#genes3alt2 <- c("SERPINA3", "GPX3", "VWF")

genesUse <- genes4alt
#genesUse <- genes3alt.vwf.serpina3
#genesUse <- genes3alt.vwf.fcgbp
#genesUse <- genes6
grpselect <- {geneList %in% genesUse}

# # # # # # # # # #  #
# read in data
# # # # # # # # # # #

setwd("C:\\Users\\mooredf\\Box\\gp10\\Res\\MariaGennaroMISC\\manuscript\\RCode\\data")
log2out.orig <- read.csv(file="mis-c_dataSet2.csv")
log2out.test <- log2out.orig[,-1]  # just the data
rownames(log2out.test) <- log2out.orig[,1]
# same data as log2out.c
names(log2out.test)
# strip off the period and numerical suffix to just leave the group indicator name
grpIndUse <- matrix(unlist(strsplit(x=names(log2out.test), split="[.]")),
                    ncol=2, byrow=TRUE)[,1]

grpIndNames <- let.sample.f   # this needs fixing
grpIndUse <- grpIndNames
grpInd <- as.numeric(grpIndNames == "M")
table(grpInd)

# log2out.c <- log2out.test   # (if reading in the data from a file)
# # # # # # # # # # # # # # #
# prediction for a gene signature
# # # # # # # # # # # # # # #

geneList <- matrix(unlist(strsplit(x=rownames(log2out.c), split ="_")),
                   ncol=2, byrow=TRUE)[,1]
# gene signatures
genes4alt <- c("FCGR3A", "LCP1", "SERPINA3", "BCHE")
genes3alt <- c("VWF", "FCGBP", "SERPINA3")

set.seed(823390)
predictSignature <-
  function(genesUse=genes4alt, geneNamesA=geneNames,
           log2out=log2out.c, cv=TRUE) {

  grpselect <- {geneNamesA %in% genesUse}

  log2protData_sig <- t(log2out[grpselect,])  # just the gene signature

  dim(log2protData_sig)
  # [1] 59  3
  head(log2protData_sig)
  numbMiss <- function(xxx) {sum(is.na(xxx))}
  #print(apply(log2protData_sig,2,numbMiss))
  if (ncol(log2protData_sig) > 1) {
    numMissing <- apply(log2protData_sig,2,numbMiss)
  }
  if (ncol(log2protData_sig) == 1) {
    numMissing <- numbMiss(log2protData_sig)
  }

  # 0  1
  #40 19    # 19 MIS-C vs 40 others

# create indicators for columns with complete data; just keep those
  sampleMissInd <- complete.cases(log2protData_sig)
  log2protData_sigUse <- log2protData_sig[sampleMissInd,]
  grpIndUse <- grpInd[sampleMissInd]  # group indicators for groups being used

  if (is.matrix(log2protData_sigUse)) {
    numMissingAfter <- apply(log2protData_sigUse,2,numbMiss)
  }
  if (!is.matrix(log2protData_sigUse)) {
    numMissingAfter <- numbMiss(log2protData_sigUse)
  }

  #numMissingAfter <- apply(log2protData_sigUse,2,numbMiss)
# all should be complete now


  xxAll <- log2protData_sigUse

  table(grpIndUse)
  names(log2out.c) <- grpIndUse
  #grpIndUseM <- as.numeric({grpIndUse == "M"})
  #table(grpIndUseM)
  #yy <- as.factor(grpIndUseM)
  yy <- as.factor(grpIndUse)


  resultM_Other <- svmProt(xx=xxAll, yy=yy)
  library(crossval)
  cv.out <- NULL
  out.cv <- NULL

  if (cv) cv.out = crossval(svmProt2, X=xxAll, Y=yy, K=5, B=5)

  #print(resultM_Other)
  #print(cv.out)
  out.svm <- resultM_Other  # SVM prediction results
  out2all <- rbind(cv.out$stat, cv.out$stat.se)
  if (cv) out.cv <- out2all[,1:4]   # cross-validated SVM results
  result <- list(numMissing=numMissing, numMissingAfter=numMissingAfter,
                 out.svm=out.svm, out.cv=out.cv)
  result
  }

# gene signatures
genes4Nygaard2024 <- c("FCGR3A", "LCP1", "SERPINA3", "BCHE") #Nygaard2024
genes3alt <- c("VWF", "FCGBP", "SERPINA3")  # Jeisac/Gennaro
genesDiorio2021 <- "PLA2G2A"
genesPatel2024  <- c("MRPL58", "LTA4H", "BTLA")  # only LTA4H found
genesGruber  <- c("CXCL5", "CXCL11", "CXCL1", "CXCL6",
                  "CD40")
genesYeoh <- c("CD163", "PCSK9")

pred.3 <- predictSignature(genesUse=genes3alt)
pred.3$out.svm
pred.3$out.cv
pred.4 <- predictSignature(genesUse=genes4Nygaard2024)

predDiorio <- predictSignature(genesUse=genesDiorio2021, cv=FALSE)
PredPatel <- predictSignature(genesUse=genesPatel2024, cv=FALSE)
predGruber <- predictSignature(genesUse=genesGruber, cv=FALSE)
predYeoh <- predictSignature(genesUse=genesYeoh, cv=TRUE)
#plot(rocAll)
#sens.ci <- ci.se(rocAll)
#plot(sens.ci, type="shape", col="lightblue")

#resultReturn

#dev.off()
# # # # # # # # # # # # # # # # # # # #
#



# UpSetR plot
library(UpSetR)
temp <- split(geneHolmList[,2], geneHolmList[,1])
temp2 <- split(geneHolmList[,1], geneHolmList[,2])
head(temp)
pdf(file="upsetPlotA.pdf")
upset(fromList(temp), nsets=6, order.by="degree", text.scale=1.5)
dev.off()

pdf(file="upsetPlot.pdf")
upset(fromList(temp2), nsets=104, order.by="degree",
      point.size=1, line.size=0.4, text.scale=0.8, mb.ratio=c(0.1,0.9),
      mainbar.y.label="n")
dev.off()

# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # #
# prediction
head(log2out.c)
dim(log2out.c)
# [1] 269  59

genesNoMissing <- rownames(dataOutAllSQLite.ord)[dataOutAllSQLite.ord$sumMiss == 0]
numPredictors <- rep(NA, 15)
rawAccuracy <- rep(NA, 15)
cvAccuracy <- rep(NA, 15)

set.seed(5829021)  # for reproducibility
for (jj in 1:15) {
  #jj=3
  genesUse <- genesNoMissing[1:jj]
  log2outX <- log2out.c[rownames(log2out.c) %in% genesUse,]
  #dim(log2outX)
  #[1]  3 59
  xx <- t(log2outX)
  table(let.sample.f)
  # M  K  P  A
  #19 13 17 10


  yy <- let.sample.f

  model.orig <- svm(yy ~ ., probability=TRUE, cross=5, data=xx)
  accuracy <- model.orig$tot.accuracy
  pred <- predict(model.orig, newdata=xx, decision.values = FALSE,
                probability = TRUE)
  predictedValues <- (levels(pred))[pred]
  tableOut <- table(as.character(yy), predictedValues)
  rawAccuracy.jj <- 100*sum(diag(tableOut))/length(y.orig.f)
  rawAccuracy[jj] <- rawAccuracy.jj
  accuracy.jj <- model.orig$tot.accuracy
  cvAccuracy[jj] <- accuracy.jj
  numPredictors[jj] <- jj
  #rawAccuracy
  # [1] 79.66102
  #accuracy
  # [1] 69.49153  # corrected by cross-validation
}
resultCV <- data.frame(numPredictors, rawAccuracy, cvAccuracy)
# windows(height=5, width=7)
# par(mar=c(5,3,5,1))
plot(rawAccuracy ~ numPredictors, type="b",
     xlab="Number of proteins", ylab="Accuracy (percent)",
     ylim=c(60, 100),
     main="Predictive accuracy by number of proteins")
lines(cvAccuracy ~ numPredictors, type="b", col="red")
legend(x="topleft", col=c("black", "red"),
       lty=1, legend=c("raw accuracy", "cv accuracy"))

tempOrderedProteinList <- data.frame(genesNoMissing[1:15])
names(tempOrderedProteinList) <- "Protein"

  # # # # # # # # # # # # # # # # # #
  # # # # # # # # # # # # # # # # # #

all.pca <- prcomp(xxAll, center=TRUE, scale=TRUE)
#all.pca <- prcomp(t(log2outN), center=TRUE, scale=TRUE)
summary(all.pca)
pca.x <- all.pca$x
sd.pca <- all.pca[[1]]
var.pca <- sd.pca^2
prop.var.pca <- var.pca/sum(var.pca)
percentExplained <- round(prop.var.pca*100, digits=1)

windows()
par(mar=c(5,4,2,2))
plot(all.pca)

windows(width=7, height=7)
#par(mfrow=c(2,2))
plot(pca.x[,2] ~ pca.x[,1],
     xlab=paste("PC 1: ", percentExplained[1], "% variance", sep=''),
     ylab=paste("PC 2: ", percentExplained[2], "% variance", sep=''),
     type="n")
#    xlim=c(-7,6), ylim=c(-7, 5.5))
#points(pca.x[,2] ~ pca.x[,1], col=col.sample.b, pch=16, cex=0.9)
#points(pca.x[,2] ~ pca.x[,1], col=col.sample.b, pch=pch.sample.b, cex=0.9)
text(pca.x[,2] ~ pca.x[,1], col=col.sample.b, labels=let.sample.b , cex=0.9)
#text(pca.x[,2] ~ pca.x[,1], col=(ABind+1), labels=ABind , cex=0.9)
#title("K vs P")

windows()
biplot(x=all.pca$x[,1:2], all.pca$rotation[,1:2], xlabs=let.sample.b,
       cex=0.8, xlim=c(-4.8,3.6))








# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # #

# following not up to date as of June 21, 2024


#log2outPCAx <- log2out.b[pValHolm < 0.0001,]
log2outPCAx <- log2out.c[pValHolm < 0.000003,]
rownames(log2outPCAx)
# [1] "SERPINA3_P01011" "VWF_P04275"      "PFN1_P07737"     "FCGBP_Q9Y6R7"
# remove immunoglobins numbers 1, 2, 4
#log2outPCA <- log2outPCAx[-c(1,2,4),]
log2outPCA <- log2outPCAx
rownames(log2outPCA)
# [1] "SERPINA3_P01011" "VWF_P04275"      "PFN1_P07737"     "FCGBP_Q9Y6R7"

# # # # # # # # # # #
# for three genes only
# # # # # # # # # # #
xmat <- as.matrix(log2outPCA)
windows(height=10,width=7)
par(mfrow=c(2,2))
plot(xmat[2,] ~ xmat[1,],
     type="n")
#    xlim=c(-7,6), ylim=c(-7, 5.5))
points(xmat[2,] ~ xmat[1,], col=col.sample, pch=16, cex=0.9)

plot(xmat[3,] ~ xmat[1,],
     type="n")
#    xlim=c(-7,6), ylim=c(-7, 5.5))
points(xmat[3,] ~ xmat[1,], col=col.sample, pch=16, cex=0.9)

plot(xmat[2,] ~ xmat[3,],
     type="n")
#    xlim=c(-7,6), ylim=c(-7, 5.5))
points(xmat[2,] ~ xmat[3,], col=col.sample, pch=16, cex=0.9)

# # # # # # # # # # #
# end for three genes only
# # # # # # # # # # #

log2outPCA <- log2out.c[complete.cases(log2out.c),]
dim(log2outPCA)
# [1] 176  59
# "SERPINA3_P01011" "PGLYRP2_Q96PD5"  "ATRN_O75882-2"   "THBS1_P07996"    "LBP_P18428"      "BCHE_P06276"
# "VWF_P04275"      "FCGBP_Q9Y6R7"

all.pca <- prcomp(t(log2outPCA), center=TRUE, scale=TRUE)
#all.pca <- prcomp(t(log2outN), center=TRUE, scale=TRUE)
summary(all.pca)
pca.x <- all.pca$x
sd.pca <- all.pca[[1]]
var.pca <- sd.pca^2
prop.var.pca <- var.pca/sum(var.pca)
percentExplained <- round(prop.var.pca*100, digits=1)

windows()
par(mar=c(5,4,2,2))
plot(all.pca)

windows(width=7, height=7)
par(mfrow=c(2,2))
plot(pca.x[,2] ~ pca.x[,1],
     xlab=paste("PC 1: ", percentExplained[1], "% variance", sep=''),
     ylab=paste("PC 2: ", percentExplained[2], "% variance", sep=''),
     type="n")
#    xlim=c(-7,6), ylim=c(-7, 5.5))
#points(pca.x[,2] ~ pca.x[,1], col=col.sample.b, pch=16, cex=0.9)
#points(pca.x[,2] ~ pca.x[,1], col=col.sample.b, pch=pch.sample.b, cex=0.9)
text(pca.x[,2] ~ pca.x[,1], col=col.sample.b, labels=let.sample.b , cex=0.9)
#text(pca.x[,2] ~ pca.x[,1], col=(ABind+1), labels=ABind , cex=0.9)
#title("K vs P")

windows()
biplot(princomp(t(log2outPCA)) )


windows(width=7, height=11)

plot(pca.x[,3] ~ pca.x[,1],
     xlab=paste("PC 1: ", percentExplained[1], "% variance", sep=''),
     ylab=paste("PC 3: ", percentExplained[3], "% variance", sep=''),
     type="n")
#    xlim=c(-7,6), ylim=c(-7, 5.5))
#points(pca.x[,2] ~ pca.x[,1], col=col.sample.b, pch=16, cex=0.9)
#points(pca.x[,3] ~ pca.x[,1], col=col.sample.b, pch=pch.sample.b, cex=0.9)
text(pca.x[,3] ~ pca.x[,1], col=col.sample.b, labels=let.sample.b , cex=0.9)
#text(pca.x[,3] ~ pca.x[,1], col=(ABind+1), labels=ABind , cex=0.9)
#title("K vs P")

windows(width=7, height=7)
plot(pca.x[,3] ~ pca.x[,2],
     xlab=paste("PC 1: ", percentExplained[1], "% variance", sep=''),
     ylab=paste("PC 3: ", percentExplained[3], "% variance", sep=''),
     type="n")
#    xlim=c(-7,6), ylim=c(-7, 5.5))
#points(pca.x[,2] ~ pca.x[,1], col=col.sample.b, pch=16, cex=0.9)
#points(pca.x[,3] ~ pca.x[,1], col=col.sample.b, pch=pch.sample.b, cex=0.9)
text(pca.x[,3] ~ pca.x[,2], col=col.sample.b, labels=let.sample.b , cex=0.9)
#text(pca.x[,3] ~ pca.x[,2], col=(ABind+1), labels=ABind , cex=0.9)
#title("K vs P")

rownames(log2outPCA)


# # # # # # # # # # # # # # # # # # #
# compare misc to pneumonia
# # # # # # # # # # # # # # # # # # #

ncol(kawasakiData)
#  13
ncol(mildAsympData)
# 10
ncol(miscData)
# 19
ncol(pneumoniaData)
# 17

grpInd <- let.sample.b
table(grpInd)
# grpInd
# A  K  M  P
# 10 13 19 17

# A: Asymptomatic/mild
# K: Kawasaki
# M: MIS-C
# P: Pneumonia

dim(log2out.b)
#[1] 185  59
dim(log2outPCA)
# [1] 8 59
# the first 13 are kawasaki and the next 10 are mildAsymp
# remove those
#miscPneumon <- log2outPCA[,-c(1:23)]
#dim(miscPneumon)
# [1]  8 36

dataSubsetAvsB <- function(log2out, grpInd, A = "M", B = "P") {
  compareA <- (1:length(grpInd))[grpInd == A]
  compareB <- (1:length(grpInd))[grpInd == B]
  #compareInd <- (1:length(grpInd))[{{grpInd == A} | {grpInd == B}}]
  compareInd <- c(compareA, compareB)
  grpAgrpB <- log2outPCA[,compareInd]
  ABind <-c(rep(1, sum(grpInd == A)), rep(0, sum(grpInd == B)))

# first 19 are misc, next 17 are pneumonia
  # miscInd <- c(rep(1, 19), rep(0,17))


#log2outAllN <- miscPneumon
  #log2outAllN <- grpAgrpB
  result <- list(grpAgrpB=grpAgrpB, ABind=ABind)
  result
}

result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "M", B = "P")
log2outN <- result$grpAgrpB
ABind <- result$ABind

rocAnal <- function(log2outN, ABind) {
# now set up for svm:
library(e1071)
library(pROC)
library(crossval)  # need for "confusionMatrix" and "diagnosticErrors"
log2outAllN <- as.matrix(log2outN)
x=t((log2outAllN))


#y=as.factor(miscInd)  # misc is 1, pneumonia is 0
y=as.factor(ABind)  # misc is 1, pneumonia is 0
# x.orig <- x
# y.orig <- y
# log2outOrig <- log2outAllN

nAll <- length(y)  # number of observations  # 29

#model <- svm(y ~ x, probability=TRUE)
model.orig <- svm(y ~ ., probability=TRUE, data=x, cross=5)

pred <- predict(model.orig, newdata=x, decision.values = FALSE, probability = TRUE)

# never do this!!
# predictedValues <- as.numeric(pred[1:nAll])

# see https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-do-I-convert-factors-to-numeric_003f
predictedValues <- as.numeric(levels(pred))[pred]
# these functions are from the "crossval" library:
cm <- confusionMatrix(a=y, p=predictedValues, negative=0)
#print(diagnosticErrors(cm) )
temp <- diagnosticErrors(cm)
sens <- temp[2]
spec <- temp[3]

table(y, predictedValues)


rocAll <- roc(response=y,
              predictor=attr(pred, "probabilities")[,"1"], ci=TRUE,
              smoothed=TRUE)
#pick off correct column, the one with column name "1";
# sometimes it is the first column, sometimes the second

AUC <- as.numeric(auc(rocAll))
result <- c(sens, spec, AUC)
result
}

ncol(kawasakiData)
#  13
ncol(mildAsympData)
# 10
ncol(miscData)
# 19
ncol(pneumoniaData)
# 17a
# A: Asymptomatic/mild
# K: Kawasaki
# M: MIS-C
# P: Pneumonia
result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "M", B = "P")
log2outN <- result$grpAgrpB
ABind <- result$ABind
rocAnal(log2outN, ABind)

# # # # # # # # # #
result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "P", B = "M")
log2outN <- result$grpAgrpB
ABind <- result$ABind
rocAnal(log2outN, ABind)

result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "K", B = "M")
log2outN <- result$grpAgrpB
ABind <- result$ABind
rocAnal(log2outN, ABind)

result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "K", B = "P")
log2outN <- result$grpAgrpB
ABind <- result$ABind
rocAnal(log2outN, ABind)


result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "K", B = "A")
log2outN <- result$grpAgrpB
ABind <- result$ABind
rocAnal(log2outN, ABind)

result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "P", B = "A")
log2outN <- result$grpAgrpB
ABind <- result$ABind
rocAnal(log2outN, ABind)

result <- dataSubsetAvsB(log2out=log2outPCA, grpInd=grpInd, A = "A", B = "M")
log2outN <- result$grpAgrpB
ABind <- result$ABind
rocAnal(log2outN, ABind)





plot(rocAll)
sens.ci <- ci.se(rocAll)
plot(sens.ci, type="shape", col="lightblue")

# # # # # # # # # # # # # # # # # # #
# compare misc to kawasaki
# # # # # # # # # # # # # # # # # # #

ncol(kawasakiData)   # 1 to 13
#  13
ncol(mildAsympData)  # 14 to 23
# 10
ncol(miscData)       # 24 to 42
# 19
ncol(pneumoniaData)  # 43 to 59
# 17

dim(log2out.b)
#[1] 185  59
dim(log2outPCA)
# [1] 8 59
# the first 13 are kawasaki and the next 10 are mildAsymp
# remove those
kawasakiMisc <- log2outPCA[,c(1:13, 24:42)]
dim(kawasakiMisc)
# [1]  8 32
# first 8 are kawasaki, next 19 are misc
miscInd2 <- c(rep(1, 13), rep(0,19))


log2outAllN <- kawasakiMisc



# now set up for svm:
library(e1071)
library(pROC)
library(crossval)  # need for "confusionMatrix" and "diagnosticErrors"

x=t((log2outAllN))


y=as.factor(miscInd2)  # kawasaki is 1, misc is 0
# x.orig <- x
# y.orig <- y
# log2outOrig <- log2outAllN

nAll <- length(y)  # number of observations  # 29

#model <- svm(y ~ x, probability=TRUE)
model.orig <- svm(y ~ ., probability=TRUE, data=x, cross=5)

pred <- predict(model.orig, newdata=x, decision.values = FALSE, probability = TRUE)

# never do this!!
# predictedValues <- as.numeric(pred[1:nAll])

# see https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-do-I-convert-factors-to-numeric_003f
predictedValues <- as.numeric(levels(pred))[pred]
# these functions are from the "crossval" library:
cm <- confusionMatrix(a=y, p=predictedValues, negative=0)
diagnosticErrors(cm)

table(y, predictedValues)








rocAll <- roc(response=y,
              predictor=attr(pred, "probabilities")[,"1"], ci=TRUE,
              smoothed=TRUE)
#pick off correct column, the one with column name "1";
# sometimes it is the first column, sometimes the second

auc(rocAll)
plot(rocAll)
sens.ci <- ci.se(rocAll)
plot(sens.ci, type="shape", col="lightblue")



