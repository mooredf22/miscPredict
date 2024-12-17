
# date:  December 15, 2024
#create the "logMISCdata.rda") data set from "AllSQLite.RData"
#Put it here (using the "save" command) for the "miscPredict" package
  "C:\\Users\\mooredf\\Desktop\\Rprojects\\miscPredict\\data"

# This code is derived from "AllSQLitePCA.R" in this folder:
#  C:\Users\mooredf\Box\gp9\NJacts\investigators\MariaGennaro\targetedMSexpand\progs

  # Here are the types of samples (all children) and sample sizes:
  #  M  K  P  A
  # 19 13 17 10

  # M: MIS-C
  # K: Kawasaki
  # P: Pneumonia
  # A: Asymptomatic/mild

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
if (FALSE) {
  #k.ind <- completeX(datax=kawasakiData, X=4)
  #table(k.ind)
  k.ind <- rep(TRUE, nrow(kawasakiData))
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
}


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
logMISCdata <- log2out.c
setwd("C:\\Users\\mooredf\\Desktop\\Rprojects\\miscPredict\\data")
#save(logMISCdata, file="logMISCdata.rda")
