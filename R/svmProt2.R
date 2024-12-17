# # # # # # # # #
# cross-validation of svm
# # # # # # # # #

#' create svm model for crossvalidation; for internal use
#' @param XTrain training data
#' @param YTrain response for training data
#' @param XTest test data
#' @param YTest response for test data
#' @return result data structure with properties
#' @import e1071
#' @import pROC
#' @import crossval
#' @export

svmProt2 <- function(XTrain, YTrain, XTest, YTest) {
  library(crossval)

  XTrain <- as.data.frame(XTrain)
  XTest <- as.data.frame(XTest)
  # xx is a matrix of samples with columns indicating protein signature names
  # now set up for svm:
  library(e1071)
  library(pROC)
  library(crossval)  # need for "confusionMatrix" and "diagnosticErrors"
  #browser()
  #x.orig <- xx
  #y.orig <- as.factor(yy)
  #nAll <- length(y.orig)  # number of observations

  #model.orig <- svm(y.orig.f ~ ., probability=TRUE, data=x.orig)
  model.orig <- svm(YTrain ~ ., probability=TRUE, cross=5, data=XTrain)
  accuracy <- model.orig$tot.accuracy
  names(accuracy) <- "accuracy"
  pred <- predict(model.orig, newdata=XTest, decision.values = FALSE, probability = TRUE)

  # never do this!!
  # predictedValues <- as.numeric(pred[1:nAll])

  # see https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-do-I-convert-factors-to-numeric_003f
  if (length(levels(pred)) == 2) {
    predictedValues <- as.numeric(levels(pred))[pred]
  }
  if (length(levels(pred)) > 2) {
    predictedValues <- (levels(pred))[pred]
  }


  # these functions are from the "crossval" library:
  cm <- confusionMatrix(a=YTest, p=predictedValues, negative=0)
  tempDiagErr <- diagnosticErrors(cm)
  sens <- tempDiagErr[2]
  spec <- tempDiagErr[3]

  tableOut <- table(YTest, predictedValues)


  rocAll <- roc(response=YTest,
                predictor=attr(pred, "probabilities")[,"1"], ci=TRUE,
                smoothed=TRUE)
  #pick off correct column, the one with column name "1";
  # sometimes it is the first column, sometimes the second

  aucOut <- as.numeric(auc(rocAll))
  names(aucOut) <- "AUC"
  aucCI <- as.numeric(rocAll$ci)[c(1,3)]

  #result <- list(sens=sens, spec=spec, tableOut=tableOut,
  #aucOut=aucOut, aucCI=aucCI, rocAll=rocAll)
  tableOutVec <- as.numeric(tableOut)  # make matrix into a vector
  #names(tableOutVec) <- c("n11", "n21", "n12", "n22")
  result <- c(accuracy, sens, spec, aucOut, tableOutVec)
  result
}
