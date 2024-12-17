
#' create svm model to predict yy from xx and report the properties
#' @param xx data frame of data
#' @param yy actual reponse variable to be predicted
#' @return result data structure with properties
#' @import e1071
#' @import pROC
#' @import crossval
#' @export
svmProt <- function(xx, yy) {
  # xx is a matrix of samples with columns indicating protein signature names
  # yy is a factor variable of 1's and 0's
  # now set up for svm:
  #library(e1071)
  #library(pROC)
  #library(crossval)  # need for "confusionMatrix" and "diagnosticErrors"
  x.orig <- as.data.frame(xx)
  y.orig <- as.factor(yy)
  nAll <- length(y.orig)  # number of observations

  #model.orig <- svm(y.orig.f ~ ., probability=TRUE, data=x.orig)
  #browser()
  model.orig <- svm(y.orig ~ ., probability=TRUE, cross=5, data=x.orig)
  accuracy <- model.orig$tot.accuracy
  pred <- predict(model.orig, newdata=x.orig, decision.values = FALSE, probability = TRUE)

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
  cm <- confusionMatrix(a=y.orig, p=predictedValues, negative=0)
  tempDiagErr <- diagnosticErrors(cm)
  sens <- tempDiagErr[2]
  spec <- tempDiagErr[3]

  tableOut <- table(y.orig, predictedValues)

  predMatrix <- attr(pred, "probabilities")
  # column names of predMatrix should be 1 and 0
  # otherwise predMatrix[,"1"] will not work!!
  rocAll <- roc(response=y.orig,
                predictor=predMatrix[,"1"], ci=TRUE,
                smoothed=TRUE)
  #pick off correct column, the one with column name "1";
  # sometimes it is the first column, sometimes the second

  aucOut <- as.numeric(auc(rocAll))
  aucCI <- as.numeric(rocAll$ci)[c(1,3)]

  result <- list(sens=sens, spec=spec, tableOut=tableOut,
                 aucOut=aucOut, aucCI=aucCI, rocAll=rocAll)
  result
}
