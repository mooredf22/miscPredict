#' create svm model to predict yy from xx and report the properties
#' @param genesUse vector of gene names comprising the signature
#' @param geneNames vector of all gene names correspondint to the rows of data
#' @param grpInd 0/1 indicator for group (control / MIS-C)
#' @param log2out log2 transformed data frame of protein abundances
#' @param cv if this is true, calculate cross-validation error rates
#' @return result a dataframe with properties of the predictive model
#' @import pROC
#' @import crossval
#' @export

predictSignature <-
  function(genesUse=genesUse, geneNames=geneNames, grpInd,
           log2out=log2out, cv=TRUE) {

    grpselect <- {geneNames %in% genesUse}

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

    table(grpIndUse)  # should be 0 or 1
    names(log2out) <- grpIndUse
    #grpIndUseM <- as.numeric({grpIndUse == "M"})
    #table(grpIndUseM)
    #yy <- as.factor(grpIndUseM)
    yy <- as.factor(grpIndUse)


    out.svm <- svmProt(xx=xxAll, yy=yy)
    #library(crossval)
    cv.out <- NULL
    out.cv <- NULL

    if (cv) cv.out = crossval(svmProt2, X=xxAll, Y=yy, K=5, B=5,
                              verbose=FALSE)

    #print(resultM_Other)
    #print(cv.out)
    # SVM prediction results
    out.auc <- matrix(c(out.svm$aucOut, out.svm$aucCI),
                         nrow=1, ncol=3,
                         dimnames=list("", c("AUC", "lowerCI", "upperCI")))
    out2all <- rbind(cv.out$stat, cv.out$stat.se)
    if (cv) out.cv <- out2all[,1:4]   # cross-validated SVM results
    #if (rocPlot) {
    #  rocAll <- out.svm$rocAll
    #  plot(rocAll, cex=1.5)
      #sens.ci <- ci.se(rocAll)
      #plot(sens.ci, type="shape", col="lightblue")
    #}
    result <- list(numMissing=numMissing, numMissingAfter=numMissingAfter,
                   out.auc=out.auc, out.svm=out.svm, out.cv=out.cv)
    result
  }
