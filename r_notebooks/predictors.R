
# train-test --------------------------------------------------------------

perfSVM <- function(model = 'ksvm', prefix, xdata, grp, tr2tstFolds, kf = c('linear', 'kendall', 'rbf', 'poly'), 
                    Cpara_list = 10^(-3:3), nfolds = 5, nrepeats = 1, seed = 64649601){
  ## model and prefix are characters that are stored to indicate which method is implemented but essentially do nothing
  ## perfSVM basically runs multi-class kernel SVM
  ## tr2tstFolds are training indices corresp to rows in xdata
  ## returns the average acc over tr2tstFolds with C parameter tuned by additional inner loop (seed set)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.vector(tr2tstFolds)){
    tr2tstFolds <- list(tr2tstFolds)
  }
  
  if (!is.character(kf)) stop('provide kernel function as a character')
  kf <- match.arg(kf)
  
  sigma <- sigest(xdata, scaled=F)['50%']
  
  FUN <- switch(kf, 
                linear = vanilladot(), 
                kendall = cor.fk, 
                rbf = rbfdot(sigma = sigma), 
                poly = polydot(degree = 2, scale = 1, offset = 0))
  
  kmat <- computeKernelMatrix(xdata = xdata, kf = FUN)
  
  # outter loop
  tr2tstfoldscore <- lapply(tr2tstFolds, function(indepfold){
    foldIndices <- createMultiFolds(indepfold, k=nfolds, times=nrepeats)
    
    message("Cross-validating ", appendLF = TRUE)
    foldscore <- lapply(foldIndices, function(fold){
      message(' + New cv fold + ', appendLF = FALSE)
      # inner loop
      kmcs <- centerScaleKernelMatrix(kmat[indepfold, indepfold], fold)
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = FALSE)
        pred <- classifierSVM(km=kmcs, trainidx=fold, traingrp=grp[indepfold[fold]], cpm=cpm)
        return(evaluateAcc(pred,grp[indepfold[-fold]]))
      })
      
      return(s)
    })
    
    cvacc <- apply(as.matrix(as.data.frame(foldscore)),1,function(u){mean(u,na.rm=TRUE)})
    names(cvacc) <- names(Cpara_list) # SVM cv score (varying C)
    
    message("Independent-validating ", appendLF = TRUE)
    kmcs <- centerScaleKernelMatrix(kmat, indepfold)
    
    cpm <- Cpara_list[which.max(cvacc)]
    pred <- classifierSVM(km=kmcs, trainidx=indepfold, traingrp=grp[indepfold], cpm=cpm)
    acc <- evaluateAcc(pred,grp[-indepfold])
    
    message('DONE!', appendLF = TRUE)
    return(acc)
  })
  acc <- mean(unlist(tr2tstfoldscore), na.rm = TRUE)
  
  return(list(model=model, prefix=prefix, kf=kf, Cpara_list=Cpara_list, acc=acc))
}


# classifiers -------------------------------------------------------------

classifierSVM <- function(km, trainidx, traingrp, cpm){
  
  # classifierSVM returns a ntst-vector of predicted binary classes
  # NOTE those indices not in trainidx is ready for test!!
  # NOTE that (row/sample) NAMES of kernel matrix is necessary for naming predicted vector
  
  # training
  mo <- ksvm(as.kernelMatrix(km[trainidx,trainidx]), y=traingrp, C=cpm, type="C-svc")
  # predicting
  pred <- predict(mo, as.kernelMatrix(km[-trainidx,trainidx,drop=F][,SVindex(mo),drop=F]), type = "response")
  names(pred) <- rownames(km)[-trainidx]
  
  return(pred)
}


# compute kernel matrix ---------------------------------------------------

computeKernelMatrix <- function(xdata, kf){
  ### xdata n*p data matrix
  n <- nrow(xdata); samplenames <- rownames(xdata)
  
  class(kf) <- 'kernel'
  kmat <- kernelMatrix(kf, as.matrix(xdata))
  dimnames(kmat) <- list(samplenames,samplenames)
  
  return(as.kernelMatrix(kmat))
}

centerScaleKernelMatrix <- function(kmat, trainidx){
  n <- nrow(kmat); ntr <- length(trainidx); nts <- n - ntr;
  
  #centering
  ed <- matrix(0,nrow=n,ncol=n)
  ed[trainidx,] <- 1
  kmcs <- kmat - t(ed)%*%kmat/ntr - kmat%*%ed/ntr + t(ed)%*%kmat%*%ed/(ntr^2)
  #scaling
  dsr <- sqrt(diag(kmcs))
  kmcs <- sweep(
    sweep(kmcs, 1, dsr, FUN='/'),
    2, dsr, FUN='/')
  rownames(kmcs) <- rownames(kmat)
  colnames(kmcs) <- colnames(kmat)
  
  return(as.kernelMatrix(kmcs))
}


# model evaluation --------------------------------------------------------

evaluateAcc <- function(predictions,observations){
  if(!is(predictions,"factor")) stop('Predictions are not factors!')
  if(all(levels(predictions)==levels(observations)))
    return(sum(predictions==observations,na.rm=T)/length(predictions))
  else
    stop('Predictions and observations have different levels!')
}
