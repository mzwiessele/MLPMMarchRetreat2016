
# train-test --------------------------------------------------------------

perfSVM <- function(model = 'ksvm', x_prefix = 'xdata', y_prefix = 'ydata', 
                    xdata, grp, tr2tstFolds, kf = c('linear', 'kendall', 'rbf', 'poly'), 
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

    message("\nIndependent-validating ", appendLF = TRUE)
    kmcs <- centerScaleKernelMatrix(kmat, indepfold)

    cpm <- Cpara_list[which.max(cvacc)]
    pred <- classifierSVM(km=kmcs, trainidx=indepfold, traingrp=grp[indepfold], cpm=cpm)
    acc <- evaluateAcc(pred,grp[-indepfold])
    accind <- evaluateIndivAcc(pred,grp[-indepfold])

    message('DONE!', appendLF = TRUE)
    return(c(acc = acc, unlist(accind)))
  })
  tr2tstfoldscore <- do.call('rbind', tr2tstfoldscore)
  tr2tstfoldscore <- apply(tr2tstfoldscore, 2, mean, na.rm = TRUE)
  acc <- tr2tstfoldscore[1]
  accind <- tr2tstfoldscore[-1]

  return(list(model=model, x_prefix=x_prefix, y_prefix=y_prefix, kf=kf, Cpara_list=Cpara_list, 
              acc=acc, accind = accind))
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

evaluateIndivAcc <- function(predictions,observations)
{
  classes <- sort(unique(as.character(observations)))
  accind <- list()
  for (classname in classes) {
    id <- which(as.character(observations) == classname)
    s <- sum(as.character(predictions)[id] == classname) / length(id) # similar to True Positive
    accind[[classname]] <- ifelse(is.na(s), 0, s)
  }
  return(accind)
}

# remove constants --------------------------------------------------------

removeConst <- function (xtr, xtst=NULL, tol=1e-6)
{

  numid <- which(sapply(as.data.frame(xtr), is.numeric))
  constid <- apply(xtr[,numid], 2, function(u){
    (mean(u, na.rm = TRUE) < tol && diff(range(u, na.rm = TRUE)) < tol) || (diff(range(u, na.rm = TRUE)) / mean(u, na.rm = TRUE) < tol)
  })
  xtr <- xtr[ , numid[!constid], drop = F]
  if(is.null(xtst)){
    return(xtr)
  } else{
    xtst <- xtst[ , numid[!constid], drop = F]
    return(list(xtr=xtr, xtst=xtst))
  }
}

# Normalize data ----------------------------------------------------------

normalizeData <- function (xtr, xtst=NULL, do.center=TRUE, do.scale=TRUE, ...)
{
  # normalizeData() normalizes col-wisely to mean 0 (do.center=T) and sd 1 (do.scale=T)
  # NOTE: xtst is normalized wrt xtr

  if(is.vector(xtr)){
    warning("vector features are coerced to a matrix!")
    xtr <- matrix(xtr, ncol=1)
    if(!is.null(xtst)){
      xtst <- matrix(xtst, ncol=1)
    }
  }

  numid <- sapply(as.data.frame(xtr), is.numeric)
  t <- scale(xtr[,numid], center=do.center, scale=do.scale)
  xtr[,numid] <- t

  if(is.null(xtst)){
    return(xtr)
  } else{
    if(do.center){
      xtst[,numid] <- sweep(xtst[,numid], 2, attr(t,"scaled:center"), FUN="-")
    }
    if(do.scale){
      xtst[,numid] <- sweep(xtst[,numid], 2, attr(t,"scaled:scale"), FUN="/")
    }
    return(list(xtr=xtr, xtst=xtst))
  }
}


# Sparse SVM --------------------------------------------------------------

SparseSVM <- function(xtr,
                      xtst,
                      ytr,
                      cost = 1,
                      do.normalize = TRUE,
                      ...){
  # L1-regularized L2-loss support vector classification

  library(LiblineaR)

  if(do.normalize){
    d <- normalizeData(xtr, xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }

  o <- order(ytr, decreasing=TRUE)
  model <- LiblineaR(xtr[o,], ytr[o], type=5, cost=cost, ...)

  pred <- predict(model, xtst, decisionValues = TRUE)
  return(pred$predictions)
}


# Perform Sparse SVM ------------------------------------------------------

perfSparseSVM <- function(model = 'kPasrseSVM',
                          prefix,
                          xdata,
                          grp,
                          tr2tstFolds,
                          # kf = c('linear', 'kendall', 'rbf', 'poly'),
                          Cpara_list = 10^(-3:3),
                          nfolds = 5,
                          nrepeats = 1,
                          seed = 64649601){
  ## model and prefix are characters that are stored to indicate which method is implemented but essentially do nothing
  ## perfSVM basically runs multi-class kernel SVM
  ## tr2tstFolds are training indices corresp to rows in xdata
  ## returns the average acc over tr2tstFolds with C parameter tuned by additional inner loop (seed set)

  if (!is.null(seed)) set.seed(seed)

  if (is.vector(tr2tstFolds)){
    tr2tstFolds <- list(tr2tstFolds)
  }

#   if (!is.character(kf)) stop('provide kernel function as a character')
#   kf <- match.arg(kf)
#
#   sigma <- sigest(xdata, scaled=F)['50%']
#
#   FUN <- switch(kf,
#                 linear = vanilladot(),
#                 kendall = cor.fk,
#                 rbf = rbfdot(sigma = sigma),
#                 poly = polydot(degree = 2, scale = 1, offset = 0))
#
#   kmat <- computeKernelMatrix(xdata = xdata, kf = FUN)

  # outter loop
  tr2tstfoldscore <- lapply(tr2tstFolds, function(indepfold){
    foldIndices <- caret::createMultiFolds(indepfold, k=nfolds, times=nrepeats)
    message("Cross-validating ", appendLF = TRUE)
    foldscore <- lapply(foldIndices, function(fold){
      message(' + New cv fold + ', appendLF = FALSE)
      xtr.xtst <- removeConst(xtr = xdata[indepfold[fold], ], xtst = xdata[indepfold[-fold], ])

      ### parse SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = FALSE)
        pred <- SparseSVM(xtr = xtr.xtst$xtr,
                          xtst = xtr.xtst$xtst,
                          ytr = grp[indepfold[fold]],
                          cost = cpm)
        return(evaluateAcc(pred, grp[indepfold[-fold]]))
      })
      return(s)
    })

    cvacc <- apply(as.matrix(as.data.frame(foldscore)),1,function(u){mean(u,na.rm=TRUE)})
    names(cvacc) <- names(Cpara_list) # SVM cv score (varying C)

    message("\nIndependent-validating ", appendLF = TRUE)
    cpm <- Cpara_list[which.max(cvacc)]
    xtr.xtst <- removeConst(xtr = xdata[indepfold, ], xtst = xdata[-indepfold, ])
    pred <- SparseSVM(xtr = xtr.xtst$xtr,
                      xtst = xtr.xtst$xtst,
                      ytr = grp[indepfold],
                      cost = cpm)
    acc <- evaluateAcc(pred, grp[-indepfold])

    message('DONE!', appendLF = TRUE)
    return(acc)
  })
  acc <- mean(unlist(tr2tstfoldscore), na.rm = TRUE)

  return(list(model=model, prefix=prefix, acc=acc))
}

