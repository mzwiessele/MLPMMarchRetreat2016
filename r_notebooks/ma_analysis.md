# Notebook : Microarray analysis


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, fig.width = 12, fig.height = 8, dev = "pdf", fig.keep = "high", fig.path = "ma_figure/", cache.path = "ma_cache/")
set.seed(35875954)

library(caret)
library(pcaPP)
library(kernlab)
source('predictors.R')

datapath <- '../Data/MicroArray/'
```

# Load data


```r
# read in datasets
fnames <- list.files(path = datapath, pattern = '[.]csv$', full.names = FALSE)
objnames <- gsub('[.]csv$', '', fnames)
for (obj in objnames) {
  x <- read.csv(paste0(datapath, obj, '.csv'), sep = ',', header = TRUE, row.names = 1)
  rownames(x) <- paste0('X_', gsub('[^[:alnum:]_]', '_', rownames(x)))
  # a bit further modification
  if (grepl('design', obj)) {
    colnames(x) <- gsub('NA[.]', 'Train.Test', colnames(x))
  } else {
    x <- t(x) # observations in rows and features in col
  }
  assign(obj, x)
  cat('==============> ', obj, '<==============\n')
  print(dim(get(obj)))
}
```

```
## ==============>  microarray_all_genes <==============
## [1]   105 13613
## ==============>  microarray_design <==============
## [1] 105   4
## ==============>  microarray_effector_path_vals <==============
## [1]  105 1044
## ==============>  microarray_metabolic_mod_activities <==============
## [1] 105  89
## ==============>  microarray_metabolic_mod_genevalues <==============
## [1] 105 201
## ==============>  microarray_metabolic_mod_nodevalues <==============
## [1] 105 462
## ==============>  microarray_path_vals <==============
## [1]  105 4457
## ==============>  microarray_signaling_genes <==============
## [1]  105 2184
```

```r
# just some checking on samplelist
xlist <- objnames[!grepl('design', objnames)]
designobj <- setdiff(objnames, xlist)
stopifnot(length(designobj) == 1)
for (xname in xlist) {
  stopifnot(isTRUE(all.equal(rownames(get(xname)), 
                             as.character(get(designobj)[,grep('sample', colnames(get(designobj)), ignore.case = TRUE)]))))
}

# generate train2test
trainidx <- grep('train', get(designobj)[,'Train.Test'], ignore.case = TRUE)
table(get(designobj)[,'Train.Test'])
```

```
## 
##     Test set Training Set 
##           42           63
```

```r
# generate labels
effect.grps <- as.factor(get(designobj)[,grep('effect', colnames(get(designobj)), ignore.case = TRUE)])
table(effect.grps[trainidx])
```

```
## 
##  Control Effect_1 Effect_2 Effect_3 Effect_4 Effect_5 Effect_6 Effect_7 
##       18        9        9        9        9        9        0        0
```

```r
table(effect.grps[-trainidx])
```

```
## 
##  Control Effect_1 Effect_2 Effect_3 Effect_4 Effect_5 Effect_6 Effect_7 
##        6        9        9        0        0        0        9        9
```

```r
# groups with unknown labels
ratio <- 0.25
key <- 'unknown'
new.grps <- as.character(effect.grps)
train.char <- unique(new.grps[trainidx])
test.char <- unique(new.grps[-trainidx])
idx <- split(trainidx, factor(new.grps[trainidx]))
for (i in seq(length(idx))) {
  new.grps[sample(idx[[i]], ceiling(length(idx[[i]]) * ratio))] <- key # train
}
new.grps[which(new.grps %in% setdiff(test.char, train.char))] <- key # test
new.grps <- as.factor(new.grps)
table(new.grps[trainidx])
```

```
## 
##  Control Effect_1 Effect_2 Effect_3 Effect_4 Effect_5  unknown 
##       13        6        6        6        6        6       20
```

```r
table(new.grps[-trainidx])
```

```
## 
##  Control Effect_1 Effect_2 Effect_3 Effect_4 Effect_5  unknown 
##        6        9        9        0        0        0       18
```

```r
ylist <- c('new.grps')
```

# SVM


```r
kflist <- list('rbf', 'kendall', 'linear', 
               c('rbf','kendall'), c('linear', 'kendall'), c('linear', 'rbf'), 
               c('linear', 'rbf', 'kendall'))
res <- list()
i <- 1
for (xname in xlist) {
  for (yname in ylist) {
    for (kf in kflist) {
      message('\n==========================> ',i,' <==========================\n')
      res[[i]] <- perfSVM(model = 'ksvm', x_prefix = xname, y_prefix = yname, 
                          xdata = get(xname), grp = get(yname), kf = kf, tr2tstFolds = trainidx)
      i <- i + 1
    }
  }
}
```


```r
# score table
dd <- data.frame()
for (i in seq(length(res))) {
  dd <- rbind(dd, 
              data.frame(model = res[[i]]$model, kf = res[[i]]$kf, 
                         x_prefix = res[[i]]$x_prefix, y_prefix = res[[i]]$y_prefix, 
                         scores = names(res[[i]]$scores),
                         values = unname(res[[i]]$scores)))
}
lapply(split(dd, dd$scores), function(u) split(u[,setdiff(colnames(u), c('scores','kf'))], u$kf))
```

```
## $acc.Control
## $acc.Control$rbf
##     model                            x_prefix y_prefix    values
## 2    ksvm                microarray_all_genes new.grps 0.3333333
## 65   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 128  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 191  ksvm microarray_metabolic_mod_genevalues new.grps 0.6666667
## 254  ksvm microarray_metabolic_mod_nodevalues new.grps 0.5000000
## 317  ksvm                microarray_path_vals new.grps 0.0000000
## 380  ksvm          microarray_signaling_genes new.grps 0.6666667
## 
## $acc.Control$kendall
##     model                            x_prefix y_prefix    values
## 11   ksvm                microarray_all_genes new.grps 0.0000000
## 74   ksvm       microarray_effector_path_vals new.grps 0.1666667
## 137  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 200  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 263  ksvm microarray_metabolic_mod_nodevalues new.grps 0.8333333
## 326  ksvm                microarray_path_vals new.grps 0.1666667
## 389  ksvm          microarray_signaling_genes new.grps 0.1666667
## 
## $acc.Control$linear
##     model                            x_prefix y_prefix values
## 20   ksvm                microarray_all_genes new.grps    0.0
## 83   ksvm       microarray_effector_path_vals new.grps    0.0
## 146  ksvm microarray_metabolic_mod_activities new.grps    0.5
## 209  ksvm microarray_metabolic_mod_genevalues new.grps    0.0
## 272  ksvm microarray_metabolic_mod_nodevalues new.grps    0.0
## 335  ksvm                microarray_path_vals new.grps    0.5
## 398  ksvm          microarray_signaling_genes new.grps    0.5
## 
## $acc.Control$`rbf+kendall`
##     model                            x_prefix y_prefix    values
## 29   ksvm                microarray_all_genes new.grps 0.1666667
## 92   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 155  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 218  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 281  ksvm microarray_metabolic_mod_nodevalues new.grps 0.5000000
## 344  ksvm                microarray_path_vals new.grps 0.0000000
## 407  ksvm          microarray_signaling_genes new.grps 0.5000000
## 
## $acc.Control$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 38   ksvm                microarray_all_genes new.grps 0.0000000
## 101  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 164  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 227  ksvm microarray_metabolic_mod_genevalues new.grps 0.5000000
## 290  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 353  ksvm                microarray_path_vals new.grps 0.0000000
## 416  ksvm          microarray_signaling_genes new.grps 0.1666667
## 
## $acc.Control$`linear+rbf`
##     model                            x_prefix y_prefix    values
## 47   ksvm                microarray_all_genes new.grps 0.1666667
## 110  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 173  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 236  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 299  ksvm microarray_metabolic_mod_nodevalues new.grps 0.5000000
## 362  ksvm                microarray_path_vals new.grps 0.6666667
## 425  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $acc.Control$`linear+rbf+kendall`
##     model                            x_prefix y_prefix    values
## 56   ksvm                microarray_all_genes new.grps 0.1666667
## 119  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 182  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 245  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 308  ksvm microarray_metabolic_mod_nodevalues new.grps 0.5000000
## 371  ksvm                microarray_path_vals new.grps 0.5000000
## 434  ksvm          microarray_signaling_genes new.grps 0.3333333
## 
## 
## $acc.Effect_1
## $acc.Effect_1$rbf
##     model                            x_prefix y_prefix    values
## 3    ksvm                microarray_all_genes new.grps 0.0000000
## 66   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 129  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 192  ksvm microarray_metabolic_mod_genevalues new.grps 0.2222222
## 255  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 318  ksvm                microarray_path_vals new.grps 0.0000000
## 381  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_1$kendall
##     model                            x_prefix y_prefix values
## 12   ksvm                microarray_all_genes new.grps      0
## 75   ksvm       microarray_effector_path_vals new.grps      0
## 138  ksvm microarray_metabolic_mod_activities new.grps      0
## 201  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 264  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 327  ksvm                microarray_path_vals new.grps      0
## 390  ksvm          microarray_signaling_genes new.grps      0
## 
## $acc.Effect_1$linear
##     model                            x_prefix y_prefix    values
## 21   ksvm                microarray_all_genes new.grps 0.0000000
## 84   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 147  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 210  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 273  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 336  ksvm                microarray_path_vals new.grps 0.1111111
## 399  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_1$`rbf+kendall`
##     model                            x_prefix y_prefix values
## 30   ksvm                microarray_all_genes new.grps      0
## 93   ksvm       microarray_effector_path_vals new.grps      0
## 156  ksvm microarray_metabolic_mod_activities new.grps      0
## 219  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 282  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 345  ksvm                microarray_path_vals new.grps      0
## 408  ksvm          microarray_signaling_genes new.grps      0
## 
## $acc.Effect_1$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 39   ksvm                microarray_all_genes new.grps 0.0000000
## 102  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 165  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 228  ksvm microarray_metabolic_mod_genevalues new.grps 0.1111111
## 291  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 354  ksvm                microarray_path_vals new.grps 0.0000000
## 417  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_1$`linear+rbf`
##     model                            x_prefix y_prefix values
## 48   ksvm                microarray_all_genes new.grps      0
## 111  ksvm       microarray_effector_path_vals new.grps      0
## 174  ksvm microarray_metabolic_mod_activities new.grps      0
## 237  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 300  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 363  ksvm                microarray_path_vals new.grps      0
## 426  ksvm          microarray_signaling_genes new.grps      0
## 
## $acc.Effect_1$`linear+rbf+kendall`
##     model                            x_prefix y_prefix values
## 57   ksvm                microarray_all_genes new.grps      0
## 120  ksvm       microarray_effector_path_vals new.grps      0
## 183  ksvm microarray_metabolic_mod_activities new.grps      0
## 246  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 309  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 372  ksvm                microarray_path_vals new.grps      0
## 435  ksvm          microarray_signaling_genes new.grps      0
## 
## 
## $acc.Effect_2
## $acc.Effect_2$rbf
##     model                            x_prefix y_prefix    values
## 4    ksvm                microarray_all_genes new.grps 0.0000000
## 67   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 130  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 193  ksvm microarray_metabolic_mod_genevalues new.grps 0.1111111
## 256  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 319  ksvm                microarray_path_vals new.grps 0.0000000
## 382  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_2$kendall
##     model                            x_prefix y_prefix values
## 13   ksvm                microarray_all_genes new.grps      0
## 76   ksvm       microarray_effector_path_vals new.grps      0
## 139  ksvm microarray_metabolic_mod_activities new.grps      0
## 202  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 265  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 328  ksvm                microarray_path_vals new.grps      0
## 391  ksvm          microarray_signaling_genes new.grps      0
## 
## $acc.Effect_2$linear
##     model                            x_prefix y_prefix    values
## 22   ksvm                microarray_all_genes new.grps 0.0000000
## 85   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 148  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 211  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 274  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 337  ksvm                microarray_path_vals new.grps 0.2222222
## 400  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_2$`rbf+kendall`
##     model                            x_prefix y_prefix values
## 31   ksvm                microarray_all_genes new.grps      0
## 94   ksvm       microarray_effector_path_vals new.grps      0
## 157  ksvm microarray_metabolic_mod_activities new.grps      0
## 220  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 283  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 346  ksvm                microarray_path_vals new.grps      0
## 409  ksvm          microarray_signaling_genes new.grps      0
## 
## $acc.Effect_2$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 40   ksvm                microarray_all_genes new.grps 0.0000000
## 103  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 166  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 229  ksvm microarray_metabolic_mod_genevalues new.grps 0.2222222
## 292  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 355  ksvm                microarray_path_vals new.grps 0.0000000
## 418  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_2$`linear+rbf`
##     model                            x_prefix y_prefix values
## 49   ksvm                microarray_all_genes new.grps      0
## 112  ksvm       microarray_effector_path_vals new.grps      0
## 175  ksvm microarray_metabolic_mod_activities new.grps      0
## 238  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 301  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 364  ksvm                microarray_path_vals new.grps      0
## 427  ksvm          microarray_signaling_genes new.grps      0
## 
## $acc.Effect_2$`linear+rbf+kendall`
##     model                            x_prefix y_prefix values
## 58   ksvm                microarray_all_genes new.grps      0
## 121  ksvm       microarray_effector_path_vals new.grps      0
## 184  ksvm microarray_metabolic_mod_activities new.grps      0
## 247  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 310  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 373  ksvm                microarray_path_vals new.grps      0
## 436  ksvm          microarray_signaling_genes new.grps      0
## 
## 
## $acc.overall
## $acc.overall$rbf
##     model                            x_prefix y_prefix    values
## 1    ksvm                microarray_all_genes new.grps 0.4285714
## 64   ksvm       microarray_effector_path_vals new.grps 0.4285714
## 127  ksvm microarray_metabolic_mod_activities new.grps 0.4761905
## 190  ksvm microarray_metabolic_mod_genevalues new.grps 0.2619048
## 253  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4523810
## 316  ksvm                microarray_path_vals new.grps 0.4285714
## 379  ksvm          microarray_signaling_genes new.grps 0.4761905
## 
## $acc.overall$kendall
##     model                            x_prefix y_prefix    values
## 10   ksvm                microarray_all_genes new.grps 0.4285714
## 73   ksvm       microarray_effector_path_vals new.grps 0.4523810
## 136  ksvm microarray_metabolic_mod_activities new.grps 0.4285714
## 199  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 262  ksvm microarray_metabolic_mod_nodevalues new.grps 0.1904762
## 325  ksvm                microarray_path_vals new.grps 0.4523810
## 388  ksvm          microarray_signaling_genes new.grps 0.4523810
## 
## $acc.overall$linear
##     model                            x_prefix y_prefix    values
## 19   ksvm                microarray_all_genes new.grps 0.4285714
## 82   ksvm       microarray_effector_path_vals new.grps 0.4285714
## 145  ksvm microarray_metabolic_mod_activities new.grps 0.3333333
## 208  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 271  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4285714
## 334  ksvm                microarray_path_vals new.grps 0.2142857
## 397  ksvm          microarray_signaling_genes new.grps 0.3333333
## 
## $acc.overall$`rbf+kendall`
##     model                            x_prefix y_prefix    values
## 28   ksvm                microarray_all_genes new.grps 0.4047619
## 91   ksvm       microarray_effector_path_vals new.grps 0.4285714
## 154  ksvm microarray_metabolic_mod_activities new.grps 0.4761905
## 217  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 280  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4761905
## 343  ksvm                microarray_path_vals new.grps 0.4285714
## 406  ksvm          microarray_signaling_genes new.grps 0.4761905
## 
## $acc.overall$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 37   ksvm                microarray_all_genes new.grps 0.4285714
## 100  ksvm       microarray_effector_path_vals new.grps 0.4285714
## 163  ksvm microarray_metabolic_mod_activities new.grps 0.4285714
## 226  ksvm microarray_metabolic_mod_genevalues new.grps 0.3095238
## 289  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4285714
## 352  ksvm                microarray_path_vals new.grps 0.4285714
## 415  ksvm          microarray_signaling_genes new.grps 0.3809524
## 
## $acc.overall$`linear+rbf`
##     model                            x_prefix y_prefix    values
## 46   ksvm                microarray_all_genes new.grps 0.3333333
## 109  ksvm       microarray_effector_path_vals new.grps 0.4285714
## 172  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 235  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 298  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4285714
## 361  ksvm                microarray_path_vals new.grps 0.4285714
## 424  ksvm          microarray_signaling_genes new.grps 0.4285714
## 
## $acc.overall$`linear+rbf+kendall`
##     model                            x_prefix y_prefix    values
## 55   ksvm                microarray_all_genes new.grps 0.4047619
## 118  ksvm       microarray_effector_path_vals new.grps 0.4285714
## 181  ksvm microarray_metabolic_mod_activities new.grps 0.4761905
## 244  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 307  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4285714
## 370  ksvm                microarray_path_vals new.grps 0.4047619
## 433  ksvm          microarray_signaling_genes new.grps 0.4285714
## 
## 
## $acc.unknown
## $acc.unknown$rbf
##     model                            x_prefix y_prefix    values
## 5    ksvm                microarray_all_genes new.grps 0.8888889
## 68   ksvm       microarray_effector_path_vals new.grps 1.0000000
## 131  ksvm microarray_metabolic_mod_activities new.grps 0.9444444
## 194  ksvm microarray_metabolic_mod_genevalues new.grps 0.2222222
## 257  ksvm microarray_metabolic_mod_nodevalues new.grps 0.8888889
## 320  ksvm                microarray_path_vals new.grps 1.0000000
## 383  ksvm          microarray_signaling_genes new.grps 0.8888889
## 
## $acc.unknown$kendall
##     model                            x_prefix y_prefix    values
## 14   ksvm                microarray_all_genes new.grps 1.0000000
## 77   ksvm       microarray_effector_path_vals new.grps 1.0000000
## 140  ksvm microarray_metabolic_mod_activities new.grps 1.0000000
## 203  ksvm microarray_metabolic_mod_genevalues new.grps 1.0000000
## 266  ksvm microarray_metabolic_mod_nodevalues new.grps 0.1666667
## 329  ksvm                microarray_path_vals new.grps 1.0000000
## 392  ksvm          microarray_signaling_genes new.grps 1.0000000
## 
## $acc.unknown$linear
##     model                            x_prefix y_prefix    values
## 23   ksvm                microarray_all_genes new.grps 1.0000000
## 86   ksvm       microarray_effector_path_vals new.grps 1.0000000
## 149  ksvm microarray_metabolic_mod_activities new.grps 0.6111111
## 212  ksvm microarray_metabolic_mod_genevalues new.grps 1.0000000
## 275  ksvm microarray_metabolic_mod_nodevalues new.grps 1.0000000
## 338  ksvm                microarray_path_vals new.grps 0.1666667
## 401  ksvm          microarray_signaling_genes new.grps 0.6111111
## 
## $acc.unknown$`rbf+kendall`
##     model                            x_prefix y_prefix    values
## 32   ksvm                microarray_all_genes new.grps 0.8888889
## 95   ksvm       microarray_effector_path_vals new.grps 1.0000000
## 158  ksvm microarray_metabolic_mod_activities new.grps 0.9444444
## 221  ksvm microarray_metabolic_mod_genevalues new.grps 1.0000000
## 284  ksvm microarray_metabolic_mod_nodevalues new.grps 0.9444444
## 347  ksvm                microarray_path_vals new.grps 1.0000000
## 410  ksvm          microarray_signaling_genes new.grps 0.9444444
## 
## $acc.unknown$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 41   ksvm                microarray_all_genes new.grps 1.0000000
## 104  ksvm       microarray_effector_path_vals new.grps 1.0000000
## 167  ksvm microarray_metabolic_mod_activities new.grps 1.0000000
## 230  ksvm microarray_metabolic_mod_genevalues new.grps 0.3888889
## 293  ksvm microarray_metabolic_mod_nodevalues new.grps 1.0000000
## 356  ksvm                microarray_path_vals new.grps 1.0000000
## 419  ksvm          microarray_signaling_genes new.grps 0.8333333
## 
## $acc.unknown$`linear+rbf`
##     model                            x_prefix y_prefix    values
## 50   ksvm                microarray_all_genes new.grps 0.7222222
## 113  ksvm       microarray_effector_path_vals new.grps 1.0000000
## 176  ksvm microarray_metabolic_mod_activities new.grps 1.0000000
## 239  ksvm microarray_metabolic_mod_genevalues new.grps 1.0000000
## 302  ksvm microarray_metabolic_mod_nodevalues new.grps 0.8333333
## 365  ksvm                microarray_path_vals new.grps 0.7777778
## 428  ksvm          microarray_signaling_genes new.grps 1.0000000
## 
## $acc.unknown$`linear+rbf+kendall`
##     model                            x_prefix y_prefix    values
## 59   ksvm                microarray_all_genes new.grps 0.8888889
## 122  ksvm       microarray_effector_path_vals new.grps 1.0000000
## 185  ksvm microarray_metabolic_mod_activities new.grps 0.9444444
## 248  ksvm microarray_metabolic_mod_genevalues new.grps 1.0000000
## 311  ksvm microarray_metabolic_mod_nodevalues new.grps 0.8333333
## 374  ksvm                microarray_path_vals new.grps 0.7777778
## 437  ksvm          microarray_signaling_genes new.grps 0.8888889
## 
## 
## $precision.Control
## $precision.Control$rbf
##     model                            x_prefix y_prefix    values
## 6    ksvm                microarray_all_genes new.grps 0.2222222
## 69   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 132  ksvm microarray_metabolic_mod_activities new.grps 0.3750000
## 195  ksvm microarray_metabolic_mod_genevalues new.grps 0.5714286
## 258  ksvm microarray_metabolic_mod_nodevalues new.grps 0.2727273
## 321  ksvm                microarray_path_vals new.grps 0.0000000
## 384  ksvm          microarray_signaling_genes new.grps 0.3636364
## 
## $precision.Control$kendall
##     model                            x_prefix y_prefix    values
## 15   ksvm                microarray_all_genes new.grps 0.0000000
## 78   ksvm       microarray_effector_path_vals new.grps 0.5000000
## 141  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 204  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 267  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4166667
## 330  ksvm                microarray_path_vals new.grps 0.5000000
## 393  ksvm          microarray_signaling_genes new.grps 0.3333333
## 
## $precision.Control$linear
##     model                            x_prefix y_prefix    values
## 24   ksvm                microarray_all_genes new.grps 0.0000000
## 87   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 150  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 213  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 276  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 339  ksvm                microarray_path_vals new.grps 0.2500000
## 402  ksvm          microarray_signaling_genes new.grps 0.2142857
## 
## $precision.Control$`rbf+kendall`
##     model                            x_prefix y_prefix    values
## 33   ksvm                microarray_all_genes new.grps 0.2500000
## 96   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 159  ksvm microarray_metabolic_mod_activities new.grps 0.4285714
## 222  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 285  ksvm microarray_metabolic_mod_nodevalues new.grps 0.3333333
## 348  ksvm                microarray_path_vals new.grps 0.0000000
## 411  ksvm          microarray_signaling_genes new.grps 0.3333333
## 
## $precision.Control$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 42   ksvm                microarray_all_genes new.grps 0.0000000
## 105  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 168  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 231  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 294  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 357  ksvm                microarray_path_vals new.grps 0.0000000
## 420  ksvm          microarray_signaling_genes new.grps 0.1666667
## 
## $precision.Control$`linear+rbf`
##     model                            x_prefix y_prefix    values
## 51   ksvm                microarray_all_genes new.grps 0.1666667
## 114  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 177  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 240  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 303  ksvm microarray_metabolic_mod_nodevalues new.grps 0.2727273
## 366  ksvm                microarray_path_vals new.grps 0.2857143
## 429  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $precision.Control$`linear+rbf+kendall`
##     model                            x_prefix y_prefix    values
## 60   ksvm                microarray_all_genes new.grps 0.2000000
## 123  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 186  ksvm microarray_metabolic_mod_activities new.grps 0.4285714
## 249  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 312  ksvm microarray_metabolic_mod_nodevalues new.grps 0.3000000
## 375  ksvm                microarray_path_vals new.grps 0.2500000
## 438  ksvm          microarray_signaling_genes new.grps 0.2500000
## 
## 
## $precision.Effect_1
## $precision.Effect_1$rbf
##     model                            x_prefix y_prefix    values
## 7    ksvm                microarray_all_genes new.grps 0.0000000
## 70   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 133  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 196  ksvm microarray_metabolic_mod_genevalues new.grps 0.2857143
## 259  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 322  ksvm                microarray_path_vals new.grps 0.0000000
## 385  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $precision.Effect_1$kendall
##     model                            x_prefix y_prefix values
## 16   ksvm                microarray_all_genes new.grps      0
## 79   ksvm       microarray_effector_path_vals new.grps      0
## 142  ksvm microarray_metabolic_mod_activities new.grps      0
## 205  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 268  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 331  ksvm                microarray_path_vals new.grps      0
## 394  ksvm          microarray_signaling_genes new.grps      0
## 
## $precision.Effect_1$linear
##     model                            x_prefix y_prefix    values
## 25   ksvm                microarray_all_genes new.grps 0.0000000
## 88   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 151  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 214  ksvm microarray_metabolic_mod_genevalues new.grps 0.0000000
## 277  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 340  ksvm                microarray_path_vals new.grps 0.3333333
## 403  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $precision.Effect_1$`rbf+kendall`
##     model                            x_prefix y_prefix values
## 34   ksvm                microarray_all_genes new.grps      0
## 97   ksvm       microarray_effector_path_vals new.grps      0
## 160  ksvm microarray_metabolic_mod_activities new.grps      0
## 223  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 286  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 349  ksvm                microarray_path_vals new.grps      0
## 412  ksvm          microarray_signaling_genes new.grps      0
## 
## $precision.Effect_1$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 43   ksvm                microarray_all_genes new.grps 0.0000000
## 106  ksvm       microarray_effector_path_vals new.grps 0.0000000
## 169  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 232  ksvm microarray_metabolic_mod_genevalues new.grps 0.1666667
## 295  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 358  ksvm                microarray_path_vals new.grps 0.0000000
## 421  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $precision.Effect_1$`linear+rbf`
##     model                            x_prefix y_prefix values
## 52   ksvm                microarray_all_genes new.grps      0
## 115  ksvm       microarray_effector_path_vals new.grps      0
## 178  ksvm microarray_metabolic_mod_activities new.grps      0
## 241  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 304  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 367  ksvm                microarray_path_vals new.grps      0
## 430  ksvm          microarray_signaling_genes new.grps      0
## 
## $precision.Effect_1$`linear+rbf+kendall`
##     model                            x_prefix y_prefix values
## 61   ksvm                microarray_all_genes new.grps      0
## 124  ksvm       microarray_effector_path_vals new.grps      0
## 187  ksvm microarray_metabolic_mod_activities new.grps      0
## 250  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 313  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 376  ksvm                microarray_path_vals new.grps      0
## 439  ksvm          microarray_signaling_genes new.grps      0
## 
## 
## $precision.Effect_2
## $precision.Effect_2$rbf
##     model                            x_prefix y_prefix    values
## 8    ksvm                microarray_all_genes new.grps 0.0000000
## 71   ksvm       microarray_effector_path_vals new.grps 0.0000000
## 134  ksvm microarray_metabolic_mod_activities new.grps 0.0000000
## 197  ksvm microarray_metabolic_mod_genevalues new.grps 0.3333333
## 260  ksvm microarray_metabolic_mod_nodevalues new.grps 0.0000000
## 323  ksvm                microarray_path_vals new.grps 0.0000000
## 386  ksvm          microarray_signaling_genes new.grps 0.0000000
## 
## $precision.Effect_2$kendall
##     model                            x_prefix y_prefix values
## 17   ksvm                microarray_all_genes new.grps      0
## 80   ksvm       microarray_effector_path_vals new.grps      0
## 143  ksvm microarray_metabolic_mod_activities new.grps      0
## 206  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 269  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 332  ksvm                microarray_path_vals new.grps      0
## 395  ksvm          microarray_signaling_genes new.grps      0
## 
## $precision.Effect_2$linear
##     model                            x_prefix y_prefix values
## 26   ksvm                microarray_all_genes new.grps    0.0
## 89   ksvm       microarray_effector_path_vals new.grps    0.0
## 152  ksvm microarray_metabolic_mod_activities new.grps    0.0
## 215  ksvm microarray_metabolic_mod_genevalues new.grps    0.0
## 278  ksvm microarray_metabolic_mod_nodevalues new.grps    0.0
## 341  ksvm                microarray_path_vals new.grps    0.5
## 404  ksvm          microarray_signaling_genes new.grps    0.0
## 
## $precision.Effect_2$`rbf+kendall`
##     model                            x_prefix y_prefix values
## 35   ksvm                microarray_all_genes new.grps      0
## 98   ksvm       microarray_effector_path_vals new.grps      0
## 161  ksvm microarray_metabolic_mod_activities new.grps      0
## 224  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 287  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 350  ksvm                microarray_path_vals new.grps      0
## 413  ksvm          microarray_signaling_genes new.grps      0
## 
## $precision.Effect_2$`linear+kendall`
##     model                            x_prefix y_prefix values
## 44   ksvm                microarray_all_genes new.grps    0.0
## 107  ksvm       microarray_effector_path_vals new.grps    0.0
## 170  ksvm microarray_metabolic_mod_activities new.grps    0.0
## 233  ksvm microarray_metabolic_mod_genevalues new.grps    0.5
## 296  ksvm microarray_metabolic_mod_nodevalues new.grps    0.0
## 359  ksvm                microarray_path_vals new.grps    0.0
## 422  ksvm          microarray_signaling_genes new.grps    0.0
## 
## $precision.Effect_2$`linear+rbf`
##     model                            x_prefix y_prefix values
## 53   ksvm                microarray_all_genes new.grps      0
## 116  ksvm       microarray_effector_path_vals new.grps      0
## 179  ksvm microarray_metabolic_mod_activities new.grps      0
## 242  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 305  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 368  ksvm                microarray_path_vals new.grps      0
## 431  ksvm          microarray_signaling_genes new.grps      0
## 
## $precision.Effect_2$`linear+rbf+kendall`
##     model                            x_prefix y_prefix values
## 62   ksvm                microarray_all_genes new.grps      0
## 125  ksvm       microarray_effector_path_vals new.grps      0
## 188  ksvm microarray_metabolic_mod_activities new.grps      0
## 251  ksvm microarray_metabolic_mod_genevalues new.grps      0
## 314  ksvm microarray_metabolic_mod_nodevalues new.grps      0
## 377  ksvm                microarray_path_vals new.grps      0
## 440  ksvm          microarray_signaling_genes new.grps      0
## 
## 
## $precision.unknown
## $precision.unknown$rbf
##     model                            x_prefix y_prefix    values
## 9    ksvm                microarray_all_genes new.grps 0.5000000
## 72   ksvm       microarray_effector_path_vals new.grps 0.4285714
## 135  ksvm microarray_metabolic_mod_activities new.grps 0.5483871
## 198  ksvm microarray_metabolic_mod_genevalues new.grps 0.5000000
## 261  ksvm microarray_metabolic_mod_nodevalues new.grps 0.5161290
## 324  ksvm                microarray_path_vals new.grps 0.4285714
## 387  ksvm          microarray_signaling_genes new.grps 0.5161290
## 
## $precision.unknown$kendall
##     model                            x_prefix y_prefix    values
## 18   ksvm                microarray_all_genes new.grps 0.4285714
## 81   ksvm       microarray_effector_path_vals new.grps 0.4500000
## 144  ksvm microarray_metabolic_mod_activities new.grps 0.4285714
## 207  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 270  ksvm microarray_metabolic_mod_nodevalues new.grps 0.2142857
## 333  ksvm                microarray_path_vals new.grps 0.4500000
## 396  ksvm          microarray_signaling_genes new.grps 0.4615385
## 
## $precision.unknown$linear
##     model                            x_prefix y_prefix    values
## 27   ksvm                microarray_all_genes new.grps 0.4285714
## 90   ksvm       microarray_effector_path_vals new.grps 0.4285714
## 153  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 216  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 279  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4285714
## 342  ksvm                microarray_path_vals new.grps 0.7500000
## 405  ksvm          microarray_signaling_genes new.grps 0.4074074
## 
## $precision.unknown$`rbf+kendall`
##     model                            x_prefix y_prefix    values
## 36   ksvm                microarray_all_genes new.grps 0.4324324
## 99   ksvm       microarray_effector_path_vals new.grps 0.4285714
## 162  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 225  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 288  ksvm microarray_metabolic_mod_nodevalues new.grps 0.5151515
## 351  ksvm                microarray_path_vals new.grps 0.4285714
## 414  ksvm          microarray_signaling_genes new.grps 0.5151515
## 
## $precision.unknown$`linear+kendall`
##     model                            x_prefix y_prefix    values
## 45   ksvm                microarray_all_genes new.grps 0.4285714
## 108  ksvm       microarray_effector_path_vals new.grps 0.4285714
## 171  ksvm microarray_metabolic_mod_activities new.grps 0.4285714
## 234  ksvm microarray_metabolic_mod_genevalues new.grps 0.7000000
## 297  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4285714
## 360  ksvm                microarray_path_vals new.grps 0.4285714
## 423  ksvm          microarray_signaling_genes new.grps 0.4411765
## 
## $precision.unknown$`linear+rbf`
##     model                            x_prefix y_prefix    values
## 54   ksvm                microarray_all_genes new.grps 0.4062500
## 117  ksvm       microarray_effector_path_vals new.grps 0.4285714
## 180  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 243  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 306  ksvm microarray_metabolic_mod_nodevalues new.grps 0.5000000
## 369  ksvm                microarray_path_vals new.grps 0.5000000
## 432  ksvm          microarray_signaling_genes new.grps 0.4285714
## 
## $precision.unknown$`linear+rbf+kendall`
##     model                            x_prefix y_prefix    values
## 63   ksvm                microarray_all_genes new.grps 0.4444444
## 126  ksvm       microarray_effector_path_vals new.grps 0.4285714
## 189  ksvm microarray_metabolic_mod_activities new.grps 0.5000000
## 252  ksvm microarray_metabolic_mod_genevalues new.grps 0.4285714
## 315  ksvm microarray_metabolic_mod_nodevalues new.grps 0.4838710
## 378  ksvm                microarray_path_vals new.grps 0.4666667
## 441  ksvm          microarray_signaling_genes new.grps 0.4848485
```

```r
# plot
for (yname in ylist) {
  dds <- subset(dd, dd$y_prefix == yname & grepl("(?=.*acc.*)(?=.*overall.*)", dd$scores, ignore.case = TRUE, perl = TRUE))
  p1 <- ggplot(dds, aes(x = x_prefix, y = values)) + 
    geom_bar(aes(fill = x_prefix), stat="identity", position = "dodge") + 
    facet_wrap(~ kf) + ylim(0,1) + 
    ggtitle(paste0('Overall acc for ', yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
}
```

![plot of chunk svm_plot](ma_figure/svm_plot-1.pdf)

# Session info


```r
devtools::session_info()
```

```
## Session info --------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.2.3 (2015-12-10)
##  system   x86_64, darwin13.4.0        
##  ui       X11                         
##  language (EN)                        
##  collate  C                           
##  tz       Europe/Paris                
##  date     2016-03-21
```

```
## Packages ------------------------------------------------------------------
```

```
##  package      * version date       source        
##  MASS           7.3-45  2015-11-10 CRAN (R 3.2.3)
##  Matrix         1.2-3   2015-11-28 CRAN (R 3.2.3)
##  MatrixModels   0.4-1   2015-08-22 CRAN (R 3.2.0)
##  Rcpp           0.12.3  2016-01-10 CRAN (R 3.2.3)
##  SparseM        1.7     2015-08-15 CRAN (R 3.2.0)
##  car            2.1-1   2015-12-14 CRAN (R 3.2.3)
##  caret        * 6.0-64  2016-01-06 CRAN (R 3.2.3)
##  codetools      0.2-14  2015-07-15 CRAN (R 3.2.3)
##  colorspace     1.2-6   2015-03-11 CRAN (R 3.2.0)
##  devtools       1.10.0  2016-01-23 CRAN (R 3.2.3)
##  digest         0.6.9   2016-01-08 CRAN (R 3.2.3)
##  evaluate       0.8     2015-09-18 CRAN (R 3.2.0)
##  foreach        1.4.3   2015-10-13 CRAN (R 3.2.0)
##  formatR        1.2.1   2015-09-18 CRAN (R 3.2.0)
##  ggplot2      * 2.1.0   2016-03-01 CRAN (R 3.2.3)
##  gtable         0.2.0   2016-02-26 CRAN (R 3.2.3)
##  iterators      1.0.8   2015-10-13 CRAN (R 3.2.0)
##  kernlab      * 0.9-23  2016-01-26 CRAN (R 3.2.3)
##  knitr          1.12.3  2016-01-22 CRAN (R 3.2.3)
##  labeling       0.3     2014-08-23 CRAN (R 3.1.2)
##  lattice      * 0.20-33 2015-07-14 CRAN (R 3.2.3)
##  lme4           1.1-11  2016-02-12 CRAN (R 3.2.3)
##  magrittr       1.5     2014-11-22 CRAN (R 3.2.0)
##  memoise        1.0.0   2016-01-29 CRAN (R 3.2.3)
##  mgcv           1.8-11  2016-01-24 CRAN (R 3.2.3)
##  minqa          1.2.4   2014-10-09 CRAN (R 3.1.2)
##  munsell        0.4.3   2016-02-13 CRAN (R 3.2.3)
##  mvtnorm        1.0-5   2016-02-02 CRAN (R 3.2.3)
##  nlme           3.1-125 2016-02-27 CRAN (R 3.2.3)
##  nloptr         1.0.4   2014-08-04 CRAN (R 3.1.2)
##  nnet           7.3-12  2016-02-02 CRAN (R 3.2.3)
##  pbkrtest       0.4-6   2016-01-27 CRAN (R 3.2.3)
##  pcaPP        * 1.9-60  2014-10-22 CRAN (R 3.1.2)
##  plyr           1.8.3   2015-06-12 CRAN (R 3.2.0)
##  quantreg       5.21    2016-02-13 CRAN (R 3.2.3)
##  reshape2       1.4.1   2014-12-06 CRAN (R 3.1.2)
##  rstudioapi     0.5     2016-01-24 CRAN (R 3.2.3)
##  scales         0.4.0   2016-02-26 CRAN (R 3.2.3)
##  stringi        1.0-1   2015-10-22 CRAN (R 3.2.0)
##  stringr        1.0.0   2015-04-30 CRAN (R 3.2.0)
```

