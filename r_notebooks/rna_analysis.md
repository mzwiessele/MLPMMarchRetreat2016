# Notebook : RNAseq analysis


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, fig.width = 12, fig.height = 8, dev = "pdf", fig.keep = "high", fig.path = "rna_figure/", cache.path = "rna_cache/")
set.seed(35875954)

library(caret)
library(pcaPP)
library(kernlab)
source('predictors.R')

datapath <- '../Data/RNASeq/'
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
## ==============>  RNAseq_all_genes <==============
## [1]   104 11216
## ==============>  rnaseq_design <==============
## [1] 104   4
## ==============>  rnaseq_effector_path_vals <==============
## [1]  104 1044
## ==============>  rnaseq_metabolic_mod_activities <==============
## [1] 104  89
## ==============>  rnaseq_metabolic_mod_genevalues <==============
## [1] 104 172
## ==============>  rnaseq_metabolic_mod_nodevalues <==============
## [1] 104 462
## ==============>  rnaseq_path_vals <==============
## [1]  104 4457
## ==============>  rnaseq_signaling_genes <==============
## [1]  104 2184
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
##     TEST SET TRAINING SET 
##           42           62
```

```r
# generate labels
effect.grps <- as.factor(get(designobj)[,grep('effect', colnames(get(designobj)), ignore.case = TRUE)])
table(effect.grps[trainidx])
```

```
## 
##  Control Effect_1 Effect_2 Effect_3 Effect_4 Effect_5 Effect_6 Effect_7 
##       18        9        9        9        8        9        0        0
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
##       13        6        6        6        6        6       19
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
##     model                        x_prefix y_prefix    values
## 2    ksvm                RNAseq_all_genes new.grps 1.0000000
## 65   ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 128  ksvm rnaseq_metabolic_mod_activities new.grps 0.8333333
## 191  ksvm rnaseq_metabolic_mod_genevalues new.grps 1.0000000
## 254  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 317  ksvm                rnaseq_path_vals new.grps 0.8333333
## 380  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $acc.Control$kendall
##     model                        x_prefix y_prefix    values
## 11   ksvm                RNAseq_all_genes new.grps 1.0000000
## 74   ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 137  ksvm rnaseq_metabolic_mod_activities new.grps 1.0000000
## 200  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8333333
## 263  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 326  ksvm                rnaseq_path_vals new.grps 1.0000000
## 389  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $acc.Control$linear
##     model                        x_prefix y_prefix    values
## 20   ksvm                RNAseq_all_genes new.grps 1.0000000
## 83   ksvm       rnaseq_effector_path_vals new.grps 0.8333333
## 146  ksvm rnaseq_metabolic_mod_activities new.grps 0.6666667
## 209  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8333333
## 272  ksvm rnaseq_metabolic_mod_nodevalues new.grps 1.0000000
## 335  ksvm                rnaseq_path_vals new.grps 0.5000000
## 398  ksvm          rnaseq_signaling_genes new.grps 0.8333333
## 
## $acc.Control$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 29   ksvm                RNAseq_all_genes new.grps 1.0000000
## 92   ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 155  ksvm rnaseq_metabolic_mod_activities new.grps 0.8333333
## 218  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.6666667
## 281  ksvm rnaseq_metabolic_mod_nodevalues new.grps 1.0000000
## 344  ksvm                rnaseq_path_vals new.grps 1.0000000
## 407  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $acc.Control$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 38   ksvm                RNAseq_all_genes new.grps 1.0000000
## 101  ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 164  ksvm rnaseq_metabolic_mod_activities new.grps 0.8333333
## 227  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8333333
## 290  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 353  ksvm                rnaseq_path_vals new.grps 1.0000000
## 416  ksvm          rnaseq_signaling_genes new.grps 0.6666667
## 
## $acc.Control$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 47   ksvm                RNAseq_all_genes new.grps 0.8333333
## 110  ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 173  ksvm rnaseq_metabolic_mod_activities new.grps 0.5000000
## 236  ksvm rnaseq_metabolic_mod_genevalues new.grps 1.0000000
## 299  ksvm rnaseq_metabolic_mod_nodevalues new.grps 1.0000000
## 362  ksvm                rnaseq_path_vals new.grps 0.8333333
## 425  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $acc.Control$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 56   ksvm                RNAseq_all_genes new.grps 1.0000000
## 119  ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 182  ksvm rnaseq_metabolic_mod_activities new.grps 0.8333333
## 245  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8333333
## 308  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 371  ksvm                rnaseq_path_vals new.grps 1.0000000
## 434  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## 
## $acc.Effect_1
## $acc.Effect_1$rbf
##     model                        x_prefix y_prefix    values
## 3    ksvm                RNAseq_all_genes new.grps 0.7777778
## 66   ksvm       rnaseq_effector_path_vals new.grps 0.4444444
## 129  ksvm rnaseq_metabolic_mod_activities new.grps 0.1111111
## 192  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.4444444
## 255  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5555556
## 318  ksvm                rnaseq_path_vals new.grps 0.5555556
## 381  ksvm          rnaseq_signaling_genes new.grps 0.8888889
## 
## $acc.Effect_1$kendall
##     model                        x_prefix y_prefix    values
## 12   ksvm                RNAseq_all_genes new.grps 0.4444444
## 75   ksvm       rnaseq_effector_path_vals new.grps 0.7777778
## 138  ksvm rnaseq_metabolic_mod_activities new.grps 0.6666667
## 201  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.2222222
## 264  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.0000000
## 327  ksvm                rnaseq_path_vals new.grps 0.4444444
## 390  ksvm          rnaseq_signaling_genes new.grps 0.4444444
## 
## $acc.Effect_1$linear
##     model                        x_prefix y_prefix    values
## 21   ksvm                RNAseq_all_genes new.grps 0.8888889
## 84   ksvm       rnaseq_effector_path_vals new.grps 0.7777778
## 147  ksvm rnaseq_metabolic_mod_activities new.grps 0.4444444
## 210  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.2222222
## 273  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.1111111
## 336  ksvm                rnaseq_path_vals new.grps 0.8888889
## 399  ksvm          rnaseq_signaling_genes new.grps 0.8888889
## 
## $acc.Effect_1$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 30   ksvm                RNAseq_all_genes new.grps 0.7777778
## 93   ksvm       rnaseq_effector_path_vals new.grps 0.6666667
## 156  ksvm rnaseq_metabolic_mod_activities new.grps 0.3333333
## 219  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.2222222
## 282  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.0000000
## 345  ksvm                rnaseq_path_vals new.grps 0.4444444
## 408  ksvm          rnaseq_signaling_genes new.grps 0.5555556
## 
## $acc.Effect_1$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 39   ksvm                RNAseq_all_genes new.grps 0.5555556
## 102  ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 165  ksvm rnaseq_metabolic_mod_activities new.grps 0.8888889
## 228  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5555556
## 291  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5555556
## 354  ksvm                rnaseq_path_vals new.grps 0.8888889
## 417  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $acc.Effect_1$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 48   ksvm                RNAseq_all_genes new.grps 0.8888889
## 111  ksvm       rnaseq_effector_path_vals new.grps 0.5555556
## 174  ksvm rnaseq_metabolic_mod_activities new.grps 0.2222222
## 237  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 300  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.7777778
## 363  ksvm                rnaseq_path_vals new.grps 0.2222222
## 426  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $acc.Effect_1$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 57   ksvm                RNAseq_all_genes new.grps 0.6666667
## 120  ksvm       rnaseq_effector_path_vals new.grps 0.4444444
## 183  ksvm rnaseq_metabolic_mod_activities new.grps 0.2222222
## 246  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5555556
## 309  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5555556
## 372  ksvm                rnaseq_path_vals new.grps 0.5555556
## 435  ksvm          rnaseq_signaling_genes new.grps 0.8888889
## 
## 
## $acc.Effect_2
## $acc.Effect_2$rbf
##     model                        x_prefix y_prefix    values
## 4    ksvm                RNAseq_all_genes new.grps 0.1111111
## 67   ksvm       rnaseq_effector_path_vals new.grps 0.1111111
## 130  ksvm rnaseq_metabolic_mod_activities new.grps 0.1111111
## 193  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 256  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.1111111
## 319  ksvm                rnaseq_path_vals new.grps 0.0000000
## 382  ksvm          rnaseq_signaling_genes new.grps 0.1111111
## 
## $acc.Effect_2$kendall
##     model                        x_prefix y_prefix    values
## 13   ksvm                RNAseq_all_genes new.grps 0.0000000
## 76   ksvm       rnaseq_effector_path_vals new.grps 0.1111111
## 139  ksvm rnaseq_metabolic_mod_activities new.grps 0.1111111
## 202  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 265  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.2222222
## 328  ksvm                rnaseq_path_vals new.grps 0.1111111
## 391  ksvm          rnaseq_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_2$linear
##     model                        x_prefix y_prefix    values
## 22   ksvm                RNAseq_all_genes new.grps 0.2222222
## 85   ksvm       rnaseq_effector_path_vals new.grps 0.1111111
## 148  ksvm rnaseq_metabolic_mod_activities new.grps 0.1111111
## 211  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 274  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.1111111
## 337  ksvm                rnaseq_path_vals new.grps 0.1111111
## 400  ksvm          rnaseq_signaling_genes new.grps 0.1111111
## 
## $acc.Effect_2$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 31   ksvm                RNAseq_all_genes new.grps 0.0000000
## 94   ksvm       rnaseq_effector_path_vals new.grps 0.1111111
## 157  ksvm rnaseq_metabolic_mod_activities new.grps 0.1111111
## 220  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 283  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.1111111
## 346  ksvm                rnaseq_path_vals new.grps 0.1111111
## 409  ksvm          rnaseq_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_2$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 40   ksvm                RNAseq_all_genes new.grps 0.0000000
## 103  ksvm       rnaseq_effector_path_vals new.grps 0.1111111
## 166  ksvm rnaseq_metabolic_mod_activities new.grps 0.1111111
## 229  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 292  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.1111111
## 355  ksvm                rnaseq_path_vals new.grps 0.0000000
## 418  ksvm          rnaseq_signaling_genes new.grps 0.0000000
## 
## $acc.Effect_2$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 49   ksvm                RNAseq_all_genes new.grps 0.2222222
## 112  ksvm       rnaseq_effector_path_vals new.grps 0.2222222
## 175  ksvm rnaseq_metabolic_mod_activities new.grps 0.2222222
## 238  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 301  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.1111111
## 364  ksvm                rnaseq_path_vals new.grps 0.1111111
## 427  ksvm          rnaseq_signaling_genes new.grps 0.1111111
## 
## $acc.Effect_2$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 58   ksvm                RNAseq_all_genes new.grps 0.0000000
## 121  ksvm       rnaseq_effector_path_vals new.grps 0.1111111
## 184  ksvm rnaseq_metabolic_mod_activities new.grps 0.1111111
## 247  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 310  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.0000000
## 373  ksvm                rnaseq_path_vals new.grps 0.1111111
## 436  ksvm          rnaseq_signaling_genes new.grps 0.1111111
## 
## 
## $acc.overall
## $acc.overall$rbf
##     model                        x_prefix y_prefix    values
## 1    ksvm                RNAseq_all_genes new.grps 0.7142857
## 64   ksvm       rnaseq_effector_path_vals new.grps 0.5952381
## 127  ksvm rnaseq_metabolic_mod_activities new.grps 0.3571429
## 190  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.4523810
## 253  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5000000
## 316  ksvm                rnaseq_path_vals new.grps 0.5000000
## 379  ksvm          rnaseq_signaling_genes new.grps 0.7142857
## 
## $acc.overall$kendall
##     model                        x_prefix y_prefix    values
## 10   ksvm                RNAseq_all_genes new.grps 0.6190476
## 73   ksvm       rnaseq_effector_path_vals new.grps 0.6666667
## 136  ksvm rnaseq_metabolic_mod_activities new.grps 0.6190476
## 199  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5476190
## 262  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5952381
## 325  ksvm                rnaseq_path_vals new.grps 0.5714286
## 388  ksvm          rnaseq_signaling_genes new.grps 0.6190476
## 
## $acc.overall$linear
##     model                        x_prefix y_prefix    values
## 19   ksvm                RNAseq_all_genes new.grps 0.6666667
## 82   ksvm       rnaseq_effector_path_vals new.grps 0.5476190
## 145  ksvm rnaseq_metabolic_mod_activities new.grps 0.4523810
## 208  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.3571429
## 271  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5476190
## 334  ksvm                rnaseq_path_vals new.grps 0.4523810
## 397  ksvm          rnaseq_signaling_genes new.grps 0.6904762
## 
## $acc.overall$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 28   ksvm                RNAseq_all_genes new.grps 0.6666667
## 91   ksvm       rnaseq_effector_path_vals new.grps 0.6190476
## 154  ksvm rnaseq_metabolic_mod_activities new.grps 0.4047619
## 217  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5238095
## 280  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5238095
## 343  ksvm                rnaseq_path_vals new.grps 0.5714286
## 406  ksvm          rnaseq_signaling_genes new.grps 0.6428571
## 
## $acc.overall$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 37   ksvm                RNAseq_all_genes new.grps 0.6190476
## 100  ksvm       rnaseq_effector_path_vals new.grps 0.6190476
## 163  ksvm rnaseq_metabolic_mod_activities new.grps 0.4761905
## 226  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5952381
## 289  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5952381
## 352  ksvm                rnaseq_path_vals new.grps 0.5952381
## 415  ksvm          rnaseq_signaling_genes new.grps 0.6428571
## 
## $acc.overall$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 46   ksvm                RNAseq_all_genes new.grps 0.7142857
## 109  ksvm       rnaseq_effector_path_vals new.grps 0.5714286
## 172  ksvm rnaseq_metabolic_mod_activities new.grps 0.3333333
## 235  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5000000
## 298  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5952381
## 361  ksvm                rnaseq_path_vals new.grps 0.3809524
## 424  ksvm          rnaseq_signaling_genes new.grps 0.6904762
## 
## $acc.overall$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 55   ksvm                RNAseq_all_genes new.grps 0.6428571
## 118  ksvm       rnaseq_effector_path_vals new.grps 0.5952381
## 181  ksvm rnaseq_metabolic_mod_activities new.grps 0.4523810
## 244  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5714286
## 307  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.4523810
## 370  ksvm                rnaseq_path_vals new.grps 0.5000000
## 433  ksvm          rnaseq_signaling_genes new.grps 0.7142857
## 
## 
## $acc.unknown
## $acc.unknown$rbf
##     model                        x_prefix y_prefix    values
## 5    ksvm                RNAseq_all_genes new.grps 0.8888889
## 68   ksvm       rnaseq_effector_path_vals new.grps 0.7777778
## 131  ksvm rnaseq_metabolic_mod_activities new.grps 0.4444444
## 194  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5000000
## 257  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5555556
## 320  ksvm                rnaseq_path_vals new.grps 0.6111111
## 383  ksvm          rnaseq_signaling_genes new.grps 0.8333333
## 
## $acc.unknown$kendall
##     model                        x_prefix y_prefix    values
## 14   ksvm                RNAseq_all_genes new.grps 0.8888889
## 77   ksvm       rnaseq_effector_path_vals new.grps 0.7777778
## 140  ksvm rnaseq_metabolic_mod_activities new.grps 0.7222222
## 203  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8888889
## 266  ksvm rnaseq_metabolic_mod_nodevalues new.grps 1.0000000
## 329  ksvm                rnaseq_path_vals new.grps 0.7222222
## 392  ksvm          rnaseq_signaling_genes new.grps 0.8888889
## 
## $acc.unknown$linear
##     model                        x_prefix y_prefix    values
## 23   ksvm                RNAseq_all_genes new.grps 0.6666667
## 86   ksvm       rnaseq_effector_path_vals new.grps 0.5555556
## 149  ksvm rnaseq_metabolic_mod_activities new.grps 0.5555556
## 212  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.4444444
## 275  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 338  ksvm                rnaseq_path_vals new.grps 0.3888889
## 401  ksvm          rnaseq_signaling_genes new.grps 0.8333333
## 
## $acc.unknown$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 32   ksvm                RNAseq_all_genes new.grps 0.8333333
## 95   ksvm       rnaseq_effector_path_vals new.grps 0.7222222
## 158  ksvm rnaseq_metabolic_mod_activities new.grps 0.4444444
## 221  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8888889
## 284  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 347  ksvm                rnaseq_path_vals new.grps 0.7222222
## 410  ksvm          rnaseq_signaling_genes new.grps 0.8888889
## 
## $acc.unknown$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 41   ksvm                RNAseq_all_genes new.grps 0.8333333
## 104  ksvm       rnaseq_effector_path_vals new.grps 0.5555556
## 167  ksvm rnaseq_metabolic_mod_activities new.grps 0.3333333
## 230  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8333333
## 293  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.7777778
## 356  ksvm                rnaseq_path_vals new.grps 0.6111111
## 419  ksvm          rnaseq_signaling_genes new.grps 0.7777778
## 
## $acc.unknown$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 50   ksvm                RNAseq_all_genes new.grps 0.8333333
## 113  ksvm       rnaseq_effector_path_vals new.grps 0.6111111
## 176  ksvm rnaseq_metabolic_mod_activities new.grps 0.3888889
## 239  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8333333
## 302  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.6111111
## 365  ksvm                rnaseq_path_vals new.grps 0.4444444
## 428  ksvm          rnaseq_signaling_genes new.grps 0.7222222
## 
## $acc.unknown$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 59   ksvm                RNAseq_all_genes new.grps 0.8333333
## 122  ksvm       rnaseq_effector_path_vals new.grps 0.7777778
## 185  ksvm rnaseq_metabolic_mod_activities new.grps 0.6111111
## 248  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.7777778
## 311  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5000000
## 374  ksvm                rnaseq_path_vals new.grps 0.5000000
## 437  ksvm          rnaseq_signaling_genes new.grps 0.8333333
## 
## 
## $precision.Control
## $precision.Control$rbf
##     model                        x_prefix y_prefix    values
## 6    ksvm                RNAseq_all_genes new.grps 1.0000000
## 69   ksvm       rnaseq_effector_path_vals new.grps 0.5000000
## 132  ksvm rnaseq_metabolic_mod_activities new.grps 0.8333333
## 195  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.3000000
## 258  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.6250000
## 321  ksvm                rnaseq_path_vals new.grps 0.5000000
## 384  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $precision.Control$kendall
##     model                        x_prefix y_prefix    values
## 15   ksvm                RNAseq_all_genes new.grps 1.0000000
## 78   ksvm       rnaseq_effector_path_vals new.grps 0.7500000
## 141  ksvm rnaseq_metabolic_mod_activities new.grps 0.8571429
## 204  ksvm rnaseq_metabolic_mod_genevalues new.grps 1.0000000
## 267  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 330  ksvm                rnaseq_path_vals new.grps 0.6000000
## 393  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $precision.Control$linear
##     model                        x_prefix y_prefix    values
## 24   ksvm                RNAseq_all_genes new.grps 1.0000000
## 87   ksvm       rnaseq_effector_path_vals new.grps 0.8333333
## 150  ksvm rnaseq_metabolic_mod_activities new.grps 0.4000000
## 213  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.3571429
## 276  ksvm rnaseq_metabolic_mod_nodevalues new.grps 1.0000000
## 339  ksvm                rnaseq_path_vals new.grps 0.7500000
## 402  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $precision.Control$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 33   ksvm                RNAseq_all_genes new.grps 1.0000000
## 96   ksvm       rnaseq_effector_path_vals new.grps 0.6666667
## 159  ksvm rnaseq_metabolic_mod_activities new.grps 0.8333333
## 222  ksvm rnaseq_metabolic_mod_genevalues new.grps 1.0000000
## 285  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.6666667
## 348  ksvm                rnaseq_path_vals new.grps 0.5000000
## 411  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $precision.Control$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 42   ksvm                RNAseq_all_genes new.grps 0.8571429
## 105  ksvm       rnaseq_effector_path_vals new.grps 0.8571429
## 168  ksvm rnaseq_metabolic_mod_activities new.grps 0.6250000
## 231  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.8333333
## 294  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 357  ksvm                rnaseq_path_vals new.grps 0.8571429
## 420  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $precision.Control$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 51   ksvm                RNAseq_all_genes new.grps 1.0000000
## 114  ksvm       rnaseq_effector_path_vals new.grps 0.6666667
## 177  ksvm rnaseq_metabolic_mod_activities new.grps 1.0000000
## 240  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5000000
## 303  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.7500000
## 366  ksvm                rnaseq_path_vals new.grps 0.4166667
## 429  ksvm          rnaseq_signaling_genes new.grps 0.8571429
## 
## $precision.Control$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 60   ksvm                RNAseq_all_genes new.grps 0.8571429
## 123  ksvm       rnaseq_effector_path_vals new.grps 0.7500000
## 186  ksvm rnaseq_metabolic_mod_activities new.grps 0.8333333
## 249  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.7142857
## 312  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.8333333
## 375  ksvm                rnaseq_path_vals new.grps 0.6000000
## 438  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## 
## $precision.Effect_1
## $precision.Effect_1$rbf
##     model                        x_prefix y_prefix    values
## 7    ksvm                RNAseq_all_genes new.grps 0.7777778
## 70   ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 133  ksvm rnaseq_metabolic_mod_activities new.grps 1.0000000
## 196  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.6666667
## 259  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.7142857
## 322  ksvm                rnaseq_path_vals new.grps 0.8333333
## 385  ksvm          rnaseq_signaling_genes new.grps 0.7272727
## 
## $precision.Effect_1$kendall
##     model                        x_prefix y_prefix    values
## 16   ksvm                RNAseq_all_genes new.grps 0.6666667
## 79   ksvm       rnaseq_effector_path_vals new.grps 0.7777778
## 142  ksvm rnaseq_metabolic_mod_activities new.grps 0.6000000
## 205  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5000000
## 268  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.0000000
## 331  ksvm                rnaseq_path_vals new.grps 0.6666667
## 394  ksvm          rnaseq_signaling_genes new.grps 0.6666667
## 
## $precision.Effect_1$linear
##     model                        x_prefix y_prefix    values
## 25   ksvm                RNAseq_all_genes new.grps 0.6153846
## 88   ksvm       rnaseq_effector_path_vals new.grps 0.7777778
## 151  ksvm rnaseq_metabolic_mod_activities new.grps 1.0000000
## 214  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5000000
## 277  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.3333333
## 340  ksvm                rnaseq_path_vals new.grps 1.0000000
## 403  ksvm          rnaseq_signaling_genes new.grps 0.8000000
## 
## $precision.Effect_1$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 34   ksvm                RNAseq_all_genes new.grps 0.7777778
## 97   ksvm       rnaseq_effector_path_vals new.grps 0.7500000
## 160  ksvm rnaseq_metabolic_mod_activities new.grps 0.5000000
## 223  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5000000
## 286  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.0000000
## 349  ksvm                rnaseq_path_vals new.grps 0.6666667
## 412  ksvm          rnaseq_signaling_genes new.grps 0.7142857
## 
## $precision.Effect_1$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 43   ksvm                RNAseq_all_genes new.grps 0.7142857
## 106  ksvm       rnaseq_effector_path_vals new.grps 0.6923077
## 169  ksvm rnaseq_metabolic_mod_activities new.grps 0.7272727
## 232  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.7142857
## 295  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.7142857
## 358  ksvm                rnaseq_path_vals new.grps 0.8000000
## 421  ksvm          rnaseq_signaling_genes new.grps 0.6923077
## 
## $precision.Effect_1$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 52   ksvm                RNAseq_all_genes new.grps 0.8000000
## 115  ksvm       rnaseq_effector_path_vals new.grps 0.7142857
## 178  ksvm rnaseq_metabolic_mod_activities new.grps 1.0000000
## 241  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 304  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.7777778
## 367  ksvm                rnaseq_path_vals new.grps 1.0000000
## 430  ksvm          rnaseq_signaling_genes new.grps 0.6923077
## 
## $precision.Effect_1$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 61   ksvm                RNAseq_all_genes new.grps 0.6666667
## 124  ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 187  ksvm rnaseq_metabolic_mod_activities new.grps 0.6666667
## 250  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.7142857
## 313  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5555556
## 376  ksvm                rnaseq_path_vals new.grps 0.8333333
## 439  ksvm          rnaseq_signaling_genes new.grps 0.7272727
## 
## 
## $precision.Effect_2
## $precision.Effect_2$rbf
##     model                        x_prefix y_prefix values
## 8    ksvm                RNAseq_all_genes new.grps    1.0
## 71   ksvm       rnaseq_effector_path_vals new.grps    1.0
## 134  ksvm rnaseq_metabolic_mod_activities new.grps    0.5
## 197  ksvm rnaseq_metabolic_mod_genevalues new.grps    0.0
## 260  ksvm rnaseq_metabolic_mod_nodevalues new.grps    1.0
## 323  ksvm                rnaseq_path_vals new.grps    0.0
## 386  ksvm          rnaseq_signaling_genes new.grps    1.0
## 
## $precision.Effect_2$kendall
##     model                        x_prefix y_prefix values
## 17   ksvm                RNAseq_all_genes new.grps      0
## 80   ksvm       rnaseq_effector_path_vals new.grps      1
## 143  ksvm rnaseq_metabolic_mod_activities new.grps      1
## 206  ksvm rnaseq_metabolic_mod_genevalues new.grps      0
## 269  ksvm rnaseq_metabolic_mod_nodevalues new.grps      1
## 332  ksvm                rnaseq_path_vals new.grps      1
## 395  ksvm          rnaseq_signaling_genes new.grps      0
## 
## $precision.Effect_2$linear
##     model                        x_prefix y_prefix values
## 26   ksvm                RNAseq_all_genes new.grps    1.0
## 89   ksvm       rnaseq_effector_path_vals new.grps    0.5
## 152  ksvm rnaseq_metabolic_mod_activities new.grps    1.0
## 215  ksvm rnaseq_metabolic_mod_genevalues new.grps    0.0
## 278  ksvm rnaseq_metabolic_mod_nodevalues new.grps    1.0
## 341  ksvm                rnaseq_path_vals new.grps    0.2
## 404  ksvm          rnaseq_signaling_genes new.grps    1.0
## 
## $precision.Effect_2$`rbf+kendall`
##     model                        x_prefix y_prefix values
## 35   ksvm                RNAseq_all_genes new.grps      0
## 98   ksvm       rnaseq_effector_path_vals new.grps      1
## 161  ksvm rnaseq_metabolic_mod_activities new.grps      1
## 224  ksvm rnaseq_metabolic_mod_genevalues new.grps      0
## 287  ksvm rnaseq_metabolic_mod_nodevalues new.grps      1
## 350  ksvm                rnaseq_path_vals new.grps      1
## 413  ksvm          rnaseq_signaling_genes new.grps      0
## 
## $precision.Effect_2$`linear+kendall`
##     model                        x_prefix y_prefix values
## 44   ksvm                RNAseq_all_genes new.grps      0
## 107  ksvm       rnaseq_effector_path_vals new.grps      1
## 170  ksvm rnaseq_metabolic_mod_activities new.grps      1
## 233  ksvm rnaseq_metabolic_mod_genevalues new.grps      0
## 296  ksvm rnaseq_metabolic_mod_nodevalues new.grps      1
## 359  ksvm                rnaseq_path_vals new.grps      0
## 422  ksvm          rnaseq_signaling_genes new.grps      0
## 
## $precision.Effect_2$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 53   ksvm                RNAseq_all_genes new.grps 1.0000000
## 116  ksvm       rnaseq_effector_path_vals new.grps 0.6666667
## 179  ksvm rnaseq_metabolic_mod_activities new.grps 0.5000000
## 242  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 305  ksvm rnaseq_metabolic_mod_nodevalues new.grps 1.0000000
## 368  ksvm                rnaseq_path_vals new.grps 1.0000000
## 431  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## $precision.Effect_2$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 62   ksvm                RNAseq_all_genes new.grps 0.0000000
## 125  ksvm       rnaseq_effector_path_vals new.grps 1.0000000
## 188  ksvm rnaseq_metabolic_mod_activities new.grps 1.0000000
## 251  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.0000000
## 314  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.0000000
## 377  ksvm                rnaseq_path_vals new.grps 0.3333333
## 440  ksvm          rnaseq_signaling_genes new.grps 1.0000000
## 
## 
## $precision.unknown
## $precision.unknown$rbf
##     model                        x_prefix y_prefix    values
## 9    ksvm                RNAseq_all_genes new.grps 0.6153846
## 72   ksvm       rnaseq_effector_path_vals new.grps 0.5833333
## 135  ksvm rnaseq_metabolic_mod_activities new.grps 0.4705882
## 198  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.6428571
## 261  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.7142857
## 324  ksvm                rnaseq_path_vals new.grps 0.4782609
## 387  ksvm          rnaseq_signaling_genes new.grps 0.6521739
## 
## $precision.unknown$kendall
##     model                        x_prefix y_prefix    values
## 18   ksvm                RNAseq_all_genes new.grps 0.5333333
## 81   ksvm       rnaseq_effector_path_vals new.grps 0.5833333
## 144  ksvm rnaseq_metabolic_mod_activities new.grps 0.5909091
## 207  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.4848485
## 270  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5625000
## 333  ksvm                rnaseq_path_vals new.grps 0.5652174
## 396  ksvm          rnaseq_signaling_genes new.grps 0.5333333
## 
## $precision.unknown$linear
##     model                        x_prefix y_prefix    values
## 27   ksvm                RNAseq_all_genes new.grps 0.6666667
## 90   ksvm       rnaseq_effector_path_vals new.grps 0.5000000
## 153  ksvm rnaseq_metabolic_mod_activities new.grps 0.5263158
## 216  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.3809524
## 279  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.4838710
## 342  ksvm                rnaseq_path_vals new.grps 0.5000000
## 405  ksvm          rnaseq_signaling_genes new.grps 0.6250000
## 
## $precision.unknown$`rbf+kendall`
##     model                        x_prefix y_prefix    values
## 36   ksvm                RNAseq_all_genes new.grps 0.5769231
## 99   ksvm       rnaseq_effector_path_vals new.grps 0.5416667
## 162  ksvm rnaseq_metabolic_mod_activities new.grps 0.5000000
## 225  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.4705882
## 288  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.6521739
## 351  ksvm                rnaseq_path_vals new.grps 0.5652174
## 414  ksvm          rnaseq_signaling_genes new.grps 0.5517241
## 
## $precision.unknown$`linear+kendall`
##     model                        x_prefix y_prefix    values
## 45   ksvm                RNAseq_all_genes new.grps 0.6521739
## 108  ksvm       rnaseq_effector_path_vals new.grps 0.7142857
## 171  ksvm rnaseq_metabolic_mod_activities new.grps 0.6666667
## 234  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.6250000
## 297  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.6666667
## 360  ksvm                rnaseq_path_vals new.grps 0.5238095
## 423  ksvm          rnaseq_signaling_genes new.grps 0.6086957
## 
## $precision.unknown$`linear+rbf`
##     model                        x_prefix y_prefix    values
## 54   ksvm                RNAseq_all_genes new.grps 0.6250000
## 117  ksvm       rnaseq_effector_path_vals new.grps 0.5238095
## 180  ksvm rnaseq_metabolic_mod_activities new.grps 0.4375000
## 243  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5172414
## 306  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.6111111
## 369  ksvm                rnaseq_path_vals new.grps 0.5714286
## 432  ksvm          rnaseq_signaling_genes new.grps 0.6500000
## 
## $precision.unknown$`linear+rbf+kendall`
##     model                        x_prefix y_prefix    values
## 63   ksvm                RNAseq_all_genes new.grps 0.5769231
## 126  ksvm       rnaseq_effector_path_vals new.grps 0.5384615
## 189  ksvm rnaseq_metabolic_mod_activities new.grps 0.4782609
## 252  ksvm rnaseq_metabolic_mod_genevalues new.grps 0.5185185
## 315  ksvm rnaseq_metabolic_mod_nodevalues new.grps 0.5625000
## 378  ksvm                rnaseq_path_vals new.grps 0.4736842
## 441  ksvm          rnaseq_signaling_genes new.grps 0.6521739
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

![plot of chunk svm_plot](rna_figure/svm_plot-1.pdf)

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

