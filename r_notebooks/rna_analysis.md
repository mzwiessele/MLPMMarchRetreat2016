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
ylist <- c('effect.grps', 'new.grps')
```

# SVM


```r
kflist <- c('linear', 'kendall', 'rbf')
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
dd <- data.frame()
for (i in seq(length(res))) {
  dd <- rbind(dd, 
              data.frame(model = res[[i]]$model, kf = res[[i]]$kf, 
                         x_prefix = res[[i]]$x_prefix, y_prefix = res[[i]]$y_prefix, 
                         acc = res[[i]]$acc))
}
dd
```

```
##    model      kf                        x_prefix    y_prefix       acc
## 1   ksvm  linear                RNAseq_all_genes effect.grps 0.4047619
## 2   ksvm kendall                RNAseq_all_genes effect.grps 0.3333333
## 3   ksvm     rbf                RNAseq_all_genes effect.grps 0.3809524
## 4   ksvm  linear                RNAseq_all_genes    new.grps 0.6666667
## 5   ksvm kendall                RNAseq_all_genes    new.grps 0.5714286
## 6   ksvm     rbf                RNAseq_all_genes    new.grps 0.7142857
## 7   ksvm  linear       rnaseq_effector_path_vals effect.grps 0.4047619
## 8   ksvm kendall       rnaseq_effector_path_vals effect.grps 0.3809524
## 9   ksvm     rbf       rnaseq_effector_path_vals effect.grps 0.4047619
## 10  ksvm  linear       rnaseq_effector_path_vals    new.grps 0.4285714
## 11  ksvm kendall       rnaseq_effector_path_vals    new.grps 0.6666667
## 12  ksvm     rbf       rnaseq_effector_path_vals    new.grps 0.6666667
## 13  ksvm  linear rnaseq_metabolic_mod_activities effect.grps 0.2142857
## 14  ksvm kendall rnaseq_metabolic_mod_activities effect.grps 0.3571429
## 15  ksvm     rbf rnaseq_metabolic_mod_activities effect.grps 0.2142857
## 16  ksvm  linear rnaseq_metabolic_mod_activities    new.grps 0.5476190
## 17  ksvm kendall rnaseq_metabolic_mod_activities    new.grps 0.6190476
## 18  ksvm     rbf rnaseq_metabolic_mod_activities    new.grps 0.5000000
## 19  ksvm  linear rnaseq_metabolic_mod_genevalues effect.grps 0.2619048
## 20  ksvm kendall rnaseq_metabolic_mod_genevalues effect.grps 0.4047619
## 21  ksvm     rbf rnaseq_metabolic_mod_genevalues effect.grps 0.3095238
## 22  ksvm  linear rnaseq_metabolic_mod_genevalues    new.grps 0.3333333
## 23  ksvm kendall rnaseq_metabolic_mod_genevalues    new.grps 0.5238095
## 24  ksvm     rbf rnaseq_metabolic_mod_genevalues    new.grps 0.5000000
## 25  ksvm  linear rnaseq_metabolic_mod_nodevalues effect.grps 0.3809524
## 26  ksvm kendall rnaseq_metabolic_mod_nodevalues effect.grps 0.4047619
## 27  ksvm     rbf rnaseq_metabolic_mod_nodevalues effect.grps 0.4047619
## 28  ksvm  linear rnaseq_metabolic_mod_nodevalues    new.grps 0.5714286
## 29  ksvm kendall rnaseq_metabolic_mod_nodevalues    new.grps 0.6190476
## 30  ksvm     rbf rnaseq_metabolic_mod_nodevalues    new.grps 0.5000000
## 31  ksvm  linear                rnaseq_path_vals effect.grps 0.3571429
## 32  ksvm kendall                rnaseq_path_vals effect.grps 0.3571429
## 33  ksvm     rbf                rnaseq_path_vals effect.grps 0.3809524
## 34  ksvm  linear                rnaseq_path_vals    new.grps 0.6190476
## 35  ksvm kendall                rnaseq_path_vals    new.grps 0.5952381
## 36  ksvm     rbf                rnaseq_path_vals    new.grps 0.5000000
## 37  ksvm  linear          rnaseq_signaling_genes effect.grps 0.4047619
## 38  ksvm kendall          rnaseq_signaling_genes effect.grps 0.3809524
## 39  ksvm     rbf          rnaseq_signaling_genes effect.grps 0.4047619
## 40  ksvm  linear          rnaseq_signaling_genes    new.grps 0.6666667
## 41  ksvm kendall          rnaseq_signaling_genes    new.grps 0.5714286
## 42  ksvm     rbf          rnaseq_signaling_genes    new.grps 0.7142857
```

```r
for (yname in ylist) {
  p1 <- ggplot(subset(dd, y_prefix == yname), aes(x = x_prefix, y = acc)) + 
    geom_bar(aes(fill = x_prefix), stat="identity", position = "dodge") + 
    facet_wrap(~ kf) + ylim(0,1) + 
    ggtitle(paste0('predict for labels = ', yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
}
```

![plot of chunk svm_plot](rna_figure/svm_plot-1.pdf)![plot of chunk svm_plot](rna_figure/svm_plot-2.pdf)

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
##  ui       RStudio (0.99.473)          
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       Europe/Paris                
##  date     2016-03-17
```

```
## Packages ------------------------------------------------------------------
```

```
##  package      * version date       source        
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
##  MASS           7.3-45  2015-11-10 CRAN (R 3.2.3)
##  Matrix         1.2-3   2015-11-28 CRAN (R 3.2.3)
##  MatrixModels   0.4-1   2015-08-22 CRAN (R 3.2.0)
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
##  Rcpp           0.12.3  2016-01-10 CRAN (R 3.2.3)
##  reshape2       1.4.1   2014-12-06 CRAN (R 3.1.2)
##  scales         0.4.0   2016-02-26 CRAN (R 3.2.3)
##  SparseM        1.7     2015-08-15 CRAN (R 3.2.0)
##  stringi        1.0-1   2015-10-22 CRAN (R 3.2.0)
##  stringr        1.0.0   2015-04-30 CRAN (R 3.2.0)
```

