# Notebook : MLPM retreat Valencia


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, fig.width = 12, fig.height = 8, dev = "pdf", fig.keep = "high", fig.path = "figure/", cache.path = "cache/")
set.seed(35875954)

library(caret)
library(pcaPP)
library(kernlab)
source('predictors.R')
```

# Load data


```r
# read in datasets
datapath <- '../Data/RNASeq/'
fnames <- list.files(path = datapath, pattern = '[.]csv$', full.names = FALSE)
objnames <- gsub('[.]csv$', '', fnames)
for (obj in objnames) {
  x <- read.csv(paste0(datapath, obj, '.csv'), sep = ',', header = TRUE, row.names = 1)
  rownames(x) <- paste0('X_', gsub('[^[:alnum:]_]', '_', rownames(x)))
  if (!grepl('design', obj)) x <- t(x) # observations in rows and features in col
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
for (xname in xlist) {
  stopifnot(isTRUE(all.equal(rownames(get(xname)), as.character(as.character(rnaseq_design$SampleID_RNASeq)))))
}
# generate train2test
trainidx <- grep('train', rnaseq_design$NA., ignore.case = TRUE)
table(rnaseq_design$NA.)
```

```
## 
##     TEST SET TRAINING SET 
##           42           62
```

```r
# generate labels
effect.grps <- as.factor(rnaseq_design$Effect)
```

# SVM


```r
xname <- 'rnaseq_effector_path_vals'
yname <- 'effect.grps'
kf <- 'linear'
res <- perfSVM(model = paste('ksvm', kf, sep = '_AAA_'), prefix = paste(xname, yname, sep = '_AAA_'), 
               xdata = get(xname), grp = get(yname), kf = kf, tr2tstFolds = trainidx)
res$acc
```

```
## [1] 0.4047619
```

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
##  date     2016-03-14
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
##  rstudioapi     0.5     2016-01-24 CRAN (R 3.2.3)
##  scales         0.4.0   2016-02-26 CRAN (R 3.2.3)
##  SparseM        1.7     2015-08-15 CRAN (R 3.2.0)
##  stringi        1.0-1   2015-10-22 CRAN (R 3.2.0)
##  stringr        1.0.0   2015-04-30 CRAN (R 3.2.0)
```

