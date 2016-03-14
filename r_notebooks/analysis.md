#### Notebook : MLPM retreat Valencia


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, fig.width = 8, fig.height = 8, dev = "pdf", fig.keep = "high", fig.path = "figure/", cache.path = "cache/")
set.seed(35875954)
```

# Load data


```r
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

