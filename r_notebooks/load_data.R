
# read data ---------------------------------------------------------------

# RNAseq

datapath <- '../Data/RNASeq/'

fnames <- list.files(path = datapath, pattern = '[.]csv$', full.names = FALSE)
objnames <- gsub('[.]csv$', '', fnames)
for (obj in objnames) {
  x <- read.csv(paste0(datapath, obj, '.csv'), sep = ',', header = TRUE, row.names = 1)
  rownames(x) <- gsub('[^[:alnum:]_]', '_', rownames(x))
  if (!grepl('design', obj)) x <- t(x) # observations in rows and features in col
  assign(obj, x)
}

