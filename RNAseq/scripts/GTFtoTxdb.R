




library(GenomicFeatures)

args = commandArgs(trailingOnly=TRUE)


txdb <- makeTxDbFromGFF(args[1])
saveDb(txdb, file = args[2])