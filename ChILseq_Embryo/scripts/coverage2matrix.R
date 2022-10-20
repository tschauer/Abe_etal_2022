




##################################################################################################################################
##################################################################################################################################



library(ShortRead)
library(rtracklayer)
library(tsTools)


args = commandArgs(trailingOnly=TRUE)



##################################################################################################################################
##################################################################################################################################



my_gtf <- import(args[1])


my_genes <- my_gtf[my_gtf$type == "gene" & my_gtf$gene_biotype == "protein_coding"]


my_chromosomes <- seqlevels(my_genes)[1:21]
my_genes <- my_genes[seqnames(my_genes) %in% my_chromosomes]


my_TSS <- data.frame(chr = seqnames(my_genes),
                     center = ifelse(strand(my_genes) == "+", start(my_genes), end(my_genes)),
                     strand = strand(my_genes), 
                     row.names = my_genes$gene_id)


my_TTS <- data.frame(chr = seqnames(my_genes),
                     center = ifelse(strand(my_genes) == "-", start(my_genes), end(my_genes)),
                     strand = strand(my_genes), 
                     row.names = my_genes$gene_id)


##################################################################################################################################
##################################################################################################################################



norm_cov <- readRDS(args[2])


mat.TSS <- coverageWindowsCenteredStranded(centers = my_TSS,
                                           coverage = norm_cov,
                                           window.size = 10000)

saveRDS(mat.TSS, file=args[3])



mat.TTS <- coverageWindowsCenteredStranded(centers = my_TTS,
                                           coverage = norm_cov,
                                           window.size = 10000)

saveRDS(mat.TTS, file=args[4])



##################################################################################################################################
##################################################################################################################################


