




##################################################################################################################################
##################################################################################################################################



library(ShortRead)
library(rtracklayer)



args = commandArgs(trailingOnly=TRUE)


my_gtf <-  import(grep("gtf", args, value = TRUE))

my_fai <- read.delim(grep("fai", args, value = TRUE), header = FALSE)



##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ########################################################## 




my_gene_ranges <- my_gtf[my_gtf$type == "gene" & my_gtf$gene_biotype == "protein_coding"]

my_chromosomes <- seqlevels(my_gene_ranges)[1:21]
my_gene_ranges <- keepSeqlevels(my_gene_ranges, my_chromosomes, pruning.mode = "coarse")

my_gene_ranges <- my_gene_ranges[width(my_gene_ranges) > 1000]

mcols(my_gene_ranges) <- mcols(my_gene_ranges)[,c("gene_id","gene_name")]
names(my_gene_ranges) <- my_gene_ranges$gene_id




##################################################################################################################################
##################################################################################################################################




my_genebody_ranges <- my_gene_ranges

start(my_genebody_ranges) <- start(my_genebody_ranges)+500
end(my_genebody_ranges)   <- end(my_genebody_ranges)-500


export.bed(object = my_genebody_ranges, gsub("rda","bed", grep("genebody", args, value = TRUE)))
save(my_genebody_ranges, file = grep("genebody", args, value = TRUE))




##################################################################################################################################
##################################################################################################################################




my_TSS_ranges <- promoters(my_gene_ranges, upstream = 2500, downstream = 2500)


stopifnot(identical(names(my_genebody_ranges), names(my_TSS_ranges)))
stopifnot(identical(my_TSS_ranges$gene_id, names(my_TSS_ranges)))


export.bed(object = my_TSS_ranges, gsub("rda","bed",grep("TSS", args, value = TRUE)))
save(my_TSS_ranges, file = grep("TSS", args, value = TRUE))





##################################################################################################################################
##################################################################################################################################





my_anti_ranges <- my_gene_ranges

strand(my_anti_ranges)[strand(my_gene_ranges) == "+"] <- "-"
strand(my_anti_ranges)[strand(my_gene_ranges) == "-"] <- "+"


my_TTS_ranges <- promoters(my_anti_ranges, upstream = 2500, downstream = 2500)

strand(my_TTS_ranges)[strand(my_anti_ranges) == "+"] <- "-"
strand(my_TTS_ranges)[strand(my_anti_ranges) == "-"] <- "+"



stopifnot(identical(my_TSS_ranges$gene_id, my_TTS_ranges$gene_id))
stopifnot(identical(names(my_TTS_ranges), my_TTS_ranges$gene_id))


export.bed(object = my_TTS_ranges, gsub("rda","bed", grep("TTS", args, value = TRUE)))
save(my_TTS_ranges, file = grep("TTS", args, value = TRUE))





##################################################################################################################################
##################################################################################################################################





my_GB5_ranges <- my_TSS_ranges

start(my_GB5_ranges[strand(my_TSS_ranges) == "+"]) <- start(my_TSS_ranges[strand(my_TSS_ranges) == "+"])+5000
end(my_GB5_ranges[  strand(my_TSS_ranges) == "+"]) <-   end(my_TSS_ranges[strand(my_TSS_ranges) == "+"])+5000

start(my_GB5_ranges[strand(my_TSS_ranges) == "-"]) <- start(my_TSS_ranges[strand(my_TSS_ranges) == "-"])-5000
end(my_GB5_ranges[  strand(my_TSS_ranges) == "-"]) <-   end(my_TSS_ranges[strand(my_TSS_ranges) == "-"])-5000



stopifnot(identical(my_genebody_ranges$gene_id, my_GB5_ranges$gene_id))
stopifnot(identical(names(my_GB5_ranges), my_GB5_ranges$gene_id))



export.bed(object = my_GB5_ranges, gsub("rda","bed", grep("GB5", args, value = TRUE)))
save(my_GB5_ranges, file = grep("GB5", args, value = TRUE))




##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################

#################################################       Genome SeqInfo       ##################################################### 



my_seqinfo <- Seqinfo(seqnames = as.character(my_fai$V1), 
                      seqlengths = my_fai$V2, 
                      isCircular = rep(FALSE, nrow(my_fai)), 
                      genome = "GRCm38")


my_bins <- tileGenome(seqlengths = my_seqinfo, 
                      tilewidth = 1e5,
                      cut.last.tile.in.chrom = TRUE)


my_bins <- keepSeqlevels(my_bins, my_chromosomes, pruning.mode = "coarse")
my_bins$bin_id <- paste(seqnames(my_bins), start(my_bins), end(my_bins), sep = "_" )
names(my_bins) <- my_bins$bin_id


export.bed(object = my_bins, gsub("rda","bed", grep("bins", args, value = TRUE)))
save(my_bins, file = grep("bins", args, value = TRUE))



##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################

#################################################        10 kb bins          ##################################################### 




my_b10_ranges <- tileGenome(seqlengths = my_seqinfo, 
                      tilewidth = 1e4,
                      cut.last.tile.in.chrom = TRUE)


my_b10_ranges <- keepSeqlevels(my_b10_ranges, my_chromosomes, pruning.mode = "coarse")
my_b10_ranges$b10_id <- paste(seqnames(my_b10_ranges), start(my_b10_ranges), end(my_b10_ranges), sep = "_" )
names(my_b10_ranges) <- my_b10_ranges$b10_id


export.bed(object = my_b10_ranges, gsub("rda","bed", grep("b10", args, value = TRUE)))
save(my_b10_ranges, file = grep("b10", args, value = TRUE))



##################################################################################################################################
##################################################################################################################################














##################################################################################################################################
##################################################################################################################################





my_exon_ranges <- my_gtf[my_gtf$type == "exon" & my_gtf$gene_biotype == "protein_coding"]
my_exon_ranges <- my_exon_ranges[my_exon_ranges$gene_id %in% my_gene_ranges$gene_id]
my_exon_ranges <- reduce(my_exon_ranges, ignore.strand=TRUE)

names(my_exon_ranges) <- paste("exon", seq_along(my_exon_ranges), sep="_")

export.bed(object = my_exon_ranges, gsub("rda","bed", grep("exons", args, value = TRUE)))
save(my_exon_ranges, file = grep("exons", args, value = TRUE))




##################################################################################################################################
##################################################################################################################################




my_intron_ranges <- setdiff(my_gene_ranges, my_exon_ranges, ignore.strand=TRUE)
my_intron_ranges <- reduce(my_intron_ranges, ignore.strand=TRUE)

names(my_intron_ranges) <- paste("intron", seq_along(my_intron_ranges), sep="_")

export.bed(object = my_intron_ranges, gsub("rda","bed", grep("introns", args, value = TRUE)))
save(my_intron_ranges, file = grep("introns", args, value = TRUE))





##################################################################################################################################
##################################################################################################################################





my_intergenic_ranges <- gaps(reduce(c(my_TSS_ranges, my_genebody_ranges, my_TTS_ranges), ignore.strand=TRUE))
my_intergenic_ranges <- reduce(my_intergenic_ranges, ignore.strand=TRUE)

names(my_intergenic_ranges) <- paste("intergenic", seq_along(my_intergenic_ranges), sep="_")

export.bed(object = my_intergenic_ranges, gsub("rda","bed", grep("intergenic", args, value = TRUE)))
save(my_intergenic_ranges, file = grep("intergenic", args, value = TRUE))





##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################


