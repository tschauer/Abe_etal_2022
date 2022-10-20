


##################################################################################################################################
##################################################################################################################################


library(ShortRead)
library(rtracklayer)


args = commandArgs(trailingOnly=TRUE)


##################################################################################################################################
##################################################################################################################################



gap <- readGAlignmentPairs(grep(".bam$", args, value = TRUE))
grs <- granges(gap, on.discordant.seqnames="drop")

my_chromosomes <- seqlevels(grs)[!(grepl("M|J|G",seqlevels(grs)))]
#grs <- keepSeqlevels(grs, my_chromosomes, pruning.mode = "coarse")

saveRDS(grs, file=grep("ranges.rds$", args, value = TRUE))


if(as.logical(args[5])){
        
        grs$chr_start <- paste(seqnames(grs), start(grs), sep="_")
        grs$chr_end <-   paste(seqnames(grs), end(grs),   sep="_")
        
        grs <- grs[!(duplicated(grs$chr_start)) & !(duplicated(grs$chr_end))]
        
}


##################################################################################################################################
##################################################################################################################################


covs <- coverage(grs)
normcovs <- covs/length(grs)*1e6

export.bw(normcovs, con=grep("coverage.bw$", args, value = TRUE))
saveRDS(normcovs,  file=grep("coverage.rds$", args, value = TRUE))


##################################################################################################################################
##################################################################################################################################


