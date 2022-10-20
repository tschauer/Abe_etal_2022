




##################################################################################################################################
##################################################################################################################################



library(ShortRead)
library(rtracklayer)
library(tsTools)

library(org.Mm.eg.db)


args = commandArgs(trailingOnly=TRUE)


output_file <- grep("pdf",args, value=TRUE)

my_assay <- gsub(".*examples.|.pdf","",output_file)


##################################################################################################################################
##################################################################################################################################



my_gtf <- import(args[1])


my_genes <- my_gtf[my_gtf$type == "gene" & my_gtf$gene_biotype == "protein_coding"]

my_exons <- my_gtf[my_gtf$type == "exon" & my_gtf$gene_biotype == "protein_coding"]



##################################################################################################################################
##################################################################################################################################






##################################################################################################################################
##################################################################################################################################

######################################################      DBTMEE      ########################################################## 



my_DBTMEE <- read.delim("external_data//cluster_gene.tsv", stringsAsFactors = FALSE)
my_DBTMEE <- my_DBTMEE[!is.na(my_DBTMEE$Entrez),]


my_DBTMEE$Entrez <- as.character(my_DBTMEE$Entrez)

my_DBTMEE$gene_id <- as.character(mapIds(org.Mm.eg.db, keys = as.character(my_DBTMEE$Entrez), 
                                         keytype = "ENTREZID", column = "ENSEMBL", multiVals = "first"))

my_DBTMEE$gene_symbol <- as.character(mapIds(org.Mm.eg.db, keys = as.character(my_DBTMEE$Entrez), 
                                             keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"))

my_DBTMEE$Cluster <- gsub(" ", "_",my_DBTMEE$Cluster)
my_DBTMEE <- my_DBTMEE[,!(grepl("Correlation|ID", colnames(my_DBTMEE)))]
my_DBTMEE$UID <- apply(my_DBTMEE, 1, function(x){ paste(x, collapse = "_") })
my_DBTMEE <- my_DBTMEE[!(duplicated(my_DBTMEE$UID)),]


my_DBTMEE <- my_DBTMEE[!(is.na(my_DBTMEE$gene_symbol)),]
my_DBTMEE <- my_DBTMEE[!(is.na(my_DBTMEE$gene_id)),]




##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################



my_res <- read.table("Output/deseq2_genebody/results/res.S2P.2C-Zy.txt", stringsAsFactors = FALSE)


my_res_merged <- merge(my_res, my_DBTMEE, by.x = "row.names", by.y = "gene_id")


my_fav_genes <- my_res_merged$gene_name[my_res_merged$Cluster %in% c("4-Cell_Transient","2-Cell_Transient") & 
                                                my_res_merged$log2FoldChange > 1.5 & 
                                                my_res_merged$padj < 0.005 & 
                                                my_res_merged$baseMean > 10]

my_fav_genes <- my_fav_genes[-grep("Sfi1", my_fav_genes)]

my_fav_ranges <- my_genes[my_genes$gene_name %in% my_fav_genes]



start(my_fav_ranges) <- floor(start(my_fav_ranges)/1e4)*1e4
end(my_fav_ranges) <- ceiling(end(my_fav_ranges)/1e4)*1e4


my_fav_ranges <- resize(my_fav_ranges, width = 250000, fix = "center")


##################################################################################################################################
##################################################################################################################################





# my_coordinates_list <- list(c("5",  103550000, 103950000),
#                             c("12",   4750000,   4950000),
#                             c("4",   74150000,  74600000),
#                             c("2",  160500000, 160850000),
#                             c("18", 84750000,   85150000),
#                             c("8" ,111650000,  112050000))

my_coordinates_list <- lapply(1:length(my_fav_ranges), function(x){as.data.frame(my_fav_ranges)[x,1:3]})


my_samples <- unique(gsub(".*\\/|_norm.*", "", grep("coverage.*rds", args, value = TRUE)))



##################################################################################################################################
##################################################################################################################################




pdf(file = output_file, height = 6, width = 10)

j=1
i=1


for(j in seq_along(my_coordinates_list)){
        
        
        my_coordinates <- makeGRangesFromDataFrame(data.frame(chr = my_coordinates_list[[j]][1], 
                                                              start= my_coordinates_list[[j]][2], 
                                                              end =  my_coordinates_list[[j]][3]))
        
        #for(my_assay in my_assays){
        
        par(mfrow=c(6,1), oma=c(2,2,2,2), mar = c(1,4,1,2))
        
        
        
        
        my_samples_sub <- grep(my_assay, my_samples, value = TRUE)
        my_samples_sub <- my_samples_sub[order(my_samples_sub)]
        my_samples_sub <- my_samples_sub[c(grep("Zy", my_samples_sub),grep("2C|8C|ES",my_samples_sub))]
        my_samples_sub <- my_samples_sub[!grepl("DRB",my_samples_sub)]
        
        
        if(length(my_samples_sub) == 15){
                
                my_color_palette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")
                
        } else {
                my_color_palette <- c("#999999", "#D55E00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")
                
        }
        
        
        #####################################################################################################
        
        
        for(i in seq(1, length(my_samples_sub),3)){
                
                
                my_cov1 <- grep(my_samples_sub[i+0], args, value = TRUE)
                my_cov2 <- grep(my_samples_sub[i+1], args, value = TRUE)
                my_cov3 <- grep(my_samples_sub[i+2], args, value = TRUE)
                
                
                stopifnot(identical(gsub(".*\\/|_[1-3]_norm.*", "", my_cov1),
                                    gsub(".*\\/|_[1-3]_norm.*", "", my_cov2)))
                
                stopifnot(identical(gsub(".*\\/|_[1-3]_norm.*", "", my_cov1),
                                    gsub(".*\\/|_[1-3]_norm.*", "", my_cov3)))
                
                print(c(my_cov1,my_cov2,my_cov3))
                
                norm_cov1 <- readRDS(my_cov1)
                norm_cov2 <- readRDS(my_cov2)
                norm_cov3 <- readRDS(my_cov3)
                
                
                norm_cov <- (norm_cov1 + norm_cov2 + norm_cov3) / 3
                
                norm_cov_coord <- norm_cov[my_coordinates][[1]]
                
                norm_cov_coord_binned <- binMeans(y = as.numeric(norm_cov_coord),
                                                  x = seq_along(norm_cov_coord),
                                                  bx = seq(1, length(norm_cov_coord), 1000))
                
                plot(norm_cov_coord_binned, type="l", 
                     main = "",
                     xlab = "", ylab="", xaxt="n", bty = "n",
                     col = rep(my_color_palette, each = 3)[i],
                     ylim = c(0, ifelse(j %in% c(1,2,3), 1.6, 2.0)), 
                     lwd = 1.5)
                
                
                xx <- c((1:length(norm_cov_coord_binned)), (length(norm_cov_coord_binned):1))
                yy <- c(norm_cov_coord_binned, rep(-0.5, length(norm_cov_coord_binned)))
                
                polygon(xx,yy, border = NA, 
                        col = rep(my_color_palette, each = 3)[i])
                
                title(main = unique(gsub(".*\\/|_[1-3]_norm.*","",c(my_cov1,my_cov2,my_cov3))), 
                      adj = 0.05)
                
                axis(side = 1, labels = FALSE, lwd.ticks = 0, lwd = 1,)
                
                
        }
        
        #####################################################################################################
        
        
        mtext(text = "Normalized Coverage", side = 2, line = 0, outer = TRUE)
        
        mtext(text = paste0("chr",my_coordinates[1]), side = 1, line = -1,outer = TRUE)
        
        
        #####################################################################################################
        
        
        par(mar = c(4,3,0,2))
        
        plot(x = c(start(my_coordinates), end(my_coordinates)) - start(my_coordinates),
             y = c(0,1), 
             type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
        
        
        try(
                rect(xleft =   start(subsetByOverlaps(my_genes[strand(my_genes) == "+"], my_coordinates)) - start(my_coordinates),
                     xright =  end(  subsetByOverlaps(my_genes[strand(my_genes) == "+"], my_coordinates)) - start(my_coordinates),
                     ybottom = 0.65, ytop = 0.65, lwd=2, border = "#666666"), silent = TRUE)
        
        try(
                rect(xleft =   start(subsetByOverlaps(my_genes[strand(my_genes) == "-"], my_coordinates)) - start(my_coordinates),
                     xright =  end(  subsetByOverlaps(my_genes[strand(my_genes) == "-"], my_coordinates)) - start(my_coordinates),
                     ybottom = 0.35, ytop = 0.35, lwd=2, border = "#666666"), silent = TRUE)
        
        my_genes_sub <- subsetByOverlaps(my_genes[my_genes$gene_name %in% my_fav_genes], my_coordinates)
        
        try(
                text(x =  mean(c(start(my_genes_sub) - start(my_coordinates),
                                 end(my_genes_sub) - start(my_coordinates))),
                     y =  ifelse(strand(my_genes_sub) == "+",0.1,0.9), 
                     col = "#666666", font = 2,
                     labels = paste0(my_genes_sub$gene_name,
                                     " *[",
                                     my_DBTMEE$Cluster[my_DBTMEE$gene_id == my_genes_sub$gene_id],
                                     "]"))
        )
        
        
        try(
                rect(xleft =   start(subsetByOverlaps(my_exons[strand(my_exons) == "+"], my_coordinates)) - start(my_coordinates),
                     xright =  end(  subsetByOverlaps(my_exons[strand(my_exons) == "+"], my_coordinates)) - start(my_coordinates),
                     ybottom = 0.50, ytop = 0.80, lwd=2, border = "#666666", col = "#666666"), silent = TRUE)
        
        try(
                rect(xleft =   start(subsetByOverlaps(my_exons[strand(my_exons) == "-"], my_coordinates)) - start(my_coordinates),
                     xright =  end(  subsetByOverlaps(my_exons[strand(my_exons) == "-"], my_coordinates)) - start(my_coordinates),
                     ybottom = 0.20, ytop = 0.50, lwd=2, border = "#666666", col = "#666666"), silent = TRUE)
        
        
        
        axis(side = 1,
             at =    seq(0, (end(my_coordinates) - start(my_coordinates)), 1e5),
             labels = as.character(seq(start(my_coordinates), end(my_coordinates), 1e5)), 
             line = 0, lwd = 0, lwd.ticks = 1)
        
        
        #####################################################################################################
        
        
}

#}


dev.off()

##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

