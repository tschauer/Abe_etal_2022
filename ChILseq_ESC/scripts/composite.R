






##################################################################################################################################
##################################################################################################################################




library(DESeq2)
library(ShortRead)
library(rtracklayer)
library(HelpersforDESeq2)
library(HelpersforChIPSeq)
library(topGO)
library(pheatmap)
library(sva)

library(multcomp)
library(gplots)

library(reshape2)
library(ggplot2)
library(ggpubr)


library(org.Mm.eg.db)


##################################################################################################################################
##################################################################################################################################




args = commandArgs(trailingOnly=TRUE)



input_table <- grep("SampleTable", args, value = TRUE)

output_file <- grep("composite.*.pdf", args, value = TRUE)


my_group <- gsub(".*\\.", "", gsub(".pdf", "",output_file))
my_assay <- gsub(".*_","",gsub("_DRB","", my_group))
my_stage <- gsub(paste0("_",my_assay),"", my_group)


matrix_files <- grep(paste0(my_group, "_T.*matrix.rds"), args, value = TRUE)
matrix_files <- matrix_files[order(matrix_files)]


stopifnot(length(matrix_files) == 2)
stopifnot(identical(gsub("TSS", "",matrix_files[1]), gsub("TTS", "",matrix_files[2])))



##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

#################################################        SampleTable         #####################################################




SampleTable <- read.delim(input_table, header = T, stringsAsFactors = F)
SampleTable <- SampleTable[, 2, drop = FALSE]

SampleTable <- SampleTable[!duplicated(SampleTable$ForeignID), , drop = FALSE]
SampleTable <- SampleTable[order(SampleTable$ForeignID), , drop = FALSE]

SampleTable$Assay <-  gsub(".*_","",gsub("_[1-9]$|_DRB","", SampleTable$ForeignID))
SampleTable$Stages <- gsub(paste(c(paste0("_",unique(SampleTable$Assay)),"_[1-9]$"), collapse = "|"),"", SampleTable$ForeignID)

SampleTable$Conditions <- paste(SampleTable$Stages, SampleTable$Assay, sep = "_")



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

######################################################      RNA-seq      ########################################################## 




TPM <- read.table("../../Public/GSE38495_GSE45719/Output/TPM/TPM_means.txt", check.names = FALSE)

log2_TPM <- log2(TPM+1)

log2_TPM$gene_id <- rownames(log2_TPM)


##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################

######################################################      Setup    ########################################################## 





SampleMatchTable <- data.frame(Stages = c("ChIL1000C", "ChIL100C",  "ChIP"),
                               RNAseq = c("ESC", "ESC", "ESC" ),
                               Colors = c("#009E73", "#009E73", "#009E73"), 
                               stringsAsFactors = FALSE)



my_color_palette1 <- colorRampPalette(c("white",
                                        SampleMatchTable$Colors[SampleMatchTable$Stages == my_stage],
                                        "black"))(7)[c(-1,-7)]

my_color_palette3 <- brewer.pal(9,"Set1")[c(2,9,1)]



##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################

######################################################     Read Mats    ########################################################## 


i=1

for(i in seq_along(matrix_files)){
        
        my_name <- gsub(".*matrix_ave/|_matrix.*","",matrix_files[i])
        my_site <- gsub(".*_","",my_name)
        my_name <- paste0(my_site,".",gsub(paste0("_",my_site),"", my_name))
        
        
        my_mat <- readRDS(matrix_files[i])
        
        assign(my_name, my_mat)
        
        rm(list = "my_mat")
}


my_mats <- ls(pattern = paste0("T[S,T]S.", my_group))


for(i in seq_along(my_mats)){
        
        stopifnot(identical(rownames(get(my_mats[1])), rownames(get(my_mats[i]))))
        stopifnot(sum(grepl("ENS", rownames(get(my_mats[i])))) == nrow(get(my_mats[i])))
}




##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

######################################################      Quantiles     ######################################################## 





log2_TPM_stage <- log2_TPM[, c("gene_id",  SampleMatchTable$RNAseq[SampleMatchTable$Stages == my_stage])]

log2_TPM_stage <- na.omit(log2_TPM_stage)


log2_TPM_stage <- log2_TPM_stage[!(log2_TPM_stage$gene_id %in% my_DBTMEE$gene_id[my_DBTMEE$Cluster == "Maternal_RNA"]), ]
log2_TPM_stage <- log2_TPM_stage[  log2_TPM_stage$gene_id %in% rownames(get(my_mats[1])), ]



log2_TPM_stage$Groups <- ifelse(round(log2_TPM_stage[,2], 4) == 0.0000, "q1", "exp")

is_nonzero <- log2_TPM_stage$Groups == "exp"

log2_TPM_stage$Groups[is_nonzero] <- c("q2","q3","q4","q5")[cut(x = log2_TPM_stage[is_nonzero, 2],
                                                                breaks = quantile(log2_TPM_stage[is_nonzero, 2]),
                                                                include.lowest = TRUE)]




##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

##############################################    RNA-seq Cdk9in Spt5 results   ###################################################


# padj_cutoff <- 0.05
# 
# 
# `res_Cdk9in-Noinject` <- read.table("../211011_RNAseq/Output/deseq2/tables/res.Cdk9in-Noinject.txt")
# `res_Cdk9in-Noinject`$Groups <- "NS"
# `res_Cdk9in-Noinject`$Groups[`res_Cdk9in-Noinject`$log2FoldChange > 0 & `res_Cdk9in-Noinject`$padj < padj_cutoff] <- "up-reg"
# `res_Cdk9in-Noinject`$Groups[`res_Cdk9in-Noinject`$log2FoldChange < 0 & `res_Cdk9in-Noinject`$padj < padj_cutoff] <- "down-reg"
# 
# #`res_Cdk9in-Noinject` <- `res_Cdk9in-Noinject`[!(`res_Cdk9in-Noinject`$gene_id %in% my_DBTMEE$gene_id[my_DBTMEE$Cluster == "Maternal_RNA"]), ]
# `res_Cdk9in-Noinject` <- `res_Cdk9in-Noinject`[  `res_Cdk9in-Noinject`$gene_id %in% rownames(get(my_mats[1])), ]
# 
# 
# 
# `res_Spt5KD-IgG` <- read.table("../211011_RNAseq/Output/deseq2/tables/res.Spt5KD-IgG.txt")
# `res_Spt5KD-IgG`$Groups <- "NS"
# `res_Spt5KD-IgG`$Groups[`res_Spt5KD-IgG`$log2FoldChange > 0 & `res_Spt5KD-IgG`$padj < padj_cutoff] <- "up-reg"
# `res_Spt5KD-IgG`$Groups[`res_Spt5KD-IgG`$log2FoldChange < 0 & `res_Spt5KD-IgG`$padj < padj_cutoff] <- "down-reg"
# 
# #`res_Spt5KD-IgG` <- `res_Spt5KD-IgG`[!(`res_Spt5KD-IgG`$gene_id %in% my_DBTMEE$gene_id[my_DBTMEE$Cluster == "Maternal_RNA"]), ]
# `res_Spt5KD-IgG` <- `res_Spt5KD-IgG`[  `res_Spt5KD-IgG`$gene_id %in% rownames(get(my_mats[1])), ]



##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

######################################################      Plotting     ######################################################### 




# my_classes <- unique(my_DBTMEE$Cluster) 
# my_class <- my_classes[1]





my_site = "TSS"

my_grouping_names <- c("log2_TPM_stage")

i=2

######################################################      


for(i in seq_along(my_grouping_names)){
        
        my_file_name <- ifelse(my_grouping_names[i] == "log2_TPM_stage", 
                               output_file,
                               paste0(gsub(".pdf", paste0(".",gsub("res_","",my_grouping_names[i]),".pdf"), output_file)))
        
        
        pdf(file = my_file_name, width = 6.5, height = 4.5)
        par(oma=c(2,2,2,2), mar=c(4,4,4,2))
        
        
        my_grouping <- get(my_grouping_names[i])
        
        
        ######################################################      
        
        stopifnot(identical(gsub("TSS","",my_mats[1]), gsub("TTS","",my_mats[2])))
        
        my_mat1 <- get(my_mats[1])
        my_mat2 <- get(my_mats[2])
        
        stopifnot(identical(rownames(my_mat1), rownames(my_mat2)))
        
        mat.combined.all <- cbind(my_mat1, matrix(NA, nrow = nrow(my_mat2), ncol = 1), my_mat2)
        
        
        for(gi in unique(my_grouping$Groups)){
                
                assign(paste0("mat.combined.",gi), 
                       mat.combined.all[rownames(mat.combined.all) %in% my_grouping$gene_id[my_grouping$Groups == gi],])
        }
        
        
        ######################################################      
        
        my_sample_mats <- ls(pattern = paste0("^mat.combined.[q,u,N,d]"))
        
        if(my_grouping_names[i] == "log2_TPM_stage"){
                my_colors_composite <- my_color_palette1        
        } else {
                my_colors_composite <- my_color_palette3
        }
        
        plotComposite(my_sample_mats = my_sample_mats, 
                      my_title = "", 
                      site_label = my_site,
                      my_sub_range = 1:ncol(get(my_sample_mats[1])),
                      ylims = c(0,1.0),
                      my_binning = 1,
                      my_colors_composite = my_colors_composite,
                      line_lwd = 2, 
                      smoother = 251, 
                      yaxt = "s",
                      add_axis = FALSE,
                      log_scale = FALSE, 
                      add_line = FALSE, 
                      add_range = NULL)
        
        axis(side = 1,
             at = seq(1, ncol(get(my_sample_mats[1])), length.out = 5),
             labels = c("-5kb", "","","","+5kb"))
        
        axis(side = 1, line = 1, lwd = 0, lwd.ticks = 1,
             at = seq(1, ncol(get(my_sample_mats[1])), length.out = 5)[c(2,4)],
             labels = c("TSS","TTS"))
        
        abline(v =seq(1, ncol(get(my_sample_mats[1])), length.out = 5)[c(2,4)], lty = 3)
        
        
        ######################################################      
        
        
        mtext(text = "Normalized Coverage", side = 2, line = 3, cex = 1.25)
        
        
        ######################################################      
        
        
        n_rows <- sapply(my_sample_mats, function(x){  nrow(get(x)) })
        
        if(my_grouping_names[i] == "log2_TPM_stage"){
                names(n_rows) <- toupper(gsub("mat.combined.","", names(n_rows)))
        } else {
                names(n_rows) <- (gsub("mat.combined.","", names(n_rows)))
        }
        
        legend("topright", 
               legend = paste((rev(names(n_rows))), 
                              rev(n_rows), sep =" - "), 
               title = ifelse(my_grouping_names[i] == "log2_TPM_stage", 
                              paste("Exp. Levels",colnames(my_grouping)[2]), 
                              gsub("res_","DE ",my_grouping_names[i])),
               fill = rev(my_colors_composite), cex=0.7, bty="n")   
        
        if(my_grouping_names[i] == "log2_TPM_stage"){
                mtext(text = "*DBTMEE maternal genes removed", side = 1, outer = T, adj = 1, cex = 0.8)
                
        }
        
        mtext(text = paste0(my_assay," - ", my_stage), side = 3, outer = T, line = -1, font = 2, cex = 1.25)
        
        
        ######################################################      
        
        rm(list = ls(pattern = "^mat.combined"))
        
        
        
        dev.off()
        
        
}



##################################################################################################################################
##################################################################################################################################





# ##################################################################################################################################
# ##################################################################################################################################
# 
# ######################################################      Plotting     ######################################################### 
# 
# 
# 
# 
# # my_classes <- unique(my_DBTMEE$Cluster) 
# # my_class <- my_classes[1]
# 
# 
# 
# 
# my_site = "TSS"
# 
# 
# ######################################################      
# 
# 
# 
# pdf(file = output_file, width = 7.5, height = 4.5)
# 
# par(mfrow=c(1,2), oma=c(2,2,2,2), mar=c(4,4,4,2))
# 
# for(my_site in c("TSS", "TTS")){
#         
#         if(my_site == "TSS"){
#                 par(mar=c(4,4,4,1))
#         } else {
#                 par(mar=c(4,1,4,4))
#         }
#         
#         ######################################################      
#         
#         mat.combined.all <- get(grep(my_site, my_mats, value = TRUE))
#         
#         
#         for(qi in unique(log2_TPM_stage$Groups)){
#                 
#                 assign(paste0("mat.combined.",qi), 
#                        mat.combined.all[rownames(mat.combined.all) %in% log2_TPM_stage$gene_id[log2_TPM_stage$Groups == qi],])
#         }
#         
#         
#         ######################################################      
#         
#         my_sample_mats <- ls(pattern = paste0("^mat.combined.q[1-5]$"))
#         
#         
#         plotComposite(my_sample_mats = my_sample_mats, 
#                       my_title = "", 
#                       site_label = my_site,
#                       my_sub_range = 1:ncol(get(my_sample_mats[1])),
#                       ylims = c(0,1.0),
#                       my_binning = 1,
#                       my_colors_composite = my_colors_composite,
#                       line_lwd = 2, 
#                       smoother = 251, 
#                       yaxt = ifelse(my_site == "TSS", "s", "n"), 
#                       log_scale = FALSE, 
#                       add_line = FALSE, 
#                       add_range = NULL)
#         
#         abline(v = (ncol(get(my_sample_mats[1]))/2)+1, lty=2)
#         
#         ######################################################      
#         
#         if(my_site == "TSS"){
#                 mtext(text = "Normalized Coverage", side = 2, line = 3, cex = 1.25)
#         }
#         
#         ######################################################      
#         
#         if(my_site == "TTS"){
#                 
#                 n_rows <- sapply(my_sample_mats, function(x){  nrow(get(x)) })
#                 
#                 legend("topright", 
#                        legend = paste(toupper(rev(gsub("mat.combined.","", names(n_rows)))), 
#                                       rev(n_rows), sep =" - "), 
#                        title = colnames(log2_TPM_stage)[2],
#                        fill = rev(my_colors_composite), cex=0.7, bty="n")   
#                 
#                 mtext(text = "*DBTMEE maternal genes removed", side = 1, outer = T, adj = 1, cex = 0.8)
#                 
#                 mtext(text = paste0(my_assay," - ", my_stage), side = 3, outer = T, line = -1, font = 2, cex = 1.25)
#         }
#         
#         ######################################################      
#         
#         rm(list = ls(pattern = "^mat.combined"))
#         
# }
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# ##################################################################################################################################
# ##################################################################################################################################















##################################################################################################################################
##################################################################################################################################


