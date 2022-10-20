






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



matrix_files <- grep("S[2,5]P.*_T[T,S]S_matrix.rds", args, value = TRUE)
matrix_files <- matrix_files[order(matrix_files)]


stopifnot(length(matrix_files) == 18)


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





SampleMatchTable <- data.frame(Stages = c("Zy","2C","2C_DRB","8C","ES"),
                               RNAseq = c("Zygote", "Late_2C", "Late_2C", "8C",  "ESC" ),
                               Colors = c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73"), 
                               stringsAsFactors = FALSE)


SampleMatchTable$Dark <- sapply(SampleMatchTable$Colors, function(x){ colorRampPalette(c("white",x,"black"))(21)[12]})



##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################

######################################################     Read Mats    ########################################################## 


i=1

for(i in seq(1, length(matrix_files), 2)){
        
        my_name1 <- gsub(".*matrix_ave/|_matrix.*","",matrix_files[i])
        my_name2 <- gsub(".*matrix_ave/|_matrix.*","",matrix_files[i+1])
        
        stopifnot(identical(gsub("TSS","",my_name1), gsub("TTS","",my_name2)))
        
        my_mat1 <- readRDS(matrix_files[i])
        my_mat2 <- readRDS(matrix_files[i+1])
        
        my_site <- gsub(".*_","",my_name1)
        my_name <- paste0("GB",".",gsub(paste0("_",my_site),"", my_name1),"_A")
        
        stopifnot(identical(rownames(my_mat1), rownames(my_mat2)))
        
        my_mat <- cbind(my_mat1, NA, my_mat2)
        
        assign(my_name, my_mat)
        
        rm(list = "my_name")
        rm(list = "my_mat")
}


my_mats <- ls(pattern = "GB.*S[2,5]P.*_A$")


stopifnot(length(my_mats) > 0)

for(i in seq_along(my_mats)){
        
        stopifnot(identical(rownames(get(my_mats[1])), rownames(get(my_mats[i]))))    
}




##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

######################################################      Quantiles     ######################################################## 




log2_TPM <- log2_TPM[!(log2_TPM$gene_id %in% my_DBTMEE$gene_id[my_DBTMEE$Cluster == "Maternal_RNA"]), ]
log2_TPM <- log2_TPM[  log2_TPM$gene_id %in% rownames(get(my_mats[1])), ]


RNAstage <- SampleMatchTable$RNAseq[1]

for(RNAstage in unique(SampleMatchTable$RNAseq)){
        
        Qstage <- paste0("Q.", RNAstage)
        
        log2_TPM[,Qstage] <- ifelse(round(log2_TPM[,RNAstage], 4) == 0.0000, "q1", "exp")
        
        is_nonzero <- log2_TPM[,Qstage] == "exp"
        
        log2_TPM[is_nonzero, Qstage] <- c("q2","q3","q4","q5")[cut(x = log2_TPM[is_nonzero, RNAstage], 
                                                                   breaks = quantile(log2_TPM[is_nonzero, RNAstage]), 
                                                                   include.lowest = TRUE)]
}


##################################################################################################################################
##################################################################################################################################




my_name <- my_mats[1]


for(my_name in my_mats){
        
        my_mat <- get(my_name)
        
        my_stage <- gsub(paste0("GB.|_S[5,2]P|_A"),"", my_name)
        my_name <- paste0(my_name, ".Q5")
        
        RNAqid <- paste0("Q.",SampleMatchTable$RNAseq[SampleMatchTable$Stages == my_stage])
        
        my_mat_sub <- my_mat[rownames(my_mat) %in% log2_TPM$gene_id[log2_TPM[,RNAqid] == "q5"],]
        
        assign(my_name, my_mat_sub)
        
        rm(list = "my_mat_sub")
        rm(list = "my_name")
        rm(list = "my_mat")       
}





for(my_name in my_mats){
    
    my_mat <- get(my_name)
    
    my_name <- paste0(my_name, ".HK")
    
    my_mat_sub <- my_mat[rownames(my_mat) %in% log2_TPM$gene_id[log2_TPM$Q.Late_2C == "q5" & log2_TPM$Q.8C == "q5"],]
    
    assign(my_name, my_mat_sub)
    
    rm(list = "my_mat_sub")
    rm(list = "my_name")
    rm(list = "my_mat")       
}


##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################

######################################################      Plotting     ######################################################### 


# my_classes <- unique(my_DBTMEE$Cluster) 
# my_class <- my_classes[1]




my_site = "GB"



my_plot_mats <- list(c("GB.2C_S5P_A.Q5", "GB.2C_S5P_DRB_A.Q5"),
                     c("GB.Zy_S5P_A.Q5","GB.2C_S5P_A.Q5","GB.8C_S5P_A.Q5", "GB.ES_S5P_A.Q5"),
                     c("GB.Zy_S2P_A.Q5","GB.2C_S2P_A.Q5","GB.8C_S2P_A.Q5", "GB.ES_S2P_A.Q5"))



######################################################      



for(my_style in c("comb","sep")){
        
        
        if(my_style == "comb"){
                range_end <- ncol(get(my_plot_mats[[1]][1]))
        } else {
                
                range_end <- 10001
        }
        
        for(my_scaling in c("scaled","unscaled")){

                
                
                pdf(file = gsub("2.pdf", paste0("_1_",my_style,".pdf"), gsub("scaled",my_scaling,output_file)), 
                    width = ifelse(my_style == "comb", 6.5, 4.5), 
                    height = 4.5)
                

                
                
                for(my_sample_mats in my_plot_mats){
                        
                        my_site <- unique(gsub("\\..*","",my_sample_mats))
                        my_assay <- unique(gsub("_.*","",gsub(".*_S","S",my_sample_mats)))
                        
                        stopifnot(length(my_site) == 1)        
                        stopifnot(length(my_assay) == 1)        
                        
                        
                        par(mfrow=c(1,1), oma=c(2,2,2,2), mar=c(4,4,4,2))
                        par(mar=c(4,4,4,1))
                        
                        
                        ######################################################      
                        
                        
                        my_sample_names <- gsub(paste0(my_site,".|_",my_assay,"|_A.Q5"),"", my_sample_mats)
                        
                        my_colors_composite <- SampleMatchTable$Dark[match(my_sample_names, SampleMatchTable$Stages, nomatch = FALSE)]
                        
                        if(my_style == "comb"){
                                
                        plotComposite(my_sample_mats = my_sample_mats, 
                                      my_title = paste0(my_assay," at Q5 genes"), 
                                      site_label = my_site,
                                      my_sub_range = 1:range_end,
                                      ylims = c(0, ifelse(my_scaling == "scaled", 1, 0.9)),
                                      my_binning = 1,
                                      my_colors_composite = my_colors_composite,
                                      line_lwd = 2, 
                                      smoother = ifelse(my_scaling == "scaled", 501, 251), 
                                      yaxt = "s",
                                      add_axis = FALSE,
                                      min_max_scale = my_scaling == "scaled",
                                      log_scale = FALSE, 
                                      add_line = FALSE, 
                                      add_range = NULL)

                                axis(side = 1,
                                     at = seq(1, range_end, length.out = 5),
                                     labels = c("-5kb", "","","","+5kb"))
                                
                                axis(side = 1, line = 1, lwd = 0, lwd.ticks = 1,
                                     at = seq(1, range_end, length.out = 5)[c(2,4)],
                                     labels = c("TSS","TTS"))    
                                
                                abline(v =seq(1, range_end, length.out = 5)[c(2,4)], lty = 3)
                                
                                
                                ######################################################      
                                
                                if(my_scaling == "scaled"){
                                        mtext(text = "Scaled Coverage", side = 2, line = 3, cex = 1.25)
                                } else {
                                        mtext(text = "Normalized Coverage", side = 2, line = 3, cex = 1.25)
                                }
                                
                                ######################################################    
                                
                        } else {
                                
                                plotComposite(my_sample_mats = my_sample_mats, 
                                              my_title = paste0(my_assay," at Q5 genes"), 
                                              site_label = my_site,
                                              my_sub_range = 1:range_end,
                                              ylims = c(0, ifelse(my_scaling == "scaled", 1, 0.9)),
                                              my_binning = 1,
                                              my_colors_composite = my_colors_composite,
                                              line_lwd = 2, 
                                              smoother = ifelse(my_scaling == "scaled", 501, 251), 
                                              yaxt = "s",
                                              add_axis = FALSE,
                                              min_max_scale = my_scaling == "scaled",
                                              log_scale = FALSE, 
                                              add_line = FALSE, 
                                              add_range = NULL)
                                
                                axis(side = 1,
                                     at = seq(1, range_end, length.out = 3),
                                     labels = c("-5kb", "TSS","+5kb"))
                                
                                abline(v = (range_end/2)+1, lty=2)
                                
                                ######################################################      
                                
                                if(my_scaling == "scaled"){
                                        mtext(text = "Scaled Coverage", side = 2, line = 3, cex = 1.25)
                                } else {
                                        mtext(text = "Normalized Coverage", side = 2, line = 3, cex = 1.25)
                                }
                                
                                ######################################################   
                                
                                plotComposite(my_sample_mats = my_sample_mats, 
                                              my_title = paste0(my_assay," at Q5 genes"), 
                                              site_label = my_site,
                                              my_sub_range = range_end:ncol(get(my_plot_mats[[1]][1])),
                                              ylims = c(0, ifelse(my_scaling == "scaled", 1, 0.9)),
                                              my_binning = 1,
                                              my_colors_composite = my_colors_composite,
                                              line_lwd = 2, 
                                              smoother = ifelse(my_scaling == "scaled", 501, 251), 
                                              yaxt = "s",
                                              add_axis = FALSE,
                                              min_max_scale = my_scaling == "scaled",
                                              log_scale = FALSE, 
                                              add_line = FALSE, 
                                              add_range = NULL)
                                
                                axis(side = 1,
                                     at = seq(1, range_end, length.out = 3),
                                     labels = c("-5kb", "TTS","+5kb"))
                                
                                abline(v = (range_end/2)+1, lty=2)
                                
                                ######################################################   
                                
                        }
                        
                        
                        
                        ######################################################      
                        
                        
                        mtext(text = "*DBTMEE maternal genes removed", side = 1, outer = T, adj = 1, cex = 0.8)
                        
                        
                        n_rows <- sapply(my_sample_mats, function(x){  nrow(get(x)) })
                        
                        legend("topright", 
                               legend = paste0(gsub("_A.Q5"," - ", 
                                                    gsub(paste0("_",my_assay),"",
                                                         gsub(paste0(my_site,"."),"", names(n_rows)))),
                                               n_rows), 
                               title = "Exp. Q5 genes",
                               fill = my_colors_composite, cex=0.7, bty="n")   
                        
                        ######################################################
                }
                
                dev.off()
        }
}


##################################################################################################################################
##################################################################################################################################




my_color_palette2 <- brewer.pal(8,"Dark2")




my_plot_mats <- list(c("GB.Zy_S5P_A.Q5","GB.Zy_S2P_A.Q5"),
                     c("GB.2C_S5P_A.Q5","GB.2C_S2P_A.Q5"),
                     c("GB.8C_S5P_A.Q5","GB.8C_S2P_A.Q5"),
                     c("GB.ES_S5P_A.Q5","GB.ES_S2P_A.Q5"),
                     c("GB.Zy_S5P_A.HK","GB.Zy_S2P_A.HK"),
                     c("GB.2C_S5P_A.HK","GB.2C_S2P_A.HK"),
                     c("GB.8C_S5P_A.HK","GB.8C_S2P_A.HK"),
                     c("GB.ES_S5P_A.HK","GB.ES_S2P_A.HK"))



######################################################      


for(my_scaling in c("scaled","unscaled")){
        
        pdf(file = gsub("2.pdf","2.pdf", gsub("scaled",my_scaling,output_file)), width = 6.5, height = 4.5)
        
        
        
        for(my_sample_mats in my_plot_mats){
                
                
                my_site <- unique(gsub("\\..*","",my_sample_mats))
                my_assay <- unique(gsub("_.*","",gsub(".*_S","S",my_sample_mats)))
                my_stage <- unique(gsub("_.*","",gsub(paste0(my_site,"."),"", my_sample_mats)))
                
                stopifnot(length(my_site) == 1)
                stopifnot(length(my_stage) == 1)        
                stopifnot(length(my_assay) == 2)  
                
                
                
                par(mfrow=c(1,1), oma=c(2,2,2,2), mar=c(4,4,4,2))
                par(mar=c(4,4,4,1))
                
                
                ######################################################      
                
                
                my_colors_composite <- my_color_palette2[2:3]
                
                
                plotComposite(my_sample_mats = my_sample_mats, 
                              my_title = paste0(my_stage," ChIL-seq at Q5 genes"), 
                              site_label = my_site,
                              my_sub_range = 1:ncol(get(my_sample_mats[1])),
                              ylims = c(0,ifelse(my_scaling == "scaled", 1, 0.9)),
                              my_binning = 1,
                              my_colors_composite = my_colors_composite,
                              line_lwd = 2, 
                              smoother = ifelse(my_scaling == "scaled", 501, 251), 
                              yaxt = "s", 
                              add_axis = FALSE,
                              min_max_scale = my_scaling == "scaled",
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
                
                
                mtext(text = "*DBTMEE maternal genes removed", side = 1, outer = T, adj = 1, cex = 0.8)
                
                
                n_rows <- sapply(my_sample_mats, function(x){  nrow(get(x)) })
                
                if(all(grepl("HK", names(n_rows)))){
                    
                    legend("topright", 
                           legend = paste0(gsub("_A.HK"," - ", 
                                                gsub(".*_S","S",
                                                     gsub(paste0(my_site,"."),"", names(n_rows)))),
                                           n_rows), 
                           title = "Q5 shared in 2C and 8C",
                           fill = my_colors_composite, cex=0.7, bty="n")  
                    
                } else {
                    
                    legend("topright", 
                           legend = paste0(gsub("_A.Q5"," - ", 
                                                gsub("_S"," S",
                                                     gsub(paste0(my_site,"."),"", names(n_rows)))),
                                           n_rows), 
                           title = paste0(my_stage," Exp. Q5 genes"),
                           fill = my_colors_composite, cex=0.7, bty="n")  
                }  
                
                
                
                
                ######################################################      
                
                
                if(my_scaling == "scaled"){
                        mtext(text = "Scaled Coverage", side = 2, line = 3, cex = 1.25)
                } else {
                        mtext(text = "Normalized Coverage", side = 2, line = 3, cex = 1.25)
                }
                
        }
        
        
        dev.off()
        
}


##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################


