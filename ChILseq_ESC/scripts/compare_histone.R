






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

library(imputeLCMD)

library(org.Mm.eg.db)


##################################################################################################################################
##################################################################################################################################





args = commandArgs(trailingOnly=TRUE)



input_table <- grep("SampleTable", args, value = TRUE)



input_sessionInfo <- grep("results/sessionInfo", args, value = TRUE)
input_directory <- gsub("\\/sessionInfo.*","",input_sessionInfo)


output_sessionInfo <- grep("corrs.*sessionInfo", args, value = TRUE)
output_plots <- gsub("\\/sessionInfo.*","",output_sessionInfo)

my_scaling <- gsub(".*_", "", output_plots)
my_setting <- gsub(".*_", "", gsub(paste0("_", my_scaling), "", output_plots))

feature_type <- gsub("Output/deseq2_|/corrs_.*", "", output_plots)



##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################



load(paste0("genome/ranges_",feature_type,".rda"))


my_feature_ranges <- get(paste0("my_",feature_type,"_ranges"))

my_feature_width <- width(my_feature_ranges) / 1e3
names(my_feature_width) <- my_feature_ranges$gene_id


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




my_color_palette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")

my_color_palette2 <- brewer.pal(8,"Dark2")



##################################################################################################################################
##################################################################################################################################































##################################################################################################################################
##################################################################################################################################

#################################################         Read Table         #####################################################





my_lnc_files <- list.files(path = input_directory, pattern = "^lnc.*.log2.txt", full.names = TRUE)


for(my_lnc_file in my_lnc_files){
        
        my_lnc_name <- gsub(".*\\/|.txt","",my_lnc_file)
        
        my_lnc <- read.table(my_lnc_file, check.names=FALSE)
        my_lnc <- my_lnc[order(rownames(my_lnc)),]
        
        #my_lnc <- my_lnc[, c(grep("Zy", colnames(my_lnc)), grep("2C|8C|ES", colnames(my_lnc)))] 
        
        #my_lnc <- my_lnc[apply(my_lnc, 1, function(x){sum(is.na(x)) <= (length(x)/4)}),]
        
        # my_lnc <- my_lnc[apply(my_lnc, 1, function(x){
        #         all(aggregate(as.numeric(x), by = list(gsub("_[1-3]$","",colnames(my_lnc))), FUN=function(i){ sum(is.na(i)) <= 1})$x)
        # }),]
        
        
        if(my_setting %in% c("min")){
                
                set.seed(456)
                #my_lnc <- impute.MinProb(my_lnc, tune.sigma = 0.1)
                my_lnc[is.na(my_lnc)] <- rnorm(n = sum(is.na(my_lnc)), mean = min(my_lnc, na.rm = TRUE), sd = 0.1)
                
                
                
        } else if(my_setting %in% c("filt")){
                
                my_lnc <- na.omit(my_lnc)
                
                # my_lnc <- my_lnc[apply(my_lnc, 1, function(x){
                #         all(aggregate(as.numeric(x), by = list(gsub("_[1-3]$","",colnames(my_lnc))), FUN=function(i){ sum(is.na(i)) <= 1})$x)
                # }),]
        }
        
        
        assign(my_lnc_name, my_lnc)
        
        rm(list = "my_lnc")
        rm(list = "my_lnc_name")
        
}


my_assays <- unique(gsub("\\..*", "", gsub("lnc.","",ls(pattern = "^lnc\\."))))
log_types <- unique(gsub(".*\\.", "", ls(pattern = "^lnc\\.")))

my_assay <- my_assays[4]
log_type <- log_types[1]


##################################################################################################################################
##################################################################################################################################














##################################################################################################################################
##################################################################################################################################

#################################################       Histone Data         #####################################################


HistoneTable <- read.delim("../../Public/Histone_data/SampleTable_ALL.txt", header = T, stringsAsFactors = F)

HistoneTable <- HistoneTable[order(HistoneTable$SampleName),]
HistoneTable <- HistoneTable[!duplicated(HistoneTable$SampleName),]
HistoneTable$Conditions <- paste(HistoneTable$Stages, HistoneTable$Assay, sep = "_")



my_histones <- unique(HistoneTable$Assay)

my_histone_files <- list.files(path = paste0("../../Public/Histone_data//Output/deseq2_", feature_type,"/results/"), 
                               pattern = paste0("^lnc.*","log2.txt"), full.names = TRUE)



my_histone <- my_histones[2]


##################################################################################################################################
##################################################################################################################################


callback = function(hc, mat){
        sv = svd(t(mat))$v[,1]
        dend = (reorder(as.dendrogram(hc), wts = sv))
        as.hclust(dend)
}


#################################################        


my_lims <- data.frame(histone = c("ATACseq", "DNase-seq",  "H3K4me3", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3", "Pol2"),
                      len = c(0.8, 0.7, 0.85, 0.7, 0.8, 0.7, 0.9, 0.9),
                      zy = c(0.5, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.6), 
                      mean = c(0.5, 0.4, 0.7, 0.5, 0.7, 0.5, 0.7, 0.6), 
                      stringsAsFactors = FALSE)


#################################################        


for(my_histone in my_histones){
        
        hislnc_tmp <- read.table(grep(my_histone, my_histone_files, value = TRUE), check.names=FALSE)
        hislnc_tmp <- hislnc_tmp[,order(colnames(hislnc_tmp))]
        
        HistoneTable_tmp <- HistoneTable[HistoneTable$SampleName %in% colnames(hislnc_tmp),]
        HistoneTable_tmp <- HistoneTable_tmp[order(HistoneTable_tmp$SampleName),]
        HistoneTable_tmp$Stages <- factor(HistoneTable_tmp$Stages)
        
        #################################################           
        
        if(my_setting %in% c("min")){
                
                set.seed(456)
                hislnc_tmp[is.na(hislnc_tmp)] <- rnorm(n = sum(is.na(hislnc_tmp)), mean = min(hislnc_tmp, na.rm = TRUE), sd = 0.1)
                
                
        } else if(my_setting %in% c("filt")){
                
                stopifnot(identical(colnames(hislnc_tmp), HistoneTable_tmp$SampleName))
                
                hislnc_tmp <- na.omit(hislnc_tmp)
                
                # hislnc_tmp <- hislnc_tmp[apply(hislnc_tmp, 1, function(x){
                #         all(aggregate(as.numeric(x), by = list(HistoneTable_tmp$Stages), FUN=function(i){ sum(is.na(i)) <= 1})$x)
                # }),]
        }
        
        #################################################           
        
        if(identical(colnames(hislnc_tmp), HistoneTable_tmp$SampleName)){
                
                hislnc_repmeans_tmp <- t(apply(hislnc_tmp, 1, function(x){aggregate(x, by = list(HistoneTable_tmp$Stages), FUN = mean, na.rm = TRUE)$x}))
                
                colnames(hislnc_repmeans_tmp) <- levels(HistoneTable_tmp$Stages)
        }
        
        
        #################################################           
        
        
        if(my_scaling %in% c("mean")){
                
                hislnc_scaled_tmp <- hislnc_repmeans_tmp - rowMeans(hislnc_repmeans_tmp, na.rm = TRUE)  
                
        } else if(my_scaling %in% c("len")){
                
                my_feature_width_ordered <- my_feature_width[match(rownames(hislnc_repmeans_tmp), names(my_feature_width))]
                
                stopifnot(identical(rownames(hislnc_repmeans_tmp), names(my_feature_width_ordered)))
                hislnc_scaled_tmp <- hislnc_repmeans_tmp - log2(my_feature_width_ordered)
                
        }  else if(my_scaling %in% c("zy")){
                
                if("zygote" %in% colnames(hislnc_repmeans_tmp)){
                        
                        hislnc_scaled_tmp <- hislnc_repmeans_tmp - hislnc_repmeans_tmp[,grep("zygote", colnames(hislnc_repmeans_tmp))]  
                        
                } else if("E2C" %in% colnames(hislnc_repmeans_tmp)){
                        
                        hislnc_scaled_tmp <- hislnc_repmeans_tmp - hislnc_repmeans_tmp[,grep("E2C", colnames(hislnc_repmeans_tmp))]
                        
                } else if("oocyte" %in% colnames(hislnc_repmeans_tmp)){
                        
                        hislnc_scaled_tmp <- hislnc_repmeans_tmp - hislnc_repmeans_tmp[,grep("oocyte", colnames(hislnc_repmeans_tmp))]
                        
                } else {
                        hislnc_scaled_tmp <- hislnc_repmeans_tmp - rowMeans(hislnc_repmeans_tmp, na.rm = TRUE)  
                }
        }
        
        
        
        
        
        #################################################           
        
        
        for(log_type in log_types){
                
                for(my_assay in my_assays){
                        
                        if(my_assay == "all"){next()}
                        
                        
                        lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
                        
                        
                        my_conditions <- SampleTable$Conditions[match(colnames(lnc_tmp), SampleTable$ForeignID)]
                        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
                        
                        #################################################           
                        
                        
                        # lnc_tmp <- as.data.frame(t(apply(lnc_tmp, 1, function(x){aggregate(x, by= list(my_conditions), mean, na.rm = TRUE)[,2]})))
                        # colnames(lnc_tmp) <- paste(levels(my_conditions), my_assay, sep = "_")
                        
                        
                        if(my_scaling %in% c("mean")){
                                
                                lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp, na.rm = TRUE)
                                
                        } else if(my_scaling %in% c("len")){
                                
                                my_feature_width_ordered <- my_feature_width[match(rownames(lnc_tmp), names(my_feature_width))]
                                
                                stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_ordered)))
                                lnc_scaled_tmp <- lnc_tmp - log2(my_feature_width_ordered)
                                
                        }  else if(my_scaling %in% c("zy")){
                                
                                lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp[, grep("ChIP", colnames(lnc_tmp))], na.rm = TRUE)
                        }
                        
                        #################################################           
                        
                        
                        lnc_merged_tmp <- merge(hislnc_scaled_tmp, lnc_scaled_tmp, by = "row.names")
                        
                        
                        my_corr <- cor(lnc_merged_tmp[,-1], method = "spearman", use = "complete.obs")
                        my_corr <- my_corr[grepl(my_assay, rownames(my_corr)),]
                        my_corr <- my_corr[,!grepl(my_assay, colnames(my_corr))]
                        
                        my_corr[is.na(my_corr)] <- 0
                        
                        #################################################           
                        
                        
                        my_corr <- my_corr[,c(grep("oocyte", colnames(my_corr)),
                                              grep("PN3", colnames(my_corr)),
                                              grep("zygote", colnames(my_corr)),
                                              grep("PN5", colnames(my_corr)),
                                              grep("E2C", colnames(my_corr)),
                                              grep("^2C", colnames(my_corr)),
                                              grep("^Ctrl", colnames(my_corr)),
                                              grep("^DRB", colnames(my_corr)),
                                              grep("^4C", colnames(my_corr)),
                                              grep("^8C", colnames(my_corr)),
                                              grep("morula", colnames(my_corr)),
                                              grep("^ICM", colnames(my_corr)),
                                              grep("^TE", colnames(my_corr)))]
                        
                        #################################################           
                        
                        
                        pdf(paste0(output_plots,"/corr2",my_histone,".", my_assay,".", log_type,".pdf"), width = 10, height = 8, useDingbats = FALSE)
                        par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                        
                        
                        pheatmap(my_corr,
                                 #main = paste0("Spearman's R [n = ", nrow(lnc_merged_tmp),"]"),
                                 main = paste0(my_histone,"\nSpearman's R"),
                                 color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100),
                                 # breaks = seq(-max(abs(my_corr), na.rm = TRUE), 
                                 #               max(abs(my_corr), na.rm = TRUE), length.out = 101),
                                 # breaks = seq(ifelse(my_scaling == "len", -1, -0.7), 
                                 #              ifelse(my_scaling == "len",  1,  0.7), length.out = 101),
                                 breaks = seq(-my_lims[,my_scaling][my_lims$histone == my_histone],
                                              my_lims[,my_scaling][my_lims$histone == my_histone], length.out = 101),
                                 clustering_callback = callback,
                                 clustering_method = "ward.D2", 
                                 cluster_cols = FALSE, 
                                 cluster_rows = FALSE, 
                                 cellheight = ifelse(my_assay == "all", 15, 25), 
                                 cellwidth = ifelse(my_assay == "all", 15, 25))
                        
                        
                        rm(list = "my_corr")
                        rm(list = ls(pattern = "^lnc_.*_tmp"))
                        
                        dev.off()
                        
                        #################################################    
                        
                        
                        
                        #################################################    
                        
                        # lnc_repmeans_tmp <- t(apply(lnc_tmp, 1, function(x){
                        #         aggregate(x, by = list(my_conditions), FUN = mean, na.rm = TRUE)$x}))
                        # 
                        # colnames(lnc_repmeans_tmp) <- levels(my_conditions)
                        # 
                        # 
                        # lnc_scaled_tmp <- lnc_repmeans_tmp - rowMeans(lnc_repmeans_tmp, na.rm = TRUE)
                        # 
                        # lnc_merged_tmp <- merge(hislnc_scaled_tmp, lnc_scaled_tmp, by = "row.names")
                        # 
                        # 
                        # my_corr <- cor(lnc_merged_tmp[,-1], method = "spearman", use = "complete.obs")
                        # my_corr <- my_corr[grepl(my_assay, rownames(my_corr)),]
                        # my_corr <- my_corr[,!grepl(my_assay, colnames(my_corr))]
                        # 
                        # 
                        # #################################################           
                        # 
                        # my_corr <- my_corr[,c(grep("oocyte", colnames(my_corr)),
                        #                       grep("PN3", colnames(my_corr)),
                        #                       grep("zygote", colnames(my_corr)),
                        #                       grep("PN5", colnames(my_corr)),
                        #                       grep("E2C", colnames(my_corr)),
                        #                       grep("^2C", colnames(my_corr)),
                        #                       grep("^Ctrl", colnames(my_corr)),
                        #                       grep("^DRB", colnames(my_corr)),
                        #                       grep("^4C", colnames(my_corr)),
                        #                       grep("^8C", colnames(my_corr)),
                        #                       grep("morula", colnames(my_corr)),
                        #                       grep("^ICM", colnames(my_corr)),
                        #                       grep("^TE", colnames(my_corr))
                        # )]
                        # 
                        # 
                        # 
                        # #################################################           
                        # 
                        # pdf(paste0(output_plots,"/corr2",my_histone,".", my_assay,".", log_type,".rm.pdf"), width = 10, height = 8, useDingbats = FALSE)
                        # par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                        # 
                        # 
                        # pheatmap(my_corr,
                        #          #main = paste0("Spearman's R [n = ", nrow(lnc_merged_tmp),"]"),
                        #          main = paste0(my_histone,"\nSpearman's R"),
                        #          color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100),
                        #          breaks = seq(-0.4, 0.4, length.out = 101),
                        #          clustering_callback = callback,
                        #          clustering_method = "ward.D2", 
                        #          cluster_cols = FALSE,
                        #          cluster_rows = FALSE,
                        #          cellheight = ifelse(my_assay == "all", 15, 25), 
                        #          cellwidth = ifelse(my_assay == "all", 15, 25))
                        # 
                        # 
                        # rm(list = ls(pattern = "^lnc_.*tmp"))
                        # rm(list = "my_corr")
                        # 
                        # dev.off()
                        # 
                        # ################################################# 
                        
                }
        }
        
        rm(list = ls(pattern = "_.*tmp"))
        
        
}


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



my_DBTMEE_classes <- c("Maternal_RNA", "Major_ZGA", "MGA",  "Minor_ZGA",
                       "1-Cell_Transient", "2-Cell_Transient", "4-Cell_Transient")


my_DBTMEE_list <- list()

for(my_class in my_DBTMEE_classes){
        
        my_DBTMEE_list[[my_class]] <- my_DBTMEE$gene_id[my_DBTMEE$Cluster %in% my_class]
        
}


##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################        RNA-seq TPM         #####################################################




TPM <- read.table("../../Public/GSE38495_GSE45719/Output/TPM/TPM_means.txt", check.names = FALSE)
TPM <- TPM[,c("MII_oocyte","Zygote", "Early_2C", "Mid_2C", "Late_2C", "4C", "8C", "ICM","ESC")]

TPM[round(TPM, 6) == 0] <- NA

log2_TPM <- log2(TPM)




##################################################################################################################################
##################################################################################################################################

#################################################     RNA-seq Quantiles      #####################################################




log2_TPM_Quantiles <- log2_TPM[rownames(log2_TPM) %in% my_feature_ranges$gene_id, ]

log2_TPM_Q5 <- list()
log2_TPM_Q35 <- list()

for(RNAstage in colnames(log2_TPM_Quantiles)){
        
        Qstage <- paste0("Q.", RNAstage)
        
        log2_TPM_Quantiles[,Qstage] <- ifelse(is.na(log2_TPM_Quantiles[,RNAstage]), "q1", "exp")
        
        is_nonzero <- log2_TPM_Quantiles[,Qstage] == "exp"
        
        log2_TPM_Quantiles[is_nonzero, Qstage] <- c("q2","q3","q4","q5")[cut(x = log2_TPM_Quantiles[is_nonzero, RNAstage], 
                                                                             breaks = quantile(log2_TPM_Quantiles[is_nonzero, RNAstage]), 
                                                                             include.lowest = TRUE)]
        
        log2_TPM_Q5[[paste0("Q5_", RNAstage)]] <- rownames(log2_TPM_Quantiles)[log2_TPM_Quantiles[, Qstage] %in% c("q5")]
        
        log2_TPM_Q35[[paste0("Q35_", RNAstage)]] <- rownames(log2_TPM_Quantiles)[log2_TPM_Quantiles[, Qstage] %in% c("q3","q4","q5")]
        
        
}


log2_TPM_Quantiles <- log2_TPM_Quantiles[,grep("^Q", colnames(log2_TPM_Quantiles))]


##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################


if(my_setting %in% c("min")){
        
        set.seed(456)
        log2_TPM[is.na(log2_TPM)] <- rnorm(n = sum(is.na(log2_TPM)), mean = min(log2_TPM, na.rm = TRUE), sd = 0.1)
        
        
} else if(my_setting %in% c("filt")){
        
        log2_TPM <- na.omit(log2_TPM)
}


#################################################           


if(my_scaling %in% c("mean")){
        
        log2_scaled_TPM <- log2_TPM - rowMeans(log2_TPM, na.rm = TRUE)
        
} else if(my_scaling %in% c("len")){
        
        log2_scaled_TPM <- log2_TPM
        
}  else if(my_scaling %in% c("zy")){
        
        log2_scaled_TPM <- log2_TPM - log2_TPM[,"Zygote"]
}

#################################################  






##################################################################################################################################
##################################################################################################################################



my_classes_full <- c(log2_TPM_Q5 #,
                     #log2_TPM_Q35,
                     #my_DBTMEE_list
)


#################################################  

for(log_type in log_types){
        
        for(my_assay in rev(my_assays)){
                
                if(my_assay == "all"){next()}
                
                
                lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
                
                #lnc_tmp <- na.omit(lnc_tmp)
                
                my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
                my_conditions <- factor(my_conditions, levels = unique(my_conditions))
                

                #################################################           
                
                # lnc_tmp <- as.data.frame(t(apply(lnc_tmp, 1, function(x){aggregate(x, by= list(my_conditions), mean, na.rm = TRUE)[,2]})))
                # colnames(lnc_tmp) <- paste(levels(my_conditions), my_assay, sep = "_")
                
                
                if(my_scaling %in% c("mean")){
                        
                        lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp, na.rm = TRUE)
                        
                } else if(my_scaling %in% c("len")){
                        
                        my_feature_width_ordered <- my_feature_width[match(rownames(lnc_tmp), names(my_feature_width))]
                        
                        stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_ordered)))
                        lnc_scaled_tmp <- lnc_tmp - log2(my_feature_width_ordered)
                        
                }  else if(my_scaling %in% c("zy")){
                        
                        lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp[, grep("ChIP", colnames(lnc_tmp))], na.rm = TRUE)
                }
                
                
                #################################################  
                
                
                log2_TPM_merged_tmp <- merge(log2_scaled_TPM, lnc_scaled_tmp, by = "row.names")
                
                
                my_corr <- cor(log2_TPM_merged_tmp[,-1], method = "spearman", use = "complete.obs")
                my_corr <- my_corr[grepl(my_assay, rownames(my_corr)),]
                my_corr <- my_corr[,!grepl(my_assay, colnames(my_corr))]
                
                my_corr[is.na(my_corr)] <- 0
                
                #################################################           
                
                pdf(paste0(output_plots,"/corr2RNAseq.", my_assay,".", log_type,".pdf"), width = 6, height = 8, useDingbats = FALSE)
                par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                
                if(my_scaling == "len"){
                        hm_colors <- colorRampPalette(rev(brewer.pal(9,name = "RdBu")[1:5]))(100)
                        
                        hm_breaks <- seq(0.1, 0.75, length.out = 101)
                } else {
                        hm_colors <- colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100)
                        
                        hm_breaks <- seq(-0.4, 0.4, length.out = 101)
                }
                
                pheatmap(my_corr,
                         #main = paste0("Spearman's R [n = ", nrow(log2_TPM_merged_tmp),"]"),
                         main = "RNA-seq\nSpearman's R",
                         color = hm_colors,
                         breaks = hm_breaks,
                         clustering_callback = callback,
                         clustering_method = "ward.D2", 
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         cellheight = ifelse(my_assay == "all", 15, 25), 
                         cellwidth = ifelse(my_assay == "all", 15, 25))
                
                
                rm(list = "my_corr")
                
                dev.off()
                
                #################################################           
                
                
                
                
                
                my_xlims <- c(min(log2_TPM_merged_tmp[,grep(my_assay, colnames(log2_TPM_merged_tmp))], na.rm = TRUE),
                              max(log2_TPM_merged_tmp[,grep(my_assay, colnames(log2_TPM_merged_tmp))], na.rm = TRUE))
                
                
                my_ylims <- c(min(log2_TPM_merged_tmp[,c(-1,-grep(my_assay, colnames(log2_TPM_merged_tmp)))], na.rm = TRUE),
                              max(log2_TPM_merged_tmp[,c(-1,-grep(my_assay, colnames(log2_TPM_merged_tmp)))], na.rm = TRUE))
                
                
                my_class <- "4-Cell_Transient"
                
                for(my_NA in c("pairwise","complete")){
                        
                        if(my_NA == "complete"){
                                log2_TPM_merged_tmp <- na.omit(log2_TPM_merged_tmp)
                        }
                        
                        #for(my_class in names(my_classes_full)){
                        for(my_class in "Q5_ESC"){
                                
                                
                                png(paste0(output_plots,"/scatter2RNAseq.", my_assay,".",my_class,".",my_NA,".",log_type,".png"), width = 25, height = 25, units = "in", res = 200)
                                
                                par(mfrow=c(9,9), oma=c(5,5,5,5), mar=c(4,4,2,2), mgp=c(2.5,1,0))
                                
                                for(rid in paste(rep(c("ChIP","ChIL100C","ChIL1000C"), each=3),my_assay, 1:3, sep="_")){
                                        
                                        for(cid in c("MII_oocyte","Zygote","Early_2C","Mid_2C","Late_2C","4C","8C","ICM","ESC")){
                                                
                                                my_complete_cases <- complete.cases(log2_TPM_merged_tmp[,c(rid,cid)])
                                                
                                                plot(x = log2_TPM_merged_tmp[,rid],
                                                     y = log2_TPM_merged_tmp[,cid],
                                                     #main = paste("Spearman´s R =", round(my_cor,2)),
                                                     xlab = rid, ylab = cid,
                                                     xlim = my_xlims, 
                                                     ylim = my_ylims,
                                                     pch = 19, col = rgb(0.8,0.8,0.8,0.5), cex=0.5)
                                                
                                                my_subset <- log2_TPM_merged_tmp$Row.names %in% my_classes_full[[my_class]]
                                                
                                                points(x = log2_TPM_merged_tmp[my_subset,rid],
                                                       y = log2_TPM_merged_tmp[my_subset,cid],
                                                       pch = 19, col = rgb(0.8,0,1,0.25), cex=0.5)   
                                                
                                                lines(lowess(x = log2_TPM_merged_tmp[my_complete_cases,rid],
                                                             y = log2_TPM_merged_tmp[my_complete_cases,cid],
                                                             f = 1),
                                                      col = rgb(0.5,0.5,0.5,1), lwd=2)
                                                
                                                lines(lowess(x = log2_TPM_merged_tmp[(my_subset & my_complete_cases),rid],
                                                             y = log2_TPM_merged_tmp[(my_subset & my_complete_cases),cid],
                                                             f = 1),
                                                      col = rgb(0.6,0,0.8,1), lwd=2)
                                                
                                                points(x = mean(log2_TPM_merged_tmp[my_complete_cases,rid]),
                                                       y = mean(log2_TPM_merged_tmp[my_complete_cases,cid]),
                                                       pch = 19, col = rgb(0.5,0.5,0.5,1), cex=1.25) 
                                                
                                                points(x = mean(log2_TPM_merged_tmp[(my_subset & my_complete_cases),rid]),
                                                       y = mean(log2_TPM_merged_tmp[(my_subset & my_complete_cases),cid]),
                                                       pch = 19, col = rgb(0.6,0,0.8,1), cex=1.25) 
                                                
                                                
                                                my_cor <- cor(x = log2_TPM_merged_tmp[,rid],
                                                              y = log2_TPM_merged_tmp[,cid], 
                                                              method = "spearman", use = "complete.obs")
                                                
                                                my_cor_sub <- cor(x = log2_TPM_merged_tmp[my_subset,rid],
                                                                  y = log2_TPM_merged_tmp[my_subset,cid], 
                                                                  method = "spearman", use = "complete.obs")
                                                
                                                legend("topleft", horiz = FALSE, 
                                                       legend = c(paste0("Rs = ", round(my_cor,2),
                                                                         "; n = ", sum(my_complete_cases)),
                                                                  paste0("Rs = ", round(my_cor_sub,2),
                                                                         "; n = ", sum(my_complete_cases & my_subset))),
                                                       fill = c(rgb(0.7,0.7,0.7,0.5), rgb(0.7,0,0.9,0.75)), cex =0.7)
                                                
                                                rm(list = "my_cor")
                                                
                                        }
                                }
                                
                                mtext(text = "RNA-seq log2 TPM", side = 2, line = 1, font = 2,  outer = TRUE, cex = 1.2)
                                mtext(text = "log2 Normalized Counts", side = 1, line = 1, font = 2,  outer = TRUE, cex = 1.2)
                                
                                par(fig=c(0,1,0,1), mar=c(1,1,1,1), oma=c(0,0,0,0), new=TRUE)        
                                plot.new()
                                
                                legend("top", horiz = TRUE, 
                                       legend = c(paste0("all genes - ", length(my_subset)),
                                                  paste0(my_class, " - ", sum(my_subset))),
                                       fill = c(rgb(0.7,0.7,0.7,0.5), rgb(0.7,0,0.9,0.75)), cex = 1.5)
                                
                                dev.off()
                                
                        }
                }
                
                
                ################################################# 
                
                
                rm(list = ls(pattern = "_tmp"))
                
                
        }
}


##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################




















##################################################################################################################################
##################################################################################################################################







writeLines(capture.output(sessionInfo()), output_sessionInfo)






##################################################################################################################################
##################################################################################################################################


