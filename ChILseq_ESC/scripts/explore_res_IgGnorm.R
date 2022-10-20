






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


output_sessionInfo <- grep("plots.*sessionInfo", args, value = TRUE)
output_plots <- gsub("\\/sessionInfo.*","",output_sessionInfo)

my_setting <- gsub(".*_", "", output_plots)

feature_type <- gsub("Output/deseq2_|/plotsnew_.*", "", output_plots)



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




my_color_palette <- colorRampPalette(c("white","#009E73","black"))(21)[c(5,10,15)]

my_color_palette2 <- brewer.pal(8,"Dark2")



##################################################################################################################################
##################################################################################################################################
















##################################################################################################################################
##################################################################################################################################

#################################################       Read Results        #####################################################



my_res_files <- list.files(path = input_directory, pattern = "^res", full.names = TRUE)



for(my_res_file in my_res_files){
    
    my_res_name <- gsub(".*\\/|.txt","",my_res_file)
    
    my_res <- read.table(my_res_file)
    
    my_res$gene_id <- rownames(my_res)
    
    assign(my_res_name, my_res)
    
    rm(list = "my_res")
    rm(list = "my_res_name")
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

#################################################         RNA-seq TPM        #####################################################



TPM <- read.table("../../Public/GSE38495_GSE45719/Output/TPM/TPM_means.txt", check.names = FALSE)
TPM <- TPM[,c("MII_oocyte","Zygote", "Early_2C", "Mid_2C", "Late_2C", "4C", "8C", "ICM","ESC")]


log2_TPM <- log2(TPM+1)

log2_scaled_TPM <- log2_TPM - rowMeans(log2_TPM)

log2_zscore_TPM <- t(scale(t(log2_TPM)))




##################################################################################################################################
##################################################################################################################################

#################################################     RNA-seq Quantiles      #####################################################


log2_TPM_Quantiles <- log2_TPM[rownames(log2_TPM) %in% my_feature_ranges$gene_id, ]

log2_TPM_Q5 <- list()
log2_TPM_Q35 <- list()


for(RNAstage in colnames(log2_TPM_Quantiles)){
    
    Qstage <- paste0("Q.", RNAstage)
    
    log2_TPM_Quantiles[,Qstage] <- ifelse(round(log2_TPM_Quantiles[,RNAstage], 4) == 0.0000, "q1", "exp")
    
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

#################################################         Read Table         #####################################################





my_lnc_files <- list.files(path = input_directory, pattern = "^lnc.all.log2.txt", full.names = TRUE)


for(my_lnc_file in my_lnc_files){
    
    my_lnc_name <- gsub(".*\\/|.txt","",my_lnc_file)
    
    my_lnc <- read.table(my_lnc_file, check.names=FALSE)
    my_lnc <- my_lnc[order(rownames(my_lnc)),]
    
    my_lnc <- my_lnc[, c(grep("100C", colnames(my_lnc)), grep("1000C|ChIP", colnames(my_lnc)))] 
    
    assign(my_lnc_name, my_lnc)
    
    rm(list = "my_lnc")
    rm(list = "my_lnc_name")
    
}



##################################################################################################################################
##################################################################################################################################






##################################################################################################################################
##################################################################################################################################

#################################################      IgG normalization     #####################################################


lnc.S5P.log2 <- lnc.all.log2[,grep("S5P", colnames(lnc.all.log2))]
lnc.S2P.log2 <- lnc.all.log2[,grep("S2P", colnames(lnc.all.log2))]
lnc.IgG.log2 <- lnc.all.log2[,grep("IgG", colnames(lnc.all.log2))]



stopifnot(identical(rownames(lnc.S5P.log2), rownames(lnc.IgG.log2)))
stopifnot(identical(rownames(lnc.S2P.log2), rownames(lnc.IgG.log2)))


#################################################        


my_IgG_stages <- gsub("_.*","", colnames(lnc.IgG.log2))
my_IgG_stages <- factor(my_IgG_stages, levels = unique(my_IgG_stages))


#################################################        


my_S5P_stages <- gsub("_.*","", colnames(lnc.S5P.log2))
my_S5P_stages <- factor(my_S5P_stages, levels = unique(my_S5P_stages))


lnc.S5P.IgGnorm <- t(sapply(1:nrow(lnc.S5P.log2), function(i){
    
    unlist(sapply(levels(my_S5P_stages), function(j){ 
        
        lnc.S5P.log2[i, my_S5P_stages == j] - mean(as.numeric(lnc.IgG.log2[i, my_IgG_stages == j]), na.rm=TRUE) 
    }))
}))


colnames(lnc.S5P.IgGnorm) <- colnames(lnc.S5P.log2)
rownames(lnc.S5P.IgGnorm) <- rownames(lnc.S5P.log2)



#################################################        

my_S2P_stages <- gsub("_.*","", colnames(lnc.S2P.log2))
my_S2P_stages <- factor(my_S2P_stages, levels = unique(my_S2P_stages))


lnc.S2P.IgGnorm <- t(sapply(1:nrow(lnc.S2P.log2), function(i){
    
    unlist(sapply(levels(my_S2P_stages), function(j){ 
        
        lnc.S2P.log2[i, my_S2P_stages == j] - mean(as.numeric(lnc.IgG.log2[i, my_IgG_stages == j]), na.rm=TRUE) 
    }))
}))


colnames(lnc.S2P.IgGnorm) <- colnames(lnc.S2P.log2)
rownames(lnc.S2P.IgGnorm) <- rownames(lnc.S2P.log2)


#################################################        






##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

#################################################        Correlation         #####################################################




my_assays <- unique(gsub("\\..*", "", gsub("lnc.","",ls(pattern = "^lnc.*norm"))))
log_types <- unique(gsub(".*\\.", "", ls(pattern = "^lnc.*norm")))


my_assay <- "S5P"
log_type <- "IgGnorm"


callback = function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = (reorder(as.dendrogram(hc), wts = sv))
    as.hclust(dend)
}



for(my_assay in my_assays){
    
    for(log_type in log_types){
        
        
        pdf(paste0(output_plots,"/correlation.", my_assay,".", log_type,".pdf"), width = 12, height = 12, useDingbats = FALSE)
        par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
        
        
        lnc_tmp <- get(paste("lnc", my_assay, log_type, sep = "."))
        
        lnc_scaled_tmp <- lnc_tmp
        
        my_breaks <- seq(-1,1, length.out = 101)
        hm_colors <- colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100)
        
        
        
        my_corr <- cor(lnc_scaled_tmp, method = "spearman", use = "complete.obs")
        rownames(my_corr) <- SampleTable$ForeignID[match(colnames(my_corr), SampleTable$ForeignID)]
        
        pheatmap(my_corr,
                 main = paste0("Spearman's R"),
                 color = hm_colors,
                 breaks = my_breaks,
                 clustering_callback = callback,
                 clustering_method = "ward.D2", 
                 cellheight = ifelse(my_assay == "all", 15, 25), 
                 cellwidth = ifelse(my_assay == "all", 15, 25))
        
        
        rm(list = ls(pattern = "_tmp"))
        rm(list = "my_corr")     
        
        
        
        dev.off()
        
    }
}



##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################            PCA             #####################################################





for(log_type in log_types){
    
    
    pdf(paste0(output_plots,"/PCA.", log_type,".pdf"), height = 6.25, width = 9, useDingbats = F)
    par(mfcol=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))
    
    
    #########################################################        
    
    
    for(i in seq_along(my_assays)){
        
        my_assay <- my_assays[i]
        
        if(my_assay == "all"){next()}
        
        lnc_tmp <- get(paste("lnc", my_assay, log_type, sep = "."))
        
        #################################################       
        
        my_conditions <- SampleTable[,"Stages"][match(colnames(lnc_tmp), SampleTable$ForeignID)]
        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
        
        pca_colors <- my_color_palette  
        
        
        #################################################       
        
        plottingPCA(lnc_tmp,
                    xcomp = 1,
                    ycomp = 2,
                    conditions = my_conditions,
                    pca_colors = pca_colors,
                    main_title = paste("PCA -", my_assay),
                    quantiles = c(0,1),
                    show_labels = FALSE,
                    point_size = 1.1,
                    my_xlimits = c(-200,200),
                    my_ylimits = c(-200,200))
        
        #################################################       
        
        plot.new()
        
        legend("top",
               legend = levels(my_conditions),
               cex = 1, bty = "n",
               fill =  pca_colors[seq_along(levels(my_conditions))])
        
        rm(list = ls(pattern = "_tmp"))
    }
    
    
    #########################################################        

    
    rm(list = ls(pattern = "_tmp"))
    
    
    dev.off()
    
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


