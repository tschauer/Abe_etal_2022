






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




my_color_palette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")

my_color_palette2 <-  brewer.pal(9,"Set1")[c(7,4,8)]



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

padj_cutoff <- 0.05

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
    
    my_lnc <- my_lnc[, c(grep("Zy", colnames(my_lnc)), grep("2C|8C|ES", colnames(my_lnc)))] 
    
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
        
        lnc_scaled_tmp <- lnc_tmp #- rowMeans(lnc_tmp)

        my_corr <- cor(lnc_scaled_tmp, method = "spearman", use = "complete.obs")
        rownames(my_corr) <- SampleTable$ForeignID[match(colnames(my_corr), SampleTable$ForeignID)]
        
        
        pheatmap(my_corr,
                 #main = paste0("Spearman's R [n = ", nrow(lnc_scaled_tmp),"]"),
                 main = "Spearman's R",
                 color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100),
                 breaks = seq(-1,1, length.out = 101),
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
        
        if(my_assay == "S5P"){
            pca_colors <- my_color_palette  
        } else {
            pca_colors <- my_color_palette[-3]
        }
        
        #################################################       
        
        plottingPCA(lnc_tmp,
                    xcomp = 1,
                    ycomp = 2,
                    conditions = my_conditions,
                    pca_colors = pca_colors,
                    main_title = paste("PCA -", my_assay),
                    quantiles = c(0,1),
                    show_labels = FALSE,
                    point_size = 1.25,
                    my_xlimits = c(-100,100),
                    my_ylimits = c(-100,100))
        
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

#################################################       Signal Violin       #####################################################



for(log_type in log_types){
    
    for(my_assay in my_assays){
        
        if(my_assay == "all"){next()}
        
        lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
        
        my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
        
        #################################################           
        
        if(my_assay == "S5P"){
            vio_colors <- my_color_palette  
        } else {
            vio_colors <- my_color_palette[-3]
        }
        
        #################################################           
        
        
        lnc_scaled_tmp <-  lnc_tmp - rowMeans(lnc_tmp[, grep("Zy", colnames(lnc_tmp))], na.rm = TRUE)
        lnc_scaled_tmp$gene_id <- rownames(lnc_scaled_tmp)
        
        lnc_tmp_long <- melt(lnc_scaled_tmp)
        
        lnc_tmp_long$Conditions <- gsub(paste0("X|_",my_assay,"|_[1-3]"),"", lnc_tmp_long$variable)
        lnc_tmp_long$Conditions <- factor(lnc_tmp_long$Conditions, levels = unique(lnc_tmp_long$Conditions))
        
        lnc_tmp_long$Classes <- my_DBTMEE$Cluster[match(lnc_tmp_long$gene_id, my_DBTMEE$gene_id)]
        
        ggp1 <- ggplot(lnc_tmp_long, aes(x=variable, y=value, fill = Conditions)) + 
            theme_bw() +
            theme(text = element_text(size = 12), aspect.ratio = 0.75, 
                  axis.text.x = element_blank(), 
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  plot.title = element_text(hjust = 0.5, face = "bold")) +
            geom_violin(width = 1.5) +
            facet_wrap(~ Classes, nrow = 3) +
            ggtitle(label = paste0(feature_type," - ",my_assay)) +
            stat_summary(fun =  "median", geom="point", size=1) +
            scale_fill_manual(values = vio_colors) +
            xlab(label = "Conditions") +
            coord_cartesian(ylim = c(-4, 4)) +
            ylab(label = "log2 Normalized Counts (relative to Zy)")
        
        
        
        
        
        ggsave(filename = paste0(output_plots,"/vioplot.",my_assay,".", log_type,".pdf"),
               plot =  ggp1,
               width = 12, height = 8)
        
        
        #################################################           
        
        
        my_subclasses_list <- list(c("Maternal_RNA","1-Cell_Transient", "2-Cell_Transient", "4-Cell_Transient"),
                                   c("Maternal_RNA", "Minor_ZGA",  "Major_ZGA", "MGA"))
        
        
        for(i in seq_along(my_subclasses_list)){
            
            lnc_tmp_long_sub <- lnc_tmp_long[lnc_tmp_long$Classes %in% my_subclasses_list[[i]],]
            
            lnc_tmp_long_sub$group <- factor(paste(lnc_tmp_long_sub$Classes, gsub(".*_","",lnc_tmp_long_sub$variable)),
                                             levels = (paste(rep(my_subclasses_list[[i]], each=3), 1:3)))
            
            
            ggp2 <- ggplot(lnc_tmp_long_sub, 
                           aes(x=group, y=value, fill = Conditions)) + 
                theme_bw() +
                theme(text = element_text(size = 12), aspect.ratio = 1, 
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      strip.background = element_blank(),
                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                geom_violin(width = 1) +
                facet_wrap(~ Conditions, nrow = 1) +
                ggtitle(label = paste0(feature_type," - ",my_assay)) +
                stat_summary(fun =  "median", geom="point", size=1) +
                scale_fill_manual(values = vio_colors) +
                xlab(label = "DBTMEE Classes") +
                coord_cartesian(ylim = c(-4, 4)) +
                ylab(label = "log2 Normalized Counts (relative to Zy)") + 
                geom_hline(yintercept = 0,  color="black", linetype="dashed") +
                geom_vline(xintercept = c(3.5,6.5,9.5),  color="grey", linetype="dotted") 
            
            rm(list = "lnc_tmp_long_sub")
            
            ggsave(filename = paste0(output_plots,"/vioplot",i,".",my_assay,".", log_type,".pdf"),
                   plot =  ggp2,
                   width = 12, height = 6)
        }
        
        #################################################           
        
        rm(list = ls(pattern = "_tmp"))
    }
    
}



##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################

##############################################    RNA-seq Cdk9in Spt5 results   ###################################################



res_external_files <- c("../RNAseq_Cdk9_Spt5//Output/deseq2/tables/res.Cdk9in-Noinject.txt",
                        "../RNAseq_Cdk9_Spt5//Output/deseq2/tables/res.Spt5KD-IgG.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.4C-Late_2C.txt")

res_external_names <- gsub(".*\\/|.txt", "", res_external_files)


for(i in seq_along(res_external_files)){
    
    my_res <- read.table(res_external_files[i])
    my_res$Cluster <- "NS"
    my_res$Cluster[my_res$log2FoldChange > 0 & my_res$padj < padj_cutoff] <- "up-reg"
    my_res$Cluster[my_res$log2FoldChange < 0 & my_res$padj < padj_cutoff] <- "down-reg"
    
    assign(res_external_names[i], my_res)
    
    rm(list = "my_res")
}




##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

#################################################       Signal Dotplot       #####################################################



for(log_type in log_types){
    
    
    for(my_assay in my_assays){
        
        
        if(my_assay == "all"){next()}
        
        #################################################           
        
        lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
        
        my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
        
        #################################################           
        
        if(my_assay == "S5P"){
            dot_colors <- my_color_palette  
        } else {
            dot_colors <- my_color_palette[-3]
        }
        
        #################################################           
        
        
        lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp[, grep("Zy", colnames(lnc_tmp))], na.rm = TRUE)
        
        
        for(my_classification in c("my_DBTMEE", res_external_names)){
            
            
            #################################################           
            
            
            lnc_scaled_tmp_merged <- merge(lnc_scaled_tmp, get(my_classification), by.x = "row.names", by.y = "gene_id", all = FALSE)
            lnc_scaled_tmp_merged$Cluster <- factor(lnc_scaled_tmp_merged$Cluster)
            
            
            lnc_scaled_tmp_merged_mn <- t(sapply(levels(lnc_scaled_tmp_merged$Cluster),function(i){
                
                colMedians(as.matrix(lnc_scaled_tmp_merged[lnc_scaled_tmp_merged$Cluster == i, 
                                                           grep(my_assay, colnames(lnc_scaled_tmp_merged))]), na.rm = TRUE)
            }))
            
            #################################################           
            
            colnames(lnc_scaled_tmp_merged_mn) <- colnames(lnc_scaled_tmp_merged[, grep(my_assay, colnames(lnc_scaled_tmp_merged))])
            
            stopifnot(identical(colnames(lnc_scaled_tmp_merged_mn), colnames(lnc_scaled_tmp)))
            
            lnc_scaled_tmp_merged_mn <- rbind(colMedians(as.matrix(lnc_scaled_tmp), na.rm = TRUE), lnc_scaled_tmp_merged_mn)
            
            
            rownames(lnc_scaled_tmp_merged_mn)[1] <- "all genes"
            
            
            #################################################           
            
            
            lnc_scaled_tmp_merged_mn_long <- melt(lnc_scaled_tmp_merged_mn, varnames = c("Classes", "Samples"))
            
            lnc_scaled_tmp_merged_mn_long$Conditions <- gsub(paste0("_",my_assay,"|_[1-3]"),"", lnc_scaled_tmp_merged_mn_long$Samples)
            lnc_scaled_tmp_merged_mn_long$Conditions <- factor(lnc_scaled_tmp_merged_mn_long$Conditions, levels = unique(lnc_scaled_tmp_merged_mn_long$Conditions))
            
            stopifnot(identical(levels(lnc_scaled_tmp_merged_mn_long$Conditions), levels(my_conditions)))
            
            #################################################           
            
            my_class <- "up-reg"
            
            
            pdf(paste0(output_plots,"/dotplot.",gsub("my_","",my_classification),".",my_assay,".", log_type,".pdf"), height = 9, width = 9, useDingbats = F)
            
            
            for(my_class in levels(lnc_scaled_tmp_merged_mn_long$Classes)){
                
                
                par(mfrow=c(2,2), mar = c(5,10,3,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))
                
                
                #################################################           
                
                
                df_tmp <- lnc_scaled_tmp_merged_mn_long[lnc_scaled_tmp_merged_mn_long$Classes == my_class,]
                
                set.seed(123)
                
                plot(x = jitter(as.integer(df_tmp$Conditions),factor = 0.1),
                     y = df_tmp$value,
                     main =  paste0(feature_type," ",my_assay, " at \n",
                                    gsub("my_|res.","", ifelse(my_class == "all genes","",my_classification))," ",my_class),
                     xlim = c(min(as.integer(df_tmp$Conditions))-0.5, max(as.integer(df_tmp$Conditions))+0.5),
                     ylim = c(min(df_tmp$value)-0.5, max(df_tmp$value)+0.5),
                     xlab = "", xaxt = "n",
                     ylab = "Median log2 Normalized Counts \n(relative to Zy)",
                     pch = 19, col = dot_colors[df_tmp$Conditions])
                
                set.seed(123)
                
                points(x = jitter(as.integer(df_tmp$Conditions), factor = 0.1),
                       y =  df_tmp$value,
                       col = "#555555", pch=1, lwd=0.5)
                
                axis(side = 1, at = seq_along(levels(df_tmp$Conditions)), 
                     labels = levels(df_tmp$Conditions), las = 2)
                
                
                #################################################           
                
                par(mar = c(5,2,3,10))
                
                plot.new()
                
                legend("left",
                       legend = levels(df_tmp$Conditions),
                       cex = 1, bty = "n",
                       fill =  dot_colors[seq_along(levels(df_tmp$Conditions))])
                
                
                #################################################           
                
                
                
                fit <- lm(value ~ Conditions, df_tmp)
                
                set.seed(2)
                
                if(my_assay == "S5P"){
                    
                    glht_out <- summary(glht(fit, linfct = mcp(Conditions = c("`2C` - `Zy` = 0",
                                                                              "`2C_DRB` - `2C` = 0",
                                                                              "`8C` - `Zy` = 0",
                                                                              "`ES` - `Zy` = 0"))))
                } else {
                    
                    glht_out <- summary(glht(fit, linfct = mcp(Conditions = c("`2C` - `Zy` = 0",
                                                                              "`8C` - `Zy` = 0",
                                                                              "`ES` - `Zy` = 0"))))
                    
                }
                
                stats_tmp <- data.frame(`Estimate` =   round(glht_out$test$coefficients*100, 2)/100,
                                        `Std. Error` = round(glht_out$test$sigma*100, 2)/100,
                                        `t value` =  round(glht_out$test$tstat*100, 2)/100,
                                        `Pr(>|t|)` = round(glht_out$test$pvalues*100, 2)/100, 
                                        check.names = FALSE)
                
                
                textplot(stats_tmp, halign = "left", valign = "top", cex = 0.7)
                text(x = 0.5, y = 0.75, cex = 0.8,
                     labels = "Multiple Comparisons of Means: User-defined Contrasts")
                
                
                
                rm(list = "df_tmp")
                rm(list = "stats_tmp")
                
                
                #################################################           
            }
            
            dev.off()
            
            
            ################################################# 
            
            if(my_classification == "my_DBTMEE"){
                my_subclasses_list <- list(c("Maternal_RNA","1-Cell_Transient", "2-Cell_Transient", "4-Cell_Transient"),
                                           c("Maternal_RNA", "Minor_ZGA",  "Major_ZGA", "MGA"))         
            } else {
                my_subclasses_list <- list(c("down-reg", "NS",  "up-reg"))   
            }
            
            
            ################################################# 
            
            if(feature_type == "genebody"){
                
                if(my_assay == "S5P"){
                    my_xlims <- c(-0.5, 1.5)  
                } else {
                    my_xlims <- c(-0.5, 1.0) 
                }
            } else {
                if(my_assay == "S5P"){
                    my_xlims <- c(-0.5, 2.0)  
                } else {
                    my_xlims <- c(-0.5, 1.0) 
                }
                
                
            }
            
            ################################################# 
            
            for(i in seq_along(my_subclasses_list)){
                
                
                pdf(paste0(output_plots,"/dotplot.",gsub("my_","",my_classification),".",i,".",my_assay,".", log_type,".pdf"), height = 5, width = 12, useDingbats = F)
                par(mfrow=c(2,3), mar = c(3,8,3,2), oma = c(3,3,3,3), mgp = c(3,1,0))
                
                
                for(my_condition in levels(lnc_scaled_tmp_merged_mn_long$Conditions)){
                    
                    
                    df_tmp <- lnc_scaled_tmp_merged_mn_long[lnc_scaled_tmp_merged_mn_long$Conditions %in% my_condition,]
                    
                    df_tmp_sub <- df_tmp[df_tmp$Classes %in% c("all genes", my_subclasses_list[[i]]),]
                    df_tmp_sub$Classes <- factor(df_tmp_sub$Classes, levels = rev(c("all genes", my_subclasses_list[[i]])))
                    
                    set.seed(123)
                    
                    plot(y = jitter(as.integer(df_tmp_sub$Classes), factor = 0.1),
                         x = df_tmp_sub$value,
                         main =  my_condition,
                         ylim = c(min(as.integer(df_tmp_sub$Classes))-0.5, max(as.integer(df_tmp_sub$Classes))+0.5),
                         xlim = my_xlims,
                         ylab = "", xlab = "", yaxt = "n",
                         pch = 19, col = dot_colors[which(my_condition == levels(df_tmp_sub$Conditions))])
                    
                    set.seed(123)
                    
                    points(y = jitter(as.integer(df_tmp_sub$Classes), factor = 0.1),
                           x =  df_tmp_sub$value,
                           col = "#555555", pch=1, lwd=0.5)
                    
                    abline(v = mean(df_tmp_sub$value[df_tmp_sub$Classes == "all genes"]), lty=3)
                    
                    axis(side = 2, at = seq_along(levels(df_tmp_sub$Classes)), labels = levels(df_tmp_sub$Classes), las=2)
                    
                    rm(list = ls(pattern = "df_tmp"))
                    
                    
                }
                
                
                mtext(text = paste0(feature_type," - ",my_assay), side = 3, line = 0, font = 2,  outer = TRUE, cex = 1.2)
                mtext(text = "Median log2 Normalized Counts (relative to Zy)", side = 1, line = 0.5, outer = TRUE)
                mtext(text = paste0(gsub("my_|res.","",my_classification), " Classes"), 
                      side = 2, line = 0.5, outer = TRUE)
                
                
                
                dev.off()
                
            }
            
            rm(list = "lnc_scaled_tmp_merged")
        }
        
        rm(list = ls(pattern = "_tmp"))
    }
    
}



##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

#################################################     compare toRNA-seq      #####################################################




log2_TPM_DBTMEE <- merge(log2_zscore_TPM, my_DBTMEE, by.x = "row.names", by.y = "gene_id", all = FALSE)

log2_TPM_DBTMEE_long <- melt(log2_TPM_DBTMEE)


ggp3 <- ggplot(log2_TPM_DBTMEE_long, aes(x=variable, y=value)) + 
    theme_bw() +
    theme(text = element_text(size = 12), aspect.ratio = 0.75, 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_violin(width = 1, fill = "grey") +
    facet_wrap(~ Cluster, nrow = 2) +
    ggtitle(label = "RNA-seq in DBTMEE classes") +
    stat_summary(fun= "median", geom="point", size=1) +
    scale_fill_manual(values = vio_colors) +
    coord_cartesian(ylim = c(-4, 4)) +
    ylab(label = "log2 TPM (z-score)") +
    xlab(label = "Stages") 

ggsave(filename = paste0(output_plots,"/vioplot.RNAseq_DBTMEE.pdf"),
       plot =  ggp3,
       width = 10, height = 6)

#################################################        

##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################




for(log_type in log_types){
    
    for(my_assay in my_assays){
        
        if(my_assay == "all"){next()}
        
        ################################################# 
        
        if(feature_type == "genebody"){
            
            if(my_assay == "S5P"){
                my_xlims <- c(-0.5, 1.5)  
            } else {
                my_xlims <- c(-0.5, 1.0) 
            }
        } else {
            if(my_assay == "S5P"){
                my_xlims <- c(-0.5, 2.0)  
            } else {
                my_xlims <- c(-0.5, 1.0) 
            }
            
            
        }
        
        ################################################# 
        
        
        
        lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
        
        
        my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
        
        #################################################           
        
        if(my_assay == "S5P"){
            dot_colors <- my_color_palette  
        } else {
            dot_colors <- my_color_palette[-3]
        }
        
        #################################################           
        
        
        
        lnc_scaled_tmp <-  lnc_tmp - rowMeans(lnc_tmp[, grep("Zy", colnames(lnc_tmp))], na.rm = TRUE)
        
        
        ################################################# 
        
        
        
        
        
        pdf(paste0(output_plots,"/dotplot_RNAseq.",my_assay,".", log_type,".pdf"), height = 5.5, width = 12, useDingbats = F)
        
        
        for(qs in c("Q5","Q35")){
            
            
            log2_TPM_Qs <- get(paste0("log2_TPM_", qs))
            
            par(mfrow=c(2,3), mar = c(3,8,3,2), oma = c(3,3,3,3), mgp = c(3,1,0))
            
            
            for(my_condition in levels(my_conditions)){
                
                
                df_tmp <- sapply(names(log2_TPM_Qs), function(x){
                    
                    colMedians(as.matrix(lnc_scaled_tmp[rownames(lnc_scaled_tmp) %in% log2_TPM_Qs[[x]], my_condition == my_conditions]), na.rm = TRUE)
                })
                
                
                df_tmp_long <- melt(df_tmp)
                df_tmp_long$Var2 <- factor(df_tmp_long$Var2, levels = rev(levels(df_tmp_long$Var2)))
                
                set.seed(123)
                
                plot(y = jitter(as.integer(df_tmp_long$Var2), factor = 0.1),
                     x = df_tmp_long$value,
                     main =  my_condition,
                     ylim = c(min(as.integer(df_tmp_long$Var2))-0.5, max(as.integer(df_tmp_long$Var2))+0.5),
                     xlim = my_xlims,
                     ylab = "", xlab = "", yaxt = "n",
                     pch = 19, col = dot_colors[which(my_condition == levels(my_conditions))])
                
                
                set.seed(123)
                
                points(y = jitter(as.integer(df_tmp_long$Var2), factor = 0.1),
                       x =  df_tmp_long$value,
                       col = "#555555", pch=1, lwd=0.5)
                
                abline(v = mean(df_tmp_long$value[df_tmp_long$Var2 == paste0(qs,"_Zygote")]), lty=3)
                
                axis(side = 2, at = seq_along(levels(df_tmp_long$Var2)), labels = levels(df_tmp_long$Var2), las=2)
                
                rm(list = ls(pattern = "df_tmp"))
                
            }
            
            rm(list = "log2_TPM_Qs")
            
            mtext(text = paste0(feature_type," - ",my_assay), side = 3, line = 0, font = 2,  outer = TRUE, cex = 1.2)
            mtext(text = "Median log2 Normalized Counts (relative to Zy)", side = 1, line = 0.5, outer = TRUE)
            mtext(text = paste0("RNA-seq ",qs," genes"), side = 2, line = 0.5, outer = TRUE)
            
        }
        
        
        
        dev.off()
        
    }
    
    
    ################################################# 
    
    rm(list = ls(pattern = "_tmp"))
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


