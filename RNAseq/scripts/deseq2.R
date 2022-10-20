






##################################################################################################################################
##################################################################################################################################




library(DESeq2)
library(sva)


library(ShortRead)
library(rtracklayer)
library(GenomicFeatures)

library(HelpersforDESeq2)

library(org.Mm.eg.db)
library(topGO)

library(RColorBrewer)
library(pheatmap)


library(reshape2)
library(ggplot2)



##################################################################################################################################
##################################################################################################################################





args = commandArgs(trailingOnly=TRUE)


input_anno <- grep("txdb", args, value = TRUE)
input_table <- grep("SampleTable", args, value = TRUE)

output_tables <- gsub("\\/log2.*","", grep("log2_", args, value = TRUE))
output_plots <-  gsub("\\/PCA.*","", grep("PCA", args, value = TRUE))
output_GO    <-  gsub("plots/PCA.*","GO",  grep("PCA", args, value = TRUE))


output_sessionInfo <- grep("sessionInfo", args, value = TRUE)



##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################



txdb <- loadDb(input_anno)


my_chromosomes <- seqlevels(txdb)[1:21]
my_chromosomes <- gsub("chr","",my_chromosomes)


my_genes <- GenomicFeatures::genes(txdb)

my_genes$gene_symbol <- mapIds(org.Mm.eg.db, keys = my_genes$gene_id, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

my_genes$gene_symbol[grep("^ERCC", my_genes$gene_id)] <- my_genes$gene_id[grep("^ERCC", my_genes$gene_id)]
my_genes$gene_symbol[grep("^pGEMHE_mCherry_mTrim21", my_genes$gene_id)] <- "mCherry_mTrim21"


seqlevelsStyle(my_genes) <- "Ensembl"


my_genes_ids_MT <- my_genes[seqnames(my_genes) == "MT"]$gene_id




##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        SampleTable         #####################################################




SampleTable <- read.delim(input_table, header = T, stringsAsFactors = F)
SampleTable <- SampleTable[, 1:2, drop = FALSE]

SampleTable <- SampleTable[!duplicated(SampleTable$SampleID), , drop = FALSE]
SampleTable <- SampleTable[order(SampleTable$SampleID), , drop = FALSE]

SampleTable$Conditions <- gsub("_.*","",SampleTable$ForeignID)




##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        count table         #####################################################



my_count_files <-  grep("ReadsPerGene.out.tab", args, value = TRUE)
my_count_files <- my_count_files[order(my_count_files)]


if(identical(SampleTable$SampleID, gsub(".*\\/|.ReadsPerGene.out.tab|\\..*","",my_count_files))){
        
        my_count_table <- makeCountTable(count_names = SampleTable$SampleID,
                                         count_files = my_count_files, 
                                         stranded = FALSE)
}

my_count_table <- my_count_table[-1:-4,]




##################################################################################################################################
##################################################################################################################################


























##################################################################################################################################
##################################################################################################################################

#################################################      Quality  Control      ##################################################### 




######### cutoffs #########

min_genic_reads = 5e5
max_percent_ERCC = 10
max_percent_mito = 10
min_detected_genes = 1e3
min_CPM_mCherry_mTrim21 = 500


################################################# 


n_genic_reads <- colSums(my_count_table[grep("ENSMUS",rownames(my_count_table)),])
n_ERCC_reads <- colSums(my_count_table[grep("ERCC",rownames(my_count_table)),]) 


percent_mito <- colSums(my_count_table[rownames(my_count_table) %in% my_genes_ids_MT,]) / colSums(my_count_table) * 100
percent_ERCC <- colSums(my_count_table[grep("ERCC",rownames(my_count_table)),]) / colSums(my_count_table) * 100

n_detected_genes <- apply(my_count_table, 2, function(x){   sum(x > 0)   })

CPM_mCherry_mTrim21 <- my_count_table[rownames(my_count_table) %in% c("pGEMHE_mCherry_mTrim21"),] / colSums(my_count_table) * 1e6 
CPM_endo_mTrim21 <- my_count_table[rownames(my_count_table) %in% c("ENSMUSG00000030966"),] / colSums(my_count_table) * 1e6 



stopifnot(identical(colnames(my_count_table), SampleTable$SampleID))

my_outliers <- (n_genic_reads < min_genic_reads) | 
        (percent_ERCC > max_percent_ERCC) |
        (percent_mito > max_percent_mito) |
        (n_detected_genes < min_detected_genes) | 
        (!(SampleTable$Conditions %in% c("Noinject", "Cdk9in")) & CPM_mCherry_mTrim21 < min_CPM_mCherry_mTrim21) |
        (  SampleTable$Conditions %in% c("Noinject", "Cdk9in")  & CPM_mCherry_mTrim21 > min_CPM_mCherry_mTrim21)



##################################################################################################################################
##################################################################################################################################




my_color_palette <- c("#5026D9", "#8D9093", "#C199CE", "#DA0EA1", "#90E7F8", "#1971A9")

if(all(identical(colnames(my_count_table), SampleTable$SampleID),
       identical(colnames(my_count_table), names(n_genic_reads)))){
        
        my_conditions <- factor(SampleTable$Conditions)
} else {
        rm(list = "my_conditions")
}




################################################# 


pdf(paste0(output_plots,"/barplot_QC.pdf"), height = 7, width = 12, useDingbats = F)

par(mfrow=c(2,2), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))


bp <- barplot(n_genic_reads / 1e6, 
              ylim = c(0,8), ylab = "Number of Genic Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_count_table), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


abline(h = min_genic_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 


bp <- barplot(n_detected_genes,
              ylim = c(0,2.5e4), ylab = "Number of Detected Genes (Counts > 0)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_count_table), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)

abline(h = min_detected_genes, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 


bp <- barplot(n_ERCC_reads / 1e6, 
              ylim = c(0,2.5), ylab = "Number of ERCC Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_count_table), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)


################################################# 



bp <- barplot(percent_ERCC, 
              ylim = c(0,80), ylab = "Percentage of ERCC", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_count_table), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


abline(h = max_percent_ERCC, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)


################################################# 


plotLegend(conditions = my_conditions,
           legend_colors = my_color_palette, 
           legend_size = 0.8)

################################################# 


par(mfrow=c(2,2), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))


bp <- barplot(percent_mito,
              ylim = c(0,50), ylab = "Percentage of Mitochondrial Reads", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_count_table), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)

abline(h = max_percent_mito, lty=2)


my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)





################################################# 


bp <- barplot(CPM_mCherry_mTrim21, 
              ylim = c(0,8000), ylab = "mCherry-mTrim21 (CPM)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_count_table), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)



my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)

abline(h = min_CPM_mCherry_mTrim21, lty=2)

################################################# 


bp <- barplot(CPM_endo_mTrim21, 
              ylim = c(0,1000), ylab = "endogenous locus Trim21 (CPM)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_count_table), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)



my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)


################################################# 



plotLegend(conditions = my_conditions,
           legend_colors = my_color_palette, 
           legend_size = 0.85)



################################################# 


dev.off()


##################################################################################################################################
##################################################################################################################################










pdf(paste0(output_plots,"/dotplot_QC.pdf"), height = 7, width = 9, useDingbats = F)

if(all(identical(SampleTable$SampleID, names(n_genic_reads)),
       identical(SampleTable$SampleID, names(n_detected_genes)),
       identical(SampleTable$SampleID, names(n_ERCC_reads)),
       identical(SampleTable$SampleID, names(n_ERCC_reads)),
       identical(SampleTable$SampleID, names(percent_ERCC)),
       identical(SampleTable$SampleID, names(percent_mito)),
       identical(SampleTable$SampleID, names(CPM_mCherry_mTrim21)))){
        
        qc_list <- list(`Number of Genic Reads (Million)` = (n_genic_reads / 1e6),
                        `Number of Detected Genes (Counts > 0)` = (n_detected_genes),
                        `Number of ERCC Reads (Million)` = ((n_ERCC_reads / 1e6)),
                        `Percentage of ERCC` = percent_ERCC,
                        `Percentage of Mitochondrial Reads` = percent_mito,
                        `mCherry-mTrim21 (CPM)` = CPM_mCherry_mTrim21)
} else {
        rm(list = "qc_list")
}


qc_thresholds <- c(min_genic_reads/1e6, 
                   min_detected_genes,
                   NA,
                   max_percent_ERCC,
                   max_percent_mito,
                   min_CPM_mCherry_mTrim21)

par(mfcol=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))

for(i in seq_along(qc_list)){
        
        set.seed(123)
        
        plot(x = jitter(as.integer(my_conditions),factor = 0.3),
             y = qc_list[[i]],
             main =  "",
             xlim = c(min(as.integer(my_conditions))-0.5, max(as.integer(my_conditions))+0.5),
             ylim = c(0, ceiling(max(qc_list[[i]]))*1.2),
             xlab = "",
             ylab = names(qc_list)[i], xaxt = "n",
             pch = 19, col = my_color_palette[my_conditions])
        
        set.seed(123)
        
        points(x = jitter(as.integer(my_conditions), factor = 0.3),
               y =  qc_list[[i]],
               col = "#555555", pch=1, lwd=0.5)
        
        abline(h =  qc_thresholds[i], lty=2)
        
        axis(side = 1, at = seq_along(levels(my_conditions)), 
             labels = levels(my_conditions), las = 2)
        
}


plotLegend(conditions = my_conditions,
           legend_colors = my_color_palette, 
           legend_size = 0.85)

dev.off()




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

#################################################           Setup            ##################################################### 




SampleTable <- SampleTable[!(my_outliers),]

my_count_table <- my_count_table[, colnames(my_count_table) %in% SampleTable$SampleID]


SampleTable$Experiments <- factor(gsub(".*_", "",gsub("_[1-9]","",SampleTable$ForeignID)))


if(identical(colnames(my_count_table), SampleTable$Sample)){
        
        dds <- setupDDS(CountTableName = "my_count_table", 
                        SampleTableName = "SampleTable",
                        SampleIdName = "SampleID", 
                        ConditionName = "Conditions", 
                        BatchName = "Experiments",
                        n_samples_for_filtering = ceiling(ncol(my_count_table)*0.5),
                        min_number_of_reads = 1)
}





##################################################################################################################################
##################################################################################################################################

#################################################         Contrasts          ##################################################### 




contrast_list <- list(c("IgG",     "Noinject"),
                      c("Cdk9in",  "Noinject"),
                      c("Spt5KD",  "IgG"),
                      c("NelfEKD", "IgG"))


for(contrast_name in contrast_list){
        
        res <- getResults(dds = dds,
                          contrast = contrast_name,
                          lfc_cutoff = 0,
                          shrink = F,
                          annotation = "my_genes",
                          anno_symbol = "gene_symbol",
                          anno_id = "gene_id") 
        
        res$gene_id <- rownames(res)
        res$DBTMEE <- my_DBTMEE$Cluster[match(res$gene_id, my_DBTMEE$gene_id)]
        
        res_name <- paste0("res.", contrast_name[1],"-",contrast_name[2])
        assign(res_name, res)
        
        write.table(res, file = paste0(output_tables,"/",res_name, ".txt"), 
                    quote = F, sep = "\t", row.names = T, col.names = NA) 
        
        rm(list = "res")
        
}







##################################################################################################################################
##################################################################################################################################




res_names <- ls(pattern = "^res\\.")

padj_cutoff <- 0.05



##################################################################################################################################
##################################################################################################################################




res_internal_list <- list()


for(i in seq_along(res_names)){
        
        my_res <- get(res_names[i])
        my_res$Cluster <- "NS"
        my_res$Cluster[my_res$log2FoldChange > 0 & my_res$padj < padj_cutoff] <- "up-reg"
        my_res$Cluster[my_res$log2FoldChange < 0 & my_res$padj < padj_cutoff] <- "down-reg"
        
        
        res_internal_list[[paste0(gsub("res.","",res_names[i]),".", "down-reg")]] <- my_res$gene_id[my_res$Cluster == "down-reg"]
        res_internal_list[[paste0(gsub("res.","",res_names[i]),".", "up-reg")]] <-   my_res$gene_id[my_res$Cluster == "up-reg"]
        
        rm(list = "my_res")
}



##################################################################################################################################
##################################################################################################################################
























##################################################################################################################################
##################################################################################################################################

#################################################             RNA-seq            ##################################################### 



TPM <- read.table("../../Public/GSE38495_GSE45719/Output/TPM/TPM_means.txt", check.names = FALSE)
TPM <- TPM[,c("MII_oocyte","Zygote", "Early_2C", "Mid_2C", "Late_2C", "4C", "8C", "ICM","ESC")]


log2_TPM <- log2(TPM+1)

log2_scaled_TPM <- log2_TPM - rowMeans(log2_TPM)

log2_zscore_TPM <- t(scale(t(log2_TPM)))



##################################################################################################################################
##################################################################################################################################




log2_zscore_TPM_long <- melt(log2_zscore_TPM)


for(res_name in grep("Cdk9|Spt5", res_names, value = TRUE)){
        
        res <- as.data.frame(get(res_name))
        
        log2_zscore_TPM_long$Class <- "NS"
        log2_zscore_TPM_long$Class[log2_zscore_TPM_long$Var1 %in% res$gene_id[res$log2FoldChange > 0 & res$padj < padj_cutoff]] <- "up"
        log2_zscore_TPM_long$Class[log2_zscore_TPM_long$Var1 %in% res$gene_id[res$log2FoldChange < 0 & res$padj < padj_cutoff]] <- "down"
        
        
        ggp <- ggplot(log2_zscore_TPM_long, aes(x=Var2, y=value)) + 
                theme_bw() +
                theme(text = element_text(size = 12), aspect.ratio = 0.75, 
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      strip.background = element_blank(),
                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                geom_violin(width = 1, fill = "grey") +
                facet_wrap(~ Class, nrow = 1) +
                ggtitle(label = paste0("Differential in ", gsub("res.","",res_name))) +
                stat_summary(fun= median, geom="point", size=1) +
                #scale_fill_manual(values = vio_colors) +
                coord_cartesian(ylim = c(-4, 4)) +
                ylab(label = "log2 TPM (z-score)") +
                xlab(label = "Stages") 
        
        ggsave(filename = paste0(output_plots,"/vioplot.Deng.",gsub("res.","",res_name),".pdf"),
               plot =  ggp,
               width = 10, height = 6)
}


##################################################################################################################################
##################################################################################################################################






##################################################################################################################################
##################################################################################################################################

#################################################     RNA-seq Quantiles      #####################################################



log2_TPM_Quantiles <- log2_TPM[rownames(log2_TPM) %in% rownames(dds),]

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



# my_class <- "Q5_Zygote"
# 
# 
# for(qs in c("Q5","Q35")){
#         
#         log2_TPM_Qs <- get(paste0("log2_TPM_", qs))
#         
#         for(my_class in names(log2_TPM_Qs)){
#                 
#                 log2_zscore_TPM_long_sub <- log2_zscore_TPM_long[log2_zscore_TPM_long$Var1 %in% log2_TPM_Qs[[my_class]],]
#                 
#                 ggp <- ggplot(log2_zscore_TPM_long_sub, aes(x=Var2, y=value)) + 
#                         theme_bw() +
#                         theme(text = element_text(size = 12), aspect.ratio = 0.75, 
#                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#                               panel.border = element_rect(colour = "black", fill=NA, size=1),
#                               panel.grid.major = element_blank(),
#                               panel.grid.minor = element_blank(),
#                               strip.background = element_blank(),
#                               plot.title = element_text(hjust = 0.5, face = "bold")) +
#                         geom_violin(width = 1, fill = "grey") +
#                         ggtitle(label = paste0(my_class," - ", length(unique(log2_zscore_TPM_long_sub$Var1)))) +
#                         stat_summary(fun= median, geom="point", size=1) +
#                         #scale_fill_manual(values = vio_colors) +
#                         coord_cartesian(ylim = c(-4, 4)) +
#                         ylab(label = "log2 TPM (z-score)") +
#                         xlab(label = "Stages") 
#                 
#                 rm(list = "log2_zscore_TPM_long_sub")
#                 
#                 ggsave(filename = paste0(output_plots,"/vioplot.Deng.",my_class,".pdf"),
#                        plot =  ggp,
#                        width = 6, height = 4)  
#                 
#         }
#         
#         rm(list = "log2_TPM_Qs")
# }


##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################

##############################################    RNA-seq Cdk9in Spt5 results   ###################################################



res_external_files <- c("../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.4C-Late_2C.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.8C-Late_2C.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.ICM-Late_2C.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.Early_2C-Zygote.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.Mid_2C-Zygote.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.Late_2C-Zygote.txt")

res_external_names <- gsub(".*\\/|.txt", "", res_external_files)


res_external_list <- list()


for(i in seq_along(res_external_files)){
        
        my_res <- read.table(res_external_files[i])
        my_res$Cluster <- "NS"
        my_res$Cluster[my_res$log2FoldChange > 0 & my_res$padj < padj_cutoff] <- "up-reg"
        my_res$Cluster[my_res$log2FoldChange < 0 & my_res$padj < padj_cutoff] <- "down-reg"
        
        
        res_external_list[[paste0(gsub("res.","",res_external_names[i]),".", "down-reg")]] <- my_res$gene_id[my_res$Cluster == "down-reg"]
        res_external_list[[paste0(gsub("res.","",res_external_names[i]),".", "up-reg")]] <- my_res$gene_id[my_res$Cluster == "up-reg"]
        
        assign(res_external_names[i], my_res)
        
        rm(list = "my_res")
}




##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################

#################################################         Venn Classes       ##################################################### 


library(Vennerable)




i=1
j=2

pdf(paste0(output_plots,"/Venn_1.pdf"), height = 15, width = 15, useDingbats = F)

for(i in  grep("Cdk9|Spt5|Nelf", names(res_internal_list))){
        
                
                for(j in grep("Cdk9|Spt5|Nelf", names(res_internal_list))){
                        
                        
                        if(i == j){next()}
                        if(identical(gsub("\\..*","", names(res_internal_list)[i]),gsub("\\..*","", names(res_internal_list)[j]))){next()}
                        
                        plot.new()
                        my_venn <- Venn(c(res_internal_list[i],
                                          res_internal_list[j]))
                        
                        my_venn <- tryCatch({
                                Vennerable::compute.Venn(my_venn,doWeights=TRUE)
                        }, error = function(error_condition) {
                                Vennerable::compute.Venn(my_venn,doWeights=FALSE)
                        })
                        
                        gp <-  VennThemes(my_venn)
                        
                        gp$Set$Set1$lwd = 10
                        gp$Set$Set2$lwd = 10


                        plot(my_venn,  show = list(Faces = FALSE, Universe=FALSE), gp=gp)   
                        
                }
        }

dev.off()




##################################################################################################################################
##################################################################################################################################


i=1
j=1
k=1

pdf(paste0(output_plots,"/Venn_2.pdf"), height = 15, width = 15, useDingbats = F)

for(k in grep("Cdk9|Spt5", names(res_internal_list))){
        
        for(i in seq_along(my_DBTMEE_list)){
                
                for(j in seq_along(log2_TPM_Q5)){
                        
                        
                        plot.new()
                        my_venn <- Venn(c(my_DBTMEE_list[i],
                                          log2_TPM_Q5[j],
                                          res_internal_list[k]))
                        
                        my_venn <- tryCatch({
                                Vennerable::compute.Venn(my_venn,doWeights=TRUE)
                        }, error = function(error_condition) {
                                Vennerable::compute.Venn(my_venn,doWeights=FALSE)
                        })
                        
                        gp <-  VennThemes(my_venn)
                        
                        gp$Set$Set1$lwd = 10
                        gp$Set$Set2$lwd = 10
                        gp$Set$Set3$lwd = 10
                        
                        # gp$SetText$Set1$fontsize <- 30
                        # gp$SetText$Set2$fontsize <- 30
                        # 
                        # gp$FaceText$`11`$fontsize <- 30
                        # gp$FaceText$`10`$fontsize <- 30
                        # gp$FaceText$`01`$fontsize <- 30
                        
                        # my_venn@SetLabels$hjust <- c("right", "left")
                        
                        plot(my_venn,  show = list(Faces = FALSE, Universe=FALSE), gp=gp)   
                        
                }
        }
}

dev.off()


##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################          MA-plots          ##################################################### 



my_classes_full <- c(NoLabel = "",
                     ERCC = list(grep("ERCC", rownames(dds), value = TRUE)),
                     mCherry_mTrim21 = "pGEMHE_mCherry_mTrim21",
                     log2_TPM_Q5,
                     log2_TPM_Q35,
                     my_DBTMEE_list,
                     res_external_list)



for(my_label in names(my_classes_full)){
        
        
        pdf(paste0(output_plots,"/MA_",my_label,".pdf"), height = 6, width = 6, useDingbats = F)
        
        for(res_name in res_names){
                
                par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(4,4,4,4), mgp = c(2.5,1,0))
                
                plottingMA(res = get(res_name),
                           main_title = gsub("res.","",res_name),
                           main_title_size = 1.5,
                           selection_ids = my_classes_full[[my_label]],
                           selection_id_type = "gene_id",
                           selection_name = my_label,
                           selection_point_size = 0.5,
                           selection_text_label = (my_label == "mCherry_mTrim21"),
                           selection_text_size = 0.75,
                           selection_text_adj = -0.5,
                           selection_shadow = FALSE,
                           xlims = c(0, 6),
                           ylims = c(-15, 15),
                           x_axis_by = 1,
                           padj_cutoff = padj_cutoff,
                           show_legend = TRUE)
                
        }
        
        dev.off()
}


##################################################################################################################################
##################################################################################################################################














##################################################################################################################################
##################################################################################################################################

#################################################       log2FC plots        ##################################################### 





for(my_selection in names(my_classes_full)){
        
        
        if(length(my_selection) == 0){next()}
        if(length(res_names) < 2){next()}
        
        pdf(paste0(output_plots,"/log2FC_",my_selection,".pdf"), height = 6, width = 6, useDingbats = F)
        
        for(res_name1 in res_names){
                
                for(res_name2 in res_names){
                        
                        par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(4,4,4,4), mgp = c(2.5,1,0))
                        
                        if(res_name1 == res_name2){next()}
                        
                        plotLog2FC(res1 = get(res_name1), 
                                   res2 = get(res_name2), 
                                   x_label = paste0("log2FC [", gsub("res.","", res_name1),"]"), 
                                   y_label = paste0("log2FC [", gsub("res.","", res_name2),"]"), 
                                   lims = c(-15,15), 
                                   selection_text_label = (my_selection == "mCherry_mTrim21"),
                                   selection_ids = my_classes_full[[my_selection]],
                                   selection_id_type = "gene_id",
                                   selection_color = rgb(0.7,0,0.9,0.25),
                                   selection_point_size = 0.33,
                                   selection_legend = gsub("Labeled","", my_selection))         
                }
        }
        
        dev.off()
        
}




##################################################################################################################################
##################################################################################################################################



















##################################################################################################################################
##################################################################################################################################

#################################################             rld            ##################################################### 




rld <- rlog(dds, blind = FALSE)

rlog_norm_counts_raw <- assay(rld)

log2_norm_counts_raw <- log2(counts(dds, normalized = TRUE)+1)


##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################      Batch correction      ##################################################### 



if(is.null(colData(dds)$Batch)){
        
        rlog_norm_counts <- rlog_norm_counts_raw
        
        log2_norm_counts <- log2_norm_counts_raw
        
} else {
        
        batchVar <- colData(dds)$Batch
        
        modcombat <- model.matrix(~Sample, data = colData(dds))
        
        
        
        rlog_norm_counts <- ComBat(dat = rlog_norm_counts_raw,
                                   batch = batchVar, mod = modcombat,
                                   par.prior = TRUE, prior.plots = FALSE)
        
        
        log2_norm_counts <- ComBat(dat = log2_norm_counts_raw,
                                   batch = batchVar, mod = modcombat,
                                   par.prior = TRUE, prior.plots = FALSE)
}


write.table(log2_norm_counts, file = paste0(output_tables,"/log2_norm_counts.txt"), 
            quote = F, sep = "\t", row.names = T, col.names = NA) 



##################################################################################################################################
##################################################################################################################################


my_conditions <- factor(SampleTable$Conditions)


stopifnot(identical(SampleTable$Conditions, as.character(my_conditions)))
stopifnot(identical(SampleTable$SampleID, colnames(log2_norm_counts)))

annotation_col = data.frame(Conditions = my_conditions)
rownames(annotation_col) <- colnames(log2_norm_counts)

ann_colors = list(Conditions = my_color_palette[seq_along(levels(my_conditions))])
names(ann_colors$Conditions) <- levels(my_conditions)





##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################

#################################################        Correlation         ##################################################### 





callback = function(hc, mat){
        sv = svd(t(mat))$v[,1]
        dend = (reorder(as.dendrogram(hc), wts = sv))
        as.hclust(dend)
}



#################################################  


pdf(paste0(output_plots,"/correlation.pdf"), width = 12, height = 12, useDingbats = FALSE)
par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)


#################################################  

for(corr_type in c("Pearson","Spearman")){
        
        for(log2_type in c("log2", "rlog")){
                
                my_log2_data <- get(paste0(log2_type,"_norm_counts")) 
                my_log2_data <- my_log2_data - rowMeans(my_log2_data)
                
                
                my_corr <- cor(my_log2_data, method = tolower(corr_type))
                
                rownames(my_corr) <- SampleTable$ForeignID[match(rownames(my_corr), SampleTable$SampleID)]
                
                pheatmap(my_corr, 
                         main = paste0(corr_type,"´s R at genes n = ", nrow(log2_norm_counts),
                                       " on relative ", log2_type, 
                                       " counts"), 
                         annotation_col = annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100), 
                         breaks = seq(-1,1, length.out = 101),
                         cellheight = 12, cellwidth = 12,
                         clustering_callback = callback,
                         clustering_method = "ward.D2")
                
        }
        
}

#################################################  


dev.off()



##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################

#################################################            PCA             ##################################################### 





xycomps <- list(c(1,2),
                c(1,3),
                c(2,3))



pdf(paste0(output_plots,"/PCA.pdf"), height = 3.5, width = 9, useDingbats = F)


#################################################  

for(log2_type in c("log2", "rlog")){
        
        
        my_log2_data <- get(paste0(log2_type,"_norm_counts")) 
        
        
        par(mfrow=c(1,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))
        
        for(xycomp in xycomps){
                
                plottingPCA(my_log2_data,
                            xcomp = xycomp[1],
                            ycomp = xycomp[2],
                            conditions = my_conditions,
                            pca_colors = my_color_palette,
                            main_title = "PCA",
                            quantiles = c(0,1),
                            show_labels = FALSE,
                            point_size = 0.9,
                            my_xlimits = c(-200,200),
                            my_ylimits = c(-100,100))  
                
        }
        
        mtext(text = paste(log2_type, "counts"), side = 1, outer = TRUE, adj = 1, cex = 0.75)
        
        plotLegend(conditions = my_conditions,
                   legend_colors = my_color_palette, 
                   legend_size = 0.8)
        
}

#################################################  


dev.off()


##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

#################################################        Top Hits         ##################################################### 



n_hits <- 50




pdf(paste0(output_plots,"/Heatmaps.pdf"), width = 14, height = 11, useDingbats = FALSE)
par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)


#################################################  

for(res_name in c("ERCC", res_names)){
        
        for(log2_type in c("log2", "rlog")){
                
                
                my_log2_data <- get(paste0(log2_type,"_norm_counts")) 
                my_log2_data <- my_log2_data - rowMeans(my_log2_data)
                
                if(res_name == "ERCC"){
                        
                        my_hm_data <- my_log2_data[grep("^ERCC-", rownames(my_log2_data)),]
                        
                        my_title <- paste(res_name, nrow(my_hm_data))
                        
                } else {
                        
                        res <- get(res_name)
                        res <- res[order(res$pvalue),]
                        
                        my_fav_genes <- c(res$gene_id[res$log2FoldChange > 0][1:(n_hits/2)],
                                          res$gene_id[res$log2FoldChange < 0][1:(n_hits/2)])
                        
                        
                        my_hm_data <- my_log2_data[rownames(my_log2_data) %in% my_fav_genes,]
                        
                        my_title <- paste0("Top",(n_hits/2), " up, ","Top",(n_hits/2), " down"," in ", gsub("res.","",res_name))
                        
                        rm(list = "res")
                        rm(list = "my_fav_genes")
                }
                
                my_hm_data <- my_hm_data - rowMeans(my_hm_data)
                
                my_hm_rownames <- my_genes$gene_symbol[match(rownames(my_hm_data),my_genes$gene_id)]
                rownames(my_hm_data)[!is.na(my_hm_rownames)] <- my_hm_rownames[!is.na(my_hm_rownames)]
                
                
                pheatmap(my_hm_data, 
                         main = paste0(my_title,
                                       " (relative ", log2_type, " counts)"),
                         annotation_col = annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100), 
                         breaks = seq(-5,5, length.out = 101),
                         cellheight = 12, cellwidth = 12,
                         clustering_callback = callback,
                         clustering_method = "ward.D2")
                
                
        }
        
}

#################################################  



dev.off()






##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

#################################################            GO              #####################################################



system(paste("mkdir",output_GO))



for(res_name in res_names){
        
        res <- get(res_name)
        res <- res[!(is.na(res$gene_symbol)),]
        
        for(ont in c("BP","CC","MF")){
                
                
                for(reg in c("up","down")){
                        
                        
                        #################################################
                        
                        
                        if(reg == "up"){
                                
                                all_genes <- factor(as.integer(res$padj < padj_cutoff & res$log2FoldChange > 0))
                                names(all_genes) <- res$gene_symbol
                                
                        } else if(reg == "down"){
                                
                                all_genes <- factor(as.integer(res$padj < padj_cutoff & res$log2FoldChange < 0))
                                names(all_genes) <- res$gene_symbol
                        }
                        
                        
                        #################################################
                        
                        
                        if(sum(all_genes == 1) < 20){next()}
                        
                        
                        gt <- makeGOTable(all_genes = all_genes,
                                          shown_terms = 20,
                                          min_signficant = 20,
                                          select_ontology = ont,
                                          select_organism = "org.Mm.eg.db",
                                          select_ID = "SYMBOL")
                        
                        
                        #################################################
                        
                        
                        write.table(gt,
                                    file = paste0(output_GO,"/GO_",ont,".",gsub("res.","", res_name),".",reg,".Table.txt"),
                                    quote = F, sep = "\t", row.names = F)
                        
                        if(nrow(gt) < 2){next()}
                        
                        
                        #################################################
                        
                        
                        pdf(paste0(output_GO,"/GO_",ont,".",gsub("res.","", res_name),".",reg,".Bubbles.pdf"),
                            height = 12, width = 12, useDingbats = F)
                        par(mfrow=c(1,1), mar = c(2,2,2,2), oma = c(1,1,1,1), mgp = c(2.5,1,0), cex.main=2)
                        
                        plotGObBubbles(gt = gt, main_title = paste0(reg,"-regulated in ", gsub("res.","", res_name)))
                        
                        dev.off()
                        
                        
                        #################################################
                        
                        
                        pdf(paste0(output_GO,"/GO_",ont,".",gsub("res.","", res_name),".",reg,".Heatmaps.pdf"),
                            height = 9, width = 9, useDingbats = F)
                        par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(4,4,4,4), mgp = c(2.5,1,0))
                        
                        for(i in 1:nrow(gt)){
                                
                                for(log2_type in c("log2", "rlog")){
                                        
                                        my_log2_data <- get(paste0(log2_type,"_norm_counts"))
                                        my_log2_data <- my_log2_data - rowMeans(my_log2_data)
                                        
                                        res$GO <- res$gene_symbol %in%  unlist(strsplit(as.character(gt$gene_names[i]), split = "/"))
                                        
                                        my_fav_genes <- head(res$gene_id[res$GO], 100)
                                        
                                        my_log2_data_sub <- my_log2_data[rownames(my_log2_data) %in% my_fav_genes,]
                                        rownames(my_log2_data_sub) <- my_genes$gene_symbol[match(rownames(my_log2_data_sub), my_genes$gene_id)]
                                        
                                        stopifnot(identical(rownames(annotation_col), colnames(my_log2_data_sub)))
                                        
                                        pheatmap(my_log2_data_sub,
                                                 main = paste0("Top Sign. GO: ",gt$Term[i]," n=", nrow(my_log2_data_sub)," (",log2_type,")"),
                                                 annotation_col = annotation_col,
                                                 annotation_colors = ann_colors,
                                                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                                 breaks = seq(-5,5, length.out = 101),
                                                 cellheight = 5, cellwidth = 5, fontsize = 5,
                                                 cluster_cols = TRUE,
                                                 show_rownames = TRUE)
                                        
                                        res$GO <- NA
                                        
                                        rm(list = ls(pattern = "_sub$"))
                                        rm(list = "my_fav_genes")
                                        
                                }
                        }
                        
                        dev.off()
                        
                        #################################################
                }
                
        }
        
}



##################################################################################################################################
##################################################################################################################################







writeLines(capture.output(sessionInfo()), output_sessionInfo)






##################################################################################################################################
##################################################################################################################################


