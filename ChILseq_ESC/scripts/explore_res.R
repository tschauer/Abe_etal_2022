






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

feature_type <- gsub("Output/deseq2_|/plots_.*", "", output_plots)



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

#################################################           MA-plots         #####################################################




res_names <- ls(pattern = "^res\\.")

my_assays <- gsub("\\..*","",gsub("res.","",res_names))


my_classes_full <- c(NoLabel = "",
                     log2_TPM_Q5 #,
                     #log2_TPM_Q35,
                     #my_DBTMEE_list
)

padj_cutoff <- 0.05
NoLabel <- ""
my_label <- "NoLabel"
res_name <- res_names[7]


#for(my_label in names(my_classes_full)){
for(my_label in c("NoLabel","Q5_ESC")){
        
        
        for(my_assay in my_assays){
                
                pdf(paste0(output_plots,"/MA.",my_assay,".",my_label,".pdf"), height = 6, width = 6, useDingbats = F)
                
                par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(4,4,4,4), mgp = c(2.5,1,0))
                
                res_names_sub <- grep(my_assay, res_names, value = TRUE)
                
                for(res_name in res_names_sub){
                        
                        plottingMA(res = get(res_name),
                                   main_title = gsub("res|\\."," ",res_name),
                                   main_title_size = 1.5,
                                   selection_ids = my_classes_full[[my_label]],
                                   selection_id_type = "gene_id",
                                   selection_name = my_label,
                                   selection_point_size = 0.5,
                                   selection_text_label = FALSE,
                                   selection_shadow = FALSE,
                                   xlims = c( 0, 5),
                                   ylims = c(-6, 6),
                                   x_axis_by = 1,
                                   padj_cutoff = padj_cutoff,
                                   show_legend = TRUE)       
                }
                
                dev.off()
        }
}


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
        
        my_lnc <- my_lnc[, c(grep("100C", colnames(my_lnc)), grep("1000C|ChIP", colnames(my_lnc)))] 
        
        #my_lnc <- my_lnc[apply(my_lnc, 1, function(x){sum(is.na(x)) <= (length(x)/4)}),]
        
        # my_lnc <- my_lnc[apply(my_lnc, 1, function(x){
        #         all(aggregate(as.numeric(x), by = list(gsub("_[1-3]$","",colnames(my_lnc))), FUN=function(i){ sum(is.na(i)) <= 1})$x)
        # }),]
        
        
        if(my_setting %in% c("minmean","minmedian")){
                
                set.seed(456)
                #my_lnc <- impute.MinProb(my_lnc, tune.sigma = 0.1)
                my_lnc[is.na(my_lnc)] <- rnorm(n = sum(is.na(my_lnc)), mean = min(my_lnc, na.rm = TRUE), sd = 0.1)
                
                
                
        } else if(my_setting %in% c("filtmean","filtmedian")){
                
                my_lnc <- na.omit(my_lnc)
                
                # my_lnc <- my_lnc[apply(my_lnc, 1, function(x){
                #         all(aggregate(as.numeric(x), by = list(gsub("_[1-3]$","",colnames(my_lnc))), FUN=function(i){ sum(is.na(i)) <= 1})$x)
                # }),]
        }
        
        
        assign(my_lnc_name, my_lnc)
        
        rm(list = "my_lnc")
        rm(list = "my_lnc_name")
        
}







##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        Correlation         #####################################################




my_assays <- unique(gsub("\\..*", "", gsub("lnc.","",ls(pattern = "^lnc\\."))))
log_types <- unique(gsub(".*\\.", "", ls(pattern = "^lnc\\.")))


my_assay <- "S5P"
log_type <- "log2"



callback = function(hc, mat){
        sv = svd(t(mat))$v[,1]
        dend = (reorder(as.dendrogram(hc), wts = sv))
        as.hclust(dend)
}



for(my_assay in my_assays){
        
        for(log_type in log_types){
                
                
                pdf(paste0(output_plots,"/correlation.", my_assay,".", log_type,".pdf"), width = 12, height = 12, useDingbats = FALSE)
                par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                
                for(my_scaling in c("mean","len")){
                        
                        lnc_tmp <- get(paste("lnc", my_assay, log_type, sep = "."))
                        
                        if(my_scaling == "mean"){
                                
                                lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp) 
                                
                                my_breaks <- seq(-1,1, length.out = 101)
                                hm_colors <- colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100)
                                
                        } else {
                                
                                my_feature_width_ordered <- my_feature_width[match(rownames(lnc_tmp), names(my_feature_width))]
                                
                                stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_ordered)))
                                lnc_scaled_tmp <- lnc_tmp - log2(my_feature_width_ordered)
                                
                                my_breaks <- seq(0,1, length.out = 101)
                                hm_colors <- colorRampPalette(rev(brewer.pal(9,name = "RdBu"))[4:9])(100)
                        }
                        
                        my_corr <- cor(lnc_scaled_tmp, method = "spearman", use = "complete.obs")
                        rownames(my_corr) <- SampleTable$ForeignID[match(colnames(my_corr), SampleTable$ForeignID)]
                        
                        pheatmap(my_corr,
                                 main = paste0("Spearman's R [",my_scaling,"]"),
                                 color = hm_colors,
                                 breaks = my_breaks,
                                 clustering_callback = callback,
                                 clustering_method = "ward.D2", 
                                 cellheight = ifelse(my_assay == "all", 15, 25), 
                                 cellwidth = ifelse(my_assay == "all", 15, 25))
                        
                        
                        rm(list = ls(pattern = "_tmp"))
                        rm(list = "my_corr")     
                        
                }
                
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
        
        
        par(mfrow=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))
        
        
        my_assay <- "all"
        lnc_tmp <- get(paste("lnc", my_assay, log_type, sep = "."))
        
        
        for(my_color_group in c("Stages", "Assay") ){
                
                
                my_conditions <- SampleTable[,my_color_group][match(colnames(lnc_tmp), SampleTable$ForeignID)]
                my_conditions <- factor(my_conditions, levels = unique(my_conditions))
                
                #################################################       
                
                if(my_color_group == "Stages"){
                        pca_colors <- my_color_palette
                        #my_conditions <- relevel(my_conditions, ref = "Zy")
                } else {
                        pca_colors <- my_color_palette2
                        my_conditions <- relevel(my_conditions,  ref = "IgG")
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
                            point_size = 1.1,
                            my_xlimits = c(-200,200),
                            my_ylimits = c(-200,200))
                
                plottingPCA(lnc_tmp,
                            xcomp = 3,
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
                
                legend("left",
                       legend = levels(my_conditions),
                       cex = 1, bty = "n",
                       fill =  pca_colors[seq_along(levels(my_conditions))])
                
        }
        
        
        #########################################################
        
        
        rm(list = ls(pattern = "_tmp"))
        
        
        dev.off()
        
}



##################################################################################################################################
##################################################################################################################################




















##################################################################################################################################
##################################################################################################################################

#################################################       Signal Density       #####################################################






# for(log_type in log_types){
#         
#         
#         pdf(paste0(output_plots,"/Density.", log_type,".pdf"), height = 8, width = 14, useDingbats = F)
#         
#         
#         for(my_assay in my_assays){
#                 
#                 if(my_assay == "all"){next()}
#                 
#                 par(mfcol=c(3,5), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))
#                 
#                 #################################################       
#                 
#                 
#                 lnc_tmp <- get(paste("lnc", my_assay, log_type, sep = "."))
#                 
#                 my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
#                 my_conditions <- factor(my_conditions, levels = unique(my_conditions))
#                 
#                 #################################################       
#                 
#                         dens_colors <- my_color_palette  
#                 
#                 #################################################       
#                 
#                 # my_feature_width_tmp <- my_feature_width[names(my_feature_width) %in% rownames(lnc_tmp)]
#                 # my_feature_width_tmp <- my_feature_width_tmp[order(names(my_feature_width_tmp))]
#                 # 
#                 # stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_tmp)))
#                 
#                 lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp[, grep("ChIP", colnames(lnc_tmp))], na.rm = TRUE)
#                 
#                 for(i in 1:ncol(lnc_scaled_tmp)){
#                         
#                         plot(density(lnc_scaled_tmp[,i], from = -5, to = 5,  na.rm = TRUE),
#                              col = dens_colors[my_conditions][i], lwd = 2,
#                              main = colnames(lnc_scaled_tmp)[i], xlab = "")
#                         
#                         abline(v = seq(-10,10,2), lty=3, col = "grey")
#                         
#                         if(i == 1){mtext(text = paste0("log2 Normalized Counts at ",feature_type," (relative to Zy)"), side = 1,line = -1, outer = TRUE)}
#                 }
#                 
#                 plotLegend(conditions = my_conditions,
#                            legend_colors = dens_colors,
#                            legend_size = 1.2)
#                 
#                 rm(list = ls(pattern = "_tmp"))
#         }
#         
#         dev.off()
# }



##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################

#################################################       Signal Violin       #####################################################



# for(log_type in log_types){
#         
#         for(my_assay in my_assays){
#                 
#                 if(my_assay == "all"){next()}
#                 
#                 lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
#                 
#                 my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
#                 my_conditions <- factor(my_conditions, levels = unique(my_conditions))
#                 
#                 #################################################           
#                 
#                         vio_colors <- my_color_palette  
#                 
#                 #################################################           
#                 
#                 # my_feature_width_tmp <- my_feature_width[names(my_feature_width) %in% rownames(lnc_tmp)]
#                 # my_feature_width_tmp <- my_feature_width_tmp[order(names(my_feature_width_tmp))]
#                 # 
#                 # stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_tmp)))
#                 # 
#                 # lnc_tmp_wnorm <- lnc_tmp - log2(my_feature_width_tmp)
#                 # lnc_tmp_wnorm$gene_id <- rownames(lnc_tmp_wnorm)
#                 
#                 lnc_scaled_tmp <-  lnc_tmp - rowMeans(lnc_tmp[, grep("ChIP", colnames(lnc_tmp))], na.rm = TRUE)
#                 lnc_scaled_tmp$gene_id <- rownames(lnc_scaled_tmp)
#                 
#                 lnc_tmp_long <- melt(lnc_scaled_tmp)
#                 
#                 lnc_tmp_long$Conditions <- gsub(paste0("X|_",my_assay,"|_[1-3]"),"", lnc_tmp_long$variable)
#                 lnc_tmp_long$Conditions <- factor(lnc_tmp_long$Conditions, levels = unique(lnc_tmp_long$Conditions))
#                 
#                 lnc_tmp_long$Classes <- my_DBTMEE$Cluster[match(lnc_tmp_long$gene_id, my_DBTMEE$gene_id)]
#                 
#                 ggp1 <- ggplot(lnc_tmp_long, aes(x=variable, y=value, fill = Conditions)) + 
#                         theme_bw() +
#                         theme(text = element_text(size = 12), aspect.ratio = 0.75, 
#                               axis.text.x = element_blank(), 
#                               panel.border = element_rect(colour = "black", fill=NA, size=1),
#                               panel.grid.major = element_blank(),
#                               panel.grid.minor = element_blank(),
#                               strip.background = element_blank(),
#                               plot.title = element_text(hjust = 0.5, face = "bold")) +
#                         geom_violin(width = 1.5) +
#                         facet_wrap(~ Classes, nrow = 3) +
#                         ggtitle(label = paste0(feature_type," - ",my_assay)) +
#                         stat_summary(fun = get(ifelse(my_setting %in% c("NAmedian","minmedian","filtmedian") , "median", "mean")), geom="point", size=1) +
#                         scale_fill_manual(values = vio_colors) +
#                         xlab(label = "Conditions") +
#                         coord_cartesian(ylim = c(-4, 4)) +
#                         ylab(label = "log2 Normalized Counts (relative to ChIP)")
#                 
#                 
#                 
#                 
#                 
#                 ggsave(filename = paste0(output_plots,"/vioplot.",my_assay,".", log_type,".pdf"),
#                        plot =  ggp1,
#                        width = 12, height = 8)
#                 
#                 
#                 #################################################           
#                 
#                 
#                 my_subclasses_list <- list(c("Maternal_RNA","1-Cell_Transient", "2-Cell_Transient", "4-Cell_Transient"),
#                                            c("Maternal_RNA", "Minor_ZGA",  "Major_ZGA", "MGA"))
#                 
#                 
#                 for(i in seq_along(my_subclasses_list)){
#                         
#                         lnc_tmp_long_sub <- lnc_tmp_long[lnc_tmp_long$Classes %in% my_subclasses_list[[i]],]
#                         
#                         lnc_tmp_long_sub$group <- factor(paste(lnc_tmp_long_sub$Classes, gsub(".*_","",lnc_tmp_long_sub$variable)),
#                                                          levels = (paste(rep(my_subclasses_list[[i]], each=3), 1:3)))
#                         
#                         
#                         ggp2 <- ggplot(lnc_tmp_long_sub, 
#                                        aes(x=group, y=value, fill = Conditions)) + 
#                                 theme_bw() +
#                                 theme(text = element_text(size = 12), aspect.ratio = 1, 
#                                       #axis.text.x = element_blank(), 
#                                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#                                       panel.border = element_rect(colour = "black", fill=NA, size=1),
#                                       panel.grid.major = element_blank(),
#                                       panel.grid.minor = element_blank(),
#                                       strip.background = element_blank(),
#                                       plot.title = element_text(hjust = 0.5, face = "bold")) +
#                                 geom_violin(width = 1) +
#                                 facet_wrap(~ Conditions, nrow = 1) +
#                                 ggtitle(label = paste0(feature_type," - ",my_assay)) +
#                                 stat_summary(fun = get(ifelse(my_setting %in% c("NAmedian","minmedian","filtmedian") , "median", "mean")), geom="point", size=1) +
#                                 scale_fill_manual(values = vio_colors) +
#                                 xlab(label = "DBTMEE Classes") +
#                                 coord_cartesian(ylim = c(-4, 4)) +
#                                 ylab(label = "log2 Normalized Counts (relative to ChIP)") + 
#                                 #coord_flip() +
#                                 geom_hline(yintercept = 0,  color="black", linetype="dashed") +
#                                 geom_vline(xintercept = c(3.5,6.5,9.5),  color="grey", linetype="dotted") 
#                         
#                         rm(list = "lnc_tmp_long_sub")
#                         
#                         ggsave(filename = paste0(output_plots,"/vioplot",i,".",my_assay,".", log_type,".pdf"),
#                                plot =  ggp2,
#                                width = 12, height = 6)
#                 }
#                 
#                 #################################################           
#                 
#                 
#                 # lnc_scaled_tmp <-  lnc_tmp - rowMeans(lnc_tmp[, grep("ChIP", colnames(lnc_tmp))], na.rm = TRUE)
#                 # 
#                 # lnc_scaled_tmp.repmeans <- as.data.frame(t(apply(lnc_scaled_tmp, 1, function(x){aggregate(x, by= list(my_conditions), mean, na.rm = TRUE)[,2]})))
#                 # 
#                 # colnames(lnc_scaled_tmp.repmeans) <- levels(my_conditions)
#                 # lnc_scaled_tmp.repmeans <- round(lnc_scaled_tmp.repmeans, 6)
#                 # lnc_scaled_tmp.repmeans$gene_id <- rownames(lnc_scaled_tmp.repmeans)
#                 # 
#                 # 
#                 # lnc_scaled_tmp.repmeans_long <- melt(lnc_scaled_tmp.repmeans)
#                 # 
#                 # lnc_scaled_tmp.repmeans_long$Conditions <- factor(gsub(paste0("X|_",my_assay,"|_[1-3]"),"", lnc_scaled_tmp.repmeans_long$variable))
#                 # lnc_scaled_tmp.repmeans_long$Conditions <- factor(lnc_scaled_tmp.repmeans_long$Conditions, levels = unique(lnc_scaled_tmp.repmeans_long$Conditions))
#                 # lnc_scaled_tmp.repmeans_long$Classes <- my_DBTMEE$Cluster[match(lnc_scaled_tmp.repmeans_long$gene_id, my_DBTMEE$gene_id)]
#                 # 
#                 # 
#                 # 
#                 # ggf <- ggplot(lnc_scaled_tmp.repmeans_long, aes(x=variable, y=value, fill = Conditions)) + 
#                 #         theme_bw() +
#                 #         theme(text = element_text(size = 12), aspect.ratio = 0.75, 
#                 #               axis.text.x = element_blank(), 
#                 #               panel.border = element_rect(colour = "black", fill=NA, size=1),
#                 #               panel.grid.major = element_blank(),
#                 #               panel.grid.minor = element_blank(),
#                 #               strip.background = element_blank(),
#                 #               plot.title = element_text(hjust = 0.5, face = "bold")) +
#                 #         geom_violin(width = 1.5) +
#                 #         facet_wrap(~ Classes, nrow = 3) +
#                 #         ggtitle(label = paste0(feature_type," - ",my_assay)) +
#                 #         stat_summary(fun = get(ifelse(my_setting %in% c("NAmedian","minmedian","filtmedian") , "median", "mean")), geom="point", size=1) +
#                 #         scale_fill_manual(values = vio_colors) +
#                 #         xlab(label = "Conditions") +
#                 #         coord_cartesian(ylim = c(-4, 4)) +
#                 #         ylab(label = "log2 Normalized Counts (relative to Zy)")
#                 # 
#                 # 
#                 # ggsave(filename = paste0(output_plots,"/vioplot.repmeans.",my_assay,".", log_type,".pdf"),
#                 #        plot =  ggf,
#                 #        width = 12, height = 8)
#                 
#                 
#                 rm(list = ls(pattern = "_tmp"))
#         }
#         
# }



##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################

##############################################    RNA-seq Cdk9in Spt5 results   ###################################################



res_external_files <- c("../RNAseq_Cdk9_Spt5/Output/deseq2/tables/res.Cdk9in-Noinject.txt",
                        "../RNAseq_Cdk9_Spt5/Output/deseq2/tables/res.Spt5KD-IgG.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.4C-Late_2C.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.8C-Late_2C.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.ICM-Late_2C.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.Early_2C-Zygote.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.Mid_2C-Zygote.txt",
                        "../../Public/GSE38495_GSE45719/Output/deseq2/tables/res.Late_2C-Zygote.txt")

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



# for(log_type in log_types){
#         
#         
#         for(my_assay in my_assays){
#                 
#                 
#                 if(my_assay == "all"){next()}
#                 
#                 #################################################           
#                 
#                 lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
#                 
#                 #lnc_tmp <- na.omit(lnc_tmp)
#                 
#                 my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
#                 my_conditions <- factor(my_conditions, levels = unique(my_conditions))
#                 
#                 #################################################           
#                 
#                         dot_colors <- my_color_palette  
#                 
#                 #################################################           
#                 
#                 
#                 # my_feature_width_tmp <- my_feature_width[names(my_feature_width) %in% rownames(lnc_tmp)]
#                 # my_feature_width_tmp <- my_feature_width_tmp[order(names(my_feature_width_tmp))]
#                 # 
#                 # stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_tmp)))
#                 # 
#                 # lnc_tmp_wnorm <- lnc_tmp - log2(my_feature_width_tmp)
#                 
#                 lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp[, grep("ChIP", colnames(lnc_tmp))], na.rm = TRUE)
#                 
#                 
#                 for(my_classification in c("my_DBTMEE", res_external_names)){
#                         
#                         
#                         
#                         #################################################           
#                         
#                         
#                         lnc_scaled_tmp_merged <- merge(lnc_scaled_tmp, get(my_classification), by.x = "row.names", by.y = "gene_id", all = FALSE)
#                         lnc_scaled_tmp_merged$Cluster <- factor(lnc_scaled_tmp_merged$Cluster)
#                         
#                         
#                         lnc_scaled_tmp_merged_mn <- t(sapply(levels(lnc_scaled_tmp_merged$Cluster),function(i){
#                                 
#                                 if(my_setting %in% c("NAmedian","minmedian","filtmedian")){
#                                         colMedians(as.matrix(lnc_scaled_tmp_merged[lnc_scaled_tmp_merged$Cluster == i, grep(my_assay, colnames(lnc_scaled_tmp_merged))]), na.rm = TRUE)
#                                 } else {
#                                         colMeans(as.matrix(lnc_scaled_tmp_merged[lnc_scaled_tmp_merged$Cluster == i, grep(my_assay, colnames(lnc_scaled_tmp_merged))]), na.rm = TRUE)
#                                 }
#                         }))
#                         
#                         #################################################           
#                         
#                         colnames(lnc_scaled_tmp_merged_mn) <- colnames(lnc_scaled_tmp_merged[, grep(my_assay, colnames(lnc_scaled_tmp_merged))])
#                         
#                         stopifnot(identical(colnames(lnc_scaled_tmp_merged_mn), colnames(lnc_scaled_tmp)))
#                         
#                         if(my_setting %in% c("NAmedian","minmedian","filtmedian")){
#                                 lnc_scaled_tmp_merged_mn <- rbind(colMedians(as.matrix(lnc_scaled_tmp), na.rm = TRUE), lnc_scaled_tmp_merged_mn)
#                         } else {
#                                 lnc_scaled_tmp_merged_mn <- rbind(colMeans(as.matrix(lnc_scaled_tmp), na.rm = TRUE), lnc_scaled_tmp_merged_mn)
#                         }
#                         
#                         rownames(lnc_scaled_tmp_merged_mn)[1] <- "all genes"
#                         
#                         
#                         #################################################           
#                         
#                         
#                         lnc_scaled_tmp_merged_mn_long <- melt(lnc_scaled_tmp_merged_mn, varnames = c("Classes", "Samples"))
#                         
#                         lnc_scaled_tmp_merged_mn_long$Conditions <- gsub(paste0("_",my_assay,"|_[1-3]"),"", lnc_scaled_tmp_merged_mn_long$Samples)
#                         lnc_scaled_tmp_merged_mn_long$Conditions <- factor(lnc_scaled_tmp_merged_mn_long$Conditions, levels = unique(lnc_scaled_tmp_merged_mn_long$Conditions))
#                         
#                         stopifnot(identical(levels(lnc_scaled_tmp_merged_mn_long$Conditions), levels(my_conditions)))
#                         
#                         #################################################           
#                         
#                         my_class <- "up-reg"
#                         
#                         
#                         pdf(paste0(output_plots,"/dotplot.",gsub("my_","",my_classification),".",my_assay,".", log_type,".pdf"), height = 9, width = 9, useDingbats = F)
#                         
#                         
#                         for(my_class in levels(lnc_scaled_tmp_merged_mn_long$Classes)){
#                                 
#                                 
#                                 par(mfrow=c(2,2), mar = c(5,10,3,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))
#                                 
#                                 
#                                 #################################################           
#                                 
#                                 
#                                 df_tmp <- lnc_scaled_tmp_merged_mn_long[lnc_scaled_tmp_merged_mn_long$Classes == my_class,]
#                                 
#                                 set.seed(123)
#                                 
#                                 plot(x = jitter(as.integer(df_tmp$Conditions),factor = 0.1),
#                                      y = df_tmp$value,
#                                      # main =  paste0(feature_type," ",my_assay, " at ", my_class, 
#                                      #                "\n[n=",  ifelse(my_class == "all genes", nrow(lnc_scaled_tmp), 
#                                      #                                 table(lnc_scaled_tmp_merged$Cluster)[my_class]), "]"),
#                                      main =  paste0(feature_type," ",my_assay, " at \n",
#                                                     gsub("my_|res.","", ifelse(my_class == "all genes","",my_classification))," ",my_class),
#                                      xlim = c(min(as.integer(df_tmp$Conditions))-0.5, max(as.integer(df_tmp$Conditions))+0.5),
#                                      ylim = c(min(df_tmp$value)-0.5, max(df_tmp$value)+0.5),
#                                      xlab = "", xaxt = "n",
#                                      ylab = ifelse(my_setting %in% c("NAmedian","minmedian","filtmedian") ,
#                                                    "Median log2 Normalized Counts \n(relative to ChIP)",
#                                                    "Mean log2 Normalized Counts \n(relative to ChIP)"),
#                                      pch = 19, col = dot_colors[df_tmp$Conditions])
#                                 
#                                 set.seed(123)
#                                 
#                                 points(x = jitter(as.integer(df_tmp$Conditions), factor = 0.1),
#                                        y =  df_tmp$value,
#                                        col = "#555555", pch=1, lwd=0.5)
#                                 
#                                 axis(side = 1, at = seq_along(levels(df_tmp$Conditions)), 
#                                      labels = levels(df_tmp$Conditions), las = 2)
#                                 
#                                 
#                                 #################################################           
#                                 
#                                 par(mar = c(5,2,3,10))
#                                 
#                                 plot.new()
#                                 
#                                 legend("left",
#                                        legend = levels(df_tmp$Conditions),
#                                        cex = 1, bty = "n",
#                                        fill =  dot_colors[seq_along(levels(df_tmp$Conditions))])
#                                 
#                                 
#                                 #################################################           
#                                 
#                                 
#                                 
#                                 fit <- lm(value ~ Conditions, df_tmp)
#                                 
#                                 set.seed(2)
#                                 
#                                 
#                                 glht_out <- summary(glht(fit, linfct = mcp(Conditions = c("`ChIL100C` - `ChIP` = 0",
#                                                                                           "`ChIL1000C` - `ChIP` = 0",
#                                                                                           "`ChIL1000C` - `ChIL100C` = 0"))))
#                                 
#                                 
#                                 stats_tmp <- data.frame(`Estimate` =   round(glht_out$test$coefficients*100, 2)/100,
#                                                         `Std. Error` = round(glht_out$test$sigma*100, 2)/100,
#                                                         `t value` =  round(glht_out$test$tstat*100, 2)/100,
#                                                         `Pr(>|t|)` = round(glht_out$test$pvalues*100, 2)/100, 
#                                                         check.names = FALSE)
#                                 
#                                 
#                                 textplot(stats_tmp, halign = "left", valign = "top", cex = 0.7)
#                                 text(x = 0.5, y = 0.75, cex = 0.8,
#                                      labels = "Multiple Comparisons of Means: User-defined Contrasts")
#                                 
#                                 
#                                 
#                                 rm(list = "df_tmp")
#                                 rm(list = "stats_tmp")
#                                 
#                                 
#                                 #################################################           
#                         }
#                         
#                         dev.off()
#                         
#                         
#                         ################################################# 
#                         
#                         if(my_classification == "my_DBTMEE"){
#                                 my_subclasses_list <- list(c("Maternal_RNA","1-Cell_Transient", "2-Cell_Transient", "4-Cell_Transient"),
#                                                            c("Maternal_RNA", "Minor_ZGA",  "Major_ZGA", "MGA"))         
#                         } else {
#                                 my_subclasses_list <- list(c("down-reg", "NS",  "up-reg"))   
#                         }
#                         
#                         
#                         ################################################# 
#                         
#                         if(feature_type == "genebody"){
#                                 
#                                 if(my_assay == "S5P"){
#                                         my_xlims <- c(-0.5, 1.5)  
#                                 } else {
#                                         my_xlims <- c(-0.5, 1.5) 
#                                 }
#                         } else {
#                                 if(my_assay == "S5P"){
#                                         my_xlims <- c(-0.5, 2.0)  
#                                 } else {
#                                         my_xlims <- c(-0.5, 1.0) 
#                                 }
#                                 
#                                 
#                         }
#                         
#                         ################################################# 
#                         
#                         for(i in seq_along(my_subclasses_list)){
#                                 
#                                 
#                                 pdf(paste0(output_plots,"/dotplot.",gsub("my_","",my_classification),".",i,".",my_assay,".", log_type,".pdf"), height = 5, width = 12, useDingbats = F)
#                                 par(mfrow=c(2,3), mar = c(3,8,3,2), oma = c(3,3,3,3), mgp = c(3,1,0))
#                                 
#                                 
#                                 for(my_condition in levels(lnc_scaled_tmp_merged_mn_long$Conditions)){
#                                         
#                                         
#                                         df_tmp <- lnc_scaled_tmp_merged_mn_long[lnc_scaled_tmp_merged_mn_long$Conditions %in% my_condition,]
#                                         
#                                         df_tmp_sub <- df_tmp[df_tmp$Classes %in% c("all genes", my_subclasses_list[[i]]),]
#                                         df_tmp_sub$Classes <- factor(df_tmp_sub$Classes, levels = rev(c("all genes", my_subclasses_list[[i]])))
#                                         
#                                         set.seed(123)
#                                         
#                                         plot(y = jitter(as.integer(df_tmp_sub$Classes), factor = 0.1),
#                                              x = df_tmp_sub$value,
#                                              main =  my_condition,
#                                              ylim = c(min(as.integer(df_tmp_sub$Classes))-0.5, max(as.integer(df_tmp_sub$Classes))+0.5),
#                                              xlim = my_xlims,
#                                              ylab = "", xlab = "", yaxt = "n",
#                                              pch = 19, col = dot_colors[which(my_condition == levels(df_tmp_sub$Conditions))])
#                                         
#                                         set.seed(123)
#                                         
#                                         points(y = jitter(as.integer(df_tmp_sub$Classes), factor = 0.1),
#                                                x =  df_tmp_sub$value,
#                                                col = "#555555", pch=1, lwd=0.5)
#                                         
#                                         abline(v = mean(df_tmp_sub$value[df_tmp_sub$Classes == "all genes"]), lty=3)
#                                         
#                                         axis(side = 2, at = seq_along(levels(df_tmp_sub$Classes)), labels = levels(df_tmp_sub$Classes), las=2)
#                                         
#                                         rm(list = ls(pattern = "df_tmp"))
#                                         
#                                         
#                                 }
#                                 
#                                 # textplot(data.frame(`n class` = as.integer(table(lnc_scaled_tmp_merged$Cluster)),
#                                 #                     row.names = names(table(lnc_scaled_tmp_merged$Cluster))), halign = "center", valign = "center", cex = 1)
#                                 
#                                 mtext(text = paste0(feature_type," - ",my_assay), side = 3, line = 0, font = 2,  outer = TRUE, cex = 1.2)
#                                 mtext(text = ifelse(my_setting %in% c("NAmedian","minmedian","filtmedian") ,
#                                                     "Median log2 Normalized Counts (relative to ChIP)",
#                                                     "Mean log2 Normalized Counts (relative to ChIP)"), side = 1, line = 0.5, outer = TRUE)
#                                 mtext(text = paste0(gsub("my_|res.","",my_classification), " Classes"), 
#                                       side = 2, line = 0.5, outer = TRUE)
#                                 
#                                 
#                                 
#                                 dev.off()
#                                 
#                         }
#                         
#                         rm(list = "lnc_scaled_tmp_merged")
#                 }
#                 
#                 rm(list = ls(pattern = "_tmp"))
#         }
#         
# }



##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

#################################################     compare toRNA-seq      #####################################################




# log2_TPM_DBTMEE <- merge(log2_zscore_TPM, my_DBTMEE, by.x = "row.names", by.y = "gene_id", all = FALSE)
# 
# log2_TPM_DBTMEE_long <- melt(log2_TPM_DBTMEE)
# 
# 
# ggp3 <- ggplot(log2_TPM_DBTMEE_long, aes(x=variable, y=value)) + 
#         theme_bw() +
#         theme(text = element_text(size = 12), aspect.ratio = 0.75, 
#               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#               panel.border = element_rect(colour = "black", fill=NA, size=1),
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               strip.background = element_blank(),
#               plot.title = element_text(hjust = 0.5, face = "bold")) +
#         geom_violin(width = 1, fill = "grey") +
#         facet_wrap(~ Cluster, nrow = 2) +
#         ggtitle(label = "RNA-seq in DBTMEE classes") +
#         stat_summary(fun= get(ifelse(my_setting %in% c("NAmedian","minmedian","filtmedian") , "median", "mean")), geom="point", size=1) +
#         scale_fill_manual(values = vio_colors) +
#         coord_cartesian(ylim = c(-4, 4)) +
#         ylab(label = "log2 TPM (z-score)") +
#         xlab(label = "Stages") 
# 
# ggsave(filename = paste0(output_plots,"/vioplot.RNAseq_DBTMEE.pdf"),
#        plot =  ggp3,
#        width = 10, height = 6)

#################################################        

##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################




# for(log_type in log_types){
#         
#         for(my_assay in my_assays){
#                 
#                 if(my_assay == "all"){next()}
#                 
#                 ################################################# 
#                 
#                 if(feature_type == "genebody"){
#                         
#                         if(my_assay == "S5P"){
#                                 my_xlims <- c(-0.5, 2.0)  
#                         } else {
#                                 my_xlims <- c(-0.5, 1.5) 
#                         }
#                 } else {
#                         if(my_assay == "S5P"){
#                                 my_xlims <- c(-0.5, 2.0)  
#                         } else {
#                                 my_xlims <- c(-0.5, 1.0) 
#                         }
#                         
#                         
#                 }
#                 
#                 ################################################# 
#                 
#                 
#                 
#                 lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
#                 
#                 #lnc_tmp <- na.omit(lnc_tmp)
#                 
#                 my_conditions <- SampleTable$Stages[match(colnames(lnc_tmp), SampleTable$ForeignID)]
#                 my_conditions <- factor(my_conditions, levels = unique(my_conditions))
#                 
#                 #################################################           
#                 
#                         dot_colors <- my_color_palette  
#                 
#                 #################################################           
#                 
#                 # lnc_tmp <- as.data.frame(t(apply(lnc_tmp, 1, function(x){aggregate(x, by= list(my_conditions), mean, na.rm = TRUE)[,2]})))
#                 # colnames(lnc_tmp) <- paste(levels(my_conditions), my_assay, sep = "_")
#                 
#                 
#                 lnc_scaled_tmp <-  lnc_tmp - rowMeans(lnc_tmp[, grep("ChIP", colnames(lnc_tmp))], na.rm = TRUE)
#                 
#                 
#                 ################################################# 
#                 
#                 
#                 
#                 
#                 
#                 pdf(paste0(output_plots,"/dotplot_RNAseq.",my_assay,".", log_type,".pdf"), height = 5.5, width = 12, useDingbats = F)
#                 
#                 
#                 for(qs in c("Q5","Q35")){
#                         
#                         
#                         log2_TPM_Qs <- get(paste0("log2_TPM_", qs))
#                         
#                         par(mfrow=c(2,3), mar = c(3,8,3,2), oma = c(3,3,3,3), mgp = c(3,1,0))
#                         
#                         
#                         for(my_condition in levels(my_conditions)){
#                                 
#                                 
#                                 df_tmp <- sapply( names(log2_TPM_Qs), function(x){
#                                         
#                                         
#                                         if(my_setting %in% c("NAmedian","minmedian","filtmedian")){
#                                                 colMedians(as.matrix(lnc_scaled_tmp[rownames(lnc_scaled_tmp) %in% log2_TPM_Qs[[x]], my_condition == my_conditions]), na.rm = TRUE)
#                                         } else {
#                                                 colMeans(as.matrix(lnc_scaled_tmp[rownames(lnc_scaled_tmp) %in% log2_TPM_Qs[[x]], my_condition == my_conditions]), na.rm = TRUE)
#                                         }
#                                 })
#                                 
#                                 
#                                 df_tmp_long <- melt(df_tmp)
#                                 df_tmp_long$Var2 <- factor(df_tmp_long$Var2, levels = rev(levels(df_tmp_long$Var2)))
#                                 
#                                 set.seed(123)
#                                 
#                                 plot(y = jitter(as.integer(df_tmp_long$Var2), factor = 0.1),
#                                      x = df_tmp_long$value,
#                                      main =  my_condition,
#                                      ylim = c(min(as.integer(df_tmp_long$Var2))-0.5, max(as.integer(df_tmp_long$Var2))+0.5),
#                                      xlim = my_xlims,
#                                      ylab = "", xlab = "", yaxt = "n",
#                                      pch = 19, col = dot_colors[which(my_condition == levels(my_conditions))])
#                                 
#                                 
#                                 set.seed(123)
#                                 
#                                 points(y = jitter(as.integer(df_tmp_long$Var2), factor = 0.1),
#                                        x =  df_tmp_long$value,
#                                        col = "#555555", pch=1, lwd=0.5)
#                                 
#                                 abline(v = mean(df_tmp_long$value[df_tmp_long$Var2 == paste0(qs,"_Zygote")]), lty=3)
#                                 
#                                 axis(side = 2, at = seq_along(levels(df_tmp_long$Var2)), labels = levels(df_tmp_long$Var2), las=2)
#                                 
#                                 rm(list = ls(pattern = "df_tmp"))
#                                 
#                         }
#                         
#                         rm(list = "log2_TPM_Qs")
#                         
#                         mtext(text = paste0(feature_type," - ",my_assay), side = 3, line = 0, font = 2,  outer = TRUE, cex = 1.2)
#                         mtext(text = ifelse(my_setting %in% c("NAmedian","minmedian","filtmedian") ,
#                                             "Median log2 Normalized Counts (relative to ChIP)",
#                                             "Mean log2 Normalized Counts (relative to ChIP)"), side = 1, line = 0.5, outer = TRUE)
#                         mtext(text = paste0("RNA-seq ",qs," genes"), side = 2, line = 0.5, outer = TRUE)
#                         
#                 }
#                 
#                 
#                 
#                 dev.off()
#                 
#         }
#         
#         
#         ################################################# 
#         
#         rm(list = ls(pattern = "_tmp"))
# }



##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################




















##################################################################################################################################
##################################################################################################################################







writeLines(capture.output(sessionInfo()), output_sessionInfo)






##################################################################################################################################
##################################################################################################################################


