






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


my_lnc_file <- my_lnc_files[4]


for(my_lnc_file in my_lnc_files){
        
        
        my_lnc_name <- gsub(".*\\/|.txt","",my_lnc_file)
        
        my_lnc <- read.table(my_lnc_file, check.names=FALSE)
        my_lnc <- my_lnc[order(rownames(my_lnc)),]
        my_lnc <- my_lnc[, c(grep("Zy", colnames(my_lnc)), grep("2C|8C|ES", colnames(my_lnc)))] 
        
        
        if(my_setting %in% c("min")){
                
                set.seed(456)
                
                my_lnc[is.na(my_lnc)] <- rnorm(n = sum(is.na(my_lnc)), mean = min(my_lnc, na.rm = TRUE), sd = 0.1)
                
                
        } else if(my_setting %in% c("filt")){
                
                my_lnc <- na.omit(my_lnc)
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


my_histone_files <- list.files(path = paste0("../../Public/Histone_data/Output/deseq2_", feature_type,"/results/"), 
                               pattern = paste0("^lnc.*","log2.txt"), full.names = TRUE)


# H3K4me3 - TSS should be taken not genebody !!!
my_histone_files <- gsub("deseq2_genebody/results//lnc.H3K4me3","deseq2_TSS/results//lnc.H3K4me3",my_histone_files)

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
                      len = c(0.8, 0.6, 0.8, 0.7, 0.8, 0.7, 0.9, 0.8),
                      zy = c(0.5, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.6), 
                      mean = c(0.5, 0.4, 0.7, 0.5, 0.7, 0.5, 0.7, 0.6), 
                      stringsAsFactors = FALSE)


#################################################        


for(my_histone in c("ATACseq","H3K4me3", "H3K36me3", "Pol2")){
        
        
        my_histone_file <- grep(my_histone, my_histone_files, value = TRUE)
        
        hislnc_tmp <- read.table(my_histone_file, check.names=FALSE)
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
                
                if(grepl("TSS", my_histone_file)){
                        
                        hislnc_scaled_tmp <- hislnc_repmeans_tmp - log2(5)      
                        
                } else {
                        
                        my_feature_width_ordered <- my_feature_width[match(rownames(hislnc_repmeans_tmp), names(my_feature_width))]
                        
                        stopifnot(identical(rownames(hislnc_repmeans_tmp), names(my_feature_width_ordered)))
                        
                        hislnc_scaled_tmp <- hislnc_repmeans_tmp - log2(my_feature_width_ordered)      
                }
                
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
                        
                        
                        if(my_scaling %in% c("mean")){
                                
                                lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp, na.rm = TRUE)
                                
                        } else if(my_scaling %in% c("len")){
                                
                                my_feature_width_ordered <- my_feature_width[match(rownames(lnc_tmp), names(my_feature_width))]
                                
                                stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_ordered)))
                                
                                lnc_scaled_tmp <- lnc_tmp - log2(my_feature_width_ordered)
                                
                        }  else if(my_scaling %in% c("zy")){
                                
                                lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp[, grep("Zy", colnames(lnc_tmp))], na.rm = TRUE)
                        }
                        
                        
                        #################################################           
                        
                        
                        lnc_merged_tmp <- merge(hislnc_scaled_tmp, lnc_scaled_tmp, by = "row.names")
                        
                        
                        my_corr <- cor(lnc_merged_tmp[,-1], method = "spearman", use = "complete.obs")
                        my_corr <- my_corr[grepl(my_assay, rownames(my_corr)),]
                        my_corr <- my_corr[,!grepl(my_assay, colnames(my_corr))]
                        
                        my_corr[is.na(my_corr)] <- 0
                        
                        
                        #################################################           
                        
                        
                        if(my_histone == "Pol2"){
                                
                                my_corr <- my_corr[,c(grep("PN3", colnames(my_corr)),
                                                      #grep("zygote", colnames(my_corr)),
                                                      #grep("PN5", colnames(my_corr)),
                                                      #grep("E2C", colnames(my_corr)),
                                                      grep("^2C", colnames(my_corr)),
                                                      #grep("^Ctrl", colnames(my_corr)),
                                                      grep("^DRB", colnames(my_corr)),
                                                      #grep("^4C", colnames(my_corr)),
                                                      grep("^8C", colnames(my_corr))#,
                                                      #grep("morula", colnames(my_corr)),
                                                      #grep("^ICM", colnames(my_corr))#,
                                                      #grep("^TE", colnames(my_corr))
                                )]
                                
                                my_corr <- my_corr[!grepl("ES", rownames(my_corr)),]
                                
                        } else {
                                
                                my_corr <- my_corr[,c(grep("oocyte", colnames(my_corr)),
                                                      #grep("PN3", colnames(my_corr)),
                                                      grep("zygote", colnames(my_corr)),
                                                      #grep("PN5", colnames(my_corr)),
                                                      grep("E2C", colnames(my_corr)),
                                                      grep("^2C", colnames(my_corr)),
                                                      #grep("^Ctrl", colnames(my_corr)),
                                                      #grep("^DRB", colnames(my_corr)),
                                                      grep("^4C", colnames(my_corr)),
                                                      grep("^8C", colnames(my_corr)),
                                                      #grep("morula", colnames(my_corr)),
                                                      grep("^ICM", colnames(my_corr))#,
                                                      #grep("^TE", colnames(my_corr))
                                )]
                                
                                my_corr <- my_corr[!grepl("ES|DRB", rownames(my_corr)),]
                        }
                        
                        
                        #################################################           
                        
                        
                        pdf(paste0(output_plots,"/corr2",my_histone,".", my_assay,".", log_type,".reps.pdf"), width = 10, height = 8, useDingbats = FALSE)
                        par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                        
                        
                        hm_colors <- colorRampPalette(brewer.pal(9,name = "YlGnBu"))(100)
                        
                        pheatmap(my_corr,
                                 main = paste0(my_histone,"\nSpearman's R"),
                                 color = hm_colors,
                                 cluster_cols = FALSE, 
                                 cluster_rows = FALSE, 
                                 cellheight = ifelse(my_assay == "all", 15, 25), 
                                 cellwidth = ifelse(my_assay == "all", 15, 25))
                        
                        
                        dev.off()
                        
                        
                        #################################################           
                        
                        
                        my_corr_long <- melt(my_corr)
                        
                        my_corr_long$Conditions <- gsub("_[1-3]","",my_corr_long$Var1)
                        my_corr_long$Conditions <- factor(my_corr_long$Conditions, levels = unique(my_corr_long$Conditions))
                        
                        ggp <- ggbarplot(data = my_corr_long, 
                                         x="Var2", y="value", 
                                         add = c("mean_se","jitter"),
                                         fill = "#666666",
                                         facet.by  = "Conditions", nrow = 1,
                                         palette = my_color_palette[-3]) + 
                                theme_bw() +
                                theme(text = element_text(size = 12), aspect.ratio = 1.25, 
                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                      axis.ticks.x = element_blank(),
                                      panel.border = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      strip.background = element_blank(),
                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                scale_y_continuous(breaks = seq(0,0.8,0.2)) +
                                geom_segment(data = data.frame(x=-Inf,xend=-Inf, y=0,yend=0.8), aes(x=x,xend=xend,y=y,yend=yend), inherit.aes = FALSE) +
                                #ggtitle(label = my_assay) +
                                ylab(label = "Spearman's R\n") +
                                xlab(label = paste(my_histone,"Stages")) 
                        
                        ggsave(filename = paste0(output_plots,"/corr2",my_histone,".", my_assay,".", log_type,".bar.pdf"),
                               plot =  ggp,
                               width = 8, height = 4)
                        
                        
                        ################################################# 
                        
                        lnc_scaled_tmp <- lnc_scaled_tmp[,grep("Zy|2C_S[2,5]P_[1-3]|2C_IgG_[1-3]|8C", colnames(lnc_scaled_tmp))]
                        
                        vio_colors <- my_color_palette[-3]
                        
                        
                        hislnc_Qs_tmp <- as.data.frame(na.omit(hislnc_scaled_tmp))
                        
                        hisStage <- colnames(hislnc_Qs_tmp)[1]
                        
                        for(hisStage in  colnames(hislnc_Qs_tmp)){
                            
                            Qstage <- paste0("Q.", hisStage)
                            
                            hislnc_Qs_tmp[, Qstage] <- c("q1","q2","q3","q4","q5")[cut(x = hislnc_Qs_tmp[, hisStage], 
                                                                                       breaks = quantile(hislnc_Qs_tmp[, hisStage], seq(0,1,0.2)), 
                                                                                       include.lowest = TRUE)]
                        }
                        
                        hislnc_Qs_tmp <- hislnc_Qs_tmp[, grep("^Q.", colnames(hislnc_Qs_tmp)), drop = FALSE]
                        
                        hisStage <- colnames(hislnc_Qs_tmp)[1]
                        
                        for(hisStage in  grep("zygote|2C|8C", colnames(hislnc_Qs_tmp), value = TRUE)){
                            
                            
                            lnc_merged_Qs_tmp <- merge(hislnc_Qs_tmp[, hisStage, drop=FALSE], 
                                                       lnc_scaled_tmp, by = "row.names")
                            
                            lnc_merged_Qs_tmp_long <- melt(lnc_merged_Qs_tmp)
                            
                            lnc_merged_Qs_tmp_long$Stage <- gsub("_S5P|_IgG|_S2P|_[1-3]$","",lnc_merged_Qs_tmp_long$variable)
                            lnc_merged_Qs_tmp_long$Stage <- factor(lnc_merged_Qs_tmp_long$Stage, levels = unique(lnc_merged_Qs_tmp_long$Stage))
                            
                            
                            ggp <- ggplot(lnc_merged_Qs_tmp_long, 
                                          aes_string(x= hisStage, y="value", fill = "Stage")) + 
                                theme_bw() +
                                theme(text = element_text(size = 12), aspect.ratio = 1, 
                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      strip.background = element_blank(),
                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                geom_violin(width = 1) +
                                facet_wrap(~ variable, nrow = 3, dir = "v") +
                                ggtitle(label = paste0(feature_type," - ",my_assay)) +
                                stat_summary(fun =  "median", geom="point", size=1) +
                                scale_fill_manual(values = vio_colors) +
                                xlab(label = paste(my_histone, hisStage)) +
                                coord_cartesian(ylim = c(-6, 4)) +
                                scale_y_continuous(breaks = seq(-6,4,2)) +
                                ylab(label = paste(my_assay, "\nlog2 Normalized Counts"))
                            
                            ggsave(filename = paste0(output_plots,"/vioplot.",my_histone,".",hisStage,".",my_assay,".", log_type,".reps.pdf"),
                                   plot =  ggp,
                                   width = 12, height = 6)
                            
                            rm(list = "ggp")
                        }
                        
                        ################################################# 
                        
                        
                        rm(list = "my_corr")
                        rm(list = "lnc_scaled_tmp")
                        rm(list = "lnc_merged_tmp")
                        
                        
                        
                        #################################################    
                        
                        
                        if(identical(gsub(paste0(my_assay,"|_[1-3]|_"),"", colnames(lnc_tmp)),
                                     gsub(paste0(my_assay,"|_[1-3]|_"),"", as.character(my_conditions)))){
                                
                                
                                lnc_repmeans_tmp <- t(apply(lnc_tmp, 1, function(x){
                                        aggregate(x, by = list(my_conditions), FUN = mean, na.rm = TRUE)$x}))
                                
                                colnames(lnc_repmeans_tmp) <- levels(my_conditions) 
                        }
                        
                        
                        #################################################    
                        
                        
                        if(my_scaling %in% c("mean")){
                                
                                lnc_scaled_tmp <- lnc_repmeans_tmp - rowMeans(lnc_repmeans_tmp, na.rm = TRUE)
                                
                        } else if(my_scaling %in% c("len")){
                                
                                my_feature_width_ordered <- my_feature_width[match(rownames(lnc_repmeans_tmp), names(my_feature_width))]
                                
                                stopifnot(identical(rownames(lnc_repmeans_tmp), names(my_feature_width_ordered)))
                                
                                lnc_scaled_tmp <- lnc_repmeans_tmp - log2(my_feature_width_ordered)
                                
                        }  else if(my_scaling %in% c("zy")){
                                
                                lnc_scaled_tmp <- lnc_repmeans_tmp - rowMeans(lnc_repmeans_tmp[, grep("Zy", colnames(lnc_repmeans_tmp))], na.rm = TRUE)
                        }
                        
                        #################################################
                        
                        
                        
                        lnc_merged_tmp <- merge(hislnc_scaled_tmp, lnc_scaled_tmp, by = "row.names")
                        
                        
                        my_corr <- cor(lnc_merged_tmp[,-1], method = "spearman", use = "complete.obs")
                        my_corr <- my_corr[grepl(my_assay, rownames(my_corr)),]
                        my_corr <- my_corr[,!grepl(my_assay, colnames(my_corr))]
                        
                        
                        #################################################           
                        
                        
                        if(my_histone == "Pol2"){
                                
                                my_corr <- my_corr[,c(grep("PN3", colnames(my_corr)),
                                                      #grep("zygote", colnames(my_corr)),
                                                      #grep("PN5", colnames(my_corr)),
                                                      #grep("E2C", colnames(my_corr)),
                                                      grep("^2C", colnames(my_corr)),
                                                      #grep("^Ctrl", colnames(my_corr)),
                                                      grep("^DRB", colnames(my_corr)),
                                                      #grep("^4C", colnames(my_corr)),
                                                      grep("^8C", colnames(my_corr))#,
                                                      #grep("morula", colnames(my_corr)),
                                                      #grep("^ICM", colnames(my_corr))#,
                                                      #grep("^TE", colnames(my_corr))
                                )]
                                
                                my_corr <- my_corr[!grepl("ES", rownames(my_corr)),]
                                
                        } else {
                                
                                my_corr <- my_corr[,c(grep("oocyte", colnames(my_corr)),
                                                      #grep("PN3", colnames(my_corr)),
                                                      grep("zygote", colnames(my_corr)),
                                                      #grep("PN5", colnames(my_corr)),
                                                      grep("E2C", colnames(my_corr)),
                                                      grep("^2C", colnames(my_corr)),
                                                      #grep("^Ctrl", colnames(my_corr)),
                                                      #grep("^DRB", colnames(my_corr)),
                                                      grep("^4C", colnames(my_corr)),
                                                      grep("^8C", colnames(my_corr)),
                                                      #grep("morula", colnames(my_corr)),
                                                      grep("^ICM", colnames(my_corr))#,
                                                      #grep("^TE", colnames(my_corr))
                                )]
                                
                                my_corr <- my_corr[!grepl("ES|DRB", rownames(my_corr)),]
                        }
                        
                        
                        ################################################# 
                        
                        pdf(paste0(output_plots,"/corr2",my_histone,".", my_assay,".", log_type,".pool.pdf"), width = 10, height = 8, useDingbats = FALSE)
                        par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                        
                        
                        pheatmap(my_corr,
                                 main = paste0(my_histone,"\nSpearman's R"),
                                 color = hm_colors,
                                 cluster_cols = FALSE,
                                 cluster_rows = FALSE,
                                 cellheight = ifelse(my_assay == "all", 15, 25),
                                 cellwidth = ifelse(my_assay == "all", 15, 25))
                        
                        
                        rm(list = "my_corr")
                        
                        dev.off()
                        
                        #################################################
                        
                        
                        lnc_scaled_tmp <- lnc_scaled_tmp[,grep("Zy|2C_S|2C_I|8C", colnames(lnc_scaled_tmp))]
                        
                        vio_colors <- my_color_palette[-3]
                        
                        
                        hislnc_Qs_tmp <- as.data.frame(na.omit(hislnc_scaled_tmp))
                        
                        hisStage <- colnames(hislnc_Qs_tmp)[1]
                        
                        for(hisStage in  colnames(hislnc_Qs_tmp)){
                            
                            Qstage <- paste0("Q.", hisStage)
                            
                            hislnc_Qs_tmp[, Qstage] <- c("q1","q2","q3","q4","q5")[cut(x = hislnc_Qs_tmp[, hisStage], 
                                                                                       breaks = quantile(hislnc_Qs_tmp[, hisStage], seq(0,1,0.2)), 
                                                                                       include.lowest = TRUE)]
                        }
                        
                        hislnc_Qs_tmp <- hislnc_Qs_tmp[, grep("^Q.", colnames(hislnc_Qs_tmp)), drop = FALSE]
                        
                        hisStage <-  colnames(hislnc_Qs_tmp)[1]
                        
                        for(hisStage in  grep("zygote|2C|8C", colnames(hislnc_Qs_tmp), value = TRUE)){
                            
                            
                            lnc_merged_Qs_tmp <- merge(hislnc_Qs_tmp[, hisStage, drop=FALSE], 
                                                       lnc_scaled_tmp, by = "row.names")
                            
                            lnc_merged_Qs_tmp_long <- melt(lnc_merged_Qs_tmp)
                            
                            lnc_merged_Qs_tmp_long$Stage <- gsub("_S5P|_S2P|_[1-3]$","",lnc_merged_Qs_tmp_long$variable)
                            lnc_merged_Qs_tmp_long$Stage <- factor(lnc_merged_Qs_tmp_long$Stage, levels = unique(lnc_merged_Qs_tmp_long$Stage))
                            
                            
                            ggp <- ggplot(lnc_merged_Qs_tmp_long, 
                                          aes_string(x= hisStage, y="value", fill = "Stage")) + 
                                theme_bw() +
                                theme(text = element_text(size = 12), aspect.ratio = 1, 
                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      strip.background = element_blank(),
                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                geom_violin(width = 1) +
                                facet_wrap(~ variable, nrow = 1, dir = "v") +
                                ggtitle(label = paste0(feature_type," - ",my_assay)) +
                                stat_summary(fun =  "median", geom="point", size=1) +
                                scale_fill_manual(values = vio_colors) +
                                xlab(label = paste(my_histone, hisStage)) +
                                coord_cartesian(ylim = c(-6, 4)) +
                                scale_y_continuous(breaks = seq(-6,4,2)) +
                                ylab(label = paste(my_assay, "\nlog2 Normalized Counts"))
                            
                            ggsave(filename = paste0(output_plots,"/vioplot.",my_histone,".",hisStage,".",my_assay,".", log_type,".mean.pdf"),
                                   plot =  ggp,
                                   width = 8, height = 4)
                            
                            rm(list = "ggp")
                        }
                        
                        
                        
                        ################################################# 
                        
                        lnc_scaled_tmp <- lnc_scaled_tmp[,grep("Zy|2C_S|2C_I|8C", colnames(lnc_scaled_tmp))]
                        
                        if(my_histone == "ATACseq"){
                            hisStages <- c("E2C","2C","8C")
                        } else {
                            hisStages <- c("zygote","2C","8C")
                        }
                        
                        vio_colors <- my_color_palette[-3]
                        
                        
                        lnc_Qs_tmp <- as.data.frame(na.omit(lnc_scaled_tmp))
                        
                        myStage <- colnames(lnc_Qs_tmp)[1]
                        
                        for(myStage in  colnames(lnc_Qs_tmp)){
                            
                            Qstage <- paste0("Q.", myStage)
                            
                            lnc_Qs_tmp[, Qstage] <- c("q1","q2","q3","q4","q5")[cut(x = lnc_Qs_tmp[, myStage], 
                                                                                    breaks = quantile(lnc_Qs_tmp[, myStage], seq(0,1,0.2)), 
                                                                                    include.lowest = TRUE)]
                        }
                        
                        lnc_Qs_tmp <- lnc_Qs_tmp[, grep("^Q.", colnames(lnc_Qs_tmp)), drop = FALSE]
                        
                        myStage <-  colnames(lnc_Qs_tmp)[1]
                        
                        for(myStage in  colnames(lnc_Qs_tmp)){
                            
                            
                            lnc_merged_Qs_tmp <- merge(lnc_Qs_tmp[, myStage, drop=FALSE], 
                                                       hislnc_scaled_tmp[, hisStages], by = "row.names")
                            
                            lnc_merged_Qs_tmp_long <- melt(lnc_merged_Qs_tmp)
                            
                            lnc_merged_Qs_tmp_long$Stage <- factor(lnc_merged_Qs_tmp_long$variable, levels = unique(lnc_merged_Qs_tmp_long$variable))
                            
                            
                            ggp <- ggplot(lnc_merged_Qs_tmp_long, 
                                          aes_string(x= myStage, y="value", fill = "Stage")) + 
                                theme_bw() +
                                theme(text = element_text(size = 12), aspect.ratio = 1, 
                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      strip.background = element_blank(),
                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                geom_violin(width = 1) +
                                facet_wrap(~ variable, nrow = 1, dir = "v") +
                                ggtitle(label = paste0(feature_type," - ",my_assay)) +
                                stat_summary(fun =  "median", geom="point", size=1) +
                                scale_fill_manual(values = vio_colors) +
                                xlab(label = gsub("_",": ",myStage)) +
                                coord_cartesian(ylim = c(-4, 8)) +
                                scale_y_continuous(breaks = seq(-4,8,2)) +
                                ylab(label = my_histone)
                            
                            ggsave(filename = paste0(output_plots,"/vioplot.",my_assay,".", log_type, ".",myStage,".",my_histone,".mean.pdf"),
                                   plot =  ggp,
                                   width = 8, height = 4)
                            
                            rm(list = "ggp")
                        }
                        
                        #################################################
                        
                        
                        rm(list = ls(pattern = "^lnc_.*tmp"))
                        
                        
                        #################################################
                        
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



my_classes_full <- c(log2_TPM_Q5,
                     log2_TPM_Q35,
                     my_DBTMEE_list)


#################################################  

for(log_type in log_types){
        
        for(my_assay in rev(my_assays)){
                
                if(my_assay == "all"){next()}
                
                
                lnc_tmp <- as.data.frame(get(paste("lnc", my_assay, log_type, sep = ".")))
                
                my_conditions <- SampleTable$Conditions[match(colnames(lnc_tmp), SampleTable$ForeignID)]
                my_conditions <- factor(my_conditions, levels = unique(my_conditions))
                
                
                #################################################           
                
                
                if(my_scaling %in% c("mean")){
                        
                        lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp, na.rm = TRUE)
                        
                } else if(my_scaling %in% c("len")){
                        
                        my_feature_width_ordered <- my_feature_width[match(rownames(lnc_tmp), names(my_feature_width))]
                        
                        stopifnot(identical(rownames(lnc_tmp), names(my_feature_width_ordered)))
                        
                        lnc_scaled_tmp <- lnc_tmp - log2(my_feature_width_ordered)
                        
                }  else if(my_scaling %in% c("zy")){
                        
                        lnc_scaled_tmp <- lnc_tmp - rowMeans(lnc_tmp[, grep("Zy", colnames(lnc_tmp))], na.rm = TRUE)
                }
                
                #################################################  
                
                
                log2_TPM_merged_tmp <- merge(log2_scaled_TPM, lnc_scaled_tmp, by = "row.names")
                
                
                my_corr <- cor(log2_TPM_merged_tmp[,-1], method = "spearman", use = "complete.obs")
                my_corr <- my_corr[grepl(my_assay, rownames(my_corr)),]
                my_corr <- my_corr[,!grepl(my_assay, colnames(my_corr))]
                
                my_corr[is.na(my_corr)] <- 0
                
                
                my_corr <- my_corr[,!grepl("ES|DRB", colnames(my_corr))]
                my_corr <- my_corr[!grepl("ES|DRB", rownames(my_corr)),]
                
                #################################################    
                
                
                pdf(paste0(output_plots,"/corr2RNAseq.", my_assay,".", log_type,".reps.pdf"), width = 6, height = 8, useDingbats = FALSE)
                par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                
                hm_colors <- colorRampPalette(brewer.pal(9,name = "YlGnBu"))(100)
                
                pheatmap(my_corr,
                         main = "RNA-seq\nSpearman's R",
                         color = hm_colors,
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         cellheight = ifelse(my_assay == "all", 15, 25), 
                         cellwidth = ifelse(my_assay == "all", 15, 25))
                
                
                dev.off()
                
                
                #################################################           
                
                
                my_corr_long <- melt(my_corr)
                
                my_corr_long$Conditions <- gsub("_[1-3]","",my_corr_long$Var1)
                my_corr_long$Conditions <- factor(my_corr_long$Conditions, levels = unique(my_corr_long$Conditions))
                
                ggp <- ggbarplot(data = my_corr_long, 
                                 x="Var2", y="value", 
                                 add = c("mean_se","jitter"),
                                 fill = "#666666",
                                 facet.by  = "Conditions",
                                 palette = my_color_palette[-3]) + 
                        theme_bw() +
                        theme(text = element_text(size = 12), aspect.ratio = 1.25, 
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                              axis.ticks.x = element_blank(),
                              panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              strip.background = element_blank(),
                              plot.title = element_text(hjust = 0.5, face = "bold")) +
                        scale_y_continuous(breaks = seq(0,0.8,0.2)) +
                        geom_segment(data = data.frame(x=-Inf,xend=-Inf, y=0,yend=0.8), aes(x=x,xend=xend,y=y,yend=yend), inherit.aes = FALSE) +
                        #ggtitle(label = my_assay) +
                        ylab(label = "Spearman's R\n") +
                        xlab(label = "RNA-seq Stages") 
                
                ggsave(filename = paste0(output_plots,"/corr2RNAseq.", my_assay,".", log_type,".bar.pdf"),
                       plot =  ggp,
                       width = 8, height = 4)
                
                
                ################################################# 
                
                
                rm(list = "my_corr")
                
                
                ################################################# 
                
                lnc_scaled_tmp <- lnc_scaled_tmp[,grep("Zy|2C_S[2,5]P_[1-3]|2C_IgG_[1-3]|8C", colnames(lnc_scaled_tmp))]
                
                
                vio_colors <- my_color_palette[-3]
                
                
                
                rnaStage <- colnames(log2_TPM_Quantiles)[5]
                
                
                for(rnaStage in  c("Q.Zygote","Q.Late_2C","Q.8C")){
                    
                    
                    lnc_merged_Qs_tmp <- merge(log2_TPM_Quantiles[, rnaStage, drop=FALSE], 
                                               lnc_scaled_tmp, by = "row.names")
                    
                    lnc_merged_Qs_tmp_long <- melt(lnc_merged_Qs_tmp)
                    
                    lnc_merged_Qs_tmp_long$Stage <- gsub("_S5P|_IgG|_S2P|_[1-3]$","",lnc_merged_Qs_tmp_long$variable)
                    lnc_merged_Qs_tmp_long$Stage <- factor(lnc_merged_Qs_tmp_long$Stage, levels = unique(lnc_merged_Qs_tmp_long$Stage))
                    
                    
                    ggp <- ggplot(lnc_merged_Qs_tmp_long, 
                                  aes_string(x= rnaStage, y="value", fill = "Stage")) + 
                        theme_bw() +
                        theme(text = element_text(size = 12), aspect.ratio = 1, 
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              strip.background = element_blank(),
                              plot.title = element_text(hjust = 0.5, face = "bold")) +
                        geom_violin(width = 1) +
                        facet_wrap(~ variable, nrow = 3, dir = "v") +
                        ggtitle(label = paste0(feature_type," - ",my_assay)) +
                        stat_summary(fun =  "median", geom="point", size=1) +
                        scale_fill_manual(values = vio_colors) +
                        xlab(label = paste("RNA-seq", rnaStage)) +
                        coord_cartesian(ylim = c(-6, 4)) +
                        scale_y_continuous(breaks = seq(-6,4,2)) +
                        ylab(label = paste(my_assay, "\nlog2 Normalized Counts"))
                    
                    ggsave(filename = paste0(output_plots,"/vioplot.RNAseq.",rnaStage,".",my_assay,".", log_type,".reps.pdf"),
                           plot =  ggp,
                           width = 12, height = 6)
                    
                    rm(list = "ggp")
                }
                
                ################################################# 
                
                
                rm(list = "lnc_scaled_tmp")
                rm(list = "log2_TPM_merged_tmp")
                
                
                #################################################           
                
                
                if(identical(gsub(paste0(my_assay,"|_[1-3]|_"),"", colnames(lnc_tmp)),
                             gsub(paste0(my_assay,"|_[1-3]|_"),"", as.character(my_conditions)))){
                        
                        lnc_repmeans_tmp <- t(apply(lnc_tmp, 1, function(x){
                                aggregate(x, by = list(my_conditions), FUN = mean, na.rm = TRUE)$x}))
                        
                        colnames(lnc_repmeans_tmp) <- levels(my_conditions)
                        
                }
                
                
                if(my_scaling %in% c("mean")){
                        
                        lnc_scaled_tmp <- lnc_repmeans_tmp - rowMeans(lnc_repmeans_tmp, na.rm = TRUE)
                        
                } else if(my_scaling %in% c("len")){
                        
                        my_feature_width_ordered <- my_feature_width[match(rownames(lnc_repmeans_tmp), names(my_feature_width))]
                        
                        stopifnot(identical(rownames(lnc_repmeans_tmp), names(my_feature_width_ordered)))
                        
                        lnc_scaled_tmp <- lnc_repmeans_tmp - log2(my_feature_width_ordered)
                        
                }  else if(my_scaling %in% c("zy")){
                        
                        lnc_scaled_tmp <- lnc_repmeans_tmp - rowMeans(lnc_repmeans_tmp[, grep("Zy", colnames(lnc_repmeans_tmp))], na.rm = TRUE)
                }
                
                
                #################################################
                
                
                log2_TPM_merged_tmp <- merge(log2_scaled_TPM, lnc_scaled_tmp, by = "row.names")
                
                
                my_corr <- cor(log2_TPM_merged_tmp[,-1], method = "spearman", use = "complete.obs")
                my_corr <- my_corr[grepl(my_assay, rownames(my_corr)),]
                my_corr <- my_corr[,!grepl(my_assay, colnames(my_corr))]
                
                my_corr[is.na(my_corr)] <- 0
                
                my_corr <- my_corr[,!grepl("ES|DRB", colnames(my_corr))]
                my_corr <- my_corr[!grepl("ES|DRB", rownames(my_corr)),]
                
                
                #################################################    
                
                pdf(paste0(output_plots,"/corr2RNAseq.", my_assay,".", log_type,".pool.pdf"), width = 6, height = 8, useDingbats = FALSE)
                par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                
                
                pheatmap(my_corr,
                         main = "RNA-seq\nSpearman's R",
                         color = hm_colors,
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         cellheight = ifelse(my_assay == "all", 15, 25), 
                         cellwidth = ifelse(my_assay == "all", 15, 25))
                
                
                rm(list = "my_corr")

                
                dev.off()
                
                
                ################################################# 
                
                lnc_scaled_tmp <- lnc_scaled_tmp[,grep("Zy|2C_S|2C_I|8C", colnames(lnc_scaled_tmp))]
                
                
                vio_colors <- my_color_palette[-3]
                
                
                
                rnaStage <- colnames(log2_TPM_Quantiles)[5]
                
                
                for(rnaStage in  c("Q.Zygote","Q.Late_2C","Q.8C")){
                    
                    
                    lnc_merged_Qs_tmp <- merge(log2_TPM_Quantiles[, rnaStage, drop=FALSE], 
                                               lnc_scaled_tmp, by = "row.names")
                    
                    lnc_merged_Qs_tmp_long <- melt(lnc_merged_Qs_tmp)
                    
                    lnc_merged_Qs_tmp_long$Stage <- gsub("_S5P|_IgG|_S2P|_[1-3]$","",lnc_merged_Qs_tmp_long$variable)
                    lnc_merged_Qs_tmp_long$Stage <- factor(lnc_merged_Qs_tmp_long$Stage, levels = unique(lnc_merged_Qs_tmp_long$Stage))
                    
                    
                    ggp <- ggplot(lnc_merged_Qs_tmp_long, 
                                  aes_string(x= rnaStage, y="value", fill = "Stage")) + 
                        theme_bw() +
                        theme(text = element_text(size = 12), aspect.ratio = 1, 
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              strip.background = element_blank(),
                              plot.title = element_text(hjust = 0.5, face = "bold")) +
                        geom_violin(width = 1) +
                        facet_wrap(~ variable, nrow = 1, dir = "v") +
                        ggtitle(label = paste0(feature_type," - ",my_assay)) +
                        stat_summary(fun =  "median", geom="point", size=1) +
                        scale_fill_manual(values = vio_colors) +
                        xlab(label = paste("RNA-seq", rnaStage)) +
                        coord_cartesian(ylim = c(-6, 4)) +
                        scale_y_continuous(breaks = seq(-6,4,2)) +
                        ylab(label = paste(my_assay, "\nlog2 Normalized Counts"))
                    
                    ggsave(filename = paste0(output_plots,"/vioplot.RNAseq.",rnaStage,".",my_assay,".", log_type,".mean.pdf"),
                           plot =  ggp,
                           width = 8, height = 4)
                    
                    rm(list = "ggp")
                }
                
                
                ################################################# 
                
                
                lnc_scaled_tmp <- lnc_scaled_tmp[,grep("Zy|2C_S|2C_I|8C", colnames(lnc_scaled_tmp))]
                
                
                if(sum(colnames(log2_TPM) %in% c("Zygote","Late_2C","8C")) != 3){next()}
                
                
                vio_colors <- my_color_palette[-3]
                
                
                lnc_Qs_tmp <- as.data.frame(na.omit(lnc_scaled_tmp))
                
                myStage <- colnames(lnc_Qs_tmp)[1]
                
                if(feature_type %in% c("TSS","GB5")){
                        my_qs_incr <- 0.25
                } else {
                        my_qs_incr <- 0.2
                }
                
                for(myStage in  colnames(lnc_Qs_tmp)){
                    
                    Qstage <- paste0("Q.", myStage)
                    
                    lnc_Qs_tmp[, Qstage] <- c("q1","q2","q3","q4","q5")[cut(x = lnc_Qs_tmp[, myStage], 
                                                                            breaks = quantile(lnc_Qs_tmp[, myStage], seq(0,1,my_qs_incr)), 
                                                                            include.lowest = TRUE)]
                }
                
                lnc_Qs_tmp <- lnc_Qs_tmp[, grep("^Q.", colnames(lnc_Qs_tmp)), drop = FALSE]
                
                myStage <-  colnames(lnc_Qs_tmp)[2]
                
                for(myStage in  colnames(lnc_Qs_tmp)){
                    
                    
                    lnc_merged_Qs_tmp <- merge(lnc_Qs_tmp[, myStage, drop=FALSE], 
                                               log2_TPM[,c("Zygote","Late_2C","8C")], by = "row.names")
                    
                    lnc_merged_Qs_tmp_long <- melt(lnc_merged_Qs_tmp)
                    
                    lnc_merged_Qs_tmp_long$Stage <- factor(lnc_merged_Qs_tmp_long$variable, levels = unique(lnc_merged_Qs_tmp_long$variable))
                    
                    
                    ggp <- ggplot(lnc_merged_Qs_tmp_long, 
                                  aes_string(x= myStage, y="value", fill = "Stage")) + 
                        theme_bw() +
                        theme(text = element_text(size = 12), aspect.ratio = 1, 
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              strip.background = element_blank(),
                              plot.title = element_text(hjust = 0.5, face = "bold")) +
                        geom_violin(width = 1) +
                        facet_wrap(~ variable, nrow = 1, dir = "v") +
                        ggtitle(label = paste0(feature_type," - ",my_assay)) +
                        stat_summary(fun =  "median", geom="point", size=1) +
                        scale_fill_manual(values = vio_colors) +
                        xlab(label = gsub("_",": ",myStage)) +
                        coord_cartesian(ylim = c(-15, 15)) +
                        scale_y_continuous(breaks = seq(-15,15,5)) +
                        ylab(label = "log2 TPM")
                    
                    ggsave(filename = paste0(output_plots,"/vioplot.",my_assay,".", log_type, ".",myStage,".RNAseq.TPM.pdf"),
                           plot =  ggp,
                           width = 8, height = 4)
                    
                    rm(list = "ggp")
                }
                
                #################################################                
                
                
                
                #################################################           
                
                
                # my_xlims <- c(min(log2_TPM_merged_tmp[,grep(my_assay, colnames(log2_TPM_merged_tmp))], na.rm = TRUE),
                #               max(log2_TPM_merged_tmp[,grep(my_assay, colnames(log2_TPM_merged_tmp))], na.rm = TRUE))
                # 
                # 
                # my_ylims <- c(min(log2_TPM_merged_tmp[,c(-1,-grep(my_assay, colnames(log2_TPM_merged_tmp)))], na.rm = TRUE),
                #               max(log2_TPM_merged_tmp[,c(-1,-grep(my_assay, colnames(log2_TPM_merged_tmp)))], na.rm = TRUE))
                # 
                # 
                # my_class <- "4-Cell_Transient"
                # 
                # # for(my_NA in c("pairwise","complete")){
                #         
                #         if(my_NA == "complete"){
                #                 log2_TPM_merged_tmp <- na.omit(log2_TPM_merged_tmp)
                #         }
                #         
                #         for(my_class in names(my_classes_full)){
                #                 
                #                 
                #                 png(paste0(output_plots,"/scatter2RNAseq.", my_assay,".",my_class,".",my_NA,".",log_type,".png"), width = 25, height = 25, units = "in", res = 200)
                #                 
                #                 par(mfrow=c(9,9), oma=c(5,5,5,5), mar=c(4,4,2,2), mgp=c(2.5,1,0))
                #                 
                #                 for(rid in paste(rep(c("Zy","2C","8C"), each=3),my_assay, 1:3, sep="_")){
                #                         
                #                         for(cid in c("MII_oocyte","Zygote","Early_2C","Mid_2C","Late_2C","4C","8C","ICM","ESC")){
                #                                 
                #                                 my_complete_cases <- complete.cases(log2_TPM_merged_tmp[,c(rid,cid)])
                #                                 
                #                                 plot(x = log2_TPM_merged_tmp[,rid],
                #                                      y = log2_TPM_merged_tmp[,cid],
                #                                      #main = paste("Spearman´s R =", round(my_cor,2)),
                #                                      xlab = rid, ylab = cid,
                #                                      xlim = my_xlims, 
                #                                      ylim = my_ylims,
                #                                      pch = 19, col = rgb(0.8,0.8,0.8,0.5), cex=0.5)
                #                                 
                #                                 my_subset <- log2_TPM_merged_tmp$Row.names %in% my_classes_full[[my_class]]
                #                                 
                #                                 points(x = log2_TPM_merged_tmp[my_subset,rid],
                #                                        y = log2_TPM_merged_tmp[my_subset,cid],
                #                                        pch = 19, col = rgb(0.8,0,1,0.25), cex=0.5)   
                #                                 
                #                                 lines(lowess(x = log2_TPM_merged_tmp[my_complete_cases,rid],
                #                                              y = log2_TPM_merged_tmp[my_complete_cases,cid],
                #                                              f = 1),
                #                                       col = rgb(0.5,0.5,0.5,1), lwd=2)
                #                                 
                #                                 lines(lowess(x = log2_TPM_merged_tmp[(my_subset & my_complete_cases),rid],
                #                                              y = log2_TPM_merged_tmp[(my_subset & my_complete_cases),cid],
                #                                              f = 1),
                #                                       col = rgb(0.6,0,0.8,1), lwd=2)
                #                                 
                #                                 points(x = mean(log2_TPM_merged_tmp[my_complete_cases,rid]),
                #                                        y = mean(log2_TPM_merged_tmp[my_complete_cases,cid]),
                #                                        pch = 19, col = rgb(0.5,0.5,0.5,1), cex=1.25) 
                #                                 
                #                                 points(x = mean(log2_TPM_merged_tmp[(my_subset & my_complete_cases),rid]),
                #                                        y = mean(log2_TPM_merged_tmp[(my_subset & my_complete_cases),cid]),
                #                                        pch = 19, col = rgb(0.6,0,0.8,1), cex=1.25) 
                #                                 
                #                                 
                #                                 my_cor <- cor(x = log2_TPM_merged_tmp[,rid],
                #                                               y = log2_TPM_merged_tmp[,cid], 
                #                                               method = "spearman", use = "complete.obs")
                #                                 
                #                                 my_cor_sub <- cor(x = log2_TPM_merged_tmp[my_subset,rid],
                #                                                   y = log2_TPM_merged_tmp[my_subset,cid], 
                #                                                   method = "spearman", use = "complete.obs")
                #                                 
                #                                 legend("topleft", horiz = FALSE, 
                #                                        legend = c(paste0("Rs = ", round(my_cor,2),
                #                                                          "; n = ", sum(my_complete_cases)),
                #                                                   paste0("Rs = ", round(my_cor_sub,2),
                #                                                          "; n = ", sum(my_complete_cases & my_subset))),
                #                                        fill = c(rgb(0.7,0.7,0.7,0.5), rgb(0.7,0,0.9,0.75)), cex =0.7)
                #                                 
                #                                 rm(list = "my_cor")
                #                                 
                #                         }
                #                 }
                #                 
                #                 mtext(text = "RNA-seq log2 TPM", side = 2, line = 1, font = 2,  outer = TRUE, cex = 1.2)
                #                 mtext(text = "log2 Normalized Counts", side = 1, line = 1, font = 2,  outer = TRUE, cex = 1.2)
                #                 
                #                 par(fig=c(0,1,0,1), mar=c(1,1,1,1), oma=c(0,0,0,0), new=TRUE)        
                #                 plot.new()
                #                 
                #                 legend("top", horiz = TRUE, 
                #                        legend = c(paste0("all genes - ", length(my_subset)),
                #                                   paste0(my_class, " - ", sum(my_subset))),
                #                        fill = c(rgb(0.7,0.7,0.7,0.5), rgb(0.7,0,0.9,0.75)), cex = 1.5)
                #                 
                #                 dev.off()
                #                 
                #         }
                # }
                
                
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


