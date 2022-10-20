






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


input_ranges_TSS <-      grep("genome/ranges_TSS", args, value = TRUE)
input_ranges_TTS <-      grep("genome/ranges_TTS", args, value = TRUE)
input_ranges_genebody <- grep("genome/ranges_genebody", args, value = TRUE)


input_table <- grep("SampleTable", args, value = TRUE)


output_plots <-  "Output/traveling/"
output_sessionInfo <- grep("Output/traveling/sessionInfo", args, value = TRUE)




##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################




load(input_ranges_TSS)
load(input_ranges_TTS)
load(input_ranges_genebody)






##################################################################################################################################
##################################################################################################################################

#################################################         Gene Width         #####################################################




my_genebody_width <- width(my_genebody_ranges) /1e3
names(my_genebody_width) <- my_genebody_ranges$gene_id


my_TSS_width <- width(my_TSS_ranges) / 1e3
names(my_TSS_width) <- my_TSS_ranges$gene_id


my_TTS_width <- width(my_TTS_ranges) / 1e3
names(my_TTS_width) <- my_TTS_ranges$gene_id


##################################################################################################################################
##################################################################################################################################





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

#################################################        RNA-seq TPM         #####################################################


TPM <- read.table("../../Public/GSE38495_GSE45719/Output/TPM/TPM_means.txt", check.names = FALSE)
TPM <- TPM[,c("MII_oocyte","Zygote", "Early_2C", "Mid_2C", "Late_2C", "4C", "8C", "ICM","ESC")]


log2_TPM <- log2(TPM+1)

#log2_TPM <- t(scale(t(log2_TPM)))
log2_scaled_TPM <- log2_TPM - rowMeans(log2_TPM)


##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################     RNA-seq Quantiles      #####################################################




log2_TPM_Quantiles <- log2_TPM[rownames(log2_TPM) %in% my_genebody_ranges$gene_id, ]

log2_TPM_Q5 <- list()

for(RNAstage in colnames(log2_TPM_Quantiles)){
        
        Qstage <- paste0("Q.", RNAstage)
        
        log2_TPM_Quantiles[,Qstage] <- ifelse(round(log2_TPM_Quantiles[,RNAstage], 4) == 0.0000, "q1", "exp")
        
        is_nonzero <- log2_TPM_Quantiles[,Qstage] == "exp"
        
        log2_TPM_Quantiles[is_nonzero, Qstage] <- c("q2","q3","q4","q5")[cut(x = log2_TPM_Quantiles[is_nonzero, RNAstage], 
                                                                             breaks = quantile(log2_TPM_Quantiles[is_nonzero, RNAstage]), 
                                                                             include.lowest = TRUE)]
        
        log2_TPM_Q5[[paste0("Q5_", RNAstage)]] <- rownames(log2_TPM_Quantiles)[log2_TPM_Quantiles[, Qstage] %in% c("q5")]
        
}


log2_TPM_Quantiles <- log2_TPM_Quantiles[,grep("^Q", colnames(log2_TPM_Quantiles))]



##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################         Settings         #####################################################



my_assay <- "S5P"

log_type <- "log2"



for(my_assay in c("S5P","S2P","IgG")){
        
        for(log_type in c("log2")){
                
                
                
                if(my_assay == "S5P"){
                        
                        my_color_palette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")
                        
                } else {
                        my_color_palette <- c("#999999", "#D55E00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")
                }
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                
                
                
                
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                #################################################         Read Table         #####################################################
                
                
                
                ### old definition based on log2 normalized counts ###
                
                # lnc.TSS.all <- read.table(paste0("Output/deseq2_TSS/results/lnc.all.",log_type,".txt"))
                # lnc.TSS <- lnc.TSS.all[,grep(my_assay, colnames(lnc.TSS.all))]              
                
                
                # lnc.TSS <- read.table(paste0("Output/deseq2_TSS/results//lnc.", my_assay,".",log_type,".txt"))
                # lnc.TSS <- lnc.TSS[order(rownames(lnc.TSS)),]
                # lnc.TSS <- lnc.TSS[, c(grep("Zy", colnames(lnc.TSS)), grep("2C|8C|ES", colnames(lnc.TSS)))]
                
                
                
                # lnc.TTS.all <- read.table(paste0("Output/deseq2_TTS/results/lnc.all.",log_type,".txt"))
                # lnc.TTS <- lnc.TTS.all[,grep(my_assay, colnames(lnc.TTS.all))]        
                
                # lnc.TTS <- read.table(paste0("Output/deseq2_TTS/results/lnc.", my_assay,".",log_type,".txt"))
                # lnc.TTS <- lnc.TTS[order(rownames(lnc.TTS)),]
                # lnc.TTS <- lnc.TTS[, c(grep("Zy", colnames(lnc.TTS)), grep("2C|8C|ES", colnames(lnc.TTS)))]
                
                
                # lnc.genebody.all <- read.table(paste0("Output/deseq2_genebody/results/lnc.all.",log_type,".txt"))
                # lnc.genebody <- lnc.genebody.all[,grep(my_assay, colnames(lnc.genebody.all))]   
                
                
                # lnc.genebody <- read.table(paste0("Output/deseq2_genebody/results/lnc.", my_assay,".",log_type,".txt"))
                # lnc.genebody <- lnc.genebody[order(rownames(lnc.genebody)),]
                # lnc.genebody <- lnc.genebody[, c(grep("Zy", colnames(lnc.genebody)), grep("2C|8C|ES", colnames(lnc.genebody)))]
                
                
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                #################################################          load dds          #####################################################
                
                
                load(paste0("Output/deseq2_genebody/results/dds.",my_assay,".rda"))
                dds.genebody <- get(paste0("dds.", my_assay))
                rm(list = paste0("dds.", my_assay))
                
                sizeF.genebody <- sizeFactors(dds.genebody)
                sizeF.genebody <- sizeF.genebody[c(grep("Zy", names(sizeF.genebody)), grep("2C|8C|ES", names(sizeF.genebody)))]
                
                counts.genebody <- counts(dds.genebody)
                counts.genebody <- counts.genebody[order(rownames(counts.genebody)),]
                counts.genebody <- counts.genebody[, c(grep("Zy", colnames(counts.genebody)), grep("2C|8C|ES", colnames(counts.genebody)))]
                
                
                #################################################          
                
                
                load(paste0("Output/deseq2_TSS/results/dds.",my_assay,".rda"))
                dds.TSS <- get(paste0("dds.", my_assay))
                rm(list = paste0("dds.", my_assay))
                
                sizeF.TSS <- sizeFactors(dds.TSS)
                sizeF.TSS <- sizeF.TSS[c(grep("Zy", names(sizeF.TSS)), grep("2C|8C|ES", names(sizeF.TSS)))]
                
                counts.TSS <- counts(dds.TSS)
                counts.TSS <- counts.TSS[order(rownames(counts.TSS)),]
                counts.TSS <- counts.TSS[, c(grep("Zy", colnames(counts.TSS)), grep("2C|8C|ES", colnames(counts.TSS)))]
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                #################################################       Traveling Ratio      #####################################################
                
                
                
                
                my_common_genes <- intersect(rownames(counts.TSS), rownames(counts.genebody))
                my_common_genes <- intersect(my_common_genes, names(my_genebody_width))
                my_common_genes <- intersect(my_common_genes, names(my_TSS_width))
                
                
                my_genebody_width_tmp <- my_genebody_width[names(my_genebody_width) %in% my_common_genes]
                my_genebody_width_tmp <- my_genebody_width_tmp[order(names(my_genebody_width_tmp))]
                
                
                my_TSS_width_tmp <- my_TSS_width[names(my_TSS_width) %in% my_common_genes]
                my_TSS_width_tmp <- my_TSS_width_tmp[order(names(my_TSS_width_tmp))]
                
                
                
                
                counts.TSS <-      counts.TSS[     rownames(counts.TSS)      %in% my_common_genes,]
                counts.genebody <- counts.genebody[rownames(counts.genebody) %in% my_common_genes,]
                
                
                
                #################################################        
                
                
                
                stopifnot(identical(colnames(counts.TSS), colnames(counts.genebody)))
                stopifnot(identical(rownames(counts.TSS), rownames(counts.genebody)))
                
                
                stopifnot(identical(rownames(counts.TSS),      names(my_TSS_width_tmp)))
                stopifnot(identical(rownames(counts.genebody), names(my_genebody_width_tmp)))
                
                stopifnot(identical(names(my_TSS_width_tmp), names(my_genebody_width_tmp)))
                
                # sizefactors are the same (bin based) so no need to use them
                stopifnot(identical(sizeF.TSS, sizeF.genebody))
                
                stopifnot(identical(colnames(counts.TSS),      names(sizeF.TSS)))
                stopifnot(identical(colnames(counts.genebody), names(sizeF.genebody)))
                
                
                #################################################        
                
                
                ### old definition ###
                # lnc.TSS.wnorm <-      lnc.TSS      - log2(my_TSS_width_tmp)
                # lnc.genebody.wnorm <- lnc.genebody - log2(my_genebody_width_tmp)
                # log2.TR <- lnc.TSS.wnorm - lnc.genebody.wnorm
                
                
                
                # log2.TSS.wnorm      <- log2(counts.TSS+1)      - log2(my_TSS_width_tmp)
                # log2.genebody.wnorm <- log2(counts.genebody+1) - log2(my_genebody_width_tmp)
                
                
                # log2.TSS.wnorm      <- log2(t(t( (counts.TSS      / my_TSS_width_tmp)      ) / sizeF.TSS      ) + 0.1)
                # log2.genebody.wnorm <- log2(t(t( (counts.genebody / my_genebody_width_tmp) ) / sizeF.genebody ) + 0.1)
                
                
                ### final definition ###
                
                counts.TSS[     round(counts.TSS,      6) == 0] <- NA
                counts.genebody[round(counts.genebody, 6) == 0] <- NA
                
                nc.TSS.wnorm      <- (  ( t(t(counts.TSS)      / sizeF.TSS)      / my_TSS_width_tmp ) )
                nc.genebody.wnorm <- (  ( t(t(counts.genebody) / sizeF.genebody) / my_genebody_width_tmp ) )
                
                
                TR <- nc.TSS.wnorm / nc.genebody.wnorm
                
                TR[is.nan(TR)] <- NA
                TR[is.infinite(TR)] <- NA
                TR[round(TR, 6) == 0] <- NA
                
                #TR <- na.omit(TR)
                # TR <- TR[apply(TR, 1, function(x){
                #         all(aggregate(as.numeric(x), by = list(gsub("_[1-3]$","",colnames(TR))), FUN=function(i){sum(is.na(i)) <= 1})$x)
                # }),]
                
                log2.TR <- log2(TR)
                
                
                my_conditions <- factor(SampleTable[,"Stages"][match(gsub("X","", colnames(log2.TR)), SampleTable$ForeignID)])
                my_conditions <- factor(my_conditions, levels = unique(my_conditions))
                
                
                stopifnot(identical(gsub(paste0("_",my_assay,"|_[1-9]|X"),"",colnames(log2.TR)),  as.character(my_conditions)))
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                
                
                
                
                
                
                
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                #################################################         Boxplots           #####################################################
                
                
                
                pdf(paste0(output_plots,"/boxplot.",my_assay,".", log_type,".pdf"), height = 8, width = 8, useDingbats = F)
                par(mfrow=c(2,2), mar = c(6,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))


                for(my_feature in c("TSS","TTS","genebody","GB5")){
                        
                        my_lnc_file <- paste0("Output/deseq2_",my_feature,"/results/lnc.",my_assay,".",log_type,".txt")
                        my_lnc <- read.table(my_lnc_file, check.names=FALSE)
                        my_lnc <- my_lnc[order(rownames(my_lnc)),]
                        
                        my_lnc <- my_lnc[, c(grep("Zy", colnames(my_lnc)), grep("2C|8C|ES", colnames(my_lnc)))] 
                        
                        my_conditions <- SampleTable[,"Stages"][match(colnames(my_lnc), SampleTable$ForeignID)]
                        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
                        
                        boxplot(my_lnc,las=2, cex.axis=0.5, outline=F,
                                main= my_feature, col = my_color_palette[my_conditions])
                        
                        rm(list = ls(pattern = "my_lnc"))
                        
                }
                


                dev.off()

                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                #################################################         Heatmaps           #####################################################
                
                
                
               
                
                # for(my_feature in c("TSS","TTS","genebody","GB5")){
                #         
                #         png(paste0(output_plots,"/heatmap.",my_assay,".",my_feature,".", log_type,".png"), height = 10, width = 6, units = "in", res = 200 )
                # 
                #         my_lnc_file <- paste0("Output/deseq2_",my_feature,"/results/lnc.",my_assay,".",log_type,".txt")
                #         my_lnc <- read.table(my_lnc_file, check.names=FALSE)
                #         my_lnc <- my_lnc[order(rownames(my_lnc)),]
                #         
                #         my_lnc <- my_lnc[, c(grep("Zy", colnames(my_lnc)), grep("2C|8C|ES", colnames(my_lnc)))] 
                #         
                #         my_lnc <- my_lnc[order(rowMeans(my_lnc, na.rm = TRUE), decreasing = TRUE),]
                #         
                #         #my_lnc[is.na(my_lnc)] <- -3
                #         
                #         
                #         pheatmap(my_lnc,
                #                  main = paste0(my_feature," ", my_assay),
                #                  color = (brewer.pal(9,name = "Blues")),
                #                  breaks = seq(-1, 5, length.out = 9),
                #                  na_col = "black",
                #                  cluster_cols = FALSE,
                #                  cluster_rows = FALSE, 
                #                  show_rownames = FALSE)
                #         
                #         rm(list = ls(pattern = "my_lnc"))
                #         
                #         dev.off()
                #         
                # }
                
                
                
                
                
                
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                
                
                
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                #################################################         Dotplots           #####################################################
                
                
                
                
                
                
                
                log2.TR_merged <- merge(log2.TR, my_DBTMEE, by.x = "row.names", by.y = "gene_id", all = FALSE)
                log2.TR_merged$Cluster <- factor(log2.TR_merged$Cluster)
                
                
                log2.TR_merged_median <- t(sapply(levels(log2.TR_merged$Cluster),function(i){
                        
                        colMedians(as.matrix(log2.TR_merged[log2.TR_merged$Cluster == i, grep(my_assay, colnames(log2.TR_merged))]), na.rm = T)
                }))
                
                colnames(log2.TR_merged_median) <- colnames(log2.TR_merged[, grep(my_assay, colnames(log2.TR_merged))])
                
                stopifnot(identical(colnames(log2.TR_merged_median), colnames(log2.TR)))
                
                log2.TR_merged_median <- rbind(colMedians(as.matrix(log2.TR), na.rm = TRUE), log2.TR_merged_median)
                
                rownames(log2.TR_merged_median)[1] <- "all genes"
                
                
                #################################################           
                
                
                log2.TR_merged_median_long <- melt(log2.TR_merged_median, varnames = c("Classes", "Samples"))
                
                log2.TR_merged_median_long$Conditions <- gsub(paste0("_",my_assay,"|_[1-3]"),"", log2.TR_merged_median_long$Samples)
                log2.TR_merged_median_long$Conditions <- factor(log2.TR_merged_median_long$Conditions, levels = unique(log2.TR_merged_median_long$Conditions))
                
                stopifnot(identical(levels(log2.TR_merged_median_long$Conditions), levels(my_conditions)))
                
                #################################################           
                
                
                pdf(paste0(output_plots,"/dotplot.",my_assay,".", log_type,".pdf"), height = 9, width = 9, useDingbats = F)
                
                for(my_class in levels(log2.TR_merged_median_long$Classes)){
                        
                        
                        par(mfrow=c(2,2), mar = c(5,6,3,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))
                        
                        df_tmp <- log2.TR_merged_median_long[log2.TR_merged_median_long$Classes == my_class,]
                        
                        set.seed(123)
                        
                        plot(x = jitter(as.integer(df_tmp$Conditions),factor = 0.1),
                             y = df_tmp$value,
                             # main =  paste0("Traveling Ratio ",my_assay, "\nat ", my_class, 
                             #                " [n=",  ifelse(my_class == "all genes", nrow(log2.TR), 
                             #                                table(log2.TR_merged$Cluster)[my_class]), "]"),
                             main =  paste0("Traveling Ratio ",my_assay, "\nat ", my_class),
                             xlim = c(min(as.integer(df_tmp$Conditions))-0.5, max(as.integer(df_tmp$Conditions))+0.5),
                             ylim = c(min(df_tmp$value)-0.5, max(df_tmp$value)+0.5),
                             xlab = "Conditions",
                             ylab = "Median log2 TSS - Genebody", xaxt = "n",
                             pch = 19, col = my_color_palette[df_tmp$Conditions])
                        
                        set.seed(123)
                        
                        points(x = jitter(as.integer(df_tmp$Conditions), factor = 0.1),
                               y =  df_tmp$value,
                               col = "#555555", pch=1, lwd=0.5)
                        
                        axis(side = 1, at = seq_along(levels(df_tmp$Conditions)), labels = levels(df_tmp$Conditions))
                        
                        
                        #################################################           
                        
                        plot.new()
                        
                        legend("left",
                               legend = levels(df_tmp$Conditions),
                               cex = 1, bty = "n",
                               fill =  my_color_palette[seq_along(levels(df_tmp$Conditions))])
                        
                        
                        #################################################           
                        
                        
                        
                        fit <- lm(value ~ Conditions, df_tmp)
                        
                        set.seed(2)
                        
                        if(my_assay == "S5P"){
                                
                                glht_out <- summary(glht(fit, linfct = mcp(Conditions = c("`Zy` - `2C` = 0",
                                                                                          "`2C_DRB` - `2C` = 0",
                                                                                          "`8C` - `2C` = 0",
                                                                                          "`ES` - `2C` = 0"))))
                        } else {
                                
                                glht_out <- summary(glht(fit, linfct = mcp(Conditions = c("`Zy` - `2C` = 0",
                                                                                          "`8C` - `2C` = 0",
                                                                                          "`ES` - `2C` = 0"))))
                                
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
                
                
                my_subclasses_list <- list(c("Maternal_RNA","1-Cell_Transient", "2-Cell_Transient", "4-Cell_Transient"),
                                           c("Maternal_RNA", "Minor_ZGA",  "Major_ZGA", "MGA"))
                
                
                for(i in seq_along(my_subclasses_list)){
                        
                        
                        pdf(paste0(output_plots,"/dotplot",i,".",my_assay,".", log_type,".pdf"), height = 5, width = 12, useDingbats = F)
                        par(mfrow=c(2,3), mar = c(3,8,3,2), oma = c(3,3,3,3), mgp = c(3,1,0))
                        
                        
                        for(my_condition in levels(log2.TR_merged_median_long$Conditions)){
                                
                                
                                df_tmp <- log2.TR_merged_median_long[log2.TR_merged_median_long$Conditions %in% my_condition,]
                                
                                df_tmp_sub <- df_tmp[df_tmp$Classes %in% c("all genes", my_subclasses_list[[i]]),]
                                df_tmp_sub$Classes <- factor(df_tmp_sub$Classes, levels = rev(c("all genes", my_subclasses_list[[i]])))
                                
                                set.seed(123)
                                
                                plot(y = jitter(as.integer(df_tmp_sub$Classes), factor = 0.1),
                                     x = df_tmp_sub$value,
                                     main =  my_condition,
                                     ylim = c(min(as.integer(df_tmp_sub$Classes))-0.5, max(as.integer(df_tmp_sub$Classes))+0.5),
                                     xlim = c(-0.5, 2),
                                     ylab = "", xlab = "", yaxt = "n",
                                     pch = 19, col = my_color_palette[which(my_condition == levels(df_tmp_sub$Conditions))])
                                
                                set.seed(123)
                                
                                points(y = jitter(as.integer(df_tmp_sub$Classes), factor = 0.1),
                                       x =  df_tmp_sub$value,
                                       col = "#555555", pch=1, lwd=0.5)
                                
                                abline(v = mean(df_tmp_sub$value[df_tmp_sub$Classes == "all genes"]), lty=3)
                                
                                axis(side = 2, at = seq_along(levels(df_tmp_sub$Classes)), labels = levels(df_tmp_sub$Classes), las=2)
                                
                                rm(list = ls(pattern = "df_tmp"))
                                
                        }
                        
                        
                        mtext(text = paste0("Traveling Ratio - ",my_assay), side = 3, line = 0, font = 2,  outer = TRUE, cex = 1.2)
                        mtext(text =  "Median log2 TSS - Genebody", side = 1, line = 0.5, outer = TRUE)
                        
                        dev.off()
                        
                }
                
                
                
                
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                
                pdf(paste0(output_plots,"/dotplot_RNAseq.",my_assay,".", log_type,".pdf"), height = 5.5, width = 12, useDingbats = F)
                par(mfrow=c(2,3), mar = c(3,8,3,2), oma = c(3,3,3,3), mgp = c(3,1,0))
                
                
                for(my_condition in levels(my_conditions)){
                        
                        
                        stopifnot(identical(gsub(paste0("_", my_assay,"|_[1-3]"),"",colnames(log2.TR)), as.character(my_conditions)))
                        
                        df_tmp <- sapply( names(log2_TPM_Q5), function(x){
                                colMedians(as.matrix(log2.TR[rownames(log2.TR) %in% log2_TPM_Q5[[x]], my_condition == my_conditions]), na.rm = TRUE)
                        })
                        
                        
                        df_tmp_long <- melt(df_tmp)
                        df_tmp_long$Var2 <- factor(df_tmp_long$Var2, levels = rev(levels(df_tmp_long$Var2)))
                        
                        set.seed(123)
                        
                        plot(y = jitter(as.integer(df_tmp_long$Var2), factor = 0.1),
                             x = df_tmp_long$value,
                             main =  my_condition,
                             ylim = c(min(as.integer(df_tmp_long$Var2))-0.5, max(as.integer(df_tmp_long$Var2))+0.5),
                             xlim = c(-0.5, 2),
                             ylab = "", xlab = "", yaxt = "n",
                             pch = 19, col = my_color_palette[which(my_condition == levels(my_conditions))])
                        
                        
                        set.seed(123)
                        
                        points(y = jitter(as.integer(df_tmp_long$Var2), factor = 0.1),
                               x =  df_tmp_long$value,
                               col = "#555555", pch=1, lwd=0.5)
                        
                        abline(v = mean(df_tmp_long$value[df_tmp_long$Var2 == "Q5_Zygote"]), lty=3)
                        
                        axis(side = 2, at = seq_along(levels(df_tmp_long$Var2)), labels = levels(df_tmp_long$Var2), las=2)
                        
                        rm(list = ls(pattern = "df_tmp"))
                        
                }
                
                mtext(text = paste0("Traveling Ratio - ",my_assay), side = 3, line = 0, font = 2,  outer = TRUE, cex = 1.2)
                mtext(text =  "Median log2 TSS - Genebody", side = 1, line = 0.5, outer = TRUE)
                mtext(text ="RNA-seq Q5 genes", side = 2, line = 0.5, outer = TRUE)
                
                dev.off()
                
                
                
                
                
                
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                #################################################         Vioplots           #####################################################
                
                
                log2.TR_tmp <- as.data.frame(log2.TR)
                log2.TR_tmp$gene_id <- rownames(log2.TR)
                
                log2.TR_long <- melt(log2.TR_tmp)
                log2.TR_long$Conditions <- factor(gsub(paste0("X|_",my_assay,"|_[1-3]"),"", log2.TR_long$variable))
                log2.TR_long$Conditions <- factor(log2.TR_long$Conditions, levels = unique(log2.TR_long$Conditions))
                
                log2.TR_long$Classes <- my_DBTMEE$Cluster[match(log2.TR_long$gene_id, my_DBTMEE$gene_id)]
                
                
                #################################################           
                
                
                ggf <- ggplot(log2.TR_long, aes(x=variable, y=value, fill = Conditions)) + 
                        theme_bw() +
                        theme(text = element_text(size = 12), aspect.ratio = 0.75, 
                              axis.text.x = element_blank(), 
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              strip.background = element_blank(),
                              plot.title = element_text(hjust = 0.5, face = "bold")) +
                        geom_violin(width = 1) +
                        facet_wrap(~ Classes, nrow = 3) +
                        ggtitle(label = paste0("Traveling Ratio [", my_assay,"]")) +
                        stat_summary(fun=median, geom="point", size=1) +
                        scale_fill_manual(values = my_color_palette) +
                        xlab(label = "Conditions") +
                        coord_cartesian(ylim = c(-6, 6)) +
                        ylab(label = "log2 TSS - Genebody")
                
                
                ggsave(filename = paste0(output_plots,"/vioplot.",my_assay,".", log_type,".pdf"),
                       plot =  ggf,
                       width = 12, height = 8)
                
                
                #################################################           
                
                
                my_subclasses_list <- list(c("Maternal_RNA","1-Cell_Transient", "2-Cell_Transient", "4-Cell_Transient"),
                                           c("Maternal_RNA", "Minor_ZGA",  "Major_ZGA", "MGA"))
                
                
                for(i in seq_along(my_subclasses_list)){
                        
                        log2.TR_long_sub <- log2.TR_long[log2.TR_long$Classes %in% my_subclasses_list[[i]],]
                        
                        log2.TR_long_sub$group <- factor(paste(log2.TR_long_sub$Classes, gsub(".*_","",log2.TR_long_sub$variable)),
                                                         levels = (paste(rep(my_subclasses_list[[i]], each=3), 1:3)))
                        
                        
                        ggp2 <- ggplot(log2.TR_long_sub, 
                                       aes(x=group, y=value, fill = Conditions)) + 
                                theme_bw() +
                                theme(text = element_text(size = 12), aspect.ratio = 1, 
                                      #axis.text.x = element_blank(), 
                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      strip.background = element_blank(),
                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                geom_violin(width = 1) +
                                facet_wrap(~ Conditions, nrow = 1) +
                                ggtitle(label = paste0("Traveling Ratio [", my_assay,"]")) +
                                stat_summary(fun=median, geom="point", size=1) +
                                scale_fill_manual(values = my_color_palette) +
                                xlab(label = "DBTMEE Classes") +
                                coord_cartesian(ylim = c(-4, 4)) +
                                ylab(label = "log2 TSS - Genebody") + 
                                #coord_flip() +
                                geom_hline(yintercept = 0,  color="black", linetype="dashed") +
                                geom_vline(xintercept = c(3.5,6.5,9.5),  color="grey", linetype="dotted") 
                        
                        rm(list = "log2.TR_long_sub")
                        
                        ggsave(filename = paste0(output_plots,"/vioplot",i,".",my_assay,".", log_type,".pdf"),
                               plot =  ggp2,
                               width = 12, height = 6)
                }
                
                #################################################  
                
                
                #################################################           
                
                # log2.TR.repmeans <- as.data.frame(t(apply(log2.TR, 1, function(x){aggregate(x, by= list(my_conditions), mean)[,2]})))
                # colnames(log2.TR.repmeans) <- levels(my_conditions)
                # log2.TR.repmeans$gene_id <- rownames(log2.TR.repmeans)
                # 
                # log2.TR.repmean_long <- melt(log2.TR.repmeans)
                # 
                # log2.TR.repmean_long$Conditions <- factor(gsub(paste0("X|_",my_assay,"|_[1-3]"),"", log2.TR.repmean_long$variable))
                # log2.TR.repmean_long$Conditions <- factor(log2.TR.repmean_long$Conditions, levels = unique(log2.TR.repmean_long$Conditions))
                # log2.TR.repmean_long$Classes <- my_DBTMEE$Cluster[match(log2.TR.repmean_long$gene_id, my_DBTMEE$gene_id)]
                # 
                # 
                # ggf <- ggplot(log2.TR.repmean_long, aes(x=variable, y=value, fill = Conditions)) + 
                #         theme_bw() +
                #         theme(text = element_text(size = 12), aspect.ratio = 0.75, 
                #               axis.text.x = element_blank(), 
                #               panel.border = element_rect(colour = "black", fill=NA, size=1),
                #               panel.grid.major = element_blank(),
                #               panel.grid.minor = element_blank(),
                #               strip.background = element_blank(),
                #               plot.title = element_text(hjust = 0.5, face = "bold")) +
                #         geom_violin(width = 1.5) +
                #         facet_wrap(~ Classes, nrow = 3) +
                #         ggtitle(label = paste0("Traveling Ratio [", my_assay,"]")) +
                #         stat_summary(fun=median, geom="point", size=1) +
                #         scale_fill_manual(values = my_color_palette) +
                #         xlab(label = "Conditions") +
                #         coord_cartesian(ylim = c(-4, 4)) +
                #         ylab(label = "log2 TSS - Genebody")
                # 
                # ggsave(filename = paste0(output_plots,"/vioplot.repmeans.",my_assay,".", log_type,".pdf"),
                #        plot =  ggf,
                #        width = 12, height = 8)
                
                ##################################################################################################################################
                ##################################################################################################################################
                
                
                
                rm(list = ls(pattern = "^lnc.|^log2.TR|tmp"))
                
                
        }
}


##################################################################################################################################
##################################################################################################################################




























##################################################################################################################################
##################################################################################################################################







writeLines(capture.output(sessionInfo()), output_sessionInfo)






##################################################################################################################################
##################################################################################################################################





