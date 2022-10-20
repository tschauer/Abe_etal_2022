






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



matrix_files <- grep("S[5]P.*_TSS_matrix.rds", args, value = TRUE)
matrix_files <- matrix_files[order(matrix_files)]




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

######################################################     Read Mats    ########################################################## 


i=1

for(i in seq_along(matrix_files)){
    
    my_name <- gsub(".*matrix\\/|.*matrix_ave/|_matrix.*","",matrix_files[i])
    my_site <- gsub(".*_","",my_name)
    my_name <- paste0(my_site,".",gsub(paste0("_",my_site),"", my_name))
    
    
    my_mat <- readRDS(matrix_files[i])
    
    assign(my_name, my_mat)
    
    rm(list = "my_mat")
}


my_mats <- ls(pattern = paste0("^TSS."))


for(i in seq_along(my_mats)){
    
    stopifnot(identical(rownames(get(my_mats[1])), rownames(get(my_mats[i]))))
    stopifnot(sum(grepl("ENS", rownames(get(my_mats[i])))) == nrow(get(my_mats[i])))
}




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






##################################################################################################################################
##################################################################################################################################




my_name <- my_mats[1]


for(my_name in my_mats){
    
    my_mat <- get(my_name)
    
    my_stage <- gsub(paste0("TSS.|_S[5,2]P|_[1-3]"),"", my_name)
    my_name <- paste0(my_name, ".Q5")
    
    RNAqid <- paste0("Q.",SampleMatchTable$RNAseq[SampleMatchTable$Stages == my_stage])
    
    my_mat_sub <- my_mat[rownames(my_mat) %in% log2_TPM$gene_id[log2_TPM[,RNAqid] == "q5"],]
    
    assign(my_name, my_mat_sub)
    
    rm(list = "my_mat_sub")
    rm(list = "my_name")
    rm(list = "my_mat")       
}




##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################



my_mats_q <- ls(pattern = "TSS.*Q5")

my_mats_q <- my_mats_q[c(grep("Zy", my_mats_q), grep("2C|8C|ES",my_mats_q))]

df <- data.frame(Object = my_mats_q,
                 Samples = gsub("TSS.|.Q5", "", my_mats_q),
                 Conditions = gsub("TSS.|_S5P|.Q5|_[1-3]", "", my_mats_q), 
                 stringsAsFactors = FALSE)

df$Conditions <- factor(df$Conditions, levels = unique(df$Conditions))

for(i in 1:nrow(df)){
    
    my_mat <- get(df$Object[i])
    
    my_comp <-  colMeans(my_mat, na.rm = TRUE)
    
    my_comp <- (my_comp - min(my_comp, na.rm = TRUE)) / (max(my_comp, na.rm = TRUE) - min(my_comp, na.rm = TRUE))
    
    df$TSS[i] <- mean(my_comp[4800:5100], na.rm=TRUE)
    df$GB[i]  <- mean(my_comp[6500:10000], na.rm=TRUE)
    
    # df$TSS[i] <- mean(my_mat[,4800:5100], na.rm=TRUE)
    # df$GB[i]  <- mean(my_mat[,6500:10000], na.rm=TRUE)
    
    df$nsites[i] <- nrow(my_mat)
    
    rm(list = "my_mat")
    
}


df$Ratio <- df$GB / df$TSS


##################################################################################################################################
##################################################################################################################################






my_color_palette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")






pdf(output_file, height = 9, width = 9, useDingbats = F)
par(mfrow=c(2,2), mar = c(5,10,3,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))


#################################################           


set.seed(123)

plot(x = jitter(as.integer(df$Conditions),factor = 0.1),
     y = df$Ratio ,
     main =  "S5P at Q5 genes",
     xlim = c(min(as.integer(df$Conditions))-0.5, max(as.integer(df$Conditions))+0.5),
     ylim = c(0, ceiling(max(df$Ratio)*10+1)/10),
     xlab = "", xaxt = "n",
     ylab = "Scaled Coverage Ratio (Genebody / TSS)",
     pch = 19, col = my_color_palette[df$Conditions])

set.seed(123)

points(x = jitter(as.integer(df$Conditions), factor = 0.1),
       y =  df$Ratio,
       col = "#555555", pch=1, lwd=0.5)

axis(side = 1, at = seq_along(levels(df$Conditions)), 
     labels = levels(df$Conditions), las = 2)


#################################################           

par(mar = c(5,2,3,10))

plot.new()

legend("left",
       legend = levels(df$Conditions),
       cex = 1, bty = "n",
       fill =  my_color_palette[seq_along(levels(df$Conditions))])


#################################################           



fit <- lm(Ratio ~ Conditions, df)

set.seed(2)


glht_out <- summary(glht(fit, linfct = mcp(Conditions = c("`Zy` - `2C` = 0",
                                                          "`2C_DRB` - `2C` = 0",
                                                          "`8C` - `2C` = 0",
                                                          "`ES` - `2C` = 0"))))

stats_tmp <- data.frame(`Estimate` =   round(glht_out$test$coefficients*100, 2)/100,
                        `Std. Error` = round(glht_out$test$sigma*100, 2)/100,
                        `t value` =  round(glht_out$test$tstat*100, 2)/100,
                        `Pr(>|t|)` = round(glht_out$test$pvalues*100, 2)/100, 
                        check.names = FALSE)


textplot(stats_tmp, halign = "left", valign = "top", cex = 0.6)
text(x = 0.5, y = 0.75, cex = 0.8,
     labels = "Multiple Comparisons of Means: User-defined Contrasts")


#################################################        

dev.off()





##################################################################################################################################
##################################################################################################################################




ggp <- ggbarplot(data = df, 
                 x="Conditions", y="Ratio", 
                 add = c("mean_se","jitter"),
                 fill = "Conditions",
                 palette = my_color_palette) + 
    theme_bw() +
    theme(text = element_text(size = 12), aspect.ratio = 1.5, 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_y_continuous(breaks = seq(0,0.8,0.2)) +
    geom_segment(data = data.frame(x=-Inf,xend=-Inf, y=0,yend=0.8), aes(x=x,xend=xend,y=y,yend=yend), inherit.aes = FALSE) +
    ggtitle(label = "S5P at Q5 genes") +
    ylab(label = "Scaled Coverage Ratio (Genebody / TSS)\n") +
    xlab(label = "") 

ggsave(filename = gsub(".pdf",".v2.pdf",output_file),
       plot =  ggp,
       width = 5, height = 5)


##################################################################################################################################
##################################################################################################################################


