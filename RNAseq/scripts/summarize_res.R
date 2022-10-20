






##################################################################################################################################
##################################################################################################################################






library(org.Mm.eg.db)

library(pheatmap)
library(RColorBrewer)


##################################################################################################################################
##################################################################################################################################





args = commandArgs(trailingOnly=TRUE)



input_table <- grep("SampleTable", args, value = TRUE)


input_sessionInfo <- grep("sessionInfo", args, value = TRUE)
input_directory <- gsub("\\/sessionInfo.*","/tables/",input_sessionInfo)


output_file <- grep("sum_table.*.txt", args, value = TRUE)

my_direction <- gsub(".*\\.","",gsub(".txt", "", output_file))



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

#################################################        Summary Tables        #####################################################




res_names <- ls(pattern = "^res\\.")

padj_cutoff <- 0.05



my_classes_list <- c(list(All_Genes = rownames(my_res)),
                     my_DBTMEE_list)


test.table <- data.frame(Assay = NA,
                         Comparison = NA,
                         Class = NA,
                         Sign.inClass = NA,
                         NS.inClass = NA,
                         Sign.notClass = NA,
                         NS.notClass = NA,
                         OddsRatio = NA,
                         ConfIntLow = NA,
                         ConfIntHigh = NA)
ridx=1
cidx=1


for(ridx in seq_along(res_names)){
        
        for(cidx in seq_along(my_classes_list)){
                
                
                my_res <- get(res_names[ridx])
                
                my_class <- my_classes_list[[cidx]]
                
                Assay <- "RNA-seq"
                Comparison <- gsub(".*\\.","",res_names[ridx])
                Class <- gsub("class.","", names(my_classes_list)[cidx])
                
                
                if(my_direction == "All"){
                        Sign.inClass <-  nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj <  padj_cutoff & my_res$log2FoldChange > -100,])
                        NS.inClass <-    nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj >= padj_cutoff & my_res$log2FoldChange > -100,])
                        Sign.notClass <- nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj <  padj_cutoff & my_res$log2FoldChange > -100,])
                        NS.notClass <-   nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj >= padj_cutoff & my_res$log2FoldChange > -100,])
                } else if(my_direction == "Up"){
                        Sign.inClass <-  nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj <  padj_cutoff & my_res$log2FoldChange > 0,])
                        NS.inClass <-    nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj >= padj_cutoff & my_res$log2FoldChange > 0,])
                        Sign.notClass <- nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj <  padj_cutoff & my_res$log2FoldChange > 0,])
                        NS.notClass <-   nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj >= padj_cutoff & my_res$log2FoldChange > 0,])
                } else if(my_direction == "Down"){
                        Sign.inClass <-  nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj <  padj_cutoff & my_res$log2FoldChange < 0,])
                        NS.inClass <-    nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj >= padj_cutoff & my_res$log2FoldChange < 0,])
                        Sign.notClass <- nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj <  padj_cutoff & my_res$log2FoldChange < 0,])
                        NS.notClass <-   nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj >= padj_cutoff & my_res$log2FoldChange < 0,])
                }
                
                my_fisher_test <- fisher.test(matrix(c(Sign.inClass,Sign.notClass,NS.inClass,NS.notClass),nrow = 2, byrow = T))
                
                OddsRatio <- round(my_fisher_test$estimate, 4)
                ConfIntLow <- round(my_fisher_test$conf.int[1], 4)
                ConfIntHigh <- round(my_fisher_test$conf.int[2], 4)
                
                test.table <- rbind(test.table, 
                                    c(Assay,Comparison,Class,
                                      Sign.inClass,
                                      NS.inClass,
                                      Sign.notClass,
                                      NS.notClass,
                                      OddsRatio,
                                      ConfIntLow,
                                      ConfIntHigh))
        }
        
        
        
}


test.table <- test.table[complete.cases(test.table),]

write.table(test.table, file = output_file, sep="\t", quote = F, row.names = F)




##################################################################################################################################
##################################################################################################################################


OR <- matrix(as.numeric(test.table$OddsRatio), ncol = length(unique(test.table$Comparison)))[-1,]

colnames(OR) <- unique(test.table$Comparison)
rownames(OR) <- unique(test.table$Class)[-1]

pdf(gsub(".txt",".pdf",output_file), width = 6, height = 6, useDingbats = FALSE)
par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)

pheatmap(log(OR[,c(-2,-3)]+0.000001), 
         main = paste0("logOR [",tolower(my_direction), "-reg genes]"), 
         color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100), 
         breaks = seq(-1.5,1.5, length.out = 101),
         cellheight = 25, cellwidth = 25,
         cluster_rows = FALSE,
         cluster_cols = FALSE)


dev.off()



##################################################################################################################################
##################################################################################################################################






