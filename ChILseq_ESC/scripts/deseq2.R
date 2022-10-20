






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

output_sessionInfo <- grep("sessionInfo", args, value = TRUE)
output_results <- gsub("\\/sessionInfo.*","",output_sessionInfo)

feature_type <- gsub("Output/deseq2_|/results", "", output_results)




##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################



load(paste0("genome/ranges_",feature_type,".rda"))


if(feature_type == "intergenic"){
        
        my_intergenic_ranges$gene_id <- names(my_intergenic_ranges)
        my_intergenic_ranges$gene_name <- names(my_intergenic_ranges)
}


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

#################################################        count table         #####################################################



my_count_files <-  grep(paste0("counts/",feature_type), args, value = TRUE)


for(i in seq_along(my_count_files)){
        
        load(my_count_files[i])
        
        stopifnot(identical(names(get(gsub(".*\\/|.rda","",my_count_files[i]))),
                            names(get(gsub(".*\\/|.rda","",my_count_files[1])))))
}


my_countss <- ls(pattern = paste0("^",feature_type))


my_count_table <- lapply(my_countss, get)
my_count_table <- as.data.frame(Reduce(cbind, my_count_table))

colnames(my_count_table) <- gsub(paste0(feature_type,"."),"",my_countss)



rm(list = ls(pattern = paste0("^",feature_type,"\\.")))



##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################

#################################################        bin table         #####################################################



my_bin_files <-  grep("counts/bins", args, value = TRUE)


for(i in seq_along(my_bin_files)){
        
        load(my_bin_files[i])
        
        stopifnot(identical(names(get(gsub(".*\\/|.rda","",my_bin_files[i]))),
                            names(get(gsub(".*\\/|.rda","",my_bin_files[1])))))
}


my_binss <- ls(pattern = "^bins")


my_bin_table <- lapply(my_binss, get)
my_bin_table <- as.data.frame(Reduce(cbind, my_bin_table))

colnames(my_bin_table) <- gsub("bins.","",my_binss)



rm(list = ls(pattern = "^bins\\."))



##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################           Setup            #####################################################



my_assay <- "all"



for(my_assay in c("all", unique(SampleTable$Assay)[-4])){
        
        
        #########################################################
        
        my_count_tmp <- my_count_table
        
        if(my_assay == "all"){
                my_subset <- rep(TRUE, ncol(my_count_table))
        } else {
                my_subset <- colnames(my_count_tmp) %in% SampleTable$ForeignID[SampleTable$Assay == my_assay]
        }
        
        
        if(sum(my_subset) <= 2){next()}
        
        my_count_sub_tmp <- my_count_tmp[,my_subset]
        SampleTable_tmp  <- SampleTable[my_subset,]
        
        #########################################################
        
        my_bin_tmp <- my_bin_table
        
        if(identical(colnames(my_count_tmp), colnames(my_bin_tmp))){
                
                my_bin_sub_tmp <- my_bin_tmp[,my_subset]
        }
        
        #########################################################
        
        if(all(identical(colnames(my_count_sub_tmp), SampleTable_tmp$ForeignID),
               identical(colnames(my_count_sub_tmp), colnames(my_bin_sub_tmp)))){
                
                dds_tmp <- setupDDS(CountTableName = "my_count_sub_tmp",
                                    SampleTableName = "SampleTable_tmp",
                                    SampleIdName = "ForeignID",
                                    ConditionName = ifelse(my_assay == "all", "Conditions", "Stages"),
                                    BatchName = NULL,
                                    n_samples_for_filtering = ncol(my_count_sub_tmp)*0,
                                    min_number_of_reads = -1)
                
                sizeFactors(dds_tmp) <- estimateSizeFactorsForMatrix(my_bin_sub_tmp)
                
                dds_tmp <- DESeq(dds_tmp)
                
                assign(paste("dds", my_assay, sep = "."), dds_tmp)
                
                save(list =  paste("dds", my_assay, sep = "."),
                     file = paste0(output_results,"/",paste("dds", my_assay, sep = ".") , ".rda"))
                
        }
        
        #########################################################
        
        rm(list = ls(pattern = "_tmp"))
}


##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

#################################################           contrasts         #####################################################



my_assays <- gsub("dds\\.","",ls(pattern = "^dds.[S,I][5,2,g]"))


for(my_assay in my_assays){
        
        
        dds_tmp <- get(paste("dds", my_assay, sep = "."))
        
        contrast_list <- list(c("ChIL100C","ChIP"),
                              c("ChIL1000C", "ChIP"),
                              c("ChIL1000C","ChIL100C"))        
        
        
        
        for(contrast_name in contrast_list){
                
                res <- getResults(dds = dds_tmp,
                                  contrast = contrast_name,
                                  lfc_cutoff = 0,
                                  shrink = FALSE,
                                  annotation = paste0("my_",feature_type,"_ranges"),
                                  anno_symbol = "gene_name",
                                  anno_id = "gene_id")
                
                res_name <- paste("res", my_assay,
                                  paste0(contrast_name[1],"-",contrast_name[2]), sep = ".")
                
                assign(res_name, res)
                
                write.table(res, file = paste0(output_results,"/",res_name, ".txt"),
                            quote = F, sep = "\t", row.names = T, col.names = NA)
                
                rm(list = "res")
        }
        
        rm(list = ls(pattern = "_tmp"))
}




##################################################################################################################################
##################################################################################################################################

















##################################################################################################################################
##################################################################################################################################

#################################################             rld            #####################################################



my_assays <- gsub("dds\\.","",ls(pattern = "^dds\\."))

log_type <- "log2"


for(my_assay in my_assays){
        
        for(log_type in c("rlog", "log2")){
                
                dds_tmp <- get(paste("dds", my_assay, sep = "."))
                
                # modcombat <- model.matrix(~Sample, data = colData(dds_tmp))
                # batchVar <- colData(dds_tmp)$Batch
                
                #########################################################
                
                if(log_type == "rlog"){
                        
                        rld_tmp <- rlog(dds_tmp, blind = FALSE)
                        lnc_tmp <- assay(rld_tmp)
                        
                        # lnc_tmp <- ComBat(dat = lnc_tmp,
                        #                   batch = batchVar, mod = modcombat,
                        #                   par.prior = TRUE, prior.plots = FALSE)
                        
                        assign(paste("lnc", my_assay, log_type, sep = "."), lnc_tmp)
                        
                } else if(log_type == "log2"){
                        
                        nc_tmp <- counts(dds_tmp, normalized = TRUE)
                        
                        #lnc_tmp <- log2(nc_tmp+1)
                        
                        nc_tmp[round(nc_tmp, 6) == 0] <- NA      
                        lnc_tmp <- log2(nc_tmp)
                        
                        
                        # lnc_tmp <- ComBat(dat = lnc_tmp,
                        #                   batch = batchVar, mod = modcombat,
                        #                   par.prior = TRUE, prior.plots = FALSE)
                        
                        assign(paste("lnc", my_assay, log_type, sep = "."), lnc_tmp)
                        
                }
                
                write.table(lnc_tmp, file = paste0(output_results,"/lnc.", my_assay,".",log_type,".txt"),
                            quote = F, sep = "\t", row.names = T, col.names = NA)
                
                #########################################################
                
                rm(list = ls(pattern = "_tmp"))
        }
}


##################################################################################################################################
##################################################################################################################################


















##################################################################################################################################
##################################################################################################################################







writeLines(capture.output(sessionInfo()), output_sessionInfo)






##################################################################################################################################
##################################################################################################################################


