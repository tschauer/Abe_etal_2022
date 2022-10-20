




args = commandArgs(trailingOnly=TRUE)



##################################################################################################################################
##################################################################################################################################




SampleTable <- read.delim(args[1], stringsAsFactors = FALSE)
SampleTable$Conditions <- factor(gsub("_.*","",SampleTable$ForeignID))




##################################################################################################################################
##################################################################################################################################




TPM_files <- args[c(-1,-length(args),-(length(args)-1))]


TPM <- read.delim(TPM_files[1], stringsAsFactors = FALSE)

TPM_table <- matrix(nrow = nrow(TPM), ncol = length(TPM_files))
rownames(TPM_table) <- TPM$gene_id
colnames(TPM_table) <- gsub(".*\\/|_rsem.*","", TPM_files)


for(i in 1:length(TPM_files)){
        
        TPM <- read.delim(TPM_files[i], stringsAsFactors = FALSE)
        
        if(all(identical(rownames(TPM_table), TPM$gene_id),
               identical(colnames(TPM_table)[i], gsub(".*\\/|_rsem.*","", TPM_files)[i]))){
                
                TPM_table[,i] <- TPM$TPM
        }
}



write.table(TPM_table, file = args[(length(args)-1)], 
            quote = F, sep = "\t", row.names = T, col.names = NA) 




##################################################################################################################################
##################################################################################################################################




if(identical(colnames(TPM_table), SampleTable$SampleID)){
        
        TPM_means <- t(apply(TPM_table, 1, function(x){aggregate(x, by = list(SampleTable$Conditions), FUN = mean)$x}))
        
        colnames(TPM_means) <- levels(SampleTable$Conditions)
}

write.table(TPM_means, file = args[length(args)], 
            quote = F, sep = "\t", row.names = T, col.names = NA) 



##################################################################################################################################
##################################################################################################################################


