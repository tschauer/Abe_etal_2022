


##################################################################################################################################
##################################################################################################################################




args = commandArgs(trailingOnly=TRUE)


output_file <- grep("read_summary.txt", args, value = TRUE)


input_files <- args[!(grepl(output_file, args))]


##################################################################################################################################
##################################################################################################################################



my_ids <- unique(gsub(".*\\/|_[t,m,r,f].*","",input_files))


my_df <- data.frame(id = my_ids,
                    raw = NA,
                    trimmed = NA,
                    mapped = NA,
                    mapped2 = NA,
                    filtered = NA)

for(i in 1:nrow(my_df)){
    
    my_df$raw[i] <-      read.table(grep(my_df$id[i], grep("raw.txt",     input_files, value = TRUE), value = TRUE))$V1
    my_df$trimmed[i] <-  read.table(grep(my_df$id[i], grep("trimmed.txt", input_files, value = TRUE), value = TRUE))$V1
    my_df$mapped[i] <-   read.table(grep(my_df$id[i], grep("mapped.txt",  input_files, value = TRUE), value = TRUE))$V1
    my_df$mapped2[i] <-  read.table(grep(my_df$id[i], grep("mapped2.txt", input_files, value = TRUE), value = TRUE))$V1
    my_df$filtered[i] <- read.table(grep(my_df$id[i], grep("filtered.txt",input_files, value = TRUE), value = TRUE))$V1
    
}

#stopifnot(identical(my_df$mapped, my_df$mapped2))


#my_df <- my_df[, !(grepl("mapped2", colnames(my_df)))]



##################################################################################################################################
##################################################################################################################################



write.table(my_df, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################


