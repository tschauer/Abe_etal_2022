






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


input_ranges <- unique(gsub(".*\\/|\\..*", "", grep("counts", args, value = TRUE)))

input_table <- grep("SampleTable", args, value = TRUE)

input_counts <- grep("counts", args, value = TRUE)

output_sessionInfo <- grep("sessionInfo", args, value = TRUE)
output_plots <-  gsub("\\/session.*","", output_sessionInfo)




##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################



for(i in seq_along(input_ranges)){
        
        load(paste0("genome/ranges_", input_ranges[i], ".rda"))
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

#################################################       Features Width         #####################################################



feature_ranges <- gsub("my_|_ranges", "", ls(pattern = "my.*_ranges"))


my_feature_widths <- sapply(feature_ranges, function(x){ sum(width(get(paste0("my_", x,"_ranges"))))/1e3   })

names(my_feature_widths) <- gsub("n$","ns",names(my_feature_widths))



##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################

#################################################        Load Counts          #####################################################



for(i in seq_along(input_counts)){
        
        load(input_counts[i])
        
}


##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################        Sum Counts          #####################################################



feature_ranges <- gsub("n$","ns",feature_ranges)

my_feature_sums <- t(sapply(feature_ranges, function(x){
        sapply(SampleTable$ForeignID, function(i){
                
                sum(get(paste0(x,".", i)))
        })      
}))


if(identical(names(my_feature_widths), rownames(my_feature_sums))){
        
        
        my_feature_rates <- my_feature_sums / my_feature_widths
        
        SF <- colSums(my_feature_rates)
        #SF <- my_feature_rates['intergenic',]
        
        my_feature_rates <- t(t(my_feature_rates) / SF)
}


my_feature_rates <- my_feature_rates[, c(grep("100C", colnames(my_feature_rates)), grep("1000C|ChIP", colnames(my_feature_rates)))]


##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################

#################################################        Barpot Rates        #####################################################





bp_colors <- brewer.pal(nrow(my_feature_rates), "Pastel1")

my_assays <- unique(SampleTable$Assay)

my_assay <- "S5P"



#################################################        



for(my_assay in my_assays){
        
        pdf(paste0(output_plots,"/barplot.",my_assay,".pdf"), height = 9, width = 9, useDingbats = F)
        par(mfrow=c(2,2), mar = c(5,4,3,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))
        
        
        my_rates_tmp <- my_feature_rates[,grep(my_assay, colnames(my_feature_rates))]
        #my_rates_tmp <- my_rates_tmp[, c(grep("Zy", colnames(my_rates_tmp)), grep("2C|8C|ES", colnames(my_rates_tmp)))]
        
        my_order <- rev(c("TSS","exons","introns","TTS", "intergenic"))
        
        bp <- barplot( my_rates_tmp[my_order,], 
                       main = my_assay,
                       ylab = "Fraction of length normalized counts", las= 3, xaxt = "n",
                       col = bp_colors)
        
        axis(side = 1, at = bp, labels =  colnames(my_rates_tmp), 
             cex.axis=0.6, tick = FALSE, line = -0.5, las=2)
        
        
        plot.new()
        legend("center", 
               legend = rev(rownames(my_rates_tmp[my_order,])),
               fill =   rev(bp_colors))
        
        dev.off()
        
        
        rm(list = ls(pattern = "tmp"))
        
}

#################################################        



##################################################################################################################################
##################################################################################################################################




















##################################################################################################################################
##################################################################################################################################

#################################################        Barpot Rates        #####################################################




my_color_palette <- colorRampPalette(c("white","#009E73","black"))(21)[c(5,10,15)]

my_feature_name <- "genebody"


#################################################        



for(my_assay in my_assays[-4]){
        
        
        my_rates_tmp <- my_feature_rates[,grep(my_assay, colnames(my_feature_rates))]
        #my_rates_tmp <- my_rates_tmp[, c(grep("Zy", colnames(my_rates_tmp)), grep("2C|8C|ES", colnames(my_rates_tmp)))]
        
        
        dot_colors <- my_color_palette  
        
        pdf(paste0(output_plots,"/dotplot.",my_assay,".pdf"), height = 9, width = 9, useDingbats = F)
        
        for(my_feature_name in rownames(my_rates_tmp)){
                
                
                par(mfrow=c(2,2), mar = c(5,10,3,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))
                
                df_tmp <- data.frame(value = my_rates_tmp[my_feature_name, ],
                                     Conditions = gsub(paste0("X|_",my_assay,"|_[1-3]"),"",colnames(my_rates_tmp)))
                
                df_tmp$Conditions <- factor(df_tmp$Conditions, levels = unique(df_tmp$Conditions))
                
                
                set.seed(123)
                
                plot(x = jitter(as.integer(df_tmp$Conditions),factor = 0.1),
                     y = df_tmp$value,
                     main =  paste0(my_assay, " at ", my_feature_name),
                     xlim = c(min(as.integer(df_tmp$Conditions))-0.5, max(as.integer(df_tmp$Conditions))+0.5),
                     ylim = c(min(df_tmp$value)-0.05, max(df_tmp$value)+0.05),
                     xlab = "",
                     ylab = "Fraction of length normalized counts", xaxt = "n",
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
                
                
                
                glht_out <- summary(glht(fit, linfct = mcp(Conditions = c("`ChIL100C` - `ChIP` = 0",
                                                                          "`ChIL1000C` - `ChIP` = 0",
                                                                          "`ChIL1000C` - `ChIL100C` = 0"))))
                
                
                
                stats_tmp <- data.frame(`Estimate` =   round(glht_out$test$coefficients*100, 2)/100,
                                        `Std. Error` = round(glht_out$test$sigma*100, 2)/100,
                                        `t value` =  round(glht_out$test$tstat*100, 2)/100,
                                        `Pr(>|t|)` = round(glht_out$test$pvalues*100, 2)/100,
                                        check.names = FALSE)
                
                
                textplot(stats_tmp, halign = "left", valign = "top", cex = 0.6)
                text(x = 0.5, y = 0.75, cex = 0.8,
                     labels = "Multiple Comparisons of Means: User-defined Contrasts")
                
                
                
                rm(list = "df_tmp")
                #rm(list = "stats_tmp")
                
                
                #################################################           
                
        }
        
        
        
        dev.off()
        
}

##################################################################################################################################
##################################################################################################################################


























##################################################################################################################################
##################################################################################################################################







writeLines(capture.output(sessionInfo()), output_sessionInfo)






##################################################################################################################################
##################################################################################################################################




