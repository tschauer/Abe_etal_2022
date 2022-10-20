






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




load("genome/ranges_TSS.rda")
load("genome/ranges_TTS.rda")
load("genome/ranges_genebody.rda")
load("genome/ranges_intergenic.rda")
load("genome/ranges_b10.rda")


my_TSS_bins <- subsetByOverlaps(my_b10_ranges, my_TSS_ranges, minoverlap = 2501L)
export.bed(object = my_TSS_bins, con = "genome/bins_TSS.bed")

my_TTS_bins <- subsetByOverlaps(my_b10_ranges, my_TTS_ranges, minoverlap = 2501L)
export.bed(object = my_TTS_bins, con = "genome/bins_TTS.bed")


my_genebody_bins <- subsetByOverlaps(my_b10_ranges, my_genebody_ranges, type = "within")
export.bed(object = my_genebody_bins, con = "genome/bins_genebody.bed")


my_intergenic_bins <- subsetByOverlaps(my_b10_ranges, my_intergenic_ranges, type = "within")
export.bed(object = my_intergenic_bins, con = "genome/bins_intergenic.bed")

        




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



feature_bins <- gsub("my_|_bins", "", ls(pattern = "my.*_bins"))


my_feature_widths <- sapply(feature_bins, function(bid){ sum(width(get(paste0("my_", bid,"_bins"))))/1e3   })

names(my_feature_widths) <- gsub("n$","ns",names(my_feature_widths))



##################################################################################################################################
##################################################################################################################################





##################################################################################################################################
##################################################################################################################################

#################################################        Sum Counts          #####################################################



my_feature_rates <- t(sapply(feature_bins, function(bid){
        
        sapply(gsub(".*\\/b10.|.rda","",input_counts), function(fid){
                
                my_counts <- get(paste0("b10.", fid)) 
                #my_counts[round(my_counts,6) == 0] <- NA
                
                my_feature_bins <- get(paste0("my_",bid,"_bins"))

                my_counts <- my_counts[names(my_counts) %in% my_feature_bins$b10_id]
                
                rate <- mean(my_counts, na.rm = TRUE)
        })      
}))



my_feature_rates <- t(t(my_feature_rates) / colSums(my_feature_rates))



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
        my_rates_tmp <- my_rates_tmp[, c(grep("Zy", colnames(my_rates_tmp)), grep("2C|8C|ES", colnames(my_rates_tmp)))]
        
        my_order <- rev(c("TSS","genebody","TTS", "intergenic"))
        
        bp <- barplot( my_rates_tmp[my_order,], 
                       main = my_assay,
                       ylab = "Fraction", las= 3, xaxt = "n",
                       col = bp_colors)
        
        axis(side = 1, at = bp, labels =  colnames(my_rates_tmp), 
             cex.axis=0.6, tick = FALSE, line = -0.5, las=2)
        
        
        plot.new()
        legend("center", 
               legend = rev(rownames( my_rates_tmp[my_order,])),
               fill = rev(bp_colors))
        
        dev.off()
        
        
        rm(list = ls(pattern = "tmp"))
        
}

#################################################        



##################################################################################################################################
##################################################################################################################################




















##################################################################################################################################
##################################################################################################################################

#################################################        Barpot Rates        #####################################################




my_color_palette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#F0E442")

my_feature_name <- "genebody"


#################################################        



for(my_assay in my_assays){
        
        
        my_rates_tmp <- my_feature_rates[,grep(my_assay, colnames(my_feature_rates))]
        my_rates_tmp <- my_rates_tmp[, c(grep("Zy", colnames(my_rates_tmp)), grep("2C|8C|ES", colnames(my_rates_tmp)))]
        
        
        if(my_assay == "S5P"){
                dot_colors <- my_color_palette  
        } else {
                dot_colors <- my_color_palette[-3]
        }
        
        
        
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
                     ylab = "Fraction", xaxt = "n",
                     pch = 19, col = dot_colors[df_tmp$Conditions])
                
                set.seed(123)
                
                points(x = jitter(as.integer(df_tmp$Conditions), factor = 0.1),
                       y =  df_tmp$value, 
                       col = "#555555", pch=1, lwd=0.5)
                
                axis(side = 1, at = seq_along(levels(df_tmp$Conditions)), 
                     labels = levels(df_tmp$Conditions), las = 2,)
                
                
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

                par(mar = c(0,0,0,0))
                
                textplot(stats_tmp, halign = "left", valign = "top", cex = 0.7)
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




