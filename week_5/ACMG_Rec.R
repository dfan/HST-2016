#James Diao
#July 11, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
load(file = paste(getwd(), "ACMG.files", sep = "/HST-2016/week_5/"))

#Left-join
freq.join <- lapply(1:num.genes, function(i) {
  shared <- intersect(ACMG.1000g.id[[i]][,1], ACMG.exac.id[[i]][,1])
  data.frame("Gene"     = I(gene.list[i]),
             "ID"       = I(shared),
             "AF_1000G" = as.numeric(ACMG.1000g.id[[i]] [ ACMG.1000g.id[[i]][,1] %in% shared, 2]),
             "AF_ExAC"  = as.numeric(ACMG.exac.id[[i]]  [ ACMG.exac.id[[i]][,1]  %in% shared, 2]))
})
af_cor <- sapply(1:num.genes, function(i) cor(freq.join[[i]]['AF_1000G'], freq.join[[i]]['AF_ExAC']))
plot(af_cor, ylim = c(0,1))

#full.genes <- table(ACMG.table$Gene_Name)
#full.genes <- full.genes[full.genes > 1]

#Setup Data Table
output <- lapply(1:num.genes, function(i) {
  temp <- separate(freq.join[[i]], "ID", c("Chrom","Position","Ref","Alt"))
  temp$Position <- as.numeric(temp$Position)
  data.frame(temp,
    "Disease" = "Unknown",
    "Prevalence" = 0,
    "Citation" = "Unknown"
  )
})
output <- do.call(rbind, output)
output <- data.frame("NumKey" = 1:nrow(output), output)

setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
write.table(output, file = paste(getwd(),"ACMG_output.txt",sep="/HST-2016/week_5/"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



#output$Disease <- "Breast/Ovarian Cancer"
#output$Disease_Prev <- 0.123
#output$Citation <- "http://seer.cancer.gov/statfacts/html/breast.html"
#output$Chrom <- 17
#output$Variant_Name <- as.character(output$ID)
#output$NumKey <- 1:nrow(output)
#output <- separate(output, "ID", c("Position","Ref","Alt"))
#output <- output[,c("NumKey","Variant_Name","Gene","Chrom","Position","Ref","Alt","Disease","1000G","Exac","Disease_Prev","Citation")]
#glimpse(output)
#
#
#
#
