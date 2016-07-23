

library(dplyr)
library(tidyr)

library(scrapeR)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
#AGMC.genes <- read.table(paste(getwd(),"Rel_Incidental_Genes",sep="/"), header = FALSE, stringsAsFactors = FALSE)$V1

library(scrapeR)
ACMG.page <- scrape(url ="http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/")
ACMG.table <- sapply(1:3, function(x) as.character(readHTMLTable(ACMG.page[[1]])[[1]][,x]))
mismatch <- which(ACMG.table[,1] == "MedGen")
ACMG.table[mismatch,3] <- ACMG.table[mismatch,2]
ACMG.table <- ACMG.table[-1,-2]
ACMG.table[50,2] <- ACMG.table[49,2]
ACMG.table[41:43,1] <- ACMG.table[40,1]
ACMG.table[45,1] <- ACMG.table[44,1]
ACMG.table[19,2] <- ACMG.table[19,1]
ACMG.table[19,1] <- ACMG.table[18,1]
colnames(ACMG.table) <- c("Disease","Gene")
ACMG.table <- tbl_df(ACMG.table)
full.table <- separate(ACMG.table, Disease, sep = " \\(", into = c("Disease_Name","Disease_MIM")) %>% 
  separate(Gene, sep = " \\(", into = c("Gene_Name","Gene_MIM"))
ACMG.table <- full.table[,c(3,1)]
  
#doc <- getURL("http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/") %>% 
#  htmlParse %>% xpathSApply("//a/@href")
#sapply(links[grep("clinsig_pathogenic", links)], function(link) {
#  paste("http://www.ncbi.nlm.nih.gov", link, sep = "")
#})

import.file.exac <- function(gene) {
  name <- paste(gene,"exac.csv", sep = "_")
  output <- read.csv(paste(getwd(),name,sep="/ExAC_ACMG/"), stringsAsFactors = FALSE)
  output[nchar(paste(output$Alternate,output$Reference))==3,]
}
ACMG.exac <- lapply(ACMG.table$Gene_Name, import.file.exac)
names(ACMG.exac) <- ACMG.table$Gene_Name


#Take relevant subset (position, ref, alt, allele.freq)
ACMG.exac.id <- lapply(data.exac.full, function(x) {
  columns <- select(x, Position, Reference, Alternate, Allele.Frequency) %>% 
    arrange(Position, Reference, Alternate) %>%
      unite(sep= "-")
    #sapply(c("Position","Reference","Alternate","Allele.Frequency"), function(y) which(colnames(x) %in% y))
  #x[order(x$Position),columns] 
})
#Combine position-ref-alt for left-join
ACMG.exac.id <- lapply(data.exac, function(x) cbind(apply(x[,-4],1, function(y) paste(y,collapse = "-")),x[,4]) )

x %>% unite(ab,a,b, sep = "-")





brca1.1000g.full <- import.file.1000g("BRCA1", FALSE)
brca1.exac.full <- import.file.exac("BRCA1")

#bring back the lapply analysis from HCM.panel --> ACMG.panel
#try to at least scrape the site itself??

#Take relevant subset (position, ref, alt, allele.freq)
columns <- sapply(c("POS","REF","ALT","AF"), function(y) which(colnames(brca1.1000g.full) %in% y))
brca1.1000g <- brca1.1000g.full[order(brca1.1000g.full$POS),columns]

columns <- sapply(c("Position","Reference","Alternate","Allele.Frequency"), function(y) which(colnames(brca1.exac.full) %in% y))
brca1.exac <- brca1.exac.full[order(brca1.exac.full$Position),columns] 

#Combine position-ref-alt for left-join
brca1.1000g.id <- cbind(apply(brca1.1000g[,-4],1, function(y) paste(y,collapse = "-")),brca1.1000g[,4])
brca1.exac.id <- cbind(apply(brca1.exac[,-4],1, function(y) paste(y,collapse = "-")),brca1.exac[,4]) 

#Left-join
shared <- intersect(brca1.1000g.id[,1], brca1.exac.id[,1])
output <- data.frame(shared,
                     as.numeric(brca1.1000g.id[ brca1.1000g.id[,1] %in% shared, 2]),
                     as.numeric(brca1.exac.id[ brca1.exac.id[,1] %in% shared, 2]))
colnames(output) <- c("ID","1000G","Exac")
output
cor(output[,2], output[,3])

output$Gene <- "BRCA1"
output$Disease <- "Breast/Ovarian Cancer"
output$Disease_Prev <- 0.123
output$Citation <- "http://seer.cancer.gov/statfacts/html/breast.html"
output$Chrom <- 17
output$Variant_Name <- as.character(output$ID)
output$NumKey <- 1:nrow(output)
output <- separate(output, "ID", c("Position","Ref","Alt"))
output <- output[,c("NumKey","Variant_Name","Gene","Chrom","Position","Ref","Alt","Disease","1000G","Exac","Disease_Prev","Citation")]
glimpse(output)

write.csv(output, file = paste(getwd(),"BRCA1.csv",sep="/"), row.names = F)



