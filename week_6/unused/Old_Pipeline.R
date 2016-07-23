#James Diao
#July 19, 2016
#Kohane Lab
#Pipeline

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
library(ggbiplot)

#######################################
###  Input: Relevant list of genes  ###
#######################################
gene.list <- c("MYBPC3","MYH7","BRCA1","BRCA2","RB1")
num.genes <- length(gene.list)

####################################################################################
###  Connecting to UCSC Genome Browser to extract gene info: chrom, start, stop  ###
####################################################################################
for (con in dbListConnections(MySQL())) dbDisconnect(con)
con <- dbConnect(MySQL(), user = 'genome',
                 dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
query <- function (input) { suppressWarnings(dbGetQuery(con, input)) }

########################################
###  Modify: Sets working directory  ###
########################################
system("mkdir test")
setwd(paste(getwd(),"test",sep = "/"))

####################################################
###  Function for downloading 1000 genomes data  ###
####################################################
download_1000g <- function(gene) {
  # Progress info
  gene %>% paste(which(gene.list==gene)) %>% paste(length(gene.list), sep = "/") %>% print
  #geneReviews <- sprintf("select * from geneReviews where name = \"%s\" limit 3", gene) %>% query
  #UCSC <- select(geneReviews, chrom, start = chromStart, end = chromEnd)
  refGene <- sprintf("select * from refGene where name2 = \"%s\" limit 20", gene) %>% query
  UCSC <- select(refGene, name, chrom, start = txStart, end = txEnd)
  if (nrow(UCSC) == 0) #No hit on geneReviews/refGene
    print("NOT FOUND")
  else {
    if (nrow(UCSC) > 1) { #Multiple hits: take the widest range
      UCSC <- UCSC[which.max(UCSC$end-UCSC$start),]
      #print("MULTIPLE")
    }
    # gets [n] from chr[n]
    chrom.num <- strsplit(UCSC$chrom, split = "chr")[[1]][2]
    # different version for chromosomes X and Y
    version <- switch(chrom.num, "X" = "shapeit2_mvncall_integrated_v1b",
                      "Y" = "integrated_v2a", "shapeit2_mvncall_integrated_v5a")
    # Executes something like:
    # tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
    # ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 2:39967768-39967768 |
    # vcf-subset -c HG00099 > test.vcf

    command <- paste("tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.%s.",
                     "phase3_%s.20130502.genotypes.vcf.gz %s:%s-%s > %s_genotypes.vcf", sep = "")
    sprintf(command, UCSC$chrom, version, chrom.num, UCSC$start, UCSC$end, gene) %>% system

    # Checks whether the file exists
    exists <- grepl(paste(gene,"_genotypes.vcf",sep =""), system("ls", intern = T)) %>% sum > 0
    file.size <- strsplit(paste("stat ","_genotypes.vcf", sep = gene) %>% system(intern = T), " ")[[1]][8]
    if (exists & file.size > 0)
      print(UCSC$name)
    else print("UNKNOWN FAILURE")
  }
}

#################
###  Download ###
#################
download <- sapply(gene.list, download_1000g)
system("rm *.genotypes.vcf.gz.tbi; ls")
download

#####################
###  Dependencies ###
#####################
system("wget -O phase3map.txt ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
map <- read.table(file = paste(getwd(),"phase3map.txt",sep = "/"), header = T) %>% tbl_df
header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", map$sample %>% as.character)

import.file.1000g <- function(gene, cut) {
  name <- paste(gene,"genotypes.vcf", sep = "_")
  output <- read.table(paste(getwd(),name,sep="/"), stringsAsFactors = FALSE)
  #Add header
  if(cut) output <- output[,1:length(header)]
  names(output)[1:length(header)] <- header
  #Remove all single alt indels
  output <- output[nchar(output$REF)==1,] #deletions
  output <- output[grepl(",",output$ALT) | nchar(output$ALT)==1,] #insertions
  #Add AF Column
  af <- sapply(output$INFO, function(x) strsplit(x,";")[[1]][2] )
  af <- sapply(af, function(x) substring(x,4,nchar(x)))
  names(af) <- NULL
  output <- cbind("AF"=af, output)
  output$AF <- as.character(output$AF)
  #Separate paired and unpaired
  alt <- strsplit(output$ALT,",")
  paired <- sapply(alt,length)!=1
  if (sum(paired)!=0) {
    sub.out <- output[paired,]
    output <- output[!paired,]
    #View all rows with multiple alts
    #sub.out[,1:6]
    #Remove all paired indels
    master <- lapply(lapply(alt[paired],nchar), function(x) x==1)
    new.af <- lapply(sub.out$AF, function(x) strsplit(x,",")[[1]])
    new.alt <- lapply(sub.out$ALT, function(x) strsplit(x,",")[[1]])
    sub.out$AF <- lapply(1:nrow(sub.out), function(x) new.af[[x]][master[[x]]])
    sub.out$ALT <- lapply(1:nrow(sub.out), function(x) new.alt[[x]][master[[x]]])
    #Add back the rows that are now single missense
    add.back <- sub.out[sapply(sub.out$AF, length)==1,]
    add.back$AF <- unlist(add.back$AF)
    add.back$ALT <- unlist(add.back$ALT)
    output <- rbind(output,add.back)
    #Extract (or actually, just sum together) double missenses
    sub.out <- sub.out[sapply(sub.out$AF, length) >1,]
    if (nrow(sub.out)!=0) {
      #how many cycles
      max.alts <- max(sapply(sub.out$ALT,length))
      for (i in 1:max.alts) {
        #Each temp is used holds the ith REF-ALT match
        temp <- sub.out
        temp$ALT <- sapply(1:nrow(sub.out), function(x) sub.out$ALT[[x]][i])
        temp$AF <- sapply(1:nrow(sub.out), function(x) sub.out$AF[[x]][i])
        output <- rbind(output,temp[!is.na(temp$AF),])
      }
    }
  }
  if (!cut)
    output[,11:ncol(output)] <- apply(output[,11:ncol(output)], 2, function(y) grepl("1",y))
  output
}

# Import 1000G data
data.1000g <- lapply(gene.list, import.file.1000g, cut = F)
#Take relevant subset (position, ref, alt, allele.freq)
data.1000g.id <- lapply(data.1000g, function(x) {
  select(x, CHROM, POS, REF, ALT, AF) %>%
    arrange(CHROM, POS, REF, ALT) %>%
    unite("ID", c(CHROM, POS, REF, ALT), sep= "-")
})
names(data.1000g) = names(data.1000g.id) = gene.list

### Contains information about what populations are found where-
pop.table <- (scrape(url ="http://www.1000genomes.org/category/population/")[[1]] %>% readHTMLTable)[[1]] %>%
  tbl_df %>% select(contains("Population"))
# Super and specific population codes
super <- pop.table$`Super Population Code` %>% as.character
names(super) <- pop.table$`Population Code`

### Values
sapply(levels(map$pop), function(pop) {
  temp <- sapply(data.1000g, function(data) {
    data[,10+which(map$pop == pop)] %>% colSums
  }) %>% rowSums
  c(mean(temp), sd(temp))
}) %>% t %>% tbl_df -> values

colnames(values) <- c("Mean","SD")
rownames(values) <- levels(map$pop)
values$Population <- levels(map$pop)
values$Superpopulation <- super[levels(map$pop)]
ord <- super[levels(map$pop)] %>% order
values <- values[ord,]
values$Population <- values$Population %>% factor(levels = values$Population)

p <- ggplot(values, aes(x=Population, y=Mean, fill = Superpopulation)) +
  geom_bar(stat = "identity") + ylim(0,500) +
  geom_errorbar(aes(ymin=Mean - SD, ymax=Mean + SD, width = 0.5))
p + theme_minimal()










### (0|1)  (1|1)  (1|0)  Sites per Patient
sapply(data.1000g, function(data) {
  data[,-c(1:10)] %>% colSums
}) %>% rowSums -> var.sites

## Vector of populations in order AFR, AMR, EAS, EUR, SAS
sorted_pop <- levels(map$pop)[super[levels(map$pop)] %>% order]
# Names = populations, values = super populations
super[sorted_pop]


# Set of 26 colors for 26 populations
pop.col <- sorted_pop %>% length %>% rainbow %>% setNames(sorted_pop)
# Set of 5 colors for 5 super populations
super.col <- unique(super) %>% length %>% rainbow %>% setNames(unique(super) %>% sort)

# Puts population levels in order by superpopulation
numkey <- 1:length(unique(map$pop)) %>% setNames(names(super[sorted_pop]))
pop <- numkey[map$pop %>% as.character]
pop.set <- unique(names(pop))
# Orders populations by updated level
ord <- pop %>% order(var.sites) ### var.sites parameter makes all pops increasing
categories <- table(super[names(pop)])
pop <- names(pop)

plot(var.sites[ord], xlab = "Patient Number", ylab = "(0|1)  (1|1)  (1|0)  Sites per Patient",
     xlim = c(-150,2504), main = NULL, col= pop.col[pop[ord]], type = 'h')
legend("topleft", legend = names(pop.col), pch = 20, cex = 0.9, col = pop.col)
abline(v = c(0,categories %>% cumsum))
text(categories %>% cumsum, max(var.sites), names(categories), pos = 2)

plot(var.sites[ord], xlab = "Patient Number", ylab = "(0|1)  (1|1)  (1|0)  Sites per Patient",
     main = NULL, col = super.col[super[pop[ord]]], type = 'h')
abline(v = c(0,categories %>% cumsum))
text(categories %>% cumsum, max(var.sites), names(categories), pos = 2)


