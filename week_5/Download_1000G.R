#James Diao
#July 11, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")

####################################################################
###  Scrape ACMG website for gene names and associated diseases  ###
####################################################################

ACMG.page <- scrape(url ="http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/")
ACMG.table <- sapply(1:3, function(x) as.character(readHTMLTable(ACMG.page[[1]])[[1]][,x]))
### Formatting corrections: NAs
badrow <- which(is.na(ACMG.table[,3]))
ACMG.table[badrow,3] <- ACMG.table[badrow-1,3]
### Formatting corrections: Sliding
mismatch <- 0
while(sum(ACMG.table[,3] == "ClinVar")>0) {
  mismatch <- which(ACMG.table[,2]!="MedGen")
  ACMG.table[mismatch,2:3] <- ACMG.table[mismatch,1:2]
  for (row in mismatch) { ACMG.table[row,1] <- ACMG.table[row-1, 1] }
}
### Take 2 relevant columns: gene names and associated diseases
ACMG.table <- ACMG.table[-1,-2]
colnames(ACMG.table) <- c("Disease","Gene")
ACMG.table <- separate(tbl_df(ACMG.table), Disease, sep = " \\(", into = c("Disease_Name","Disease_MIM")) %>%
  separate(Gene, sep = " \\(", into = c("Gene_Name","Gene_MIM"))
ACMG.table$Disease_MIM <- sapply(ACMG.table$Disease_MIM, function(x) strsplit(x, "[ )]")[[1]][2] )
ACMG.table$Gene_MIM <- sapply(ACMG.table$Gene_MIM, function(x) strsplit(x, "[ )]")[[1]][2] )
gene.list <- unique(ACMG.table$Gene_Name)
num.genes <- length(gene.list)
### Cleanup
rm(ACMG.page, badrow, mismatch, row)

#################################################################################
###  Connecting to UCSC Genome Browser to extract gene info: chrom, start, stop  ###
#################################################################################
for (con in dbListConnections(MySQL())) dbDisconnect(con)
con <- dbConnect(MySQL(), user = 'genome',
                 dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
query <- function (input) { suppressWarnings(dbGetQuery(con, input)) }
setwd("/Users/jamesdiao/Documents/Kohane_Lab/1000G")

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

### DOWNLOAD ###
download <- sapply(gene.list, download_1000g)
system("rm *.genotypes.vcf.gz.tbi; ls")
download

header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
#write.table(as.matrix(download), quote = FALSE, sep = "\t", col.names = FALSE,
#            paste(getwd(),"refGene_sources.txt",sep="/1000G/"))

#save.image(file = paste(getwd(), "ACMG.files", sep = "/HST-2016/week_5/"))

