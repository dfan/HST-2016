#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1)
  stop("Requires 1 argument1 (num.files, genes.per.file)", call.=FALSE)
inputs <- as.vector(unlist(read.table(args[1], header = FALSE)))
trials <- inputs[1]
total <- inputs[2]


suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(scrapeR))
suppressMessages(library(RMySQL))
suppressMessages(library(ggbiplot))

print("Connecting to UCSC Genome Browser", quote = F)
for (con in dbListConnections(MySQL())) dbDisconnect(con)
con <- dbConnect(MySQL(), user = 'genome',dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
query <- function (input) { suppressWarnings(dbGetQuery(con, input)) }
total.genes <- "select name2 from refGene" %>% query %>% unlist %>% unique

if(!exists("trials"))
trials <- 10
total <- 50

for (i in 1:trials) {
  print(sprintf("TRIAL %s", i), quote = F)
  use.genes <- NULL
  while(length(use.genes) < total) {
    gene <- sample(total.genes,1)
    refGene <- sprintf("select * from refGene where name2 = \"%s\" limit 20", gene) %>% query
    chrom <- (select(refGene, chrom) %>% unlist)[1] %>% setNames(NULL)
    if ((chrom %>% strsplit("chr") %>% unlist)[c(F,T)] %in% c(as.character(1:23), "X","Y"))
      use.genes <- c(use.genes, gene)
    else
      print(sprinf("Removed: %s", gene), quote = F)
    print(chrom)
  }
  write.table(use.genes, file = sprintf("genes_%s.txt",i), quote = F, sep = "\t", row.names = F, col.names = F)
}

