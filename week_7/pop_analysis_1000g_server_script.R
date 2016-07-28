#!/usr/bin/env Rscript

#James Diao
#July 23, 2016
#Kohane Lab | HST
#Pipeline Server

# Rscript --vanilla pop_analysis_1000g_server_script.R genes.txt

args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1)
  stop("Requires 1 argument (table of gene names)", call.=FALSE)
gene.list <- as.vector(unlist(read.table(args[1], header = FALSE)))
num.genes <- length(gene.list)
print(paste("Imported genes:", paste(gene.list, collapse = ", ")), quote = F)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(scrapeR))
suppressMessages(library(RMySQL))
suppressMessages(library(ggbiplot))

to_min <- function(timing) {
  timing <- as.numeric(timing)
  minutes <- floor(timing/60)
  seconds <- round(timing - 60*minutes, 0)
  minute.in <- ifelse(minutes == 1, "1 minute, ",ifelse(minutes == 0, "", paste(minutes,"minutes, ")))
  second.in <- ifelse(seconds == 1, "1 second",ifelse(seconds == 0, "0 seconds", paste(seconds,"seconds")))
  paste(minute.in, second.in, sep = "")
}


###################################
####  Setting up dependencies  ####
###################################

print("Setting Up Dependencies", quote = F)

### Connect to UCSC to download all refGene names
print("Connecting to UCSC Genome Browser", quote = F)
for (con in dbListConnections(MySQL())) dbDisconnect(con)
con <- dbConnect(MySQL(), user = 'genome',dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
query <- function (input) { suppressWarnings(dbGetQuery(con, input)) }
total.genes <- "select name2 from refGene" %>% query %>% unlist %>% unique

print("Downloading Phase 3 Populations Map", quote = F)
#Download populations-superpoplations
pop.table <- (scrape(url ="http://www.1000genomes.org/category/population/")[[1]] %>% readHTMLTable)[[1]] %>% tbl_df %>% select(contains("Population"))  # Super and specific population codes
super <- pop.table$`Super Population Code` %>% as.character
names(super) <- pop.table$`Population Code`
### Download phase 3 map from 1000 genomes FTP
download.file(url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
              destfile = paste(getwd(),"phase3map.txt",sep = "/"), method = "internal")
map <- read.table(file = paste(getwd(),"phase3map.txt",sep = "/"), header = T) %>% tbl_df
ord <- super[levels(map$pop)] %>% order
map$pop <- factor(as.character(map$pop), levels = levels(map$pop)[ord])
unlink("phase3map.txt")
header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", map$sample %>% as.character)

print("Scraping LMM and Clinvar websites for HCM- and ACMG- relevant genes", quote = F)
### Scrape from LMM and Clinvar websites for HCM- and ACMG- relevant genes
LMM.page <- scrape(url="http://personalizedmedicine.partners.org/Laboratory-For-Molecular-Medicine/Tests/Cardiomyopathy/HCM-Panel.aspx", headers=FALSE, parse=TRUE)
HCM.panel <- levels(readHTMLTable(LMM.page[[1]])[[1]]$Gene)
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
ACMG.panel <- ACMG.table[-1,3] %>% strsplit(" \\(") %>% sapply(function(x) x[1]) %>% unique

m <- ceiling(c(18,23) * num.genes/60)
if (diff(m)<1) {
  print("Estimated Runtime: <1 minute", quote = F)
} else {
  print(sprintf("Estimated Runtime: %s-%s minutes", m[1], m[2]), quote = F)
}

#Set Directory
dir <- sprintf("%s/1000genomes_%s",getwd(),Sys.Date())
system(paste("mkdir",dir))
setwd(dir)
print(paste("Directory set to:",dir), quote = F)

#if(!is.null(gene.list)) {

print("Starting Downloads", quote = F)

ptm <- proc.time()
present.files <- system("ls", intern = T)

### FILE INPUT
#dt <- input$data
#if (!is.null(dt)) {
#  title <- dt$name
#  load(dt$datapath, verbose = F)
#  gene.list <- names(data.1000g)
#  num.genes <- length(gene.list)
#} else {
#  gene.list <- rvalues$genes.all
#  num.genes <- rvalues$genes.len


####################################################################################
###  Connecting to UCSC Genome Browser to extract gene info: chrom, start, stop  ###
####################################################################################

for (con in dbListConnections(MySQL())) dbDisconnect(con)
con <- dbConnect(MySQL(), user = 'genome',dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")

####################################################
###  Function for downloading 1000 genomes data  ###
####################################################

download_1000g <- function(gene) {
  #gene %>% paste(which(gene.list==gene)) %>% paste(length(gene.list), sep = "/") %>% print
  refGene <- sprintf("select * from refGene where name2 = \"%s\" limit 20", gene) %>% query
  UCSC <- select(refGene, name, chrom, start = txStart, end = txEnd)
  if (nrow(UCSC) == 0) #No hit on refGene
    print("NOT FOUND", quote = F)
  else {
    if (nrow(UCSC) > 1) { #Multiple hits: take the widest range
      UCSC <- UCSC[which.max(UCSC$end-UCSC$start),]
    }
    # gets [n] from chr[n]
    chrom.num <- strsplit(UCSC$chrom, split = "chr")[[1]][2]
    # different version for chromosomes X and Y
    version <- switch(chrom.num, "X" = "shapeit2_mvncall_integrated_v1b",
                    "Y" = "integrated_v2a", "shapeit2_mvncall_integrated_v5a")
    command <- paste("tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.%s.",
                   "phase3_%s.20130502.genotypes.vcf.gz %s:%s-%s > %s_genotypes.vcf", sep = "")
    suppressMessages(sprintf(command, UCSC$chrom, version, chrom.num, UCSC$start, UCSC$end, gene) %>% system)

    # Checks whether the file exists
    exists <- grepl(paste(gene,"_genotypes.vcf",sep =""), system("ls", intern = T)) %>% sum > 0
    file.size <- strsplit(paste("stat ","_genotypes.vcf", sep = gene) %>% system(intern = T), " ")[[1]][8]
    if (exists & file.size > 0) {
      print("SUCCESS", quote = F)
      unlist(UCSC)
    } else {
      print("UNKNOWN FAILURE", quote = F)
    }
  }
}

###################
###   Download  ###
###################

print("Step 1 of 3: Downloading VCF files from 1000 Genomes", quote = F)
drop <- NULL
#present <- sapply(gene.list, function(x) paste(x,"genotypes.vcf",sep = "_") %in% present.files)
download <- lapply(1:num.genes, function(i) {
  timing <- (proc.time()-ptm)['elapsed'] %>% to_min
  print(sprintf("Downloading VCF from 1000 Genomes - %s [%s/%s]", gene.list[[i]], i, num.genes), quote = F)
  print(paste("Elapsed Time: ", timing, sep = ""), quote = F)
  print("", quote = F)
  download_1000g(gene.list[i])
  if (!grepl("HG000",system(sprintf("tail -n 1 %s_genotypes.vcf", gene.list[i]), intern = T))) {
    drop <- c(drop,i)
    unlink(sprintf("%s_genotypes.vcf", gene.list[i]))
  }
})

failed <- NULL
for (i in length(download):1) {
  if (length(download[[i]]) != 4 | (i %in% drop)) {
    failed <- c(failed, gene.list[i])
    gene.list <- gene.list[-i]
    num.genes <- num.genes -1
    download[[i]] <- NULL
  }
}

download <- do.call("rbind",download) %>% tbl_df
download$start <- as.numeric(download$start)
download$end <- as.numeric(download$end)
system("rm *.genotypes.vcf.gz.tbi")

########################################################
###  Import VCF into R                               ###
###  Process (0|0) (0|1) (1|0) (1|1) --> FALSE/TRUE  ###
########################################################

import.file.1000g <- function(gene) {
  name <- paste(gene,"genotypes.vcf", sep = "_")
  colclass <- rep("character", 2513)
  colclass[c(2,6)] <- rep("integer",2)
  output <- read.table(paste(getwd(),name,sep="/"), stringsAsFactors = FALSE, colClasses = colclass)
  #Add header
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
  # Convert (0|0) (0|1) (1|0) (1|1) --> FALSE/TRUE
  output[,11:ncol(output)] <- apply(output[,11:ncol(output)], 2, function(y) grepl("1",y))
  output
}

##############################
###  Call Import Function  ###
##############################

data.1000g <- lapply(1:num.genes, function(i){
  timing <- (proc.time()-ptm)['elapsed'] %>% to_min
  print(sprintf("Importing %s VCF [%s/%s]",gene.list[[i]], i, num.genes), quote = F)
  print(paste("Elapsed Time: ", timing, sep = ""), quote = F)
  print("", quote = F)
  import.file.1000g(gene.list[[i]])
})
names(data.1000g) <- gene.list

#} # end of (if manual input)


### Contains information about what populations are found where-
print("Step 3 of 3: Finalizing Plots", quote = F)

### Values to be plotted
sapply(levels(map$pop), function(pop) {
  temp <- sapply(data.1000g, function(dt) {
    dt[,10+which(map$pop == pop)] %>% colSums
  }) %>% rowSums
  c(mean(temp), sd(temp))
}) %>% t %>% tbl_df -> values

colnames(values) <- c("Mean","SD")
values$Population <- factor(levels(map$pop), levels = levels(map$pop))
values$Superpopulation <- super[levels(map$pop)]

title <- paste(gene.list, collapse = "-")

tx.length = download$end - download$start
var.num = sapply(data.1000g, nrow)
var.data <- data.frame(gene = gene.list, tx.length,var.num)

plot.pop <- ggplot(values, aes(x=Population, y=Mean, fill = Superpopulation)) +
  geom_bar(stat = "identity") + ylim(0,1.1*max(values$SD+values$Mean)) +
  geom_errorbar(aes(ymin=Mean - SD, ymax=Mean + SD, width = 0.5)) + theme_minimal() +
  ggtitle(title) + xlab("Population") + ylab("Mean No. of Alternate Variants") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0.4))

label.1 <- data.frame(x=min(tx.length), y=max(var.num), label = paste("Slope =", mean(var.num/tx.length) %>% round(3), "variants per nucleotide"))

plot.line <- ggplot(var.data, aes(x = tx.length, y = var.num, label = gene.list)) +
  geom_text(check_overlap = T, vjust = 0, nudge_y = max(var.num)/70) +
  ggtitle(title) + xlab("Gene Length (tx region)") + ylab("Number of Variant Positions") +
  geom_point(size = 2) + geom_text(data = label.1, aes(x=x, y=y, label = label, hjust = 0), size = 6)

if (num.genes > 4) {
  label.2 <- data.frame(x=min(tx.length), y=(max(var.num)-min(var.num))*0.9, label = paste("Correlation =", cor(tx.length, var.num) %>% round(3)))
  plot.line <- plot.line + geom_text(data = label.2, aes(x=x, y=y, label = label, hjust = 0), size = 6) +
    geom_abline(intercept = 0, slope = mean(var.num/tx.length))
}

population.data <- sapply(data.1000g, function(data) {
  sapply(levels(map$pop), function(pop) {
    (subset(data, select = c(rep(FALSE,10), map$pop == pop)) %>% colSums > 0) %>% mean
  })
}) %>% tbl_df


### Gene-specific population fraction plotting
###
#if (input$select %in% gene.list) {
#  pop.data <- data.frame("Population" = factor(levels(map$pop), levels = levels(map$pop)),
#                         "Fraction" = population.data[,input$select] %>% unlist %>% setNames(NULL),
#                         "Superpopulation" = super[levels(map$pop)])
#  plot.frac <- ggplot(pop.data, aes(x=Population, y=Fraction, fill = Superpopulation)) + geom_bar(stat = "identity") +
#    theme_minimal() + xlab("Population") + ylab("Fraction with at least 1 alternate variant") +
#    ggtitle(label = input$select) + theme(axis.text.x = element_text(angle = -45, hjust = 0.4))
#  output$frac_plot <- renderPlot({ plot.frac })
#}

pdf(file = sprintf("Plots_%s.pdf", title), onefile = TRUE, width = 8, height = 4)
print(plot.pop)
print(plot.line)
suppressMessages(dev.off())

save(download,data.1000g, file = sprintf("Data_%s.RData", title))

write.table(values, file = sprintf("Table_%s.txt", title), col.names = T, row.names = F, quote = F, sep = "\t")

timing <- (proc.time()-ptm)['elapsed'] %>% to_min

print(paste("Runtime: ", timing, sep = ""), quote = F)

if (is.null(failed))
  failed <- "None"
print(paste("Failed Downloads:", failed), quote = F)

print(


#} #Closes (is gene.list empty)