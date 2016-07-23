#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A

library(RMySQL)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
all_cons <- dbListConnections(MySQL())
for (con in all_cons)
  dbDisconnect(con)

con <- dbConnect(MySQL(), user = 'genome',
                 dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
query <- function (input) {
  suppressWarnings(dbGetQuery(con, input))
}

#knownGene.mybpc3 <- rbind(query("select * from knownGene where name = \"uc021qir.1\""),
#                          query("select * from knownGene where name = \"uc021qis.1\""),
#                          query("select * from knownGene where name = \"uc010rhl.2\""))
#MYH7 = \"uc001wjx.3\"


#http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000

gene <- "PRKAG2"
geneReviews <- sprintf("select * from geneReviews where name = \"%s\" limit 20", gene) %>% query
geneReviews
#geneReviewsDetail.mybpc3 <- query("select * from geneReviewsDetail where geneSymbol = \"MYBPC3\" limit 3")
#geneReviewsDetail.mybpc3
#omimId_input.mybpc3 <- query("select * from omimGeneSymbol where geneSymbol = \"MYBPC3\" limit 3")
#omimId.mybpc3 <- as.character(omimId_input.mybpc3[,1])
#omimGeneMap.mybpc3 <- query(paste("select * from omimGeneMap where omimId =",omimId.mybpc3))
#t(omimGeneMap.mybpc3)

refGene <- sprintf("select * from refGene where name2 = \"%s\" limit 20", gene) %>% query
refGene[,-c(9:16)]

#starts <- as.numeric(strsplit(refGene.mybpc3$exonStarts,",")[[1]])
#ends <- as.numeric(strsplit(refGene.mybpc3$exonEnds,",")[[1]])
#exon <- sapply(1:34, function(x) seq(starts[x],ends[x]))
#variants <- as.numeric(rownames(set))
#exon.variants <- set[(variants %in% unlist(exon)),]
#mean(colSums(exon.variants))
#hist(colSums(exon.variants))


