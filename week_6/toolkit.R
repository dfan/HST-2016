#James Diao
#July 12, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
library(ggbiplot)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
HST.folder <- function(file.name) { paste(getwd(), file.name, sep = "/HST-2016/week_5/") }
load(file = HST.folder("ACMG.files"))
load(file = HST.folder("ExAC_1000G_HCM.RData"))

### Information about which individuals are from which populations
map <- read.table(file = paste(getwd(),"phase3map.txt",sep = "/1000G/"), header = T) %>% tbl_df

### Contains information about what populations are found where-
pop.table <- (scrape(url ="http://www.1000genomes.org/category/population/")[[1]] %>% readHTMLTable)[[1]] %>%
  tbl_df %>% select(contains("Population"))
# Super and specific population codes
super <- pop.table$`Super Population Code` %>% as.character
names(super) <- pop.table$`Population Code`

### Subsetted PCA
run.pca <- function(data, size) {
  if (dim(data) %>% is.null)
    data <- do.call("rbind", data)
  sub.pca <- data %>% tbl_df %>% select(contains("V")) %>%
    sample_n(size) %>% t %>% prcomp(center = TRUE)
  plot(sub.pca, type = "l")
  ggbiplot(sub.pca, obsx.scale = 1, var.scale = 1, choices = c(1,2),
    varname.abbrev = FALSE, var.axes = FALSE, varname.size = 3, groups = super[map$pop]) +
    scale_color_discrete(name = '') + theme(legend.position = 'right')
    #ellipse = TRUE, circle = FALSE, labels = tissue
}

### Subsetted PCA
#ACMG.total <- do.call("rbind", ACMG.1000g)
#HCM.total <- do.call("rbind", HCM.1000g.full)
run.pca(ACMG.total, 1000)
run.pca(HCM.total, 1000)

### Avg.var.by.pop
sapply(levels(map$pop), function(pop) {
  sapply(ACMG.1000g, function(acmg) {
     acmg[,10+which(map$pop == pop)] %>% colSums
  }) %>% rowSums %>% mean
}) %>% sort(decreasing = T) -> avg.var.by.pop

### Number of Variants maps closely with Gene Length
ACMG_tx <- read.table(file = HST.folder("ACMG_gene_info_tab_delim.txt"), header = T)
pop.var.status <- read.table(file = HST.folder("Population_Variant_Status.txt"), header = T)
tx.length <- (ACMG_tx$txend - ACMG_tx$txstart) %>% setNames(gene.list)
var.num <- sapply(ACMG.1000g, nrow)
barplot(var.num, main = "Number of Variants", las = 2)
barplot(tx.length, main = "Gene Length (tx region)", las = 2)
cor(var.num, tx.length)
plot(var.num, tx.length, ylab = "Gene Length (tx region)", xlab = "Number of Variants")
abline(a = 0, b = mean(tx.length/var.num))
text(0,max(tx.length), pos = 4, paste("Slope =", mean(tx.length/var.num) %>% round(1), "nucleotides per variant"))


### (0|1)  (1|1)  (1|0)  Sites per Patient
sapply(ACMG.1000g, function(acmg) {
  acmg[,-c(1:10)] %>% colSums
}) %>% rowSums -> var.sites

## Vector of populations in order AFR, AMR, EAS, EUR, SAS
sorted_pop <- levels(map$pop)[super[levels(map$pop)] %>% order]
# Names = populations, values = super populations
super[sorted_pop]

# Set of 26 colors for 26 populations
pop.col <- sorted_pop %>% length %>% rainbow %>% setNames(sorted_pop)
# Set of 5 colors for 5 super populations
super.col <- unique(super) %>% length %>% rainbow %>% setNames(unique(super) %>% sort)
plot(1:26, rep(0,26), col = pop.col, pch = 20, cex = 2)
text(1:26, rep(0.1,26), names(pop.col), cex = 0.6)
points(1:5, rep(0.5,5), col = super.col, pch = 20, cex = 2)
text(1:5, rep(0.6,5), names(super.col), cex = 0.6)

# Puts population levels in order by superpopulation
numkey <- 1:26 %>% setNames(names(super[sorted_pop]))
pop <- numkey[map$pop %>% as.character]
pop.set <- unique(names(pop))
# Orders populations by updated level
ord <- pop %>% order(var.sites) ### var.sites parameter makes all pops increasing
categories <- table(super[names(pop)])
categories
### Sanity checks
sum(names(pop) %in% pop.set[which(super[pop.set] == "AFR")])
sum(super[names(pop)] == "AFR")
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

### Average Number of Variant Sites Positives
col.set <- rainbow(length(unique(super)))
col.use <- col.set[super[names(avg.var.by.pop)] %>% as.numeric]
barplot(avg.var.by.pop, las = 2, col = col.use, ylab = "Average Number of Variant Sites (1|0) (0|1) (1|1)")
legend("bottomright", legend = levels(super), col = col.set, pch = 20, cex = 1)

### 0|0 Sites per Patient
sapply(ACMG.1000g, function(acmg) {
  nrow(acmg) - colSums(acmg[,-c(1:10)])
}) %>% rowSums -> null.sites
hist(null.sites, xlab = "0|0 Sites per Patient", main = NULL, breaks = 30)





