#James Diao
#July 23, 2016
#Kohane Lab

library(ggbiplot)
library(scrapeR)
library(dplyr)
library(tidyr)

########################
####  Dependencies  ####
########################

setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
HST.folder.6 <- function(file.name) { paste(getwd(), file.name, sep = "/HST-2016/week_6/") }
merged <- read.csv(file = HST.folder.6("merged_output.csv"), stringsAsFactors = F)
pop.table <- (scrape(url ="http://www.1000genomes.org/category/population/")[[1]] %>% readHTMLTable)[[1]] %>% tbl_df %>% select(contains("Population"))
# Super and specific population codes
super <- pop.table$`Super Population Code` %>% as.character
names(super) <- pop.table$`Population Code`
download.file(url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
              destfile = paste(getwd(),"phase3map.txt",sep = "/"), method = "internal")
map <- read.table(file = paste(getwd(),"phase3map.txt",sep = "/"), header = T) %>% tbl_df
unlink("phase3map.txt")

################################
####  Process clinvar data  ####
################################

process.clinvar <- function(file) {
  ### Sanity checks
  # Confirmed: start = stop
  mean(file$Start==file$Stop)
  # Confirmed: all are pathogenic
  #mean(file$Assertion=="pathogenic")
  # Confirmed: all chromosomes are present
  #sort(levels(file$Chromosome))
  # Collect all relevant rows
  file$Stop <- NULL #remove
  file$Assertion <- NULL #remove
  file <- dplyr::rename(file, Position = Start) #rename
  file <- unite(file, col = "VarID", Chromosome, Position, Ref, Alt, sep = "_", remove = FALSE) #unite
  # Remove indels
  to.rm <- (file$Ref=="-") | (file$Alt=="-")
  mean(to.rm) # == 0.05531253, not removing too many
  file <- file[-to.rm,]
  # How many diseases
  length(unique(file$Disease)) == 2858 # wow that's a lot
  file$VarID %>% length() # == 10286
  file$VarID %>% unique %>% length() # == 8441
  ### There are nonunique variants associated with different diseases
  return(file)
}

merged <- process.clinvar(merged)

#load("~/Documents/Kohane_Lab/Eigenstrat/data.1000g")
#data.1000g <- do.call("rbind",ACMG.1000g)
#save.image(file = "join_tools")

load("join_tools")
if (!("VarID" %in% colnames(data.1000g)))
  data.1000g <- data.1000g %>% unite(col = "VarID", CHROM, POS, REF, ALT, sep = "_", remove = FALSE)

### UHHH Difference between full and segmented by 66 vars

# All variants are unique
#mean(table(data.1000g$VarID)) == 1
join <- intersect(data.1000g$VarID, merged$VarID)
ind.1000g <- which(data.1000g$VarID %in% join)
ind.clinvar <- which(merged$VarID %in% join)


gene.counts <- sapply(ACMG.1000g, nrow)
key <- cumsum(gene.counts)
gene.loc <- sapply(ind.1000g, function(i) {
  gene.list[sum(i>key)+1]
}) #which genes are the joined variants in?

zeros <- gene.list[!(gene.list %in% gene.loc)]
rel.genes <- c(table(gene.loc), rep(0, length(zeros)) %>% setNames(zeros)) %>% sort(decreasing = T)
barplot(rel.genes, las =2, ylab = "Number of Clinvar Pathogenic Variants", main = "Distribution of 137 Pathogenic Variants by Gene") #

### isolate from ACMG.1000g: gives you 1000g AF, how many in each population have it.
sub.1000g <- data.frame(GENE = gene.loc, data.1000g[ind.1000g,])
sub.clinvar <- merged[ind.clinvar,]

values <- sapply(unique(gene.loc), function(gene) {
  mean((sub.1000g %>% subset(GENE==gene) %>% select(contains("HG"), contains("NA")) %>% colSums)>0)
})
### Fraction with at least 1 mutation in each
barplot(sort(values, decreasing = T), las = 2, ylim = c(0, 0.5), ylab = "Have at least 1 mutation", main = "Proportion of Patients from 1000G with at Least 1 Variant in Different Genes")
sort(values, decreasing = T) #%>% round(2)
frac.population.with.variant <- values[unique(gene.loc)]
num.variants <- rel.genes[unique(gene.loc)]
cor(frac.population.with.variant, num.variants)

###### NOT USEFUL
values <- sapply(unique(gene.loc), function(gene) {
  temp <- sub.1000g %>% subset(GENE==gene) %>% select(contains("HG"), contains("NA")) %>% colSums
  #vector of total variant positions for each patient, in 1 gene.
  c(mean(temp), sd(temp))
}) %>% t %>% data.frame(unique(gene.loc))  # 37 elements for each gene
colnames(values) <- c("Mean","SD","Gene")
values$Gene <- values$Gene %>% as.character
values <- arrange(values, desc(Mean))
values$Gene <- values$Gene %>% factor(levels = values$Gene)

plot.pop <- ggplot(values, aes(x=Gene, y=Mean)) +
  geom_bar(stat = "identity") + ylim(-0.20,max(values$SD+values$Mean)) +
  ggtitle("137 ClinVar Variants in ACMG genes") + xlab("Population") +
  ylab("Mean Number of (0|1) (1|0) (1|1) Variant Sites Across All Patients") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD, width = 0.5)) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plot.pop
###### End

#data.1000g[,11:ncol(data.1000g)] <- apply(data.1000g[,11:ncol(data.1000g)], 2, function(y) grepl("1",y))

# For every population, sum up the columns in that population
# Print the mean and sd
sapply(levels(map$pop), function(pop) {
    temp <- data.1000g[,10+which(map$pop == pop)] %>% colSums
    c(mean(temp), sd(temp))
  }) %>% t %>% tbl_df -> values

colnames(values) <- c("Mean","SD")
rownames(values) <- levels(map$pop)
values$Population <- levels(map$pop)
values$Superpopulation <- super[levels(map$pop)]
ord <- super[levels(map$pop)] %>% order
values <- values[ord,]
values$Population <- values$Population %>% factor(levels = values$Population)

plot.pop <- ggplot(values, aes(x=Population, y=Mean, fill = Superpopulation)) +
  geom_bar(stat = "identity") + ylim(-0.26,1.1*max(values$SD+values$Mean)) +
  ggtitle("Only ClinVar-Pathogenic Variants (137)") + xlab("Population") + ylab("Mean Number of (0|1) (1|0) (1|1) Variant Sites") +
  geom_errorbar(aes(ymin=Mean - SD, ymax=Mean + SD, width = 0.5)) + theme_minimal() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0.4))
plot.pop


# For every population, sum up the columns in that population ONLY in the
# 137 pathogenic variants that overlap with Clinvar
# Print the proportion of individuals in that population that have at least 1 variant
sapply(levels(map$pop), function(pop) {
  temp <- data.1000g[ind.1000g,10+which(map$pop == pop)] %>% colSums
  mean(temp>0)
})-> values

#sapply(levels(map$pop), function(pop) {
#  temp <- data.1000g[ind.1000g,10+which(map$pop == pop)] %>% colSums
#  c(mean(temp), sd(temp))
#}) %>% t %>% tbl_df -> values.137

#Calculating test statistics (F-values)
 groups <- ifelse(super[levels(map$pop)]=="EUR","EUR","Other")
 data <- data.frame(y = values, group = factor(groups))
 fit <- lm(y ~ group, data)
 plot(y ~ group, data, las = 2)
 anova(fit)

# desired plotting order of populations (alphabetal WITHIN each superpopulation)
ord <- super[levels(map$pop)] %>% order
# levels of superpopulation ("AFR","AMR","EAS","EUR","SAS")
super.levels <- sort(unique(super))
# levels of population ("ACB","ASW","ESN"..."CEU","FIN","GBR"...)
pop.levels <- levels(map$pop)[ord]

# set of all colors used (by superpopulation)
col.set <- length(super.levels) %>% rainbow %>% setNames(super.levels)
# turn into a vector for plotting
col.use <- col.set[super[pop.levels]]
barplot(values[ord], las = 2, ylim = c(0,1), col = col.use, ylab = "Proportion with at least 1 variant", main = "Carrier Fraction Across All 137 ClinVar Variants")
legend("topright", legend = sort(unique(super)), pch = 20, cex = 0.9, col = rainbow(5))

values.all <- data.frame("Proportion" = values,
                         "Population" = factor(names(values), levels = pop.levels),
                         "Superpopulation" = factor(super[names(values)], levels = super.levels))
plot.pop <- ggplot(values.all, aes(x=Population, y=Proportion, fill = Superpopulation)) +
  geom_bar(stat = "identity") + ylim(0,1) +
  ggtitle("Carrier Fraction Across All 137 ClinVar Variants") + xlab("Population") + ylab("Proportion with at least 1 variant") +
  theme_minimal() + theme(axis.text.x = element_text(angle = -45, hjust = 0.4))
plot.pop



plot.pop <- ggplot(values, aes(x=Gene, y=Mean)) +
  geom_bar(stat = "identity") + ylim(-0.20,max(values$SD+values$Mean)) +
  ggtitle("137 ClinVar Variants in ACMG genes") + xlab("Population") +
  ylab("Mean Number of (0|1) (1|0) (1|1) Variant Sites Across All Patients") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD, width = 0.5)) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plot.pop

#%>% t %>% tbl_df -> values

### Correlate allele frequencies
### First: account for different lengths by nonuniqueness
af.1000g <- data.1000g[ind.1000g,'AF'] %>% as.numeric
af.esp <- merged[ind.clinvar,'ESP_Overall_Frequency']
af.exac <- merged[ind.clinvar,'ExAc_Overall_Frequency']
cor(af.1000g, af.exac)
sum(!is.na(af.esp))



### Once you have the actual data.
setwd("/Users/jamesdiao/Documents/Kohane_Lab/HST-2016/week_7")
#save.image(file = "join_tools_7-27")
load("join_tools_7-27")

ACMG_Lit <- read.csv(file = "ACMG_Lit_Small.csv", header = TRUE, stringsAsFactors = F) #some table
prevalence <- ACMG_Lit$Inverse.Prevalence #some vector
disease <- ACMG_Lit$Disease

tags <- c("adenomatous", rep("aneurysm",3), rep("arrhythmogenic;dreifuss",5),
          rep("breast;ovarian",2),"brugada;gardner","tachycardia","dilated","ehler",
          "fabry","hypercholesterolemia",rep("hypertrophic",8),
          "medullary","hypercholesterolemia","noncompaction",
          "Fraumeni",rep("Loeys;Dietz",5),rep("QT",3),"lynch;endometrial",
          "hyperthermia","Marfan",rep("neoplasia;men2a",3),"MYH;colon","neurofibromatosis",
          rep("paraganglioma;pheochromocytoma",4),"peutz;jeghers","pilomatrixoma",
          "Cowden;PTEN;hamartoma;Merkel","retinoblastoma",rep("tuberous",2),"Hippel;Lindau","Wilms")

getExACAlleleFreq <- function(input, tags, data) {
  temp <- sapply(unique(tags), function(tag) {
    tag.vec <- strsplit(tag,";") %>% unlist
    loc <- rep(0,nrow(input))
    for(tag in tag.vec)
      loc <- loc | grepl(tag,input$Disease, ignore.case = T)
    freq.data <- input %>% subset(loc) %>% select(contains(data))
    #freq.esp <- input %>% subset(loc) %>% select(freq.use)
    missing <- is.na(freq.data)
    #freq.exac[missing] <- freq.esp[missing]
    hits <- sum(loc)
    final <- 1-prod(1-freq.data[!missing])
    c(final, hits) %>% setNames(c("AF","Hits"))
  }) %>% t %>% tbl_df
  data.frame("Tags" = unique(tags), temp)
}

get1000GAlleleFreq <- function(input, tags) {
  temp <- sapply(unique(tags), function(tag) {
    tag.vec <- strsplit(tag,";") %>% unlist
    loc <- rep(0,nrow(input))
    for(tag in tag.vec)
      loc <- loc | grepl(tag,input$Disease, ignore.case = T)
    final <- mean(colSums( input[loc,11:ncol(input)] )>0)
    hits <- sum(loc)
    c(final, hits) %>% setNames(c("AF","Hits"))
  }) %>% t %>% tbl_df
  data.frame("Tags" = unique(tags), temp)
}

sub.clinvar$G1000_Overall_Frequency <- sapply(sub.clinvar$VarID, function(id) {
  sub.1000g %>% subset(sub.1000g$VarID == id) %>% select(AF) %>% unlist %>% as.numeric
})

#Write the 1000G info into clinvar
temp <- sapply(sub.clinvar$VarID, function(id) {
  sub.1000g %>% subset(sub.1000g$VarID == id) %>% select(13:ncol(sub.1000g)) %>% unlist
}) %>% t %>% tbl_df
sub.clinvar <- cbind(sub.clinvar,temp)


freq <- getExACAlleleFreq(merged,tags, "ExAc")
freq <- get1000GAlleleFreq(sub.clinvar,tags)

bp <- freq$AF %>% setNames(freq$Tags) %>% sort(decreasing = T)
par(mar=c(5, 15, 5, 2)) #changes plotting window to have greater left-margins
barplot(bp, las = 2, pch = 'h', xlab ="P(having a variant)", main = "Carrier frequency by disease",
        horiz = T, xlim = c(0,max(bp))*1.1, las = 1)
par(mar=c(5, 4, 4, 2)+0.1) #resets margins

# 1000G frequency, ExAC frequency, prevalence estimate, penetrance estimate
#

# Map of disease name to disease tags
cbind("Disease" = disease,"PATTERN" = tags %>% unique)
named.freqs <- bp[tags] %>% setNames(disease)
named.prev <- 1/ACMG_Lit$Inverse.Prevalence %>% setNames(disease)
# Repeats allow for correct quartile calculations
allelic.het <- c(0.001,0.001,0.02,0.5,0.5)
# Point values for penetrance
penetrance <- named.prev[unique(disease)]/named.freqs[unique(disease)] * median(allelic.het)

# Matrix of penetrance values for allelic het range, capped at 1
pen.unlist <- sapply(named.prev[unique(disease)]/named.freqs[unique(disease)], function(x) x*allelic.het) %>% as.vector
inf.loc <- which(pen.unlist %>% is.infinite)[c(T,F,F,F,F)]
for(i in inf.loc)
  pen.unlist[seq(i,i+4)] <- c(0,0,1,1,1)
pen.unlist[pen.unlist>1] <- 1

# replicate each element n times to create labels
data <- data.frame("penetrance" = pen.unlist, "disease" = sapply(disease, function(x) rep(x,length(allelic.het))) %>% as.vector)
par(mar=c(5, 22, 5, 2))
boxplot(penetrance ~ disease, data, horizontal = TRUE, las = 1, xlab = "Penetrance", main = "Penetrance Range Estimates for P(V|D) = 0.001, 0.02, 0.5, AF from 1000G")
par(mar=c(5, 4, 4, 2)+0.1)



sub.clinvar$Matched <- c("PTEN hamartoma tumor syndrome",
          "PTEN hamartoma tumor syndrome",
          "MYH-associated polyposis",
          "Unknown",
          "Unknown",
          "Unknown",
          "MYH-associated polyposis",
          "MYH-associated polyposis",
          "Unknown",
          "Unknown",
          "MYH-associated polyposis",
          "MYH-associated polyposis",
          "MYH-associated polyposis",
          "MYH-associated polyposis",
          "Familial medullary thyroid carcinoma",
          "Paragangliomas 1",
          "Multiple endocrine neoplasia type 2a",
          "Unknown",
          "PTEN hamartoma tumor syndrome",
          "Paragangliomas 1",
          "Unknown",
          "Paragangliomas 1",
          "Unknown",
          "Paragangliomas 1",
          "Unknown",
          "PTEN hamartoma tumor syndrome",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 10",
          "Breast-ovarian cancer familial 2",
          "Unknown",
          "Familial hypertrophic cardiomyopathy 1",
          "Familial hypertrophic cardiomyopathy 1",
          "Breast-ovarian cancer familial 1",
          "Marfan's syndrome",
          "Arrhythmogenic right ventricular cardiomyopathy type 10",
          "Familial hypertrophic cardiomyopathy 1",
          "Unknown",
          "Unknown",
          "Unknown",
          "Familial hypertrophic cardiomyopathy 7",
          "Familial hypercholesterolemia",
          "Lynch syndrome",
          "Lynch syndrome",
          "Lynch syndrome",
          "Neurofibromatosis type 2",
          "Unknown",
          "Arrhythmogenic right ventricular cardiomyopathy type 5",
          "Lynch syndrome",
          "Lynch syndrome",
          "Long QT syndrome 1",
          "Unknown",
          "Dilated cardiomyopathy 1A",
          "Long QT syndrome 1",
          "Brugada syndrome 1",
          "Long QT syndrome 1",
          "Brugada syndrome 1",
          "Unknown",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Unknown",
          "Brugada syndrome 1",
          "Brugada syndrome 1",
          "Brugada syndrome 1",
          "Unknown",
          "Dilated cardiomyopathy 1A",
          "Brugada syndrome 1",
          "Brugada syndrome 1",
          "Brugada syndrome 1",
          "Adenomatous polyposis coli",
          "Adenomatous polyposis coli",
          "Adenomatous polyposis coli",
          "Unknown",
          "Unknown",
          "Arrhythmogenic right ventricular cardiomyopathy type 8",
          "Arrhythmogenic right ventricular cardiomyopathy type 8",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Tuberous sclerosis 1",
          "Tuberous sclerosis 1",
          "Tuberous sclerosis 1",
          "Unknown",
          "Paragangliomas 4",
          "Paragangliomas 4",
          "Paragangliomas 4",
          "Paragangliomas 4",
          "Familial hypertrophic cardiomyopathy 1",
          "Familial hypertrophic cardiomyopathy 1",
          "Familial hypertrophic cardiomyopathy 1",
          "Dilated cardiomyopathy 1A",
          "Left ventricular noncompaction 6",
          "Dilated cardiomyopathy 1A",
          "Familial hypertrophic cardiomyopathy 1",
          "Familial hypertrophic cardiomyopathy 1",
          "MYH-associated polyposis",
          "MYH-associated polyposis",
          "MYH-associated polyposis",
          "MYH-associated polyposis",
          "Unknown",
          "Multiple endocrine neoplasia type 2a",
          "Familial medullary thyroid carcinoma",
          "Multiple endocrine neoplasia type 2a",
          "Multiple endocrine neoplasia type 2a",
          "Familial medullary thyroid carcinoma",
          "Multiple endocrine neoplasia type 2a",
          "Multiple endocrine neoplasia type 2a",
          "Unknown",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Unknown",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Left ventricular noncompaction 6",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Arrhythmogenic right ventricular cardiomyopathy type 5",
          "Breast-ovarian cancer familial 2",
          "Breast-ovarian cancer familial 2",
          "Retinoblastoma",
          "Unknown",
          "Familial hypertrophic cardiomyopathy 1",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Familial hypertrophic cardiomyopathy 4",
          "Breast-ovarian cancer familial 1",
          "Breast-ovarian cancer familial 1",
          "Breast-ovarian cancer familial 1",
          "Breast-ovarian cancer familial 1",
          "Breast-ovarian cancer familial 1",
          "Breast-ovarian cancer familial 1",
          "Breast-ovarian cancer familial 1",
          "Breast-ovarian cancer familial 1",
          "Unknown",
          "Li-Fraumeni syndrome 1",
          "Unknown",
          "Unknown",
          "Unknown",
          "Unknown",
          "Unknown",
          "Unknown",
          "Ehlers-Danlos syndrome, type 4",
          "Lynch syndrome",
          "Lynch syndrome",
          "Lynch syndrome",
          "Lynch syndrome",
          "Lynch syndrome",
          "Long QT syndrome 1",
          "Brugada syndrome 1",
          "Brugada syndrome 1",
          "Long QT syndrome 1",
          "Brugada syndrome 1",
          "Long QT syndrome 1",
          "Brugada syndrome 1",
          "Brugada syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Brugada syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Brugada syndrome 1",
          "Long QT syndrome 1",
          "Familial hypertrophic cardiomyopathy 8",
          "Familial hypertrophic cardiomyopathy 8",
          "Familial hypertrophic cardiomyopathy 8",
          "Adenomatous polyposis coli",
          "Long QT syndrome 1",
          "Unknown",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Long QT syndrome 1",
          "Unknown",
          "Fabry's disease",
          "Fabry's disease"
)
### Look for the diseases on our list







