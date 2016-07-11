load("~/Documents/Kohane_Lab/week_4/ExAC_1000G_HCM.RData")
library(dplyr)
library(tidyr)

af_threshold <- c(seq(0,      0.00009, 0.00001),
                  seq(0.0001, 0.0009,  0.0001),
                  seq(0.001,  0.009,   0.001),
                  seq(0.01,   0.09,    0.01),
                  seq(0.1,    1,      0.1))
plot(log10(af_threshold))

#Uses Overall Population Allele Frequency
var.pass.full.1 <- sapply(af_threshold, function(af) {
  unlist(lapply(1:num.genes, function(i)
    data.1000g.id[[i]][,1] %in% 
      data.exac.id[[i]][(as.numeric(data.exac.id[[i]][,2]) < af),1]
  ))
})

#Uses Population Max Allele Frequency
var.pass.full.2 <- sapply(af_threshold, function(af) {
  unlist(lapply(1:num.genes, function(i) {
    max_pop <- apply((select(data.exac.full[[i]], contains("Allele.Count."))) / 
                       (select(data.exac.full[[i]], contains("Allele.Number."))),1,max)
    data.1000g.id[[i]][,1] %in% 
      data.exac.id[[i]][(max_pop < af),1]
  }))
})

var.pass.full <- var.pass.full.2

intersect.var.num <- unlist(lapply(1:num.genes, function(i)
  data.1000g.id[[i]][,1] %in% data.exac.id[[i]][,1]
))

rownames(var.pass.full) <- std_var[,1]
colnames(var.pass.full) <- af_threshold
var.pass.full[1:5,1:5] #rows = variants, cols = threshold, box = passed or not
apply(var.pass.full,2,sum)
plot(af_threshold, .Last.value, ylab = "Number Variants Passed")

std_var_pass <- std_var[,2]

switch.pt <- ncol(var.pass.full)-apply(var.pass.full,1,sum)
switch.pt <- switch.pt[switch.pt!=101]
af_threshold[switch.pt]
hist(af_threshold[switch.pt], xlab = "Threshold", ylab = "How many pass")

tpr <- colSums(var.pass.full & std_var_pass) / sum(std_var_pass)
fpr <- colSums(var.pass.full & !std_var_pass) / sum(!std_var_pass)
spec.fpr <- colSums(var.pass.full & !std_var_pass) / sum(intersect.var.num & !std_var_pass)
roc <- rbind(tpr, spec.fpr)
plot(spec.fpr,tpr)

#IDEA: make shiny app where you can move along the curve as a red dot 
#and give positional information
