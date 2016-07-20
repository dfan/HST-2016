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