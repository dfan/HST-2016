library(scrapeR)
u <- "http://browser.1000genomes.org/Homo_sapiens/Component/Gene/Web/VariationTable?db=core;g=ENSG00000134571;r=11:47352957-47374253;sub_table=ALL;update_panel=1"
pageSource <- scrape(url=u, headers=FALSE, parse=TRUE)
tables <- readHTMLTable(pageSource[[1]])
tables[[1]]$`Variant identifierID`
