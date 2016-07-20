#James Diao
#July 20, 2016
#Kohane Lab | HST Summer
#Pipeline Server

library(shiny)
library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
library(ggbiplot)

#################################
####    Main Shiny Server    ####
#################################

shinyServer(function(input, output) {

  ###################################
  ####  Setting up dependencies  ####
  ###################################

  withProgress(message = "--- Setting Up Dependencies ---", value = 0, {

    ### Connect to UCSC to download all refGene names (used for checking gene validity)
    setProgress(0.2, detail = "Connecting to UCSC Genome Browser")
    for (con in dbListConnections(MySQL())) dbDisconnect(con)
    con <- dbConnect(MySQL(), user = 'genome',dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                     unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
    query <- function (input) { suppressWarnings(dbGetQuery(con, input)) }
    total.genes <- "select name2 from refGene" %>% query %>% unlist %>% unique
    present.files <- system("ls", intern = T)

    setProgress(0.6, detail ="Downloading Phase 3 Populations Map")
    ### Download phase 3 map from 1000 genomes FTP
    if (!("phase3map.txt" %in% present.files))
      system("wget -O phase3map.txt ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
    map <- read.table(file = paste(getwd(),"phase3map.txt",sep = "/"), header = T) %>% tbl_df
    header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", map$sample %>% as.character)

    setProgress(0.9, detail ="Scraping LMM and Clinvar websites for HCM- and ACMG- relevant genes")
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

  })

  ### Three possibilities, depending on input$add (Radiobuttons)
  ### (1) manual input: split by "," and " "; (2) HCM.panel; (3) ACMG.panel
  genes.all <- reactive({
    switch(input$add,
      Manual = strsplit(input$genes," ") %>% unlist %>% strsplit(",") %>% unlist,
      HCM = HCM.panel,
      ACMG = ACMG.panel
    )
  })
  ### How many genes are we working with
  genes.len <- reactive({ length(genes.all()) })
  ### Outputs selected genes in tabulated format
  output$selected <- renderTable({
    paste(genes.all(), collapse = ", ")
    if (!is.null(genes.all())) {
      # Aim for squarish table
      row <- ceiling(sqrt(genes.len()))
      # Fill in remaining squares with "---"
      c(genes.all(),rep("---", ceiling(genes.len()/row)*row - genes.len())) %>% matrix(nrow = row, byrow = F)
    }
  } #, include.colnames=FALSE
  )

  ### Sets up new directory for relevant files
  observeEvent(input$set, {
    dir <- sprintf("%s1000_genomes_populations_analysis_%s",input$wd,Sys.Date())
    system(paste("mkdir",dir))
    setwd(dir)
  })

  ### Tracks whether genes are present in UCSC RefGene data table
  valid <- reactive({ genes.all() %in% total.genes %>% mean == 1 })
  output$update <- renderText({
      paste("All Genes Valid:",valid())
  })


  ### When you click "MAKE PLOT: "
  observeEvent(input$run, {

  if(valid()) {
    gene.list <- genes.all()
    num.genes <- genes.len()

    ####################################################################################
    ###  Connecting to UCSC Genome Browser to extract gene info: chrom, start, stop  ###
    ####################################################################################
    withProgress(message = "Step 1 of 4: Connecting to UCSC Genome Browser", value = 0, {
    for (con in dbListConnections(MySQL())) dbDisconnect(con)
    con <- dbConnect(MySQL(), user = 'genome',dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                     unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
    query <- function (input) { suppressWarnings(dbGetQuery(con, input)) }



    ####################################################
    ###  Function for downloading 1000 genomes data  ###
    ####################################################
    download_1000g <- function(gene) {
      gene %>% paste(which(gene.list==gene)) %>% paste(length(gene.list), sep = "/") %>% print
      refGene <- sprintf("select * from refGene where name2 = \"%s\" limit 20", gene) %>% query
      UCSC <- select(refGene, name, chrom, start = txStart, end = txEnd)
      if (nrow(UCSC) == 0) #No hit on refGene
        print("NOT FOUND")
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
        sprintf(command, UCSC$chrom, version, chrom.num, UCSC$start, UCSC$end, gene) %>% system

        # Checks whether the file exists
        exists <- grepl(paste(gene,"_genotypes.vcf",sep =""), system("ls", intern = T)) %>% sum > 0
        file.size <- strsplit(paste("stat ","_genotypes.vcf", sep = gene) %>% system(intern = T), " ")[[1]][8]
        if (exists & file.size > 0)
          print(UCSC$name)
        else print("UNKNOWN FAILURE")
      }
    }



    ###################
    ###   Download  ###
    ###################
    setProgress(0.1, message = sprintf("Step 2 of 4: Downloading VCF files from 1000 Genomes"))
    present <- sapply(gene.list, function(x) paste(x,"genotypes.vcf",sep = "_") %in% present.files)
    download <- sapply(1:num.genes, function(i) {
      incProgress(0.4/num.genes, detail = sprintf("%s [%s/%s]", gene.list[[i]], i, num.genes))
      if (!present[i])
        download_1000g(gene.list[[i]])
    })
    system("rm *.genotypes.vcf.gz.tbi")


    ########################################################
    ###  Import VCF into R                               ###
    ###  Process (0|0) (0|1) (1|0) (1|1) --> FALSE/TRUE  ###
    ########################################################

    import.file.1000g <- function(gene) {
      name <- paste(gene,"genotypes.vcf", sep = "_")
      output <- read.table(paste(getwd(),name,sep="/"), stringsAsFactors = FALSE)
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

    # Import 1000G data
    setProgress(0.5, message = "Step 3 of 4: Importing VCF", detail = " ")
    data.1000g <- lapply(1:num.genes, function(i){
      incProgress(0.4/num.genes, detail = sprintf("%s [%s/%s]", gene.list[[i]], i, num.genes))
      import.file.1000g(gene.list[[i]])
    })
    names(data.1000g) <- gene.list

    ### Contains information about what populations are found where-
    setProgress(10/10, message = "Step 4 of 4: Finalizing Plot", detail = " ")
    pop.table <- (scrape(url ="http://www.1000genomes.org/category/population/")[[1]] %>% readHTMLTable)[[1]] %>%
      tbl_df %>% select(contains("Population"))

    # Super and specific population codes
    super <- pop.table$`Super Population Code` %>% as.character
    names(super) <- pop.table$`Population Code`

    ### Values to be plotted
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

    title <- switch(input$add,
                    Manual = paste(gene.list, collapse = "-"),
                    paste(input$add,"Genes"))

    plot.out <- ggplot(values, aes(x=Population, y=Mean, fill = Superpopulation)) +
      geom_bar(stat = "identity") + ylim(0,1.1*max(values$SD+values$Mean)) +
      geom_errorbar(aes(ymin=Mean - SD, ymax=Mean + SD, width = 0.5)) + theme_minimal() +
      ggtitle(title) + theme(axis.text.x = element_text(angle = -45, hjust = 0.4))


    output$down <- downloadHandler(
      filename = function() {
        sprintf("Plot_%s.%s", title, input$pic)
      },
      content = function(file) {
        ggsave(filename = file, plot = plot.out, device = input$pic, units = "in", width = 8, height = 4)
      }
    )

    output$export <- downloadHandler(
      filename = function() {
        sprintf("Data_%s.%s", title, "RData")
      },
      content = function(file) {
        save(data.1000g, file = file)
      }
    )

    output$tbl <- downloadHandler(
      filename = function() {
        sprintf("Table_%s.%s", title, "txt")
      },
      content = function(file) {
        write.table(values, file = file, col.names = T, row.names = F, quote = F, sep = "\t")
      }
    )

    output$pop_plot <- renderPlot({ plot.out })


    })
  }
  })
})





