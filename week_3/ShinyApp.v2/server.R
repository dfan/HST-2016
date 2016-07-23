load("~/Documents/Kohane_Lab/week_3/ExAC_HCM.RData")
library(dplyr)
library(tidyr)
run <- T

shinyServer(function(input, output, session) {
  observe({
    gene <- input$gene
    var_id <- input$var_id
    var_list <- as.character(data.total.id[[gene]][,1])
    var_num <- which(var_list %in% var_id)
    plot.density <- as.numeric(input$plot.density)
    
    updateSelectInput(session, "var_id", choices = var_list, selected = var_id)
    
    if (run) {
      params <- list(penetrance = as.numeric(input$penetrance), prevalence = as.numeric(input$prevalence), allelic.het = as.numeric(input$allelic.het))
      params[[which(c("Penetrance","Prevalence","Allelic Heterogeneity") == input$change.on)]] <- seq(0,1,length.out = plot.density)
      af_threshold <- params$prevalence * params$allelic.het / params$penetrance
        
      if (F) {
        plot.density <- 21
        af_threshold <- 0.002*0.001/seq(0,1,length.out = plot.density)
        gene <- "MYBPC3"
        var_id <- "47369220-C-T"
        var_list <- as.character(data.total.id[[gene]][,1])
        var_num <- which(var_list %in% var_id)
      }
      gene.exac.id <- data.exac.id[[gene]]
      gene.1000g.full <- data.1000g.full[[gene]]
      gene.logical.1000g <- unlist(data.1000g.full[[gene]][var_num,11:ncol(data.1000g.full[[gene]])])
      #mean(gene.logical.1000g) #Allele Freq
      
      #Logical vector for which genes passed the cutoff
      genes.pass <- lapply(af_threshold, function(af) 
        apply(gene.1000g.full[,c(3,5,6)],1, function(y) paste(y,collapse = "-")) %in% 
          gene.exac.id[(as.numeric(gene.exac.id[,2]) < af),1] #allele_pass
      )
      #sapply(genes.pass, mean)
      #plot(.Last.value)
      
      #Logical vector for which patients are diseased
      empty.row <- rep(FALSE,ncol(gene.1000g.full)-10)
      patient.status <- lapply(genes.pass, function(gp) {
        sub <- gene.1000g.full[gp,11:ncol(gene.1000g.full)] #Subset only variants in allele_pass
        if (nrow(sub)==0) empty.row
        sub %>% colSums>0
      })
      #Proportion of HCM-positive
      #sapply(patient.status, mean)
      #plot(.Last.value)
      
      roc <- sapply(patient.status, function(ps) 
        rbind(sum(gene.logical.1000g & ps) / sum(ps), #tpr, #NaN means everyone is healthy (disease = 0)
              sum(gene.logical.1000g & !ps) / sum(!ps)) #fpr, #NaN means everyone is diseased (healthy = 0)
      )
      rownames(roc) <- c("True Positive Rate", "False Positive Rate")
      roc.print <- roc[, !is.na(colSums(roc))]
    
    }
    
    output$track <- renderText({ })
    output$roc_plot <- renderPlot({ 
      plot(roc.print[1,], roc.print[2,], ylab = "True Positive Rate", xlab = "False Positive Rate", main = "ROC Curve")
    })
    output$roc_table <- renderTable({ 
      roc
    })
    output$ps_plot <- renderPlot({ 
      plot(seq(0,1,length.out = plot.density), sapply(patient.status, mean), ylab = "Proportion HCM-Positive", xlab = input$change.on)
    })
    output$genes_plot <- renderPlot({ 
      plot(seq(0,1,length.out = plot.density), sapply(genes.pass, mean), ylab = "Proportion 1000G-variants Passed", xlab = input$change.on)
    })
    output$threshold_plot <- renderPlot({ 
      plot(seq(0,1,length.out = plot.density), af_threshold, ylab = "Allele Frequency Threshold", xlab = input$change.on)
      abline(h=min(as.numeric(gene.exac.id[,2])), col = 'red') #Marks minimum 
      abline(h=max(as.numeric(gene.exac.id[,2])), col = 'blue')
    })
    output$exac_1000g_freq <- renderText({ 
      paste("1000G: ",paste(data.total.id[[gene]][var_num,2], data.total.id[[gene]][var_num,3], sep = ", ExAC: "), sep = "")
    })
    
  })
})
