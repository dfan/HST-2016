load("~/Documents/Kohane_Lab/week_3/ExAC_HCM.RData")
library(dplyr)
library(tidyr)
counts <- sapply(data.1000g.full, nrow)
names(counts) <- HCM.panel
choice_names <- c("Penetrance","Prevalence","Allelic Heterogeneity")
defaults <- c(0.8, 0.002, 0.001)
plot.density <- 21
test <- F

shinyServer(function(input, output, session) {
  observe({
    gene <- input$gene
    var_id <- input$var_id
    var_list <- as.character(data.1000g.id[[gene]][,1])
    var_num <- which(var_list %in% var_id)
    change.on <- choice_names %in% input$change.on
    
    updateSelectInput(session, "var_id", choices = var_list, selected = var_id)
    updateSliderInput(session, "param1", label = choice_names[!change.on][1], value = defaults[!change.on][1])
    updateSliderInput(session, "param2", label = choice_names[!change.on][2], value = defaults[!change.on][2])
    
    if (test == F) {
    if (input$change.on == "Penetrance") {
      penetrance <- seq(0,1,length.out = 20)
      prevalence <- input$param1
      allelic.het <- input$param2
    }
    if (input$change.on == "Prevalence") {
      penetrance <- input$param1
      prevalence <- seq(0,1,length.out = 20)
      allelic.het <- input$param2
    }
    if (input$change.on == "Allelic Heterogeneity") {
      penetrance <- input$param1
      prevalence <- input$param2
      allelic.het <- seq(0,1,length.out = 20)
    }
    af_threshold <- prevalence*allelic.het/penetrance
    
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
    patient.status <- lapply(genes.pass, function(gp) 
      gene.1000g.full[gp,11:ncol(gene.1000g.full)] %>% rbind(FALSE) %>% colSums>0 #Subset only variants in allele_pass
    )
    #Proportion of HCM-positive
    #sapply(patient.status, mean)
    #plot(.Last.value)
    
    roc <- sapply(patient.status, function(ps) 
      rbind(sum(gene.logical.1000g & ps) / sum(ps), #tpr, #NaN means everyone is healthy (disease = 0)
            sum(gene.logical.1000g & !ps) / sum(!ps)) #fpr, #NaN means everyone is diseased (healthy = 0)
    )
    rownames(roc) <- c("True Positive Rate", "False Positive Rate")
    roc.print <- roc[, !is.na(colSums(roc))]
    
    #output$track1 <- renderText({ penetrance })
    #output$track2 <- renderText({ prevalence })
    #output$track3 <- renderText({ allelic.het })
    }
    
    x <- 1:100
    output$pen_plot <- renderPlot({ 
      plot(roc.print[1,], roc.print[2,], ylab = "True Positive Rate", xlab = "False Positive Rate", main = "ROC Curve")
    })
    output$ps_plot <- renderPlot({ 
      plot(sapply(patient.status, mean), ylab = "Proportion HCM-Positive")
    })
    output$pen_table <- renderTable({ 
      roc
    })
    output$pen_text <- renderText({ 
      paste("1000G: ",paste(data.total.id[[gene]][var_num,2], data.total.id[[gene]][var_num,3], sep = ", ExAC: "), sep = "")
    })
    
    output$prev_plot <- renderPlot({ 
      plot(roc.print[1,], roc.print[2,], ylab = "True Positive Rate", xlab = "False Positive Rate", main = "ROC Curve")
    })
    output$prev_table <- renderTable({ 
      roc
    })
    
    output$ah_plot <- renderPlot({ 
      plot(roc.print[1,], roc.print[2,], ylab = "True Positive Rate", xlab = "False Positive Rate", main = "ROC Curve")
    })
    output$ah_table <- renderTable({ 
      roc
    })
    
    updateTabsetPanel(session, "change.on", selected = paste("panel",which(change.on)))
    
  })
})
