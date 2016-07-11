load("~/Documents/Kohane_Lab/week_4/appdata.RData")
library(dplyr)

#Prints the log values in their proper form as: decimal = fraction
antilog.print <- function(x) {
  paste(signif(10^x,3),paste("1",signif(1/(10^x),3),sep = "/"),sep = " = ")
}

shinyServer(function(input, output, session) {
  observe({
    #Calculates the allele frequency threshold, based on threshold OR parameters
    if (input$use.af == "Threshold") af <- 10^input$threshold
    else af <- 10^input$prevalence * 10^input$allelic.het / 10^input$penetrance
    
    #Logical denoting allele frequency from overall population vs. population max 
    overall <- input$af.calc == "Overall Population"
    if (overall) {
      #Take all 1000g IDs that are in the exac IDs that pass the threshold, unlists to all 26,263
      var.pass <- unlist(lapply(1:num.genes, function(i) #num.genes is the number of HCM panel genes
          data.1000g.id[[i]][,1] %in% #Column 1 contains variant IDs, in the form "position-ref-alt"
            data.exac.id[[i]][(as.numeric(data.exac.id[[i]][,2]) < af),1] 
            #All ExAC variant IDs such that the allele frequencies (in column 2) are greater than the cutoff (af)
        )) 
    } else {
      var.pass <- unlist(lapply(1:num.genes, function(i) {
        max_pop <- apply((select(data.exac.full[[i]], contains("Allele.Count."))) / 
                           (select(data.exac.full[[i]], contains("Allele.Number."))),1,max)
        data.1000g.id[[i]][,1] %in% 
          data.exac.id[[i]][(max_pop < af),1]
      })) #Same thing, but take the maximum allele frequency from each population as the ExAC allele frequencies
    }
    #Get the gold standard variant status. 
    #std_var.1 was calculated using allele frequencies from overall population, std_var.2 used population max.
    if (overall) std_var_pass <- std_var.1[,2]
    else std_var_pass <- std_var.2[,2]
    
    #var.pass is a logical vector for which variants passed (using cutoff from parameteres)
    #std_var_pass is a logical vector for which variants passed (using gold-standard cutoff)
    tpr.pt <- sum(var.pass & std_var_pass) / sum(std_var_pass)
    fpr.pt <- sum(var.pass & !std_var_pass) / sum(intersect.var.num & !std_var_pass) 
    #intersect.var.num is used to only count the variants that are shared between ExAC and 1000G; otherwise, FPR peaks at 0.03.
    
    output$roc_plot <- renderPlot({ 
      tpr <- tpr.1; spec.fpr <- spec.fpr.1
      #Plots the ROC curve saved from earlier (using a wide range of values)
      plot(spec.fpr,tpr, type = 'l', xlab = "False Positive Rate", ylab = "True Positive Rate",
         main = paste(paste("FPR:",signif(fpr.pt,2)),paste("TPR:",signif(tpr.pt)), sep = ", ")) 
      #Plots the singular point calculated using parameter cutoffs
      points(fpr.pt, tpr.pt, col = 'red')
    })
    
    #Various text prints to track values
    output$allele_freq <- renderText({ 
      antilog.print(log10(af))
    })
    output$pen <- renderText({ 
      antilog.print(input$penetrance)
    })
    output$prev <- renderText({ 
      antilog.print(input$prevalence)
    })
    output$ah <- renderText({ 
      antilog.print(input$allelic.het)
    })
    output$point <- renderText({ 
      paste(paste("FPR:",signif(fpr.pt,2)),paste("TPR:",signif(tpr.pt)))
    })
    
  })
})
