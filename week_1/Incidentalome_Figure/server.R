
#server.R

library(shiny)
shinyServer(
  function(input, output) {
    test <- reactive ({
      10^input$log10test
    })
    fpr <- reactive ({
      10^input$log10fpr
    })
    
    output$distPlot <- renderPlot({
      test_num <- floor(seq(0,test(),length.out = 150))[-1]
      plot(test_num, 1-(1-fpr())^test_num, 
           main = "Figure. Percentage of Total Population with a False-Positive Test Result", 
           ylab = "Percentage of Total Population with a False-Positive Test Result", 
           xlab = "No. of Independent Tests")
    })
    output$disp_test <- renderText({
      paste("No. of Independent Tests: ",as.character(test()), sep = "")
    })
    output$disp_fpr <- renderText({
      paste("False-Positive Rate: ",as.character(fpr()), sep = "")
    })
    output$disp_equ <- renderText({
      if (fpr() == 1) "P(at least 1 false-positive) = 1"
      else paste("P(at least 1 false-positive) = 1-(",")^n",sep = as.character(1-fpr()))
    })
    output$disp_endval <- renderText({
      paste("End value = ","%", sep = as.character(100*round(1-(1-fpr())^test(),2)))
    })
  })

