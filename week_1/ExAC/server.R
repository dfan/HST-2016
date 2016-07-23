
#server.R

library(shiny)
library(RMySQL)
con <- dbConnect(MySQL(), user = 'root', 
                 password = 'root', dbname = 'kohane_lab', host = 'localhost',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
table_query <- function(gene, threshold) {
  paste("select Position, abs(1.0*Allele_Count_African/Allele_Number_African - 1.0*Allele_Count_European_Non_Finnish/Allele_Number_European_Non_Finnish) as Diff from ",
        paste(" having Diff > "," order by Diff desc", sep = as.character(threshold)), sep=gene)
}
cut_query <- function(gene, threshold) {
  paste(paste("select count(NumKey) from "," where abs(1.0*Allele_Count_African/Allele_Number_African - 1.0*Allele_Count_European_Non_Finnish/Allele_Number_European_Non_Finnish) > ", sep = gene), as.character(threshold), sep = "")
}
total_query <- function(gene) {
  paste("select count(*) from ",gene)
}

shinyServer(
  function(input, output) {
    
    db.mybpc3 <- reactive ({
      dbGetQuery(con, table_query("mybpc3",input$threshold))
    })
    db.myh7 <- reactive ({
      dbGetQuery(con, table_query("myh7",input$threshold))
    })
    
    
    #num.mybpc3 <- reactive ({
    #  dbGetQuery(con, cut_query("mybpc3",input$threshold))[1,1]
    #})
    #den.mybpc3 <- reactive ({
    #  dbGetQuery(con, total_query("mybpc3",input$threshold))[1,1]
    #})
    #num.myh7 <- reactive ({
    #  dbGetQuery(con, cut_query("myh7",input$threshold))[1,1]
    #})
    #den.myh7 <- reactive ({
    #  dbGetQuery(con, total_query("myh7",input$threshold))[1,1]
    #})
    
    output$title.mybpc3 <- renderText({"MYBPC3"})
    output$title.myh7 <- renderText({"MYH7"})
    
    output$caption.mybpc3 <- renderText({
      paste(as.character(dbGetQuery(con, cut_query("mybpc3",input$threshold))[1,1]), as.character(dbGetQuery(con, total_query("mybpc3"))[1,1]), sep = " / ")
    })
    
    output$caption.myh7 <- renderText({
      paste(as.character(dbGetQuery(con, cut_query("myh7",input$threshold))[1,1]), as.character(dbGetQuery(con, total_query("myh7"))[1,1]), sep = " / ")
    })
    
    output$db.mybpc3 <- renderTable({ db.mybpc3() })
    output$db.myh7 <- renderTable({ db.myh7() })

    output$dbcdf <- renderPlot({
      plotvar1 <- db.mybpc3()[,2]
      plotvar2 <- db.myh7()[,2]
      if (length(c(plotvar1, plotvar2))>1 ) {
        plot(cumsum(plotvar2/sum(plotvar2)), main = "Empirical CDF Plot of Allele Frequency Differences between AFR and NFE", ylab = "CDF", xlab = "Alleles Ranked by Difference")
        points(cumsum(plotvar1/sum(plotvar1)), pch = 20, col = 'red')
        legend("bottomright",c("MYBPC3","MYH7"), col = c("red","black"), pch = c(20,1))
      }
    })
    
    #output$testplot <- renderPlot({
    #  barplot(matrix(c(num.mybpc3,den.mybpc3,num.myh7,den.myh7), byrow = T, nrow = 2))
    #})
    
  })

