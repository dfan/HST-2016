
#ui.R
library(shiny)
library(RMySQL)
con <- dbConnect(MySQL(), user = 'root', unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock",
                 password = 'root', dbname = 'kohane_lab', host = 'localhost')

# Define UI for application that draws a histogram
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("No. Positions with Allele Frequency Differences Between AFR & NFE > Threshold"), 
  
  # Sidebar with slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("threshold","Allele Frequency Threshold:",
                  min = 0, max = 1, value = 0.01, step = 0.01)
    ),
    mainPanel(
      textOutput("title.mybpc3")
      ,textOutput("caption.mybpc3")
      #tableOutput("db.mybpc3")
      ,textOutput("title.myh7")
      ,textOutput("caption.myh7")
      #tableOutput("db.myh7")
      #tableOutput("stack_plot")
      ,plotOutput("dbcdf")
      ,plotOutput("testplot")
    )
  ),
  mainPanel()
))


