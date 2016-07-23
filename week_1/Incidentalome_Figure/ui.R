
#ui.R

library(shiny)

# Define UI for application that draws a histogram
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Incidentalome Figure (Kohane et al. 2006)"), 
  
  # Sidebar with slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("log10test","log10(No. of Independent Tests):",
                  min = 1, max = 6, value = 4, step = 1),
      sliderInput("log10fpr","log10(False-Positive Rate):",
                  min = -6, max = 0, value = -4, step = 1)
    ), 
    mainPanel(
      textOutput("disp_test"),
      textOutput("disp_fpr"),
      textOutput("disp_equ"),
      textOutput("disp_endval"),
      plotOutput("distPlot")
    )
  ),
  mainPanel()
))


