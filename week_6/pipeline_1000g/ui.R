#James Diao
#July 23, 2016
#Kohane Lab | HST
#Pipeline UI

library(shiny)
library(shinysky)

shinyUI(fluidPage(
  titlePanel("1000 Genomes Populations Analysis"),
  fluidRow(
    column(5, wellPanel(
      textInput("wd", "Directory:", value = getwd(), width = '100%', placeholder = getwd()),
      actionButton("set", " Set Directory", styleclass = "primary" ), #icon = icon("folder-open", lib = "font-awesome")),
      shinyalert("setdir", click.hide = FALSE, auto.close.after = 3),
      radioButtons("add", "", choices = c("Manual","ACMG","HCM"), selected = "Manual", inline = TRUE),
      textInput("genes", "List of Genes (Manual Input):", value = "", width = '100%', placeholder = "BRCA1, RB1, MYBPC3"),
      #selectInput("genes","List of Genes", choices = c("BRCA1", "BRCA2", "MYBPC3"), selected = "BRCA1", multiple = TRUE, selectize = TRUE),
      fileInput('data','Load from RData: ', accept = c('.RData')),
      actionButton("run", " Make Plots", icon = icon("bar-chart"), styleclass = "primary"),
      shinyalert("wrongfile", click.hide = FALSE, auto.close.after = 3)
      )),
    column(7,wellPanel(
           h4("Selected Genes:"),
           textOutput("track"),
           tableOutput("selected"),
           textOutput("update"),
           textOutput("est"),
           textOutput("failure"),
           shinyalert("noneselected", click.hide = FALSE, auto.close.after = 3),
           shinyalert("invalid", click.hide = FALSE, auto.close.after = 3),
           shinyalert("override", click.hide = FALSE, auto.close.after = 5)
    )),
    column(12,wellPanel(
      busyIndicator("In Progress: Please Wait", wait = 500),
      textOutput("time"),
      h2(" "),
      plotOutput("pop_plot"),
      h2(" "),
      plotOutput("line_plot"),
      h2(" "),
      selectInput("select", "Select Gene: ", choices = c("N/A"), width = '100%'),
      h2(" "),
      plotOutput("frac_plot"),
      h2(" "),
      downloadButton(outputId = "down", label = "Download Plots"),
      downloadButton(outputId = "export", label = "Export RData")
      #downloadButton(outputId = "tbl", label = "Export Table")
    ))
  )

))