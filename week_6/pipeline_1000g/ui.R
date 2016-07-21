#James Diao
#July 20, 2016
#Kohane Lab | HST Summer
#Pipeline UI

shinyUI(fluidPage(
  titlePanel("1000 Genomes Populations Analysis"),
  fluidRow(
    column(5, wellPanel(
      textInput("wd", "Directory:", value = getwd(), width = '100%', placeholder = getwd()),
      actionButton("set", " Set Directory", icon("folder-open")),
      radioButtons("add", "", choices = c("Manual","ACMG","HCM"), selected = "Manual", inline = TRUE),
      textInput("genes", "List of Genes (Manual Input):", value = "", width = '100%', placeholder = "BRCA1, RB1, MYBPC3"),
      actionButton("run", " Make Plot", length = "100%", icon("bar-chart"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
      )),
    column(7,wellPanel(
           h4("Selected Genes:"),
           textOutput("update"),
           h2(" "),
           tableOutput("selected")
    )),
    column(12,wellPanel(
      plotOutput("pop_plot"),
      h2(" "),
      plotOutput("line_plot"),
      h2(" "),
      textInput("select", "Select Gene: ", value = "", width = '100%', placeholder = "GLA"),
      actionButton("in.gene", " Make Gene-Specific Plot", length = "100%", icon("bar-chart")),
      h2(" "),
      plotOutput("frac_plot"),
      h2(" "),
      downloadButton(outputId = "down", label = "Download Plots"),
      downloadButton(outputId = "export", label = "Export RData"),
      downloadButton(outputId = "tbl", label = "Export Table")
    ))
  )

))