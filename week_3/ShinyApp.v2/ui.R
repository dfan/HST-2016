log_choices <- c(seq(0, 0.01, 0.001), seq(0.01,0.1,0.01), seq(0.1,1,0.1))

shinyUI(fluidPage(
  titlePanel("ROC Incidentalome App"),
  fluidRow(
    column(4, wellPanel(
      h4("Gene/Variant of Interest:"),
      selectInput("gene", "Select gene:", choices = HCM.panel, selected = "MYBPC3"), 
      selectInput("var_id", "Variant ID:", as.character(data.total.id[["MYBPC3"]][,1]), 
                  selected = as.character(data.total.id[["MYBPC3"]][160,1])),
      selectInput("change.on", "Change on:", choices = c("Penetrance","Prevalence","Allelic Heterogeneity")), 
      selectInput("plot.density", "Plot Density:", choices = 10:100, selected = 21),
      h4("Set Values:"),
      selectInput("penetrance", "Penetrance:", choices = log_choices, selected = 0.8),
      selectInput("prevalence", "Prevalence:", choices = log_choices, selected = 0.002),
      selectInput("allelic.het", "Allelic Heterogeneity:",  choices = log_choices, selected = 0.001)
    )),
    
    column(6,wellPanel(
           #textOutput("track"),
           h3("Allele Frequencies"),
           textOutput("exac_1000g_freq"),
           h3("ROC Plot"),
           plotOutput("roc_plot")
           #h3("ROC Table"),
           #tableOutput("roc_table"),
    )),
    
    column(6,wellPanel(
           h3("Patient Status Plot"),
           plotOutput("ps_plot"),
           h3("Variant Status Plot"),
           plotOutput("genes_plot"),
           h3("Allele Frequency Threshold Plot"),
           plotOutput("threshold_plot")
  ))
  )
  
))