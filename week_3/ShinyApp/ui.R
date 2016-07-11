

shinyUI(fluidPage(
  titlePanel("ROC Incidentalome App"),
  fluidRow(
    column(4, wellPanel(
      h4("Gene of Interest"),
      selectInput("gene", "Select gene:", choices = HCM.panel), 
      selectInput("var_id", "Variant ID:", as.character(data.1000g.id[["ACTC1"]][,1]), 
                  selected = as.character(data.1000g.id[["ACTC1"]][1,1])),
      sliderInput("param1", "Input 1:", min = 0, max = 1, value = 0.8, step = 0.001),
      sliderInput("param2", "Input 2:",  min = 0, max = 1, value = 0.002, step = 0.001)
    )),
    
    column(8,
           textOutput("track1"),
           textOutput("track2"),
           textOutput("track3"),
           tabsetPanel(id = "change.on",
                       tabPanel("Penetrance", 
                                h2("Penetrance from 0 to 1, step = 0.05"),
                                textOutput("pen_text"),
                                plotOutput("pen_plot"),
                                #tableOutput("pen_table"),
                                plotOutput("ps_plot")
                                ),
                       tabPanel("Prevalence", 
                                h2("Prevalence from 0 to 1, step = 0.05"),
                                plotOutput("prev_plot"),
                                tableOutput("prev_table")
                                ),
                       tabPanel("Allelic Heterogeneity", 
                                h2("Allelic Heterogeneity from 0 to 1, step = 0.05"),
                                plotOutput("ah_plot"),
                                tableOutput("ah_table")
                                )
           )
    )
  )
))