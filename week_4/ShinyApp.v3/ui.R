
shinyUI(fluidPage(
  titlePanel("Incidentalome ROC App"),
  fluidRow(
    column(4, wellPanel(
      h3("Threshold Settings:"),
      radioButtons("use.af", "Input Type:", choices = c("Threshold","Parameters"), selected = "Parameters", inline = FALSE),
      radioButtons("af.calc", "ExAC Variant Allele Frequency Calculation:", choices = c("Overall Population","Population Max"), selected = "Overall Population", inline = FALSE),
      sliderInput("threshold",  "Log Threshold:",  min = -5.1, max = 0, step = 0.1, value = log10(0.001)),
      sliderInput("penetrance", "Log Penetrance:", min = -5,   max = 0, step = 0.1, value = log10(0.8)),
      sliderInput("prevalence", "Log Prevalence:", min = -5, max = 0, step = 0.1, value = log10(0.002)),
      sliderInput("allelic.het", "Log Allelic Heterogeneity:",  min = -5, max = 0, step = 0.1, value = log10(0.001))
      )),
    column(6,wellPanel(
           h4("ExAC Allele Frequency Threshold"),
           textOutput("allele_freq"),
           h4("Penetrance"),
           textOutput("pen"),
           h4("Prevalence"),
           textOutput("prev"),
           h4("Allelic Heterogeneity"),
           textOutput("ah")
    )),
    column(6,wellPanel(
      plotOutput("roc_plot")
    ))
  )
  
))