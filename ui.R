library(shiny)
library(plotly)
rm(list = ls())
shinyUI(fluidPage(
  titlePanel("ANOMALOUS APPLIANCE LOCALIZATION"),
  fluidRow(
    column(3,
           helpText("Select the below mentioned items"),
           fileInput("infilepower",
                     label = "Read Energy data",
                     accept= c(
                       'text/csv',
                       'text/comma-separated-values',
                       'text/tab-separated-values',
                       'text/plain',
                       '.csv',
                       '.tsv' )
           ),
           fileInput("infileweather",
                     label = "Read Weather data",
                     accept= c(
                       'text/csv',
                       'text/comma-separated-values',
                       'text/tab-separated-values',
                       'text/plain',
                       '.csv',
                       '.tsv' )
           ),
           dateRangeInput("seldaterange",
                          label = "Date Range",
                          start = "2014-06-01",
                          end = "2014-10-30",
                          format = "yyyy-mm-dd"),
           checkboxInput('specdaterange',
                         label = "Visualize usage in above dates",
                         value = FALSE),
           dateInput("seldate",
                     label="Select Date",
                     value="2014-06-26"
                     ),
           checkboxInput('specdate',
                         label = "Visualize Usage on above Date",
                         value =FALSE),
           
           checkboxInput('localizeanomaly',
                         label = "Localize Anomalies",
                         value = FALSE),
           numericInput("anom_rate",   label = "Anomaly Detection Rate (in hours)", value = 1),
           numericInput("anom_window", label = "Anomaly window (in mins)", value = 30)
    ),
           column(9,
           fluidRow(
             column(12,
                    #textOutput("text1"),
                    plotlyOutput("lineplt1")
             )
           ), 
            fluidRow(
              column(12,
                     dataTableOutput("facetplt2")
              )
           )
           
    )
  )
)
)

