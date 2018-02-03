# KDD DEMO APPLICATION
library(shiny)
library(data.table)
library(fasttime) #for fastPosixct
library(xts)
library(ggplot2)
library(scales)
library(rwantshue) # for distint colors https://github.com/hoesler/rwantshue
library(RColorBrewer)# to increas
suppressWarnings(library(plotly))
library(cluster)
options(shiny.maxRequestSize=100*1024^2)
rm(list=ls())
Sys.setenv(TZ="Asia/Kolkata")
source("./Support_functions.R")
shinyServer(

  function(input, output) {
    
    dframe <- reactive( {
      # inputfile <- input$infilepower
       validate(
         need(input$infilepower != "","Please select a data set")
       )
      #browser()
      dframe <- fread(input$infilepower$datapath, header = TRUE, sep = ",")
      dframe$timestamp <- fastPOSIXct(dframe$localminute)-19800
     dframe_xts <- xts(dframe$use,dframe$timestamp)
      df_weather <- fread(input$infileweather$datapath, header = TRUE, sep = ",")
      df_weather$timestamp <- fastPOSIXct(df_weather$timestamp)-19800
      df_weather_xts <- xts(df_weather[,2:dim(df_weather)[2]],df_weather$timestamp)
      
      merge_start_date <- as.POSIXct(strptime('2014-06-05',format = "%Y-%m-%d"))
      merge_end_date   <- as.POSIXct(strptime('2014-08-31',format = "%Y-%m-%d"))
     # confirm_validity(data_ob, merge_start_date, merge_end_date)
    # browser()
      my_range <- paste0(merge_start_date,'/',merge_end_date)
      data_ob <- list(power_data = dframe_xts, weather_data = df_weather_xts)
      sampled_ob <- combine_energy_weather(data_ob,my_range)
      sampled_ob
   
    } )
    
    df <- reactive( {
      dframe <- dframe()
      #str(dframe)
      if (input$specdaterange| input$specdate){
        if(input$specdaterange) {
        start_date = input$seldaterange[1]
        end_date =input$seldaterange[2]
        startdate <- fastPOSIXct(paste0(start_date,' ',"00:00:00"))-19800
        enddate <- fastPOSIXct(paste0(end_date,' ',"23:59:59"))-19800
        } else {
          datex = input$seldate
          startdate <- fastPOSIXct(paste0(datex,' ',"00:00:00"))-19800
          enddate <- fastPOSIXct(paste0(datex,' ',"23:59:59"))-19800
        }
     #   browser()
      dframe <-  dframe[paste0(startdate,"/",enddate)]
        #dframe <- dframe[dframe$timestamp >= startdate & dframe$timestamp <= enddate,] #reduced
      }
      dfs <- dframe
      dfs
    } )
    
    output$lineplt1 <- renderPlotly({
      df_sub <- subset(df(),select=-c(TemperatureC,Humidity))
      df_sub <- fortify(df_sub)
      colnames(df_sub) <- c("timestamp","Power")
      g <- ggplot(df_sub, aes(timestamp, Power))
      g <- g + geom_line() + labs(x=" ",y ="power (Watts)")
      g <- g + scale_x_datetime(breaks = date_breaks("6 hour"), labels = date_format("%d-%b %H:%M",tz="Asia/Kolkata")) # use scales package
      g <- g + theme(axis.text.x = element_text(angle = 90,hjust = 1))
    ggplotly(g)
    })
  
    output$facetplt2 <- renderDataTable({
      if(input$localizeanomaly){
        anomaly_rate <- input$anom_rate
        anomaly_window <- input$anom_window
        
        sample_ob <- df()
        colnames(sample_ob) <- c("power","temperature","humidity")
        train_data <- sample_ob['2014-06-05/2014-06-21']
        test_data <- sample_ob['2014-06-22/2014-06-30']
       # browser()
        neural_result <- neuralnetwork_procedure(train_data,test_data,hourwindow = 6, daywindow = 15)

        res_neu <- find_anomalous_status(test_data,result=neural_result,anomaly_window = 1,anomalythreshold_len = 4)
        res_neu[res_neu==TRUE]
      
        applaince_features <- c(4000,2000,500,1200,250,300,700,400,100,300,70)
        names(applaince_features) <- c("Dryer","Range","Microwave","AC","KitchenApp1","KitchenApp2","Dishwasher","ClothesWasher","Fride","Furnace","Lights")
        
        
        anom_loc <- find_anomaly_appliance_from_livedata(test_data, result = neural_result, anom_status = res_neu, applaince_features,window_minutes = 5)  
        
        dtab <- as.data.table(anom_loc)
        colnames(dtab) <- c("Timestamp","Anomalous Appliance")
        
      dtab
      }
    })
    

    
  } )