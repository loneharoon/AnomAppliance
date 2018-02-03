# MAIN FILE
# if any function is found missing please mail to the the main author of paper at: haroonr@iiitd.ac.in
# THIS FILE CONTAINS THREE FUNCTION CORRESPONDING TO THREE STEPS OF THE PAPER
library(xts)
library(data.table)
library(ggplot2)
library(gtools)

prediction_part <- function() {
  # this function stores prediction results of all homes in a folder 

  path1 <- "patht to foldeer containg energy consumption (CSVs) of homes"
  fls <- mixedsort(list.files(path1,pattern = ".csv"))
  file2 <- "path to weather file"
  source("path to support_functions.R")
  
  day_pred_result <- list()
  
  for (i in 1:length(fls)) {
    file1 <- fls[i]
    data_ob <- create_weather_power_object_fromAggDataport(path1,file1,file2)
    print("DATA RANGES ARE:")
    print(paste0("Power data,","start: ",index(first(data_ob$power_data))," end: ",index(last(data_ob$power_data))))
    print(paste0("Weather data,","start: ",index(first(data_ob$weather_data))," end: ",index(last(data_ob$weather_data))))
    merge_start_date <- as.POSIXct(strptime('2014-06-05',format = "%Y-%m-%d"))
    merge_end_date   <- as.POSIXct(strptime('2014-08-31',format = "%Y-%m-%d"))
    my_range <- paste0(merge_start_date,'/',merge_end_date)
    sampled_ob <- combine_energy_weather(data_ob,my_range)
    train_data <- sampled_ob['2014-06-05/2014-06-21']
    test_data <- sampled_ob['2014-06-22/2014-08-30']
     # USING REGRESSION MODEL
    #reg_result <- regression_procedure(train_data,test_data,hourwindow = 6)
    #print("Regression done")
     # using NEURAL NETWORKS
    neural_result <- neuralnetwork_procedure(train_data,test_data,hourwindow = 6, daywindow = 15)
    #print("N-network done")
    # NOW STORE PREDICTION RESULTS
    # NOTE I STORE ONLY RESULTS FROM NEURAL RESULTS
    comb_df <- cbind(test_data$power,neural_result)
    comb_temp <- data.frame(timestamp = index(comb_df),round(coredata(comb_df),2))
    write.csv(comb_temp,paste0("prediction_results/",file1),row.names=FALSE)
  }
}

anomaly_detection_part <- function() {
  # this function detects anomalies while comparing actual consumption data with prediction results found in stage 1 [i.e, prediction]
  source("path to support_functions.R")
  pathxx <- "path to prediction results"
  pathy <- "path to store anomaly detection results"
  fls <- mixedsort(list.files(pathxx,pattern = ".csv"))
  for (i in 1:length(fls)) {
    file1 <- fls[i]
    data_ob <- fread(paste0(pathxx,file1))
    data_xts <- xts(data_ob[,c("power","fit","lwr","upr")],as.POSIXct(strptime(data_ob$timestamp,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01"))
    test_data_orig <- data_xts$power
    neu_result  <- subset(data_xts,select = -power)
    #anomaly_status <- find_anomalous_status(test_data=test_data_orig,result=neu_result,anomaly_window = 1,anomalythreshold_len = 4)
    anomaly_status <- find_anomalous_status_neg_anomalies(test_data=test_data_orig,result=neu_result,anomaly_window = 1,anomalythreshold_len = 4)
    anom_readings_online <- anomaly_status[anomaly_status == TRUE]
    write.csv(x = fortify(anom_readings_online),file = paste0(pathy,file1),row.names = FALSE)
  }
}

anomaly_localization <- function() {
# this function localizes the anomlous appliance
  # this will need three inputs
  # 1. ANomaly detection results file 
  # 2. Prediction file 
  # 3 appliance specifiction/ratings file
  AD_detection_results <- " anomaly detection results"
  fls <- mixedsort(list.files(AD_detection_results,pattern="online_*"))
  predict_data <- "path to prediction results"
  write_path <-  "path to store results"
  appliance_feature_file <- fread("path to appliance rating file")
  
  for (i in 1:length(fls)) {
    file1 <- fls[i]
    dfile <- strsplit(file1,'[_]')[[1]][2]
    df <- fread(paste0(predict_data,dfile))
    df_xts <- xts(subset(df,select=-timestamp),as.POSIXct(strptime(df$timestamp,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
    
    df_AD <- read.csv(paste0(AD_detection_results,file1))
    df_AD_xts <- xts(df_AD[,2],as.POSIXct(strptime(df_AD[,1],format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
    
    
    appliance_features <- appliance_feature_file[V1==dfile][,2:dim(appliance_feature_file)[2]]
    
    anom_location <- find_anomaly_appliance_with_stored_Results(pred_results = df_xts, anom_status = df_AD_xts,appliance_features,window_minutes = 5)
    write.csv(x = fortify(anom_location),file = paste0(write_path,dfile),row.names = FALSE)
    
  }
  
}