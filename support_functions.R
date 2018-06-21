#REPLICA OF FILE CCD_supportcode2017.R
# I.e., I have

compute_predictionband <- function(datas2_copy,predvalue,predictionband){
  # used to compute prediction band for neural networks setup
  df <- datas2_copy
  df <- df[!colnames(df)%in%c("power","temperature","humidity")]
  #mean_val <- rowMeans(df,na.rm = TRUE)
  sd_val <- apply(df,1,sd)
  dfs <- data.frame(lwr = predvalue - (predictionband * sd_val),upr = predvalue + (predictionband * sd_val))
  return(dfs)
}

compute_regression_model <- function(train_data,test_days,windowsize) {
  #browser()
  #library(MASS)
  #get train data acc to test day nature(train or test day)
  # train_data <- separate_weekend_weekday(train_data,test_days) 
  split_train <- split.xts(train_data,f="days",k=1)
  prac_day <- xts::last(split_train)
  temp <- tail(split_train,12)
  temp <- temp[1:length(temp)-1] # recent seven without last one, already used for testing
  xdat <- create_feature_matrix(temp)
  colnames(xdat) <- paste0('D',dim(xdat)[2]:1)
  #pdatsplit <- split(prac_day[[1]],lubridate::hour(index(prac_day[[1]])))
  pdatsplit <- split_hourwise(prac_day[[1]],windowsize)
  
  
  hhmodels <- list()
  for (i in 1:length(pdatsplit)) {
    #browser()
    testinterval <- pdatsplit[[i]]
    #next lines ensures that all values are not same otherwise scale function fails
    testinterval$temperature <- testinterval$temperature + rnorm(NROW(testinterval),0.01,0.01)
    testinterval$humidity <- testinterval$humidity + rnorm(NROW(testinterval),0.01,0.01)
    temp <- xdat[lubridate::hour(xdat) %in% unique(lubridate::hour(testinterval))]
    stat <- apply(temp,2,function(x) length(unique(x))==1) # remove columns which have same values -> result in NA at scale()
    temp <- temp[,!stat]
    # browser()
    print(last(index(testinterval)))
    temp_N_anom <- get_selected_nonanomalousdays(temp)
    #temp_N_anom <- temp
    datas <- cbind(coredata(temp_N_anom),coredata(testinterval))
    datas <- as.data.frame(datas)
    datas_temp <- scale(datas[,!colnames(datas)%in% c("power")]) # scaling only regressors
    datas_temp <- cbind(datas_temp,power=datas$power)# combining regressor and predictor
    modelformula <- as.formula(log(power) ~  D1 + D2 + D3 + D4 + D5 + temperature + humidity)
    # print(paste0("Rl1:",i))
    # browser()
    hhmodels[[i]] <- lm(modelformula, 
                        data = as.data.frame(datas_temp), na.action = na.exclude)
  }
  print("regression trainin done")
  #browser()
  #  NOW PREDICT FOR ACTUAL TEST DAY
  ydat <- split.xts(test_days,"days",k=1)[[1]] # data of current test day only
  temp2 <- tail(split_train,9)
  xdat2 <- create_feature_matrix(temp2)
  colnames(xdat2) <- paste0('D',dim(xdat2)[2]:1)
  
  ydatsplit <- split_hourwise(ydat,windowsize)
  pred_value <- list()
  for (i in 1:length(ydatsplit)) {
    #browser()
    testinterval2 <- ydatsplit[[i]]
    testinterval2$temperature <- testinterval2$temperature + rnorm(NROW(testinterval2),0.01,0.01)
    testinterval2$humidity <- testinterval2$humidity + rnorm(NROW(testinterval2),0.01,0.01)
    temp3 <- xdat2[lubridate::hour(xdat2) %in% unique(lubridate::hour(testinterval2))]
    stat <- apply(temp3,2,function(x) length(unique(x))==1) # remove columns which have same values -> result in NA at scale()
    temp3 <- temp3[,!stat]
    temp3_N_anom <- get_selected_nonanomalousdays(temp3) # removing anomalous days
    #temp3_N_anom <- temp3
    datas2 <- cbind(coredata(temp3_N_anom),coredata(testinterval2))
    datas2 <- as.data.frame(datas2)
    datas_temp2 <- scale(datas2[,!colnames(datas2)%in% c("power")]) # scaling only regressors
    datas_temp2 <- cbind(datas_temp2,power=datas2$power)# combining regressor and predictor
    # browser()
    #print(paste0("Rl2:",i))
    pred_value[[i]] <- predict(hhmodels[[i]], newdata = as.data.frame(datas_temp2),interval="prediction",level=0.95)
    print(paste0("prediction done",last(index(testinterval2))))
    pred_value[[i]] <- exp(pred_value[[i]])
    
  }
  pred_energy <- xts(do.call(rbind,pred_value),index(ydat))
  #pred_energy <- xts(unlist(pred_value),index(ydat))
  return(pred_energy)
}


compute_averaging_model <- function(train_data,test_days) {
  # this function predicts usage by averaging on past usage
  split_train <- split.xts(train_data,f="days",k=1)
  split_train <- tail(split_train,6)
  xdat <- create_feature_matrix(split_train)
  colnames(xdat) <- paste0('D',dim(xdat)[2]:1)
  
  ydat <- split.xts(test_days,"days",k=1)[[1]] # data of current test day only
  ydatsplit <- split(ydat,lubridate::hour(index(ydat)))
  pred_value <- list()
  for (i in 1:length(ydatsplit)) {
    # browser()
    testinterval2 <- ydatsplit[[i]]
    testinterval2$temperature <- testinterval2$temperature + rnorm(NROW(testinterval2),0.01,0.01)
    temp3 <- xdat[lubridate::hour(xdat) %in% unique(lubridate::hour(testinterval2))]
    pred_value[[i]] <- rowMeans(temp3,na.rm = TRUE)
  }
  pred_energy <- xts(unlist(pred_value),index(ydat))
  return(pred_energy)
}

compute_recurrent_Neural_model <- function(train_data,test_days,windowsize,trainingdays) {
  #browser()
  #get train data acc to test day nature(train or test day)
  # train_data <- separate_weekend_weekday(train_data,test_days) 
  library(rnn)
  split_train <- split.xts(train_data,f="days",k=1)
  prac_day <- xts::last(split_train)
  temp <- tail(split_train,trainingdays)
  temp <- temp[1:length(temp)-1] # recent seven without last one, already used for testing
  xdat <- create_feature_matrix(temp)
  colnames(xdat) <- paste0('D',dim(xdat)[2]:1)
  #pdatsplit <- split(prac_day[[1]],lubridate::hour(index(prac_day[[1]])))
  pdatsplit <- split_hourwise(prac_day[[1]],windowsize)
  
  hhmodels <- list()
  sel_columns <- c(paste0('D',1:4),'temperature','humidity')
  for (i in 1:length(pdatsplit)) {
    #browser()
    testinterval <- pdatsplit[[i]]
    #next lines ensures that all values are not same otherwise scale function fails
    testinterval$temperature <- testinterval$temperature + rnorm(NROW(testinterval),0.01,0.01)
    testinterval$humidity <- testinterval$humidity + rnorm(NROW(testinterval),0.01,0.01)
    temp <- xdat[lubridate::hour(xdat) %in% unique(lubridate::hour(testinterval))]
    # browser()
    temp_N_anom <- get_selected_nonanomalousdays(temp)
    datas <- cbind(coredata(temp_N_anom),coredata(testinterval))
    datas <- as.data.frame(datas)
    datas_temp <- datas[,!colnames(datas)%in% c("power")] # scaling only regressors
    #datas_temp <- scale(datas[,!colnames(datas)%in% c("power")]) # scaling only regressors
    
    #browser()
    datas_temp <- datas_temp[,sel_columns]
    #datas_temp <- t(datas_temp)
    #response <- t(datas$power)
    #response <- datas$power
    browser()
    tx <- array(datas_temp,dim=c(NROW(datas_temp),1,NCOL(datas_temp)))
    tx <- array(datas_temp,dim=c(dim(datas_temp),NCOL(datas_temp)))
    ty <- array(response,dim=c(NROW(response),1,NCOL(response)))
    
    hhmodels[[i]] <- trainr(X=tx,Y=ty,learningrate = 0.05, hidden_dim = 16, numepochs = 1000)
  }
  
  #  NOW PREDICT FOR ACTUAL TEST DAY
  ydat <- split.xts(test_days,"days",k=1)[[1]] # data of current test day only
  temp2 <- tail(split_train,trainingdays)
  xdat2 <- create_feature_matrix(temp2)
  colnames(xdat2) <- paste0('D',dim(xdat2)[2]:1)
  ydatsplit <- split_hourwise(ydat,windowsize)
  
  pred_value <- list()
  for (i in 1:length(ydatsplit)) {
    #browser()
    testinterval2 <- ydatsplit[[i]]
    testinterval2$temperature <- testinterval2$temperature + rnorm(NROW(testinterval2),0.01,0.01)
    testinterval2$humidity <- testinterval2$humidity + rnorm(NROW(testinterval2),0.01,0.01)
    temp3 <- xdat2[lubridate::hour(xdat2) %in% unique(lubridate::hour(testinterval2))]
    temp3_N_anom <- get_selected_nonanomalousdays(temp3) # removing anomalous days
    datas2 <- cbind(coredata(temp3_N_anom),coredata(testinterval2))
    datas2 <- as.data.frame(datas2)
    #datas_temp2 <- scale(datas2[,!colnames(datas2)%in% c("power")]) # scaling only regressors
    datas_temp2 <- datas2[,!colnames(datas2)%in% c("power")] # scaling only regressors
    temp <- datas_temp2[,sel_columns]
    temp <- t(temp)
    #response <- t(datas$power)
    tz <- array(temp,dim=c(NROW(temp),1,NCOL(temp)))
    pred_value[[i]] <- predictr(hhmodels[[i]],tz)
    
  }
  pred_energy <- xts(unlist(pred_value),index(ydat))
  return(pred_energy)
}

compute_LOF_scores <- function(train_data,test_days,windowsize,trainingdays) {
  # for each window (hour or complete day), this function calls LOF and computes anomaly score
  ydat <- split.xts(test_days,"days",k=1)[[1]] # data of current test day only
  split_train <- split.xts(train_data,f="days",k=1)
  temp2 <- tail(split_train,trainingdays)
  xdat2 <- create_feature_matrix(temp2)
  colnames(xdat2) <- paste0('D',dim(xdat2)[2]:1)
  ydatsplit <- split_hourwise(ydat,windowsize)
  
  lof_score_xts <- list()
  for (i in 1:length(ydatsplit)) {
    #browser()
    testinterval2 <- ydatsplit[[i]]
    temp3 <- xdat2[lubridate::hour(xdat2) %in% unique(lubridate::hour(testinterval2))]
    temp3_update <- cbind(temp3,coredata(testinterval2$power))
    lof_score <- return_LOF_anomaly_scores(temp3_update)
    lof_score_xts[[i]] <- xts(lof_score,last(index(testinterval2)))
  }
  combine_scores <- do.call(rbind,lof_score_xts)
  return(combine_scores)
}

compute_Neural_model <- function(train_data,test_days,windowsize,trainingdays) {
  #browser()
  #get train data acc to test day nature(train or test day)
  # train_data <- separate_weekend_weekday(train_data,test_days) 
  library(caret)
  split_train <- split.xts(train_data,f="days",k = 1)
  prac_day <- xts::last(split_train)
  temp <- tail(split_train,trainingdays)
  temp <- temp[1:length(temp)-1] # recent seven without last one, already used for testing
  xdat <- create_feature_matrix(temp)
  colnames(xdat) <- paste0('D',dim(xdat)[2]:1)
  #pdatsplit <- split(prac_day[[1]],lubridate::hour(index(prac_day[[1]])))
  pdatsplit <- split_hourwise(prac_day[[1]],windowsize)
  
  #browser()
  hhmodels <- list()
  for (i in 1:length(pdatsplit)) {
    #browser()
    testinterval <- pdatsplit[[i]]
    #next lines ensures that all values are not same otherwise scale function fails
    testinterval$temperature <- testinterval$temperature + rnorm(NROW(testinterval),0.01,0.01)
    testinterval$humidity <- testinterval$humidity + rnorm(NROW(testinterval),0.01,0.01)
    temp <- xdat[lubridate::hour(xdat) %in% unique(lubridate::hour(testinterval))]
    temp <- subset(temp,select = -D1) # temporary fix. otherwise code creates sometimes problems
    # browser()
    stat <- apply(temp,2,function(x) length(unique(x))==1) # remove columns which have same values -> result in NA at scale()
    temp <- temp[,!stat]
    temp_N_anom <- get_selected_nonanomalousdays(temp)
    temp_N_anom <- temp_N_anom[,(dim(temp_N_anom)[2]):(dim(temp_N_anom)[2] - 5)] # USING ONLY 5 DAYS FOR TRAINING
    datas <- cbind(coredata(temp_N_anom),coredata(testinterval))
    datas <- as.data.frame(datas)
    datas_temp <- scale(datas[,!colnames(datas)%in% c("power")]) # scaling only regressors
    # datas_temp <- cbind(datas_temp,power=datas$power)# combining regressor and predictor
    #modelformula <- as.formula(log(power) ~  D1 + D2 + D3 +D4 + temperature + humidity)
    hhmodels[[i]] <- avNNet(x = datas_temp, y = datas$power, size = 10, decay = 0.05, 
                            linout = TRUE, maxit = 500)
  }
  print("NN training done")
  
  #  NOW PREDICT FOR ACTUAL TEST DAY
  ydat <- split.xts(test_days,"days",k=1)[[1]] # data of current test day only
  temp2 <- tail(split_train,trainingdays)
  xdat2 <- create_feature_matrix(temp2)
  colnames(xdat2) <- paste0('D',dim(xdat2)[2]:1)
  
  #ydatsplit <- split(ydat,lubridate::hour(index(ydat)))
  ydatsplit <- split_hourwise(ydat,windowsize)
  #browser()
  pred_value <- list()
  for (i in 1:length(ydatsplit)) {
    #browser()
    testinterval2 <- ydatsplit[[i]]
    testinterval2$temperature <- testinterval2$temperature + rnorm(NROW(testinterval2),0.01,0.01)
    testinterval2$humidity <- testinterval2$humidity + rnorm(NROW(testinterval2),0.01,0.01)
    temp3 <- xdat2[lubridate::hour(xdat2) %in% unique(lubridate::hour(testinterval2))]
    stat2 <- apply(temp3,2,function(x) length(unique(x))==1) # remove columns which have same values -> result in NA at scale()
    temp3 <- temp3[,!stat2]
    temp3_N_anom <- get_selected_nonanomalousdays(temp3) # removing anomalous days
    temp3_N_anom <- temp3_N_anom[,(dim(temp3_N_anom)[2]):(dim(temp3_N_anom)[2] - 5)] #USING ONLY 5 DAYS FOR TRAINING
    datas2 <- cbind(coredata(temp3_N_anom),coredata(testinterval2))
    datas2 <- as.data.frame(datas2)
    datas2_copy <- datas2  # replica used afterwards
    datas_temp2 <- scale(datas2[,!colnames(datas2)%in% c("power")]) # scaling only regressors
    # datas_temp2 <- cbind(datas_temp2,power=datas2$power)# combining regressor and predictor
    # browser()
    print(last(index(testinterval2)))
    pred_value[[i]] <- predict(hhmodels[[i]], newdata = datas_temp2)
    pred_value[[i]] <- ifelse(pred_value[[i]]<0,100+rnorm(1,2,2),pred_value[[i]])
    bands <- compute_predictionband(datas2_copy,pred_value[[i]],2) #computing prediction bands
    pred_value[[i]] <- cbind(fit=pred_value[[i]],bands)
    
  }
  pred_energy <- xts(do.call(rbind,pred_value),index(ydat))
  #pred_energy <- xts(unlist(pred_value),index(ydat))
  return(pred_energy)
}

create_weather_power_object<- function(path1,file1,file2) {
  
  #path <- "/Volumes/MacintoshHD2/Users/haroonr/Detailed_datasets/"
  #file1 = "ap601_twomonths.csv"
  #file2 = "delhi_weatherData.csv"
  powerdata<- fread(paste0(path1,file1),header=TRUE,sep=",")
  #xx<- read.csv(paste0(path,file1),header=TRUE,sep=",")
  weatherdata <- fread(file2,header=TRUE,sep=",")
  
  power_data <- xts(powerdata$power,as.POSIXct(strptime(powerdata$timestamp,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
  names(power_data)="power"
  weather_data<- xts(data.frame(temperature=weatherdata$TemperatureC,humidity=as.numeric(weatherdata$Humidity)),order.by = as.POSIXct(weatherdata$timestamp,origin="1970-01-01"))
  
  return(list(power_data=power_data,weather_data=weather_data)) 
} 

create_weather_power_object_fromAggDataport <- function(path1,file1,file2) {
  
  powerdata<- fread(paste0(path1,file1),header=TRUE,sep=",")
  #xx<- read.csv(paste0(path,file1),header=TRUE,sep=",")
  weatherdata <- fread(file2,header=TRUE,sep=",")
  
  power_data <- xts(powerdata$use,as.POSIXct(strptime(powerdata$localminute,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
  names(power_data)="power"
  weather_data<- xts(data.frame(temperature=weatherdata$TemperatureC,humidity=as.numeric(weatherdata$Humidity)),order.by = as.POSIXct(weatherdata$timestamp,origin="1970-01-01"))
  
  return(list(power_data=power_data,weather_data=weather_data)) 
}

create_weather_power_object_from_ECO_dataset <- function(path1,file1,file2) {
  
  powerdata<- fread(paste0(path1,file1),header=TRUE,sep=",")
  #xx<- read.csv(paste0(path,file1),header=TRUE,sep=",")
  weatherdata <- fread(file2,header=TRUE,sep=",")
  
  power_data <- xts(powerdata$total,as.POSIXct(strptime(powerdata$localminute,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
  names(power_data)="power"
  weather_data<- xts(data.frame(temperature=weatherdata$TemperatureC,humidity=as.numeric(weatherdata$Humidity)),order.by = as.POSIXct(weatherdata$timestamp,origin="1970-01-01"))
  
  return(list(power_data=power_data,weather_data=weather_data)) 
}

create_weather_power_object_from_REFIT_dataset <- function(path1,file1,file2) {
  
  powerdata<- fread(paste0(path1,file1),header=TRUE,sep=",")
  #xx<- read.csv(paste0(path,file1),header=TRUE,sep=",")
  weatherdata <- fread(file2,header=TRUE,sep=",")
  
  power_data <- xts(powerdata$Aggregate,as.POSIXct(strptime(powerdata$localminute,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
  names(power_data)="power"
  
  timerange <- seq( first(index(power_data)), last(index(power_data)), by ="10 mins")
  temp = xts(rep(NA,length(timerange)),timerange)
  complete_xts = merge(power_data,temp)[,1]
  power_data_complete <- na.approx(complete_xts)
  #browser()
  weather_data<- xts(data.frame(temperature=weatherdata$TemperatureC,humidity=as.numeric(weatherdata$Humidity)),order.by = as.POSIXct(weatherdata$timestamp,origin="1970-01-01"))
  
  return(list(power_data = power_data_complete, weather_data = weather_data)) 
}


create_feature_matrix <- function(tdata) {
  #daydata <- split.xts(tdata,f="days",k=1)
  #sapply(daydata)
  matdata <- as.data.frame(sapply(tdata,function(x) coredata(x$power)))
  mat_xts <- xts(matdata,index(tdata[[1]]))
  return(mat_xts)
}


confirm_validity <- function(data_ob,merge_start_date,merge_end_date){
  if( (merge_start_date < start(data_ob$power_data)) || (merge_end_date > end(data_ob$power_data) ))
    stop("Power data within this range does't exists!!")
  if( (merge_start_date < start(data_ob$weather_data)) || (merge_end_date > end(data_ob$weather_data) ))
    stop("weather data within thisrange does't exists!!")
}

combine_energy_weather <- function(data_ob,my_range){
  # this functions simply cbinds energy and weather data
  power_data <- data_ob$power_data
  weather_data <- data_ob$weather_data
  subpower_data <- power_data[my_range]
  subweather_data <- weather_data[my_range]
  merge_data <- merge(subpower_data,subweather_data)
  no_na_data <- round(na.approx(merge_data),2) # remove NA
  return(no_na_data)
}

combine_energy_weather_ECOdata <- function(data_ob,my_range){
  # this functions simply cbinds energy and weather data
  power_data <- data_ob$power_data
  weather_data <- data_ob$weather_data
  subpower_data <- power_data[my_range]
  subweather_data <- weather_data[my_range]
  merge_data <- merge(subpower_data,subweather_data)
  merge_data <-  merge_data[!is.na(merge_data$power)]# EXPLICITLY DONE TO HANDLE MISSING POWER VALUES OF HOUSE # 5
  no_na_data <- round(na.approx(merge_data),2) # remove NA
  return(no_na_data)
}

compute_prediction_accuracy_onDatasets <- function() {
  # FOR DATAPORT
  paths <- "/Volumes/MacintoshHD2/Users/haroonr/Dropbox/R_codesDirectory/R_Codes/KDD2017/results/"
  file <- "prediction_result_dataport.csv"
  df <- fread(paste0(paths,file))
  nnet <- df[grep("NNet",df$V1),]
  regres <- df[grep("Regression",df$V1),]
  apply(regres[,2:4],2,mean)
  apply(nnet[,2:4],2,mean)
  
  # FOR AMPDS
  file <- "prediction_result_dataport.csv"
  df <- fread(paste0(paths,file))
  nnet <- df[grep("NNet",df$V1),]
  regres <- df[grep("Regression",df$V1),]
  apply(regres[,2:4],2,mean)
  apply(nnet[,2:4],2,mean)
  
  # FOR REFIT
  file <- "prediction_result_REFIT_dataset.csv"
  df <- fread(paste0(paths,file))
  nnet <- df[grep("NNet",df$V1),]
  regres <- df[grep("Regression",df$V1),]
  apply(regres[,2:4],2,mean)
  apply(nnet[,2:4],2,mean)
  
  # FOR ECO
  file <- "prediction_result_eco_dataset.csv"
  df <- fread(paste0(paths,file))
  nnet <- df[grep("NNet",df$V1),]
  regres <- df[grep("Regression",df$V1),]
  apply(regres[,2:4],2,mean,na.rm=TRUE)
  apply(nnet[,2:4],2,mean,na.rm=TRUE)
  
}

create_appliance_rating_file <- function(){
  # this stores appliance ratings of dataprot houses in a single file. Next step: Please plot each home separtely and edit this  file acc. to VI 
  setwd("/Volumes/MacintoshHD2/Users/haroonr/Dropbox/R_codesDirectory/R_Codes/KDD2017/results/")
  path1 <- "/Volumes/MacintoshHD2/Users/haroonr/Detailed_datasets/Dataport/without_car/9appliances/"
  fls <- mixedsort(list.files(path1,pattern = ".csv"))
  #file2 <- "/Volumes/MacintoshHD2/Users/haroonr/Detailed_datasets/Dataport/weather/Austin2014/10minutely_Austinweather.csv"
  source("/Volumes/MacintoshHD2/Users/haroonr/Dropbox/R_codesDirectory/R_Codes/KDD2017/CCD_supportcode2017.R")
  
  home_appliances <- list()
  
  #for (i in 1:length(fls)) {
  for (i in 1:length(fls)) {
    file1 <- fls[i]
    #data_ob <- create_weather_power_object_fromAggDataport(path1,file1,file2)
    #appliance_features <- get_appliance_features(path1, file1)
    df <- fread(paste0(path1,file1))
    df_xts <- xts(df[,2:dim(df)[2]],as.POSIXct(strptime(df$localminute,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
    appliances <- subset(df_xts,select=-use)
    app_rating <- sort(apply(appliances,2,max),decreasing = TRUE)
    home_appliances[[file1]] <- t(as.data.frame(app_rating))
    
  }
  
  home_applinace_df <- do.call(gtools::smartbind,home_appliances)
  write.csv(home_app2,file="dport_applianceRatings.csv")
  
}

metric_SMAPE <- function(object){
  numer <- sum(abs(object$fit-object$actual))
  denom <- sum(abs(object$fit)+abs(object$actual))
  smape <- numer/denom
  return (smape)
}

metric_MASE <- function(forecast_ob){
  # https://en.wikipedia.org/wiki/Mean_absolute_scaled_error 
  numerator <- abs(forecast_ob$actual - forecast_ob$fit)
  naive_error <- 0
  multiplier <- NROW(forecast_ob)/NROW(forecast_ob-1)
  for(i in 2:NROW(forecast_ob)){
    naive_difference <- abs(coredata(forecast_ob$actual[i]) - coredata(forecast_ob$actual[i-1]))
    naive_error <- naive_error + naive_difference
  }
  metric_val <- as.numeric(sum(numerator)/(multiplier*naive_error))
  return(metric_val)
}

metric_RMSE <- function(object) {
  error <- object$fit - object$actual
  return(sqrt(mean(error^2)))
}



get_selected_nonanomalousdays <- function(df_matrix) {
  
  # this function takes input matrix of energy consumption of n days and outputs selected days non anomalous days
  # browser()
  library(Rlof) # FOR LOF
  library(HighDimOut)
  df_matrix <- t(df_matrix)# transpose as it works row wise
  dist_matrix <- dist(df_matrix, method = "euclidean", diag = TRUE, upper = TRUE)
  fit <- cmdscale(dist_matrix, eig = TRUE, k = 2)
  daymat <- fit$points
  #print(row.names(daymat)[1:2])
  k <- NROW(daymat)
  dfLof <- lof(daymat, c((k-2):(k-1)), cores = 2)
  #dfLof <- lof(daymat, c(4:5), cores = 2)
  dfLofNormalized <- apply(dfLof,2,function(x) Func.trans(x,method = "FBOD"))
  FinalAnomScore <-  apply(dfLofNormalized,1,function(x) round(max(x,na.rm = TRUE),2) )#feature bagging for outlier detection
  FinalAnomScore <- data.frame(index = rownames(daymat), score = FinalAnomScore)
  FinalAnomScore$status <- FinalAnomScore$score >= 0.80
  AnomDays <- as.character(FinalAnomScore[FinalAnomScore$score >= 0.80,]$index)
  mats <- data.frame(df_matrix) # converting matrix to data frame
  
  if(length(AnomDays)!= 0){
    dbSubset <- mats[-which(rownames(mats) %in% AnomDays),] # remove anomaloous days from further operations
  }else{
    dbSubset <- mats
  }
  dbSubset <- t(dbSubset)
  colnames(dbSubset) <- paste0('D',dim(dbSubset)[2]:1)
  return(dbSubset)
}


split_hourwise <- function(tempday,windowsize) {
  #http://stackoverflow.com/a/41765212/3317829
  #tempday <- prac_day[[1]]
  .index(tempday) <- .index(tempday) - 60*30
  tempx <-split.xts(tempday,"hours",k=windowsize)
  tt <-lapply(tempx,function(x) {.index(x) <- .index(x) + 60*30;  return (x)})
  return(tt)
}

separate_weekend_weekday<- function(train_data,test_days) {
  # this function splits traindata into weekday or weekend according to the test day. if testday is weekend then only weekend training days are used 
  ydat <- split.xts(test_days,"days",k=1)[[1]]
  daytype <- unique(lubridate::wday(ydat)) #(1-> sunday, 7-> saturday)
  train_data$wday <- lubridate::wday(train_data)
  if(daytype == 1 | daytype == 7) { # saturday or sunday
    newtraindata <- train_data[train_data$wday %in% c(1,7)]
  } else {
    newtraindata <- train_data[!train_data$wday %in% c(1,7)]
  }
  return(newtraindata[,-4]) # removing wday column
}



regression_procedure <- function(train_data,test_data,hourwindow){
  days <- as.numeric( last(as.Date(index(test_data))) - first(as.Date(index(test_data))) )
  #days <- length(unique(lubridate::day(index(test_data))))
  result <- list()
  for (i in 1:days) {
    # third parameter represents the length(hours) of window
    result[[i]] <- compute_regression_model(train_data,test_data,hourwindow)
    testsplit <- split.xts(test_data,"days",k=1)
    train_data <- rbind(train_data,testsplit[[1]]) # update train data
    test_data <- do.call(rbind, testsplit[2:length(testsplit)])# update test data
  }
  finresult <- do.call(rbind,result)
  return(finresult)
}

neuralnetwork_procedure <- function(train_data,test_data,hourwindow,daywindow){
  #days <- length(unique(lubridate::day(index(test_data))))
  days <- as.numeric( last(as.Date(index(test_data))) - first(as.Date(index(test_data))) )
  result <- list()
  for (i in 1:days) {
    # browser()
    result[[i]] <- compute_Neural_model(train_data,test_data,hourwindow,daywindow)
    testsplit <- split.xts(test_data,"days",k=1)
    train_data <- rbind(train_data,testsplit[[1]]) # update train data
    test_data <- do.call(rbind, testsplit[2:length(testsplit)])# update test data
  }
  #browser()
  finresult <- do.call(rbind,result)
  return(finresult)
}

plot_graph <- function(test_data,reg_result,neural_result) {
  # takes input data from regression and neural network predictions
  library(plotly)
  res <- cbind(actual = test_data$power, regression = reg_result, neural = neural_result)
  long <- reshape2::melt(fortify(res),id.vars="Index")
  p <- ggplot(long,aes(Index,value,col=variable)) + geom_line()
  ggplotly(p)
}

plot_withhistoricaldata_graph <- function(train_data, test_data,reg_result=NULL,neural_result=NULL) {
  # input: data from regression and neural network predictions
  # input: train data(historical data), and test data of which we have predictions
  # output: ggplot with three lines
  library(plotly)
  # setting training data of some historical days for which we don't have prediction
  traindat <- split.xts(train_data,f="days",k=1)
  train_temp <- do.call(rbind,tail(traindat,5))
  if(is.null(reg_result)){
    reg_result <- test_data}
  if(is.null(neural_result)) {
    neural_result = test_data$power
  }
  train <- cbind(actual = train_temp$power, regression = rep(NA,NROW(train_temp)), neural = rep(NA,NROW(train_temp))) 
  
  # setting predicted values
  res <- cbind(actual = test_data$power, regression = reg_result, neural = neural_result)
  res <- rbind(train,res) # append res to train
  long <- reshape2::melt(fortify(res),id.vars="Index")
  p <- ggplot(long,aes(Index,value,col=variable)) + geom_line()
  ggplotly(p)
}

plot_table <- function(test_data, reg_result, neural_result, intent){
  # this function lists all the prediction results
  #1 Regression  result
  # intent shows whethe we want to print or we want to return
  resobject <- cbind(reg_result,test_data$power)
  colnames(resobject) <- c("fit","actual")
  reg_list <- list(smape = metric_SMAPE(resobject), rmse = metric_RMSE(resobject), mase = metric_MASE(resobject))
  #2 Neural network result
  resobject <- cbind(neural_result,test_data$power)
  colnames(resobject) <- c("fit","actual")
  neural_list <- list(smape = metric_SMAPE(resobject), rmse = metric_RMSE(resobject), mase = metric_MASE(resobject))
  # Combine results
  #browser()
  tab <- rbind(unlist(reg_list),unlist(neural_list))
  rownames(tab) <- c("Regression","NNet")
  if(intent=="print"){
    print(round(tab,3))
  } else{
    return(round(tab,3))
  }
}


plot_singlepredictor_graph <- function(train_data, test_data,pred_result=NULL) {
  # takes input data from regression and neural network predictions
  library(plotly)
  # setting training data of some historical days
  # browser()
  traindat <- split.xts(train_data,f="days",k=1)
  train_temp <- do.call(rbind,tail(traindat,5))
  train <- cbind(actual = train_temp$power, prediction = rep(NA,NROW(train_temp))) 
  
  # setting predicted values
  res <- cbind(actual = test_data$power, prediction = pred_result$fit)
  res <- rbind(train,res) # append res to train
  long <- reshape2::melt(fortify(res),id.vars="Index")
  p <- ggplot(long,aes(Index,value,col=variable)) + geom_line()
  p <- p + scale_x_datetime(date_breaks="1 day")
  ggplotly(p)
}


recurrent_neural_procedure <- function(train_data,test_data,hourwindow,daywindow){
  days <- length(unique(lubridate::day(index(test_data))))
  result <- list()
  for (i in 1:days) {
    result[[i]] <- compute_recurrent_Neural_model(train_data,test_data,hourwindow,daywindow)
    testsplit <- split.xts(test_data,"days",k=1)
    train_data <- rbind(train_data,testsplit[[1]]) # update train data
    test_data <- do.call(rbind, testsplit[2:length(testsplit)])# update test data
  }
  finresult <- do.call(rbind,result)
  return(finresult)
}



find_anomalous_status <- function(test_data,result,anomaly_window,anomalythreshold_len){
  # this checks if the anomaly_window (in terms of hours) contain anomaly of specific time duration (defined by anomalythreshold_len (minutes 4 mean 40))
  #browser()
  suppressWarnings(rm("status"))
  df_comp  <- cbind(result,test_data)
  df_days <- split.xts(df_comp,f="days",k=1) # daywise dataframe
  pastanomaly_vec <- as.vector("numeric") # keeps track of possible anomalies in last window
  # Loop i: Divides in terms of days
  for( i in 1:length(df_days)) {
    # Loop j: Divides in terms of anomaly window [one hour or 2 hour etc]
    day_hour <- split_hourwise(df_days[[i]],windowsize = anomaly_window)
    for( j in 1:length(day_hour)) {
      decision_vector <- as.vector(ifelse(day_hour[[j]]$power > day_hour[[j]]$upr,1,0))
      decision_vector <- c(pastanomaly_vec,decision_vector)
      length_list <- rle(decision_vector) # run length encoding
      is_anom <- any(length_list$lengths >= anomalythreshold_len & length_list$values == 1) #http://stackover
      #Logic to store only meaningful information from this window for the subsequent future window
      cnt <- 0; lcs <- 0 #last common sequence with status as 1
      last <- decision_vector[length(decision_vector) - cnt]
      while(last == 1){
        lcs <- lcs + 1
        cnt <- cnt +1
        if(lcs >= anomalythreshold_len)
          break
        last <- decision_vector[length(decision_vector)-cnt]
      }
      pastanomaly_vec <- rep(1,lcs)
      if(length(pastanomaly_vec) >= anomalythreshold_len)
        pastanomaly_vec <- as.numeric() # this sequence is already detected in current
      # create anomaly status frame
      if(!exists("status")){
        status <- xts(is_anom,index(last(day_hour[[j]])))
      } else {
        temp <- xts(is_anom,index(last(day_hour[[j]])))
        status <- rbind(status,temp)
      }
    }
  }
  return(status)
} 


get_appliance_features <- function (path1,file1) {
  # this function will return both appliane rating and the the dissagrigated energy of a home
  powerdata <- fread(paste0(path1,file1),header=TRUE,sep=",")
  powerdata_xts <- xts(powerdata[,2:dim(powerdata)[2]],as.POSIXct(strptime(powerdata$localminute,format="%Y-%m-%d %H:%M:%S"),origin="1970-01-01")) 
  
  appliance_data_xts <- subset(powerdata_xts,select = -c(use))
  app_rating <- sort(apply(appliance_data_xts,2,max),decreasing = TRUE)
  
  return(list(home_data = powerdata_xts, appliance_rating = app_rating))
  
}


find_anomaly_appliance_from_livedata <- function(test_data, result, anom_status, appliance_features, window_minutes) {
  # appliance_features <- home_details$appliance_rating
  # anom_status <- res_neu
  #  result <- neural_result
  day_dat <- split.xts(anom_status,"days",k=1)
  suppressWarnings(rm("final_result"))
  # find anomalies in different intervals
  for (i in 1:length(day_dat)) {
    df <- day_dat[[i]]
    # times of day at which anomaly was found
    anom_intervals <- index(df[df %in% c("TRUE"),] )
    if(length(anom_intervals) != 0) { # if there are anomalies
      for (j in 1:length(anom_intervals)) { # for each anomaly interval
        # finding range of anomaly detection time interval
        upr_lmt <- index(test_data[anom_intervals[j],])
        lwr_lmt <- upr_lmt - 10  * 60 * window_minutes # last one hour
        actual <- test_data[paste0(lwr_lmt,'/',upr_lmt)]$power # actual power consumed
        baseline <- result[paste0(lwr_lmt,'/',upr_lmt)]$fit   # predicted power consumed
        # we know this interval is anomalous so actual will be much higher than the baseline
        delta <- baseline - actual
        magnitude <- abs(min(delta))
        # find the applaince to which negative numbers match most
        anom_app <- which(abs(appliance_features - magnitude) == min(abs(appliance_features - magnitude)))
        if(!exists("final_result")){
          final_result <- xts(names(anom_app),upr_lmt)
        } else{
          temp <- xts(names(anom_app),upr_lmt)
          final_result <- rbind(final_result,temp)
        }
      }
    }
  }
  return(final_result)
}

find_anomaly_appliance_with_stored_Results <- function(pred_results, anom_status, appliance_features, window_minutes) {
  # this function returns the name of possible anomalou applaince
  day_dat <- split.xts(anom_status,"days",k=1)
  suppressWarnings(rm("final_result"))
  # find anomalies in different intervals
  for (i in 1:length(day_dat)) {
    df <- day_dat[[i]]
    # times of day at which anomaly was found
    anom_intervals <- index(df[df %in% c("TRUE"),] )
    if(length(anom_intervals) != 0) { # if there are anomalies
      for (j in 1:length(anom_intervals)) { # for each anomaly interval
        # finding range of anomaly detection time interval
        upr_lmt <- index(pred_results[anom_intervals[j],])
        lwr_lmt <- upr_lmt - 10  * 60 * window_minutes # last one hour
        actual <- pred_results[paste0(lwr_lmt,'/',upr_lmt)]$power # actual power consumed
        baseline <- pred_results[paste0(lwr_lmt,'/',upr_lmt)]$fit   # predicted power consumed
        # we know this interval is anomalous so actual will be much higher than the baseline
        delta <- baseline - actual
        magnitude <- abs(min(delta))
        # find the applaince to which negative numbers match most
        temp <- abs(as.numeric(appliance_features) - magnitude)
        names(temp) <- names(appliance_features)
        anom_app <- which.min(temp)
        #anom_app <- which(abs(appliance_features - magnitude) == min(abs(appliance_features - magnitude)))
        if(!exists("final_result")){
          final_result <- xts(names(anom_app),upr_lmt)
        } else{
          temp <- xts(names(anom_app),upr_lmt)
          final_result <- rbind(final_result,temp)
        }
      }
    }
  }
  return(final_result)
}

find_regrassivedays_with_BIC <- function(regression_days) {
  # I assume that I will recieve more than 21 days of data
  # Using this function we find how many historical days should we use in regression
  # Lower BIC ensures optimal estimate 
  data <- regression_days
  # browser()
  daydat <- split.xts(data,"days")
  mat <- sapply(daydat,function(x) return(coredata(x$power)))
  df_x <- as.data.frame(mat)
  optimal_hist_day <- vector("numeric") # for each different test day, this will store how may historical should we use ideally 
  index <- 1
  for(day in 21:dim(df_x)[2]){ #  for different day of input data
    
    #index <- 1
    testday <- df_x[,day] 
    bic_vector <- vector("numeric")
    
    for (i in 1:(day-1)) { # for same day with different number of historical days
      df_temp <- data.frame(df_x[,(day-1):(day-i)])
      colnames(df_temp) <- paste0("D",dim(df_temp)[2]:1)
      df_temp$y <- testday
      mod <- lm(y ~.,df_temp)
      bic_vector[i] <- BIC(mod)
    }
    optimal_daylimit <- which(bic_vector == min(bic_vector)) # How many days I should consider for minimal bic
    optimal_hist_day[index] <- optimal_daylimit
    index <- index + 1
  }
  return(optimal_hist_day)
}

incremental_lof_procedure <- function(train_data,test_data,hourwindow,daywindow){
  days <- as.numeric( last(as.Date(index(test_data))) - first(as.Date(index(test_data))) )
  #days <- length(unique(lubridate::day(index(test_data))))
  result <- list()
  for (i in 1:days) {
    # third parameter represents the length(hours) of window
    result[[i]] <- compute_LOF_scores(train_data,test_data,hourwindow,daywindow)
    testsplit <- split.xts(test_data,"days",k=1)
    train_data <- rbind(train_data,testsplit[[1]]) # update train data
    test_data <- do.call(rbind, testsplit[2:length(testsplit)])# update test data
  }
  finresult <- do.call(rbind,result)
  return(finresult)
}


return_LOF_anomaly_scores <- function(df_matrix) {
  # this function takes input matrix of energy conusmption of n days and outputs anomaly score of last day as this runs in incremental fashion
  library(Rlof) # FOR LOF
  library(HighDimOut)
  df_matrix <- t(df_matrix)# transpose as it works row wise
  dist_matrix <- dist(df_matrix, method = "euclidean", diag = TRUE, upper = TRUE)
  fit <- cmdscale(dist_matrix, eig = TRUE, k = 2)
  daymat <- fit$points
  k <- NROW(daymat)
  dfLof <- lof(daymat, c((k-2):(k-1)), cores = 2)
  #dfLof <- lof(daymat, c(4:5), cores = 2)
  dfLofNormalized <- apply(dfLof,2,function(x) Func.trans(x,method = "FBOD"))
  FinalAnomScore <-  apply(dfLofNormalized,1,function(x) round(max(x,na.rm = TRUE),2) )#feature bagging for outlier detection
  return(last(FinalAnomScore))
}

print_single_prediction_method_table <- function(test_data, method_result){
  # this function lists all the prediction results
  #method_result <- reg_result$fit
  resobject <- cbind(method_result,test_data$power)
  colnames(resobject) <- c("fit","actual")
  reg_list <- list(smape = metric_SMAPE(resobject), rmse = metric_RMSE(resobject), mase = metric_MASE(resobject))
  #print results
  tab <- unlist(reg_list)
  print(round(tab,3))
}


averaging_Case <- function() {
  # in this approach, we predict by simple averaging of the past consumption
  train_data <- sampled_ob['2014-06-05/2014-06-25']
  test_data <- sampled_ob['2014-06-26/2014-06-28']
  backupcopy <- test_data
  days <- length(unique(lubridate::day(index(test_data))))
  for (i in 1:days) {
    # third parameter represents the length(hours) of window
    result[[i]] <- compute_averaging_model(train_data,test_data)
    #plot(result[[i]])
    testsplit <- split.xts(test_data,"days",k=1)
    train_data <- rbind(train_data,testsplit[[1]]) # update train data
    test_data <- do.call(rbind, testsplit[2:length(testsplit)])# update test data
  }
  avgresult <- do.call(rbind,result)
  #res_avg <- compute_averaging_model(train_data,test_data)
  lines(avgresult,col="green")
}



daywise_anomaly_inspection <- function(home_details, anom_loc) {
  # browser()
  library("RColorBrewer")
  # VISUALIZE PORTION OF DATA ON GIVEN DATE
  #http://novyden.blogspot.in/2013/09/how-to-expand-color-palette-with-ggplot.html
  days <- unique(format(index(anom_loc),"%Y-%m-%d"))
  data_10min <- home_details
  cat(paste0("Anomalous days:: ",length(days)))
  for (i in 1:length(days)){
    #browser()
    cat("enter number")
    x = scan(what= numeric(),n =1)
    sel_date <- days[x]
    fulldata <- data_10min[paste0(as.POSIXct(paste0(sel_date," ",'00:00:00')),"/",as.POSIXct(paste0(sel_date," ",'23:59:59')))]
    dframe <- fulldata
    dframe <- data.frame(timeindex = index(dframe), coredata(dframe))
    # dframe$dataid <- NULL ; dframe$air1 <-NULL ; dframe$use<- NULL ; dframe$drye1 <- NULL
    df_long <- reshape2::melt(dframe,id.vars = "timeindex")
    colourCount = length(unique(df_long$variable))
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount) # brewer.pal(8, "Dark2") or brewer.pal(9, "Set1")
    #cat("rached g")
    g <- ggplot(df_long,aes(timeindex,value,col=variable,group=variable))
    g <- g + geom_line() + scale_colour_manual(values=getPalette) + ggtitle(as.Date(df_long$timeindex[15]))
    g <- g + scale_x_datetime(labels = date_format("%H",tz="Asia/Kolkata"), breaks = pretty_breaks(n=24))
    plot(g)
    #cat("loop done")
  }
}

similarity_of_KDD_datasets <- function() {
  # I used this function to measure how different/irregular are these datasets
  library(cluster)
  library(gridExtra)
  library(fpc)
  library(gtools)
  # DATAPORT  
  path1 <- "/Volumes/MacintoshHD2/Users/haroonr/Detailed_datasets/Dataport/without_car/9appliances/"
  fls <- mixedsort(list.files(path1,pattern = ".csv"))
  avg_cluster <- vector("numeric")
  
  for (h in 1:length(fls)){
    df <- fread(paste0(path1,fls[h]))
    dfs <-  xts(df[,df$use],fastPOSIXct(df$localminute)-19800)
    subset1 <- dfs["2014-06-11 00:00:00/2014-08-26 23:59:59"]
    
    #week_data <- split.xts(subset1,f="weeks",k=1)
    week_data <- split(subset1,lubridate::week(index(subset1)))
    clusters <- vector("numeric")
    for (i in 1:length(week_data)){
      
      day_data <- split.xts(week_data[[i]], f = "days", k = 1)
      sub_2 <- lapply(day_data,function(x) coredata(x))
      sub_mat <- do.call(cbind,sub_2)
      colnames(sub_mat) <- paste0("c",1:dim(sub_mat)[2])
      set.seed(123)
      # browser()
      clusters[[i]] <- pamk(t(sub_mat),krange=2:5)$nc
    }
    avg_cluster[h] <- ceiling(mean(clusters))
  }
  avg_cluster
  
  
  #### AMPDS####
  path1 <- "/Volumes/MacintoshHD2/Users/haroonr/Detailed_datasets/AMPds/"
  avg_cluster <- vector("numeric")
  df <- fread(paste0(path1,"ten_minutes_data.csv"))
  dfs <-  xts(df[,df$WHE],fastPOSIXct(df$localminute)-19800)
  subset1 <- dfs["2014-01-08 00:00:00/2014-03-31 23:59:59"]
  
  #week_data <- split.xts(subset1,f="weeks",k=1)
  week_data <- split(subset1,lubridate::week(index(subset1)))
  clusters <- vector("numeric")
  for (i in 1:length(week_data)){
    
    day_data <- split.xts(week_data[[i]], f = "days", k = 1)
    sub_2 <- lapply(day_data,function(x) coredata(x))
    sub_mat <- do.call(cbind,sub_2)
    colnames(sub_mat) <- paste0("c",1:dim(sub_mat)[2])
    set.seed(123)
    clusters[[i]] <- pamk(t(sub_mat),krange=2:5)$nc
  }
  clusters
  
  #### ECO####
  path1 <- "/Volumes/MacintoshHD2/Users/haroonr/Detailed_datasets/ECO_dataset/"
  fls <- mixedsort(list.files(path1,pattern = ".csv"))
  avg_cluster <- vector("numeric")
  
  for (h in 1:length(fls)){
    df <- fread(paste0(path1,fls[h]))
    dfs <-  xts(df[,df$total],fastPOSIXct(df$localminute)-19800)
    subset1 <- dfs["2012-08-12 00:00:00/2012-10-27 23:59:59"]
    
    #week_data <- split.xts(subset1,f="weeks",k=1)
    week_data <- split(subset1,lubridate::week(index(subset1)))
    clusters <- vector("numeric")
    for (i in 1:length(week_data)){
      
      day_data <- split.xts(week_data[[i]], f = "days", k = 1)
      sub_2 <- lapply(day_data,function(x) coredata(x))
      sub_mat <- do.call(cbind,sub_2)
      colnames(sub_mat) <- paste0("c",1:dim(sub_mat)[2])
      set.seed(123)
      clusters[[i]] <-pamk(t(sub_mat),krange=2:5)$nc
    }
    avg_cluster[h] <- ceiling(mean(clusters))
  }
  avg_cluster
  
  #### REFIT####
  path1 <- "/Volumes/MacintoshHD2/Users/haroonr/Detailed_datasets/REFITT/dataset_10mins/"
  fls <- mixedsort(list.files(path1,pattern = ".csv"))
  avg_cluster <- vector("numeric")
  
  for (h in 1:length(fls)){
    df <- fread(paste0(path1,fls[h]))
    dfs <-  xts(df[,df$Aggregate],fastPOSIXct(df$localminute)-19800)
    subset1 <- dfs["2014-06-10 00:00:00/2014-08-29 23:59:59"]
    
    #week_data <- split.xts(subset1,f="weeks",k=1)
    week_data <- split(subset1,lubridate::week(index(subset1)))
    clusters <- vector("numeric")
    for (i in 1:length(week_data)){
      
      day_data <- split.xts(week_data[[i]], f = "days", k = 1)
      sub_2 <- lapply(day_data,function(x) coredata(x))
      sub_mat <- do.call(cbind,sub_2)
      colnames(sub_mat) <- paste0("c",1:dim(sub_mat)[2])
      set.seed(123)
      clusters[[i]] <- pamk(t(sub_mat),krange=2:10)$nc
    }
    avg_cluster[h] <- mean(clusters)
  }
}
