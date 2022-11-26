library(ncdf4)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(parallel) # 载入parallel包
# load data
library(INLA)
library(zoo)
library(car)
library(hydroTSM)

library(rjson)
library(rray)

# filter the outlier value with the simplist methods:
# 
pipeline_filter_for_outlier <- function(array,std_n=3) { # create a function with the name my_function
  std<-sqrt(var(array,na.rm = TRUE))
  mean<-mean(array,na.rm=TRUE)
  return((!is.na(array))&((array>mean+std_n*std)|(array<mean-std_n*std)))
}

##geom_contour
visualize_frame_level<-function(value_2d,latitude=NULL,longitude=NULL,title="",nlevels=30){

  if(is.null(latitude)){
    longitude=seq(1,dim(value_2d)[1])
    latitude=seq(1,dim(value_2d)[2])


  }
  visualize_frame_raw(value_2d,latitude,longitude,title)
  contour(longitude, rev(latitude), value_2d, nlevels=nlevels,
          add = TRUE, col = "brown")

}

#plot a image
visualize_frame_raw<- function(value_2d,latitude=NULL,longitude=NULL,timestamp="",title="") { # create a function with the name my_function
  if(is.null(latitude)){
    longitude=seq(1,dim(value_2d)[1])
    latitude=seq(1,dim(value_2d)[2])


  }
  image(longitude,rev(latitude),value_2d,xlab="longitude",ylab="latitude",col = rev(brewer.pal(10,"RdBu")))
  title(main = title, font.main = 4)
}


helper_describe<-function(array){

  valid_element_number<-length(array[!is.na(array)])
  valid_percentage<-(valid_element_number/length(array)*100)
  sum_table<-summary(array[!is.na(array)],na.rm=TRUE)
  sum_table["Valid.Ent"]<-valid_element_number
  sum_table["Valid.percent"]<-valid_percentage
  sum_table["Std"]<-sqrt(var(array,na.rm=TRUE))
  userful_entry<-array[!is.nan(array)]



  hist(userful_entry,
    breaks = 100
  ) 


  print(sum_table)
}


## linear regression for 3D and reduce to 2d trend matrix
reducer_trend<-function(T_array,dataType="Chl",thres_min_data_point){
  results<-list()
  length_latitude<-dim(T_array)[2]
  length_longitude<-dim(T_array)[1]
  length_timestamp<-dim(T_array)[3]

  for(latitude_i in 1:length_latitude) {

    for (longitude_i in 1:length_longitude){
      raw_t=seq(1,length_timestamp)
      raw_y=T_array[longitude_i,latitude_i,1:length_timestamp]

      t<-raw_t[!is.na(raw_y)]
      y<-raw_y[!is.na(raw_y)]
      if (length(t)>thres_min_data_point){
        data=data.frame(t,y)
        ## specify model based on datatype
        if (dataType=="Chl"){
          model<-gam(y~t,family=Gamma(link = "log"),data=data) 
        }
        else if(dataType=="SST"){
          model<-gam(y~t,data=data) 
        }
        else{
          stopifnot(FALSE)
        }
        results<-append(results,model$coefficients["t"])
      }
      else{
        results<-append(results,NA)
      }
    }

  }
  matrix<-unlist(results)
  dim(matrix) <- c(length_longitude, length_latitude)
  return(matrix)
}


# take one 3d arrayday reduce the 3d data to 2d
# the two tensor should have the same shape
reducer2_correlation<-function(T_array1,T_array2){
    

    stopifnot(dim(T_array1)[3]==dim(T_array2)[3])
    d3_length<-dim(T_array1)[3]
    T1_mean<-rowMeans(T_array1,na.rm=TRUE,dims=2)
    T2_mean<-rowMeans(T_array2,na.rm=TRUE,dims=2)

    nominator<-rowMeans(rray_subtract(T_array1,T1_mean)*rray_subtract(T_array2,T2_mean)
        ,na.rm=TRUE,dims=2)

    var1<-rowMeans(rray_subtract(T_array1,T1_mean)**2,na.rm=TRUE,dims=2)
    var2<-rowMeans(rray_subtract(T_array2,T2_mean)**2,na.rm=TRUE,dims=2)

    denominator<-(var1*var2)**0.5

    return(nominator/denominator)

}

##deprecated since no acceration observed 
parallel_trend_linear_regression<-function(array,dataType,thres_min_data_point){
    raw_t=seq(1,length_timestamp)
    t<-raw_t[!is.na(array)]
    y<-array[!is.na(array)]
    if (length(t)>thres_min_data_point){
        if (dataType=="Chl"){
            model<-gam(y~t|season+f(month,model="RW2"),family=Gamma(link = "log"),data=data) 
        }
        else if(dataType=="SST"){
            model<-gam(y~t,data=data) 
        }
        else{
            stopifnot(FALSE)
        }
        return(results,model$coefficients["t"])
      }
    else{
        return(NA)
        
    }

    
}
##deprecated since no acceration observed 
reducer_trend_parallel<-function(T_array,dataType="Chl",thres_min_data_point){
  length_latitude<-dim(T_array)[2]
  length_longitude<-dim(T_array)[1]
  length_timestamp<-dim(T_array)[3]

  array_list<-list()
  for(latitude_i in 1:length_latitude) {

      for (longitude_i in 1:length_longitude){
          array_list<-append(array_list,T_array[longitude_i,latitude_i,1:length_timestamp])
      }
      }

  print(length(array_list))
  result<-foreach(x=1:length(array_list),.combine = ) %do% 
      parallel_trend_linear_regression(array_list[x],dataType,thres_min_data_point)
  print(class(result))
}



visualize_frame_ggplot<- function(value_2d,latitude=NULL,longitude=NULL,timestamp="",title="") { # create a function with the name my_function
  if(is.null(latitude)){
    longitude=seq(1,dim(value_2d)[1])
    latitude=seq(1,dim(value_2d)[2])


  }
  lonlat <- as.matrix(expand.grid(longitude,latitude))
  xx<-as.vector(value_2d)
  tmp_df01 <- data.frame(cbind(lonlat,xx))
  names(tmp_df01)<- c("longitude","latitude","value")
  ggplot(tmp_df01, aes(x = longitude, y = latitude, fill = value)) +
  geom_tile()+
    scale_fill_gradient(low="green", high="red")+
  ggtitle(title)
}

feature_extract_time<-function(timeIndex){
    ## depends on hydroTSM package, which is used for time-series
    if(!require("hydroTSM")){
        install.packages("hydroTSM")
    }
    month<-strtoi(format(timeIndex, "%m"),base=10)
    season<-recode(time2season(timeIndex,out.fmt = "seasons"),"'winter'=1; 'spring'=2; 'summer'=3; 'autumm'=4;")

    return (list(month,season))
}

reducer_trend_seasonal<-function(T_array,date,agg_on="Season",thres_min_data_point=25){
  ## extract the month and season information
  

  results<-list()
  length_latitude<-dim(T_array)[2]
  length_longitude<-dim(T_array)[1]
  length_timestamp<-dim(T_array)[3]
  raw_t=seq(1,length_timestamp)
  ## create result matrix to score the coefficient,. the seasonal feature for the whole T_array, size is len(time)
  if(agg_on=="Season"){
    seasonal_feature<-unlist(feature_extract_time(date)[2])
    result_matrix<-unlist(rep(list(NA),length_latitude*length_longitude*4))
    dim(result_matrix)<-c(length_latitude,length_longitude,4)
  }
  else{
    seasonal_feature<-unlist(feature_extract_time(date)[1])
    result_matrix<-unlist(rep(list(NA),length_latitude*length_longitude*12))
    dim(result_matrix)<-c(length_latitude,length_longitude,12)

  }


  ## two layer for
  for(latitude_i in 1:length_latitude) {

    for (longitude_i in 1:length_longitude){
     
      raw_y=T_array[longitude_i,latitude_i,1:length_timestamp]

      t<-raw_t[!is.na(raw_y)]
      y<-raw_y[!is.na(raw_y)]

      seasonal<-seasonal_feature[!is.na(raw_y)]


      ## for rw2, the unique value need be larger then 2
      if (length(t)>thres_min_data_point & length(unique(seasonal))>2){

        
        data=data.frame(t,y,seasonal)



        model<-inla(y~t+f(seasonal,model="rw2"),data=data,verbose=TRUE) 
        results<-append(results,model$summary.random$season$mean)


        ## store the result
        result_df<-model$summary.random$seasonal

        result_matrix[latitude_i,longitude_i,result_df$ID]=result_df$mean


      }
    }

  }
  return(result_matrix)
}
