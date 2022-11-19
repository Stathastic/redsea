library(ncdf4)
library(ggplot2)
library(lubridate)
library(RColorBrewer)

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
  visualize_frame_raw(value_2d)
  contour(longitude, rev(latitude), value_2d, nlevels=nlevels,
          add = TRUE, col = "brown")

}

#plot a image
visualize_frame_raw<- function(value_2d,lantitude=NULL,latidude=NULL,timestamp="",title="") { # create a function with the name my_function
  if(is.null(latitude)){
    longitude=seq(1,dim(value_2d)[1])
    latitude=seq(1,dim(value_2d)[2])


  }
  image(longitude,rev(latitude),value_2d,xlab="longitude",ylab="latitude",col = rev(brewer.pal(10,"RdBu")))
    
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
  xlab = "log(mg/m3)",
  breaks = 100
) 


print(sum_table)
}