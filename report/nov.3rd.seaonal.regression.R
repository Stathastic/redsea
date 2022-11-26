# load data
library(ncdf4)
library(ggplot2)
library(lubridate)
library(rjson)
library(rray)
library(parallel) # 载入parallel包
library(parallelly) # 载入parallel包library(parallel) # 载入parallel包
source("../helper.R")


database <- '/home/jianj0c/dataset/redsea/'

chl_data_folder  <- 'Chlorophyll/8-Days_composite/'
chl_dir_path <- paste(database,chl_data_folder,"Aqua_MODIS_chloA_8_Days_Composite.2003_2022.nc",sep="")
chl_nc_obj <- nc_open(chl_dir_path)
chl_T_array <- ncvar_get(chl_nc_obj,"chlorophyllA")

sst_data_folder  <- 'SST/Aqua_MODIS_8_days_L3m_4km_SST/'
sst_dir_path <- paste(database,sst_data_folder,"Aqua_MODIS_sstd_8day_Composite.2003_2022.nc",sep="")
sst_nc_obj <- nc_open(sst_dir_path)
sst_T_array <- ncvar_get(sst_nc_obj,"sstMasked")

longitude<-chl_nc_obj$dim[[1]]$vals
latitude <- chl_nc_obj$dim[[2]]$vals
timestamp_ch <- chl_nc_obj$dim[[3]]$vals
timestamp_sst <- sst_nc_obj$dim[[3]]$vals
## configuration
config<-fromJSON(file="../config.json")
task_paral_maximum<-parallel::detectCores()-5

t_base<-strsplit(ncatt_get(chl_nc_obj,"time","units")$value," ")
date<-ymd(t_base[[1]][3])+dseconds(timestamp_ch)
mean_over_2006<-rowMeans(chl_T_array[,,date>"2016-01-01"&date<"2016-12-31"],na.rm=TRUE,dims=2)
sst_t_base<-strsplit(ncatt_get(sst_nc_obj,"time","units")$value," ")
sst_date<-ymd(sst_t_base[[1]][3])+dseconds(timestamp_sst)

sst_T_array_filter_outlier<-sst_T_array
sst_T_array_filter_outlier[pipeline_filter_for_outlier(sst_T_array,5)]=NA


chl_T_array_filter_outlier<-chl_T_array
chl_T_array_filter_outlier[pipeline_filter_for_outlier(chl_T_array,5)]=NA



cl <- makeClusterPSOCK(
    workers=rep(c("10.68.171.157","10.68.170.226","10.68.170.220"),3), user = "jianj0c",
  ## Manual configuration of reverse SSH tunnelling
  revtunnel = TRUE
)
clusterEvalQ(cl,source("/home/jianj0c/project/redsea/helper.R"))


anonym_func<-function(matrix,date,shape,agg_on="Season"){
    xx<-matrix(unlist(matrix))
    dim(xx)<-shape    
    return(reducer_trend_seasonal(xx,date,agg_on,25))
}


chl_datas<-lapply(seq_len(dim(chl_T_array_filter_outlier)[1]),function(i) chl_T_array_filter_outlier[i,,])

ChlA_season<-parallel::parLapply(cl,chl_datas,anonym_func,date=date,agg_on="Season",shape=c(1,423,890))
save(ChlA_season,file = file.path("../tmp/ChlA_season.Rdata"))

ChlA_month<-parallel::parLapply(cl,chl_datas,anonym_func,date=date,agg_on="month",shape=c(1,423,890))
save(ChlA_month,file = file.path("../tmp/ChlA_month.Rdata"))

sst_datas<-lapply(seq_len(dim(sst_T_array_filter_outlier)[1]),function(i) sst_T_array_filter_outlier[i,,])


sst_month<-parallel::parLapply(cl,sst_datas,anonym_func,date=sst_date,agg_on="month",shape=c(1,423,887))
save(sst_month,file = file.path("../tmp/sst_month.Rdata"))

sst_season<-parallel::parLapply(cl,sst_datas,anonym_func,date=sst_date,agg_on="Season",shape=c(1,423,887))
save(sst_season,file = file.path("../tmp/sst_season.Rdata"))

stopCluster(cl)