#require
require(ggplot2)
require(Hmisc)
require(schoolmath)
library(raster)
library(rgdal)
library(sp)
library(rwrfhydro)
library(stringr)
require(dplyr)
install.packages("remotes")
remotes::install_github("USGS-R/EflowStats")

#create list of dates
datesWanted<-seq.Date(as.Date("2003-10-01"), as.Date("2020-09-30"), by="day")

#set output directory
outDir<-c("/Users/keirajohnson/Documents/SNODAS")

#scaling factors
nCol <- 6935
nRow <- 3351 #columns and rows number: masked version of contiguous US
dataScaleFactor  <-  1000  #multiply the data this amount, both depth and SWE

#set coordinate system
wgs84 <- CRS("+proj=longlat +datum=WGS84")

#read in HJA shapefile
HJA_shape<-readOGR("/Users/keirajohnson/Box Sync/Keira_Johnson/HJAshapefile/layers/globalwatershed.shp")

#open list
swe_HJA_matrix_mm<-list()
dateList<-list()

#constants
beg<-"us_ssmv11034tS__T0001TTNATS"
swe_files<-list.files(path=outDir, pattern = beg)

for (i in 6162:length(datesWanted)) {
  
  #this allows the loop to keep running even if there is an error
  tryCatch({
    
    #download daily snodas for one day
    #GetSnodasDepthSweDate(datesWanted[i], outputDir = outDir, overwrite = T, quiet=FALSE)
  
    #get correct URL based on swe code and date
    date<-gsub("-", "", datesWanted[i])
    file_name<-str_subset(swe_files, date)
    
    #extract SWE
    sweFile0  <- paste0(outDir, "/", file_name)
    
    #extract from dat.gz file
    sweCon  <- gzcon(file(sweFile0, "rb"))
    
    #read binary data from file above
    sweData <- readBin(sweCon, integer(), n=nRow*nCol, size=2, signed=TRUE, endian='big')/dataScaleFactor
    
    #remove error values
    sweData[which(sweData==-9.999)]<-NA
    
    #convert to matrix
    sweDataMatrix<-matrix(sweData, ncol=nCol, nrow=nRow, byrow = T)
    
    #convert to mm
    sweDataMatrix_mm<-sweDataMatrix*dataScaleFactor
    
    #georeference
    sweRaster <- raster(sweDataMatrix_mm, ymn=24.9504, ymx=52.8754, xmn=-124.7337, xmx=-66.9421, crs=wgs84)
    
    #extract data for HJA
    swe_HJA_mm<-crop(sweRaster, HJA_shape)
    swe_HJA_mm_mask<-mask(swe_HJA_mm, HJA_shape)
    
    #convert back to matrix
    swe_HJA_matrix_mm[[i]]<-as.matrix(swe_HJA_mm_mask)
    
    closeAllConnections()
    
    dateList[[i]]<-datesWanted[i]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

save(swe_HJA_matrix_mm, "swe_HJA_matrix_mm.RData")

dateList_rmNULL<-dateList[-which(sapply(dateList, is.null))]

swe_HJA_rmNULL<-swe_HJA_matrix_mm[-which(sapply(swe_HJA_matrix_mm, is.null))]

#combine matrices into array
swe_array<-array(as.numeric(unlist(swe_HJA_rmNULL)), dim=c(10, 19, length(swe_HJA_rmNULL)))

#set width and height
width<-seq(1,10,1)
height<-seq(1,19,1)

#open lists
swe_vector<-list()
swe_vector_final<-list()

#extract time series data by grid cell
for (i in 1:length(width)) {
  
  for (j in 1:length(height)) {
    
    #add vectors to list
    swe_vector[[j]]<-as.data.frame(swe_array[i,j, ,drop=FALSE])
    
  }
  
  #append lists to dataframe
  swe_dataframe<-data.frame(do.call(rbind, swe_vector))
  #switch rows and columns
  swe_df_t<-as.data.frame(t(swe_dataframe))
  
  #add dataframes to list
  swe_vector_final[[i]]<-swe_df_t
  
}

#combine dataframes
swe_df_total<-data.frame(do.call(cbind, swe_vector_final))

write.csv(swe_df_total, "HJASNODAS.csv")

#unlist date
datedf <- data.frame(matrix(unlist(dateList_rmNULL), nrow=length(dateList_rmNULL), byrow=TRUE))
#change from julian date to normal date
datedf$matrix.unlist.dateList_rmNULL...nrow...length.dateList_rmNULL...<-
  as.Date(datedf$matrix.unlist.dateList_rmNULL...nrow...length.dateList_rmNULL..., origin = "1970-01-01")
#rename columns 
names(datedf)<-"date"

#format date dataframe - add year, month, calculate water year
datedf$year<-as.numeric(format(as.Date(datedf$date), "%Y"))
datedf$month<-as.numeric(format(as.Date(datedf$date), "%m"))
datedf$julianday<-as.numeric(strftime(datedf$date, "%j"))
datedf$WY<-ifelse(datedf$month > 9, datedf$year+1, datedf$year)

#combine date df and swe df
swe_tot_date<-cbind(datedf, swe_df_total)

write.csv(swe_tot_date, "HJASNODAS_date.csv")

swe_tot_date<-read.csv("HJASNODAS_date.csv")
swe_tot_date<-swe_tot_date[,-1]
swe_tot_date$date<-as.Date(swe_tot_date$date)

#remove error values (only present in WY 2004)
swe_tot_date[swe_tot_date==-1275]<-NA

#remove columns with NA (so just shape of HJA watershed)
swe_tot_date_crop<-swe_tot_date[,colSums(is.na(swe_tot_date))<nrow(swe_tot_date)]

#prep for loop
#list all WY
WY_list<-unique(swe_tot_date_crop$WY)

#assign leap years to variable
leapYear<-c(2004, 2008, 2012, 2016, 2020)

#open lists for loop to calculate f_melt
f_melt_year_list<-list()
f_melt_final<-list()

#open loop for f_melt calculations
#skip first 5 columns - these are date columns
for (i in 6:length(swe_tot_date_crop)) {
  
  #subset all 5 date columns and one grid cell of swe data
  data<-swe_tot_date_crop[,c(1:5,i)]
  
  #subset to individual water year
  for (j in 1:length(WY_list)) {
    
    #pull one water year
    data_year<-subset(data, data$WY==WY_list[j])
    
    #rename columns
    names(data_year)<-c("date", "year", "month", "doy", "WY", "swe")
    
    #get water year day (1-365 starting on Oct 1)
    data_year$water.day<-get_waterYearDay(data_year$date)
    
    #lag data to subtract between day and day+1
    data_year$lag<-Lag(data_year$swe, shift=-1)
    
    #calculate difference between day+1 and day
    data_year$diff<-data_year$lag-data_year$swe
    data_year$diff[is.na(data_year$diff)]<-0
    
    #add column denoting whether difference between day+1 and day is negative or positive
    data_year$posneg<-ifelse(is.negative(data_year$diff), "negative", "positive")
    
    #get max SWE
    max<-max(data_year$swe, na.rm=T)
    
    #subset to row of max swe
    peakSWE<-subset(data_year, data_year$swe==max)
    
    #subset dataframe to only data BEFORE peak SWE day
    prepeakSWE<-subset(data_year, data_year$water.day < peakSWE$water.day)
    
    #further subset to only include days where SWE decreased
    prepeakSWEneg<-subset(prepeakSWE, prepeakSWE$posneg=="negative")
    
    #get all decreasing swe days for whole year (before day 300 to avoid early fall snow)
    neg<-subset(data_year, data_year$posneg=="negative" & data_year$water.day < 300)
    
    #neg$absmelt<-abs(neg$diff)
    
    #neg$cummelt<-cumsum(neg$absmelt)
    
    #ggplot()+geom_line(sage, mapping = aes(doy, swe))+
      #geom_bar(sageneg, mapping=aes(x=doy, y=absmelt), stat="identity", fill="blue")+
      #geom_line(sageneg, mapping = aes(x=doy, y=cummelt), col="red", lty=2)+
      #geom_vline(xintercept = sagepeakSWE$doy, lty=2)+labs(x="Day of Year", y="SWE")+ggtitle("Sagehen 2017")
    
    #calculate f_melt stat
    f_melt<-sum(prepeakSWEneg$diff)/sum(neg$diff)
    
    #add to list
    f_melt_year_list[[j]]<-f_melt

  }
  
  #turn f_melt_year_list into dataframe
  f_melt_dataframe<-data.frame(do.call(rbind, f_melt_year_list))
  
  #switch rows and columns
  f_melt_df_t<-as.data.frame(t(f_melt_dataframe))
  
  #add f_melt dataframes to list
  f_melt_final[[i]]<-f_melt_df_t
  
}

#combine all f_melt dataframes into one dataframe
f_melt_df_total<-data.frame(do.call(rbind, f_melt_final[-c(1:5)]))

names(f_melt_df_total)<-WY_list

write.csv(f_melt_df_total, "HJA_f_melt.csv")

mean_list<-list()
sd_list<-list()

for (i in 1:length(f_melt_df_total)) {
  
  mean_list[[i]]<-mean(f_melt_df_total[,i])
  sd_list[[i]]<-sd(f_melt_df_total[,i])
  
}

f_melt_avg_df<-data.frame(do.call(rbind, mean_list))
f_melt_sd_df<-data.frame(do.call(rbind, sd_list))
f_melt_tot_df<-cbind(WY_list, f_melt_avg_df)
f_melt_tot_df<-cbind(f_melt_tot_df, f_melt_sd_df)
names(f_melt_tot_df)<-c("WY", "f_melt", "sd")

peak_swe_year_list<-list()
peak_swe_final<-list()

#open loop for peak swe calculations
#skip first 5 columns - these are date columns
for (i in 6:length(swe_tot_date_crop)) {
  
  #subset all 5 date columns and one grid cell of swe data
  data<-swe_tot_date_crop[,c(1:5,i)]
  
  #subset to individual water year
  for (j in 1:length(WY_list)) {
    
    #pull one water year
    data_year<-subset(data, data$WY==WY_list[j])
    
    #rename columns
    names(data_year)<-c("date", "year", "month", "doy", "WY", "swe")
    
    #get max SWE
    max<-max(data_year$swe, na.rm=T)
    
    #add to list
    peak_swe_year_list[[j]]<-max
    
  }
  
  #turn f_melt_year_list into dataframe
  peak_swe_dataframe<-data.frame(do.call(rbind, peak_swe_year_list))
  
  #switch rows and columns
  peak_swe_df_t<-as.data.frame(t(peak_swe_dataframe))
  
  #add f_melt dataframes to list
  peak_swe_final[[i]]<-peak_swe_df_t
  
}

#combine all f_melt dataframes into one dataframe
peak_swe_df_total<-data.frame(do.call(rbind, peak_swe_final[-c(1:5)]))

names(peak_swe_df_total)<-WY_list

write.csv(peak_swe_df_total, "HJA_peak_swe.csv")

swe_mean_list<-list()
swe_sd_list<-list()

for (i in 1:length(peak_swe_df_total)) {
  
  swe_mean_list[[i]]<-mean(peak_swe_df_total[,i])
  swe_sd_list[[i]]<-sd(peak_swe_df_total[,i])
  
}

peak_swe_avg_df<-data.frame(do.call(rbind, swe_mean_list))
peak_swe_sd_df<-data.frame(do.call(rbind, swe_sd_list))
peak_swe_tot_df<-cbind(WY_list, peak_swe_avg_df)
peak_swe_tot_df<-cbind(peak_swe_tot_df, peak_swe_sd_df)
names(peak_swe_tot_df)<-c("WY", "peakSWE", "sd")

pdf("HJASnodasPlots.pdf")

#f_melt plot
ggplot(f_melt_tot_df)+geom_line(aes(WY, f_melt))+
  geom_pointrange(aes(x=WY, y=f_melt, ymin=f_melt-sd, ymax=f_melt+sd))

#peak swe plot
ggplot()+geom_line(peak_swe_tot_df, mapping = aes(WY, peakSWE))+
  geom_pointrange(peak_swe_tot_df, mapping = aes(x=WY, y=peakSWE, ymin=peakSWE-sd, ymax=peakSWE+sd))

#both comparison
ggplot()+geom_line(f_melt_tot_df, mapping = aes(WY, f_melt), col="red")+
  geom_line(peak_swe_tot_df, mapping = aes(WY, peakSWE/500), col="blue")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*500))

#bind together f_melt and peak swe dataframes
swe_tot<-cbind(f_melt_tot_df, peak_swe_tot_df)
swe_tot<-swe_tot[,c(1,2,5)]

#correlation
ggplot(swe_tot)+geom_point(aes(peakSWE, f_melt))

dev.off()
