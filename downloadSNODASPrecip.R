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
require(EflowStats)
#install.packages("remotes")
#remotes::install_github("USGS-R/EflowStats")

###run this for GetSnodasPrecipDate function
##modified from rwrfhydro GetSnodasDepthSweDate function

GetSnodasPrecipDate <- function(datePOSIXct, outputDir='.', overwrite=FALSE, 
                                quiet=TRUE, 
                                parallel=(foreach::getDoParWorkers()>1) ){
  
  ## This atomic function gets called below. 
  GetSnodasPrecipDate.atomic <- function(datePOSIXct, outputDir=outputDir, overwrite=FALSE, 
                                         quiet=TRUE) {
    # date parameters
    yy <- format(datePOSIXct, c("%Y")); mm <- format(datePOSIXct, c("%m"))
    mon <- format(datePOSIXct, c("%h")); dd <- format(datePOSIXct, c("%d"))
    
    # depthProdId <- '1036',  sweProdId <- '1034'
    # Can construct the filenames. Calling this "0" incase they dont match the
    # names in the tarball.
    # depthFile0<- paste0(outputDir,'/','us_ssmv1',
    #                     '1036tS__T0001TTNATS',yy,mm,dd,'05HP001.dat.gz')
    precipFile0  <- paste0(outputDir,'/','us_ssmv0',
                           '1025SlL00T0024TTNATS',yy,mm,dd,'05DP001.dat.gz')
    
    ## If either of the depth or SWE files exist, bail out unless told to overwrite.
    if( (file.exists(precipFile0)) & !overwrite) return(0)
    
    # Go to the correct directory for arciving the data
    origDir <- getwd()
    setwd(outputDir)
    
    # theFile is the tarball
    theFile <- paste0('SNODAS_',yy,mm,dd,'.tar')
    
    # if theFile (tarball) exists, then skip downloading unless told to overwrite
    if(!file.exists(theFile) | overwrite) {
      theUrl <- paste0('ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/',
                       yy,'/',mm,'_',mon,'/',theFile)
      if(!quiet) print(paste('SNODAS: ', datePOSIXct))
      try(curl::curl_download(theUrl, theFile, quiet=quiet))
      if(!file.exists(theFile)) {
        warning(paste0('Error: File not obtained via FTP: ',theFile))
        setwd(origDir)
        return(FALSE)
      }
    }
    
    # unpack tarball to depth and SWE. name the tmpDir by date so dates can be done in parallel.
    tmpDir <- paste0('tmp',yy,mm,dd)
    if(file.exists(tmpDir)) unlink(tmpDir, recursive=TRUE)
    untar(theFile, exdir=tmpDir)
    precipFile <- list.files(path = tmpDir, pattern=glob2rx('*ssmv01025SlL00*.dat.gz'))  # SWE
    file.copy(paste0(tmpDir,'/',precipFile), precipFile0)
    # depthFile <- list.files(path = tmpDir, pattern=glob2rx('*ssmv11036tS*.dat.gz'))  # depth
    # file.copy(paste0(tmpDir,'/',depthFile), depthFile0)
    unlink(c(tmpDir, theFile), recursive=TRUE) 
    Sys.chmod(c(precipFile0), mode='0777', use_umask=FALSE)
    
    setwd(origDir)
    TRUE
  }
  
  ## FormalsToDf handles vector arguments and passes collated combos 
  ## to the atomic function 
  vecDf <- FormalsToDf(GetSnodasPrecipDate)
  vecDf <- vecDf[,-which(names(vecDf) %in% c('parallel'))] ## exclude formals not passed to atomic
  ret <- plyr::mlply(vecDf, GetSnodasPrecipDate.atomic, .parallel=parallel)
  names(ret) <- datePOSIXct
  if(length(ret)==1) ret <- ret[[1]]
  ret
}


#create list of dates
datesWanted<-seq.Date(as.Date("2003-10-01"), as.Date("2020-09-30"), by="day")

#set output directory
outDir<-c("/Users/keirajohnson/Documents/SNODAS_Precip_Updated")

#scaling factors
nCol <- 6935
nRow <- 3351 #columns and rows number: masked version of contiguous US
dataScaleFactor  <-  10  #multiply the data this amount, both depth and SWE

#set coordinate system
wgs84 <- CRS("+proj=longlat +datum=WGS84")

#read in HJA shapefile
HJA_shape<-readOGR("/Users/keirajohnson/Box Sync/Keira_Johnson/HJAshapefile/layers/globalwatershed.shp")

#open list
precip_HJA_matrix_mm<-list()
dateList<-list()

#constants
beg<-"us_ssmv01025SlL00T0024TTNATS"
precip_files<-list.files(path=outDir, pattern = beg)

for (i in 4263:length(datesWanted)) {
  
  #this allows the loop to keep running even if there is an error
  tryCatch({
    
    #download daily snodas for one day
    #GetSnodasPrecipDate(datesWanted[i], outputDir = outDir, overwrite = T, quiet=FALSE)
    
    #print date
    print(datesWanted[i])
    
    #get correct URL based on swe code and date
    date<-gsub("-", "", datesWanted[i])
    file_name<-str_subset(precip_files, date)

    #extract SWE
    precipFile0  <- paste0(outDir, "/", file_name)

    #extract from dat.gz file
    precipCon  <- gzcon(file(precipFile0, "rb"))

    #read binary data from file above
    precipData <- readBin(precipCon, integer(), n=nRow*nCol, size=2, signed=TRUE, endian='big')/dataScaleFactor

    #remove error values
    precipData[which(precipData==-9.999)]<-NA

    #convert to matrix
    precipDataMatrix<-matrix(precipData, ncol=nCol, nrow=nRow, byrow = T)

    #georeference
    precipRaster <- raster(precipDataMatrix, ymn=24.9504, ymx=52.8754, xmn=-124.7337, xmx=-66.9421, crs=wgs84)

    #extract data for HJA
    precip_HJA_mm<-crop(precipRaster, HJA_shape)
    precip_HJA_mm_mask<-mask(precip_HJA_mm, HJA_shape)

    #convert back to matrix
    precip_HJA_matrix_mm[i]<-as.matrix(precip_HJA_mm_mask)

    closeAllConnections()

    dateList[[i]]<-datesWanted[i]
    # 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

save(precip_HJA_matrix_mm, file = "precip_HJA_matrix_mm.RData")
