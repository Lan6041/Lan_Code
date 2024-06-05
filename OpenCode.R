library(sp)
library(raster)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(gdm)
library(dplyr)
library(foreign)

##species data should include name, site, longitude and latitude
sppData=read.csv("sppdata.csv")

##fit raster environmental data
##sets up site-pair table
##n equals the number of environmental variables
envRast <- stack("x1.tif","x2.tif","xn.tif")

##environmental raster data
sitePairRast <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
                               siteColumn="site", predData=envRast)
##sometimes raster data returns NA in the site-pair table, these rows will have to be removed 
##before fitting gdm
sitePairRast <- na.omit(sitePairRast)

##fit raster GDM
gdmRastMod <- gdm(sitePairRast, geo=TRUE)

##input habitat condition
hc <- read.csv("hc.csv")

##prepare loop data
rtpenvrast <- as.data.frame(rasterToPoints(envRast))
rtpenvrast <- na.omit(rtpenvrast)

defaultdistance <- c(1)
defaultweight <- c(1)

sitepi <- rtpenvrast[0,1:5]
colnames(sitepi) <- c("s1.xCoord","s1.yCoord","sumsij","sumweightsij","pi")

##start loop
for (i in 1:nrow(rtpenvrast)){
  sitei <- rtpenvrast[i,]
  rtpenvrastsamp <- (rtpenvrast[-i,])
  weightstemp <- (weights[-i])
  sitepair <- cbind(rtpenvrastsamp,sitei,defaultdistance,defaultweight)[c((2*n+5),(2*n+6),(n+3),(n+4),1,2,(n+5):(2*n+4),3:(n+2))]
  colnames(sitepair) <- colnames(sitePairRast)
  dij <- predict(gdmRastMod,sitepair)
  sij <- (1-dij)
  sitepairisij <- cbind(sitepair[,c(3:6)],sij,weightstemp)
  
  sitegroup <- sitepairisij%>%
    group_by(s1.xCoord,s1.yCoord)%>%
    summarise(sumsij=sum(sij),
              sumweightsij=sum(weightstemp*sij))
  sitegroup$pi <- (sitegroup$sumweightsij/sitegroup$sumsij)^0.25
  sitepi <- rbind(sitepi,sitegroup[,c(1:5)])
}

##output results
write.csv(sitepi,file = "result.csv")
