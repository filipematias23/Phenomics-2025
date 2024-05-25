############################
### FIELDimageR pipeline ###
############################

################
### Packages ### 
################

#devtools::install_github("OpenDroneMap/FIELDimageR", dependencies=FALSE)
#devtools::install_github("filipematias23/FIELDimageR.Extra", dependencies=FALSE)

library(FIELDimageR)
library(FIELDimageR.Extra)
library(ggplot2)
library(agricolae)
library(reshape2)
library(lme4)
library(readxl)
library(terra)
library(mapview)
library(sf)
library(stars)

# Uploading one image as example and decreasing the resolution
EX3<-rast("soybean/11.jpg")
EX3<-imgLAB(EX3)
# EX3<-aggregate(EX3, fact= 4)
plotRGB(EX3)

# Select one index to identify leaves and remove the background
EX3.I1<- fieldIndex(mosaic = EX3,index = c("SI","BGI","BI"))

# Thresholding
dev.off()
par(mfrow=c(1,2))
hist(EX3.I1$BGI)
plot(EX3.I1$BGI)

# Removing the background
EX3.R<- fieldMask(mosaic = EX3, index = "BGI",
                   cropValue = 0.7,
                   cropAbove = T)

# Counting the total number of seeds
EX.P.Total<-fieldCount(mosaic = EX3.R$mask,plot = T)

dev.off()
par(mfrow=c(1,2))
hist(EX.P.Total$area)
plotRGB(EX3.R$newMosaic)

Total=EX.P.Total[EX.P.Total$area>15000,]
dim(Total) #50 seeds

# Select one index to identify green seeds
EX3.I2<- fieldIndex(mosaic = EX3.R$newMosaic,index = c("SI","BGI","BI"))

#BI index
plot(EX3.I2$BI)

# Selecting green seeds
EX3.R2<- fieldMask(mosaic = EX3.R$newMosaic, 
                   index = "BI",
                   #myIndex = "Blue",
                   cropValue = 130,
                   cropAbove = T)

# Counting the number of green seeds
EX.Green<-fieldCount(mosaic = EX3.R2$mask,plot = T)

dev.off()
par(mfrow=c(1,2))
hist(EX.Green$area)
plotRGB(EX3.R2$newMosaic)

Green=EX.Green[EX.Green$area>11000,]
dim(Green) #50 seeds

# Joying information
data.frame(Total=dim(Total)[1],
           Green=dim(Green)[1],
           Percentage=round(dim(Green)[1]/dim(Total)[1],3))

################
### Parallel ###
################

# Required packages
library(parallel)
library(foreach)
library(doParallel)

# Images names (folder directory: "./soybean/")
pics<-list.files("./soybean/")

# Number of cores
n.core<-2

# Starting parallel
cl <- makeCluster(n.core, output = "")
registerDoParallel(cl)
EX.Table.Parallel <- foreach(i = 1:length(pics), .packages = c("stars","sf","terra","FIELDimageR"), 
                             .combine = rbind) %dopar% {
                               EX3<-rast(paste("./soybean/",pics[i],sep = ""))
                               EX3<-imgLAB(EX3)
                               # Removing the background
                               EX3.R<- fieldMask(mosaic = EX3, index = "BGI",
                                                 cropValue = 0.7,
                                                 cropAbove = T)
                               # Counting the total number of seeds
                               EX.P.Total<-fieldCount(mosaic = EX3.R$mask,plot = T)
                               Total=EX.P.Total[EX.P.Total$area>15000,]
                               # Selecting green seeds
                               EX3.R2<- fieldMask(mosaic = EX3.R$newMosaic, 
                                                  index = "BI",
                                                  cropValue = 130,
                                                  cropAbove = T)
                               # Counting the number of green seeds
                               EX.Green<-fieldCount(mosaic = EX3.R2$mask,plot = T)
                               Green=EX.Green[EX.Green$area>10800,]
                               data.frame(Total=dim(Total)[1],
                                          Green=dim(Green)[1],
                                          Percentage=round(dim(Green)[1]/dim(Total)[1],2))
                             }
stopCluster(cl)
rownames(EX.Table.Parallel)<-pics
EX.Table.Parallel 

##########################
### Seeds measurements ###
##########################

# Taking individual seed measurements (Remove artifacts by changing the parameter *minArea* and observing the values on EX3.D$Dimension$area)
dev.off()

# Taking measurements:
EX.Obj.D<-fieldCount(mosaic = EX3.R$mask,
                     plot = T)

EX.Obj.D=EX.Obj.D[EX.Obj.D$area>15000,]

fieldView(EX3,EX.Obj.D,type = 2)
plotRGB(EX3)
plot(EX.Obj.D$geometry, add=T, border="red")
plot(EX.Obj.D[1,], add=T, col="yellow")
EX.Obj.I<- fieldIndex(mosaic = EX3,index = c("SI","BGI","BI"))
EX.Obj.Data<-fieldInfo_extra(mosaic = EX.Obj.I[[c("SI","BGI","BI")]], 
                             fieldShape = EX.Obj.D)
plot(EX.Obj.Data)

# Data visualization: 
library(reshape2)
Data.Obj1<-melt(EX.Obj.Data,measure.vars = c("area","perimeter","width","SI_mean","BGI_mean","BI_mean"))

library(ggplot2)
ggplot(Data.Obj1, aes(x=value, fill=variable)) +
  geom_histogram(aes(y=..density..), colour="black")+
  geom_density(alpha=.2)+
  facet_wrap(~variable, scales = "free")

###########
### END ###
###########
