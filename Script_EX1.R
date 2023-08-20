############################
### FIELDimageR pipeline ###
############################

################
### Packages ### 
################

devtools::install_github("OpenDroneMap/FIELDimageR", dependencies=FALSE)
# devtools::install_github("filipematias23/FIELDimageR.Extra", dependencies=FALSE)

library(FIELDimageR)
library(FIELDimageR.Extra)
library(raster)
library(agricolae)
library(reshape2)
library(ggplot2)
library(lme4)
library(plyr)
library(DescTools)
library(ggrepel)
library(terra)
library(mapview)
library(leafsync)

############
### Data ###
############

MOSAIC<-list.files("./MOSAIC/")
DSM<-list.files("./DSM/")

###################
### Basic steps ###
###################

# Uploading an example mosaic
Test <- rast(paste("./MOSAIC/",MOSAIC[3],sep = ""))

##################
### Shape file ###
##################

# Reading FieldData.csv
Data<-read.csv(file = "EX1_Data.csv",header = T)

# Making the field Map
Map<-fieldMap(fieldPlot = Data$Plot,
              fieldColumn = Data$Column,
              fieldRow = Data$Row,
              decreasing = T)
Map

rotate <- function(x) t(apply(x, 2, rev))

Map<-rotate(rotate(Map))
Map

# Building the plot shapefile (ncols = 14 and nrows = 10)

plotShape<-fieldShape_render(mosaic = Test, 
                      ncols = 14, 
                      nrows = 10, 
                      fieldData = Data, 
                      PlotID = "Plot", 
                      fieldMap = Map,
                      buffer = -0.05)

# Vizualizing new plot grid shape object with field data:
fieldView(mosaic = Test,
          fieldShape = plotShape,
          type = 2,
          alpha = 0.2)

# Editing plot grid shapefile:
plotShape<- fieldShape_edit(mosaic=Test,
                            fieldShape=plotShape)

# Checking the edited plot grid shapefile:
fieldView(mosaic = Test,
          fieldShape = plotShape,
          type = 2,
          alpha = 0.2)

##########################
### Vegetation indices ###
##########################

Test.Indices<- fieldIndex(mosaic = Test, 
                         Red = 1, Green = 2, Blue = 3, 
                         index = c("NGRDI","BGI", "GLI","VARI"), 
                         myIndex = c("(Red-Blue)/Green","2*Green/Blue"))

#######################################
### Removing soil and making a mask ###
#######################################

Test.kmean<-fieldKmeans(mosaic=Test.Indices,
                        clusters = 3)
fieldView(Test.kmean)

# Check which cluster is related to plants or soil based on the color. 
rgb<-fieldView(Test)
plants<-fieldView(Test.kmean==2)
sync(rgb,plants)

# Soil Mask (cluster 2) to remove soil effect from the mosaic using FIELDimageR::fieldMask :
mask<-Test.kmean==2
Test.RemSoil<-fieldMask(Test.Indices,
                        mask = mask,
                        cropValue = 1,
                        cropAbove = F) 

fieldView(Test.RemSoil$newMosaic,
          fieldShape = plotShape,
          type = 2,
          alpha_grid = 0.2)

############################
### Extracting plot data ###
############################

Test.Info<- fieldInfo_extra(mosaic = Test.RemSoil$newMosaic,
                      fieldShape = plotShape)

fieldView(Test.RemSoil$newMosaic,
          fieldShape = Test.Info,
          plotCol = "Trait",
          type = 2,
          alpha_grid = 0.7)

###############################
### Estimating plant height ###
###############################

# Uploading files from soil base and vegetative growth:
DSM0 <- rast(paste("./DSM/",DSM[1],sep = ""))
DSM1 <- rast(paste("./DSM/",DSM[3],sep = ""))

# Canopy Height Model (CHM):
DSM0.R <- resample(DSM0, DSM1)
CHM <- DSM1-DSM0.R

fieldView(CHM,
          fieldShape = Test.Info,
          colorOptions = viridisLite::inferno,
          type = 2,
          alpha_grid = 0.05)

# Removing the soil using mask from step 4:
CHM.RemSoil <- fieldMask(CHM, mask = Test.RemSoil$mask)

# Extracting the estimate plant height average (EPH):
EPH<-CHM.RemSoil$newMosaic
names(EPH)<-"EPH"
Test.Info <- fieldInfo_extra(mosaic = EPH, 
                 fieldShape = Test.Info, 
                 fun = mean) 
colnames(Test.Info)
colnames(Test.Info)[17]<-"EPH"

########################################
### Evaluating all mosaics in a loop ###
########################################

DataTotal<-NULL
for(i in 2:length(MOSAIC)){
  EX1 <- rast(paste("./MOSAIC/",MOSAIC[i],sep = ""))
  EX1.RemSoil<-fieldMask(EX1, plot = F)
  EX1.Indices<- fieldIndex(mosaic = EX1.RemSoil$newMosaic,
                            Red = 1, Green = 2, Blue = 3,
                            index = c("NGRDI","BGI", "GLI","VARI"),
                            myIndex = c("(Red-Blue)/Green","2*Green/Blue"),
                            plot = F)
  EX1.Info<- fieldInfo_extra(mosaic = EX1.Indices[[c("NGRDI","BGI", "GLI","VARI","myIndex.1","myIndex.2")]],
                        fieldShape = plotShape)
  DSM0 <- rast(paste("./DSM/",DSM[1],sep = ""))
  DSM1 <- rast(paste("./DSM/",DSM[i],sep = ""))
  DSM0.R <- resample(DSM0, DSM1)
  CHM <- DSM1-DSM0.R
  CHM.RemSoil <- fieldMask(CHM, mask = EX1.RemSoil$mask,plot = F)
  EPH<-CHM.RemSoil$newMosaic
  names(EPH)<-"EPH"
  EX1.Info <- fieldInfo_extra(EPH,
                   fieldShape = EX1.Info,
                   fun = "mean")
  DataTotal<-rbind(DataTotal,
                   data.frame(DAP=as.character(do.call(c,strsplit(MOSAIC[i],split = "_"))[2]),
                              EX1.Info))
  print(paste("### Completed: ", "Mosaic_",i," ###",sep=""))
  }

colnames(DataTotal)
colnames(DataTotal)[13]<-"EPH"
DataTotal<-DataTotal[,!colnames(DataTotal)%in%c("ID","ID.1","PlotID","geometry")] # Removing column 12 ("ID.1")
#write.csv(DataTotal,"DataTotal.csv",row.names = F,col.names = T)

DataTotal<-read.csv("DataTotal.csv",header = T)
DataTotal

################
### Graphics ###
################

DataTotal$Name<-as.factor(as.character(DataTotal$Name))
DataTotal$Row<-as.factor(as.character(DataTotal$Row))
DataTotal$Column<-as.factor(as.character(DataTotal$Column))
DataTotal$DAP<-as.numeric(as.character(DataTotal$DAP))
DataTotal$NGRDI<-as.numeric(as.character(DataTotal$NGRDI))
DataTotal$BGI<-as.numeric(as.character(DataTotal$BGI))
DataTotal$EPH<-as.numeric(as.character(DataTotal$EPH))

ggplot(DataTotal, aes(x = NGRDI,fill=as.factor(DAP))) +
  geom_density(alpha=.5,position = 'identity') +
  facet_wrap(~DAP,ncol = 1)+
  scale_fill_grey(start=1, end=0)+
  labs(y="#genotypes",x="NGRDI", fill="DAP") +
  theme_bw() 

####################
### Heritability ###
####################

DAP<-unique(DataTotal$DAP)

H2<-NULL
for(h in 1:length(DAP)){
  mod<-lmer(NGRDI~Row+Column+(1|Name),DataTotal[as.character(DataTotal$DAP)==DAP[h],])
  H2.a<-c("NGRDI",DAP[h],as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov))
  
  mod<-lmer(EPH~Row+Column+(1|Name),DataTotal[as.character(DataTotal$DAP)==DAP[h],])
  H2.b<-c("EPH",DAP[h],as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov))
  
  H2<-rbind.data.frame(H2,rbind(H2.a,H2.b))
}
colnames(H2)<-c("Trait","DAP","H2")
H2$H2<-as.numeric(as.character(H2$H2))
H2$DAP<-as.numeric(as.character(H2$DAP))

ggplot(H2,aes(x=as.factor(DAP),y=H2,fill=as.factor(DAP)))+
  geom_bar(stat="identity")+
  facet_wrap(~Trait)+
  scale_fill_grey(start=0.8, end=0.2)+
  labs(x="Days After Planting (DAP)", fill="")+
  geom_text(aes(label=round(H2,2)), vjust=1.6, color="white", size=6)+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


##################################
### Area under the curve (AUC) ###
##################################

DataTotal1<-DataTotal[as.character(DataTotal$Name)%in%c("G43","G44","G45"),]

ggplot(data=DataTotal1, aes(x=as.numeric(DAP), y= NGRDI, col= Name, group=Name)) +
  geom_point(size=6)+
  geom_line(size=1.2) +
  scale_color_grey(start=0.8, end=0.2)+
  labs(x="Days After Planting (DAP)", fill="", col="")+
  theme_linedraw()

Trait<-c("NGRDI") # c("GLI","EPH")
Plot<-as.character(unique(DataTotal$Plot))

DataAUC<-NULL
for(a1 in 1:length(Plot)){
  D1<-DataTotal[as.character(DataTotal$Plot)==Plot[a1],]
  x1<-c(0,as.numeric(D1$DAP))
  y1<-c(0,as.numeric(D1[,Trait]))
  DataAUC <- rbind(DataAUC,
                   c(NGRDI_AUC=AUC(x = x1[!is.na(y1)], y = y1[!is.na(y1)]),
                     Name=unique(as.character(D1$Name)),
                     Trait=unique(D1$Trait),
                     Row=unique(D1$Row),
                     Column=unique(D1$Column)))}

DataAUC<-as.data.frame(DataAUC)
DataAUC$NGRDI_AUC<-as.numeric(as.character(DataAUC$NGRDI_AUC))
DataAUC$Name<-as.factor(DataAUC$Name)
DataAUC$Row<-as.factor(DataAUC$Row)
DataAUC$Column<-as.factor(DataAUC$Column)
DataAUC

### AUC Heritability ###

mod<-lmer(NGRDI_AUC~Row+Column+(1|Name),DataAUC)
H2<-as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov)
H2

ggplot(DataAUC, aes(x = NGRDI_AUC)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.5,position = 'identity', fill="cadetblue") +
  labs(y="#genotypes",x=paste("Area under the curve (AUC_NGRDI: H2=",round(H2,2),")",sep="")) +
  theme_bw() 

### Linear Regression ###

DataTotal.Reg<-subset(DataTotal,DAP=="40")
DataTotal.Reg$Check<-as.character(DataTotal.Reg$Name)
DataTotal.Reg$Check[!DataTotal.Reg$Check%in%c("G43","G44","G45")]<-""

ggplot(DataTotal.Reg,aes(y=NGRDI, x=EPH)) + 
  geom_point() +
  geom_smooth(method=lm)+
  labs(y="EPH",x="NGRDI",fill="",alpha="")+
  geom_vline(aes(xintercept=0.02),col="red", linetype = 2, size=0.7) +
  theme_bw()+
  geom_text_repel(aes(label = Check),
                  size = 3.5, col="red",
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

###########
### END ###
###########
