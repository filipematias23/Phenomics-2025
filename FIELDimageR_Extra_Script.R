##################################
### FIELDimageR.Extra pipeline ###
##################################

################
### Packages ### 
################

library(FIELDimageR)
library(FIELDimageR.Extra)
library(terra)
library(mapview)
library(sf)
library(stars)
library(leafsync)
library("RColorBrewer")
library(corrplot)
library(agricolae)

############
### Data ###
############

# Uploading an example mosaic:
Test <- rast("EX1_RGB.tif")

# Visualization Option-01:
fieldView(Test)

# Visualization Option-02:
plotRGB(Test)

##################
### Shape file ###
##################

# Reading DataTable.csv:
DataTable<-read.csv("DataTable.csv",header = T)  

# Field Map:
fieldMap<-fieldMap(fieldPlot=DataTable$Plot, fieldColumn=DataTable$Row, fieldRow=DataTable$Range, decreasing=T)
fieldMap

# Making plot grid:
plotShape<-fieldShape_render(mosaic = Test,
                             ncols=14,
                             nrows=9,
                             fieldData = DataTable,
                             fieldMap = fieldMap,
                             PlotID = "Plot",
                             buffer = -0.05)

# Checking the new plot grid shapefile:
fieldView(mosaic = Test,
          fieldShape = plotShape,
          type = 2,
          alpha = 0.2)

# Editing plot grid shapefile:
editShape<- fieldShape_edit(mosaic=Test,
                            fieldShape=plotShape)

# Checking the edited plot grid shapefile:
fieldView(mosaic = Test,
          fieldShape = editShape,
          type = 2,
          alpha = 0.2)

##########################
### Vegetation indices ###
##########################

Test.Indices<- fieldIndex(mosaic = Test, 
                          Red = 1, Green = 2, Blue = 3, 
                          index = c("NGRDI","BGI", "GLI","VARI"), 
                          myIndex = c("(Red-Blue)/Green","2*Green/Blue"))

#To visualize single band layer 
single_layer<-Test.Indices$GLI

# Visualization-01:
fieldView(single_layer,
          fieldShape = editShape,
          colorOptions = viridisLite::inferno,
          type = 2,
          alpha_grid = 0.2)

# Visualization-02:
fieldView(single_layer,
          fieldShape = editShape,
          plotCol = c("Yield"),
          colorOptions = 'RdYlGn',
          type = 1,
          alpha_grid = 0.2)

#####################
### Removing soil ###
#####################

### Option_01 = Using fieldKmeans() to extract soil
Test.kmean<-fieldKmeans(mosaic=Test.Indices$GLI,
                        clusters = 2)

fieldView(Test.kmean)
#cluster 1 represents plants
#cluster 2 represents soil

# Comparing RGB and Binary:
m0<-fieldView(Test)
m1<-fieldView(Test.kmean==1)
sync(m0,m1)

# Soil Mask (cluster 2):
soil<-Test.kmean==2
Test.RemSoil<-fieldMask(Test.Indices,
                        mask = soil) 
fieldView(Test.RemSoil$newMosaic,
          fieldShape = editShape,
          type = 2,
          alpha_grid = 0.2)

### Option_02 = Using fieldSegment() - Supervised image segmentation

# Digitize soil object by drawing polygons at least 5-6 large polygon uniformly distributed
soil<-fieldView(mosaic = Test, editor = TRUE) #generate random 200 points for soil class
soil<-st_as_sf(st_sample(soil, 200))
soil$class<-'soil'

# Digitize plants object by drawing polygons. The number of polygon will depends upon the number
# of training points to be generated.
plants<-fieldView(mosaic = Test, editor = TRUE)
plants<-st_as_sf(st_sample(plants, 200)) #generate random 200 points for plants class
plants$class<-'plants'

#similarly you can digitize shadow, other objects by using draw polygon tool of editor
training_sam<-rbind(soil,plants)

# Random FOrest classification: 
classification<-fieldSegment(mosaic = Test, 
                             trainDataset = training_sam,
                             model = "rf")

# To display results of classification from randomForest
classification$sup_model
classification$pred
classification$rastPred
fieldView(classification$rastPred)
plot(classification$rastPred)

# Soil Mask (removing soil):
soil<-classification$rastPred=="soil"
Test.RemSoil<-fieldMask(Test.Indices,
                        mask = soil) 
fieldView(Test.RemSoil$newMosaic,
          fieldShape = editShape,
          type = 2,
          alpha_grid = 0.2)

### Option_03 = Using the traditional FIELDimageR::fieldMask()
Test.RemSoil<-fieldMask(Test.Indices) 
fieldView(Test.RemSoil$newMosaic,
          fieldShape = editShape,
          type = 2,
          alpha_grid = 0.2)

############################
### Extracting plot data ###
############################

DataTotal<- fieldInfo_extra(mosaic = Test.RemSoil$newMosaic,
                            fieldShape = editShape, 
                            fun = mean)

View(DataTotal)

# Visualization-01:
fieldView(mosaic = Test.RemSoil$newMosaic,
          fieldShape = DataTotal,
          plotCol = "Yield",
          col_grid = c("blue", "grey", "red"),
          type = 2,
          alpha_grid = 0.6)

# Visualization-02:
fieldView(mosaic = Test.RemSoil$newMosaic,
          fieldShape = DataTotal,
          plotCol = c("Yield"),
          type = 1)

# Visualization-03:
fieldView(mosaic = Test.RemSoil$newMosaic,
          fieldShape = DataTotal,
          plotCol = c("Yield","GLI"),
          type = 1)

# Visualization-04:
m1<-viewRGB(Test.RemSoil$newMosaic)
m2<-mapview(Test.RemSoil$newMosaic$NGRDI,layer.name="NGRDI")
m3<-mapview(Test.RemSoil$newMosaic$BGI,layer.name="BGI")
m4<-mapview(DataTotal,zcol="Yield")
sync(m1,m2,m3,m4)

#############################
### Saving Extracted Data ###
#############################

# Saving Data.csv (Removing 'geometry' info)
write.csv(as.data.frame(DataTotal)[,-dim(DataTotal)[2]],
          "DataTotal.csv",
          col.names = T, row.names = F)

# Saving grid fieldShape file
st_write(DataTotal,"grid.shp")
Saved_Grid = st_read("grid.shp")
plot(Saved_Grid)

############################################
### Cropping individual plots and saving ###
############################################

# Saving plot images format .tif according to 'Maturity': 
Field_plot_grids<- fieldCrop_grid(mosaic = rast(Test.Indices),
                                  fieldShape = editShape, 
                                  classifier = "Maturity", 
                                  plotID = "Plot",
                                  format = '.jpg',
                                  output_dir = "./")

fieldView(Field_plot_grids$'21')

fieldView(mosaic = Field_plot_grids$'21',
          fieldShape = editShape,
          plotCol ="Maturity")

# Reading .tif file and visualizing single plots with fieldView()
Field_plot_grids<- fieldCrop_grid(mosaic = rast(Test.Indices),
                                  fieldShape = editShape, 
                                  classifier = "Maturity", 
                                  plotID = "Plot",
                                  format = '.tif',
                                  output_dir = "./")
tif<-rast("./2/21.tif")
fieldView(mosaic = tif,
          fieldShape = editShape,
          plotCol ="Maturity")

###########################################
### Simple Data Extracted Visualization ###
###########################################

DataTotal<-as.data.frame(DataTotal)[,-dim(DataTotal)[2]]

### Correlation ###

colnames(DataTotal)
DataTotal.1<-scale(DataTotal[,7:20],scale = T)
r<-correlation(DataTotal.1)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.05, # Level of significance 5%
         method="color", col=col(200),  
         type="upper", 
         order="original",
         addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)

### Linear Regression - 01 ###

ggplot(DataTotal, aes(x = NGRDI,y = Yield)) +
  geom_point() +
  geom_smooth(method=lm)+
  scale_fill_grey(start=1, end=0)+
  theme_bw() 

### Linear Regression - 02 ###

Check<-rep("black",length(DataTotal$Yield))
Check[DataTotal$Yield>600]<-"blue"
Check[DataTotal$Yield<350]<-"red"
DataTotal$Check<-Check

ggplot(DataTotal, aes(x = NGRDI, y = Yield)) +
  geom_point(color=Check) +
  geom_smooth(method=lm)+
  scale_fill_grey(start=1, end=0)+
  theme_bw()+ 
  geom_vline(aes(xintercept=0.15),col="red", linetype = 2, size=0.7)

###########
### END ###
###########
