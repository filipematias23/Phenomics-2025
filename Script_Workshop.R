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
library(tidyverse)
library(gridExtra)
library(tidyr)
library(corrplot)
library(GGally)

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

# Option-01:
Test.RemSoil<-fieldMask(Test.Indices)

# Option-02:                        
Test.kmean<-fieldKmeans(mosaic=Test.Indices, clusters = 3)
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
#CHM <-fieldHeight(DSM0.R,DSM1) # The function fieldHeight() also calculates height and volume!

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
                 fun = "mean") 
colnames(Test.Info)
colnames(Test.Info)<-gsub("_mean","",colnames(Test.Info))

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
  EX1.Info<- fieldInfo_extra(mosaic = EX1.Indices,
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
colnames(DataTotal)<-gsub("_mean","",colnames(DataTotal))
DataTotal<-DataTotal[,!colnames(DataTotal)%in%c("ID","ID.1","PlotID","geometry")] # Removing column 12 ("ID.1")
#write.csv(DataTotal,"DataTotal.csv",row.names = F,col.names = T)

DataTotal<-read.csv("DataTotal.csv",header = T)
DataTotal

################
### Graphics ###
################

# Data Preparation:

DataTotal <- DataTotal %>%
  mutate(across(c(Name, Row, Column), as.factor),
         across(c(DAP, NGRDI, BGI, EPH), as.numeric))

str(DataTotal)

# Visualization:

variable<-"NGRDI"

ggplot(DataTotal, aes(x = .data[[variable]], fill = as.factor(DAP))) +
  geom_density(alpha = 0.5, position = 'identity') +
  facet_wrap(~DAP, ncol = 1) +
  scale_fill_viridis_d() +
  scale_fill_grey(start=0.8, end=0.2)+
  labs(y = "Density", x = variable, fill = "DAP") +
  theme_minimal() +
  theme(legend.position = "right")

####################
### Heritability ###
####################

# Advanced Heritability Analysis:
calculate_heritability <- function(data, traits) {
  dap_values <- unique(data$DAP)
  results <- data.frame(
    Trait = character(),
    DAP = numeric(),
    H2 = numeric(),
    Vg = numeric(),    
    Ve = numeric(),    
    stringsAsFactors = FALSE
  )
  
  for(dap in dap_values) {
    subset_data <- filter(data, DAP == dap)
    
    for(trait in traits) {
      tryCatch({
        model <- lmer(as.formula(paste(trait, "~ Row + Column + (1|Name)")), 
                      data = subset_data)
        
        vc <- VarCorr(model)
        vg <- as.numeric(vc$Name)      
        ve <- attr(vc, "sc")^2         
        
        h2 <- vg / (vg + ve)
        
        new_row <- data.frame(
          Trait = trait,
          DAP = dap,
          H2 = h2,
          Vg = vg,
          Ve = ve,
          stringsAsFactors = FALSE
        )
        
        results <- rbind(results, new_row)
      }, error = function(e) {
        warning(paste("Error processing trait:", trait, "DAP:", dap, "\n", e))
      })
    }
  }
  return(results)
}

# Execute the function
h2_results <- calculate_heritability(DataTotal, c("NGRDI", "EPH"))

# Line plot for H2 over time:
p1 <- ggplot(h2_results, aes(x = DAP, y = H2, color = Trait)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", H2)), 
            vjust = -0.5, 
            size = 3.5) +
  theme_bw() +
  labs(title = "Heritability over time",
       x = "Days after planting",
       y = "Heritability (H²)") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, max(h2_results$H2) * 1.1)) # Adjust y-axis to accommodate labels

p1 #View H2 Graphic 

# Bar plot to compare Vg and Ve:
results_long <- tidyr::pivot_longer(h2_results, 
                                    cols = c(Vg, Ve),
                                    names_to = "Variance_Component",
                                    values_to = "Value")

p2 <- ggplot(results_long, 
             aes(x = DAP, y = Value, fill = Variance_Component)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Trait) +
  scale_fill_grey(start=0.8, end=0.2)+
  theme_bw() +
  labs(title = "Variance components",
       x = "Days after planting",
       y = "Variance") +
  theme(legend.position = "bottom")

p2 #View var Graphic

# Arrange plots
grid.arrange(p1, p2, ncol = 1)

##################################
### Area under the curve (AUC) ###
##################################

#AUC Importance: 
DataTotal %>%
  filter(Name %in% c("G43", "G44", "G45")) %>%
  ggplot(aes(x = as.numeric(DAP), 
             y = NGRDI, 
             color = Name, 
             group = Name)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 6, alpha = 0.9) +
  scale_color_grey(start = 0.8, end = 0.2) +
  labs(title = "NGRDI Evolution by Genotype",
       x = "Days After Planting (DAP)", 
       y = "NGRDI",
       color = "Genotype") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = unique(as.numeric(DataTotal$DAP)))

# Calculating AUC: 

traits <- c("NGRDI", "GLI", "EPH","myIndex.1")

calculate_auc <- function(data, traits) {
  results_list <- list()
  for(trait in traits) {
    auc_data <- data.frame()
    for(plot in unique(data$Plot)) {
      d1 <- data[data$Plot == plot, ]
      x1 <- c(0, as.numeric(d1$DAP))
      y1 <- c(0, as.numeric(d1[[trait]]))
      auc_value <- AUC(x = x1[!is.na(y1)], y = y1[!is.na(y1)])
      new_row <- data.frame(
        AUC = auc_value,
        Name = unique(as.character(d1$Name)),
        Trait = trait,
        Row = unique(d1$Row),
        Column = unique(d1$Column)
      )
      auc_data <- rbind(auc_data, new_row)
    }
    auc_data$AUC <- as.numeric(auc_data$AUC)
    auc_data$Name <- as.factor(auc_data$Name)
    auc_data$Row <- as.factor(auc_data$Row)
    auc_data$Column <- as.factor(auc_data$Column)
    results_list[[trait]] <- auc_data
  }
  return(results_list)
}

auc_results <- calculate_auc(DataTotal, traits)
auc_results[["NGRDI"]] #View results

### AUC Heritability ###
calculate_heritability <- function(data, trait_name) {
  mod <- lmer(AUC ~ Row + Column + (1|Name), data = data)
  vc <- VarCorr(mod)
  vg <- as.numeric(vc$Name)      
  ve <- attr(vc, "sc")^2         
  h2 <- vg / (vg + ve)           
  return(h2)
}

# H2 Graphic:
create_density_plot <- function(data, trait_name, h2) {
  ggplot(data, aes(x = AUC)) +
    geom_histogram(aes(y = ..density..), 
                   colour = "black", 
                   fill = "white",
                   bins = 30) +
    geom_density(alpha = 0.5,
                 position = 'identity', 
                 fill = "cadetblue") +
    labs(title = paste("Distribution of", trait_name),
         y = "Density",
         x = paste("Area under the curve (AUC_", trait_name, ": H² = ", 
                   round(h2, 2), ")", sep = "")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_list <- list()
results_summary <- data.frame(
  Trait = character(),
  Heritability = numeric(),
  Mean_AUC = numeric(),
  SD_AUC = numeric(),
  stringsAsFactors = FALSE
)

for(trait in traits) {
  h2 <- calculate_heritability(auc_results[[trait]], trait)
  plot_list[[trait]] <- create_density_plot(auc_results[[trait]], trait, h2)
  results_summary <- rbind(results_summary, 
                           data.frame(
                             Trait = trait,
                             Heritability = h2,
                             Mean_AUC = mean(auc_results[[trait]]$AUC),
                             SD_AUC = sd(auc_results[[trait]]$AUC)
                           ))
}

# AUC results:
grid.arrange(grobs = plot_list, ncol = 2)
print(results_summary)

###################
### Correlation ###
###################

# DAP:
dap <- "100"

# List of traits:
selected_traits <- c("Trait","NGRDI", "BGI", "GLI", "VARI", "myIndex.1", "myIndex.2", "EPH")

# Filtering data:
DataTotal.Cor <- DataTotal %>%
  filter(DAP == dap) %>%
  select(Name, all_of(selected_traits))

# Correlation matrix using agricolae
cor_matrix <- correlation(DataTotal.Cor[,selected_traits], 
                          method = "pearson")

# Color palette for correlation plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Correlation plot
corrplot(cor_matrix$correlation, 
         p.mat = cor_matrix$pvalue,
         sig.level = 0.05, # Significance level 5%
         method = "color", 
         col = col(200),  
         type = "upper", 
         order = "alphabet",
         addCoef.col = "black", 
         tl.col = "black", 
         tl.srt = 45, 
         insig = "blank", 
         diag = FALSE)

# Scatter plots with correlation
ggpairs(DataTotal.Cor[,selected_traits],
        upper = list(continuous = wrap("cor", method = "pearson")),
        lower = list(continuous = wrap("smooth", alpha = 0.3, color = "blue")),
        diag = list(continuous = wrap("densityDiag", fill = "lightblue"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#########################
### Linear Regression ###
#########################

dap="40"

DataTotal.Reg<-subset(DataTotal,DAP==dap)
DataTotal.Reg$Check<-as.character(DataTotal.Reg$Name)
DataTotal.Reg$Check[!DataTotal.Reg$Check%in%c("G43","G44","G45")]<-""

reg_model <- lm(NGRDI ~ EPH, data = DataTotal.Reg)
r2 <- round(summary(reg_model)$r.squared, 3)
eq <- sprintf("y = %.3fx + %.3f", 
              coef(reg_model)[2], 
              coef(reg_model)[1])

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
  labs(title = paste("DAP:", dap),
       subtitle = paste("R² =", r2, "\n", eq),
       x = "EPH",
       y = "NGRDI") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

###########
### END ###
###########
