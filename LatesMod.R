#=======================================================================
#Topic:        Nile perch distribution in varrying bioclimatic conditions
#subtopic:     Species distribution modelling (sdm)
#Author:       Anthony Basooma; 
#Section:      Fish Biology and Ecology (Ecological Modelling)
#Organisation: National Fisheries Resources Research Institute
#========================================================================
#Setting working directory

setwd("E:/sdmmodel/species modelling/manuscript/sdm/data analysis")

library(sp)
library(dplyr)
library(tidyr)
library(raster)
library(mapview)
library(sdm)
library(usdm)
library(ggplot2)
library(extrafont)
library(dismo)

#Nile perch data for Uganda boundary
#-------------------------------------------------------------------------------------
Lates<-   read.csv(file = "Perch.csv", header = TRUE, strip.white = T)%>%
  filter(!is.na(lat), continent%in%c("AFRICA"))%>% dplyr::select(lon, lat)%>%  
  mutate(lates= 1)
coordinates(Lates)<- ~lon+lat

#Environmental variables data from Worldclim
#-------------------------------------------------------------------------------------
biocp <- raster::getData('worldclim', var='bio', res=2.5)
names(biocp)<- c("AMT", "MDR", "ISO", "TS", "MaxTWM", "MinTCM", "TAR", 
                 "MTWeQ","MTDQ", 
                 "MTWaQ","MTCoQ", "AP", "PWM", "PDrM","PS", "PWeQ", 
                 "PDrQ", "PWaQ",
                 "PCoQ")

#cropping according African layer 
#------------------------------------------------------------------------------------
Africashp    <- shapefile("Africa_layer/Africa.shp", warnPRJ=FALSE)
africa       <- shapefile("waterbodies_africa/waterbodies_africa.shp",
                          warnPRJ=FALSE)
#uganda       <- shapefile("ugandalakes/Uganda_Lakes.shp", warnPRJ=FALSE)
extcor       <- floor(extent(africa))
#extcug       <- floor(extent(uganda))
biod         <- crop(biocp, extcor)
bio.mask     <- mask(biod, africa)
plot(bio.mask[[1]])
plot(Africashp, add=TRUE, border="grey")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Check for multicollinearity 
#--------------------------------------------------------------------------------------
bio1<- extract(bio.mask, Lates)
bio1<- data.frame(bio1)
bio.exc<- vifstep(bio1, th=10)
biocdp<- exclude(bio.mask, bio.exc)

#-------------------------------------------------------------------------------------
#Data preparation and model fitting 
sdmpre<- sdmData(lates~., train=Lates, predictors=biocdp, bg=list(n=1000)) 

set.seed(1135)
SDModel<- sdm(lates~., data=sdmpre, methods = c("glm", "svm", "gam", "rf", 
                                                "brt","bioclim"),
              replication=c("sub"), test.percent=30, n=10)
#-------------------------------------------------------------------------------------
#Model predcitions 
PredN<- predict(SDModel, biocdp, "LatesN.img")

plot(PredN[[1]])
plot(Africashp, add=TRUE, border="grey")
#--------------------------------------------------------------------------------------
#Ensemble model predictions
ensp<- ensemble(SDModel, biocdp, filename = "Lnilo20.img",
                setting=list(method="weighted", stat="AUC"))
plot(ensp, col=col(4))
plot(Africashp, add=TRUE, border="grey")
zoom(ensp, ext=extcug)
#=======================================================================================
#Individual Model ensembles
glm1<- getModelId(sdmodel, method = c("glm"))
glme<- ensemble(sdmodel, biocdp, filename = "GLM.img",
                setting=list(method="weighted", stat="AUC", id=glm1))
rf<- getModelId(sdmodel, method = c("rf"))
rfe<- ensemble(sdmodel, biocdp, filename = "rf.img",
               setting=list(method="weighted", stat="AUC", id=rf))
gam1<- getModelId(sdmodel, method = c("gam"))
game<- ensemble(sdmodel, biocdp, filename = "gam.img",
                setting=list(method="weighted", stat="AUC", id=gam1))
svm<- getModelId(sdmodel, method = c("svm"))
svme<- ensemble(sdmodel, biocdp, filename = "svm.img",
                setting=list(method="weighted", stat="AUC", id=svm))
brt1<- getModelId(sdmodel, method = c("brt"))
brte<- ensemble(sdmodel, biocdp, filename = "brt.img",
                setting=list(method="weighted", stat="AUC", id=brt1))

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#Future predictions 
biof50<- raster::getData("CMIP5", var="bio", res=2.5, year=50, 
                         model = "AC",rcp=85)
biof70<- raster::getData("CMIP5", var="bio", res=2.5, year=70, 
                         model = "AC",rcp=85)

biof50<- crop(biof50, extent(extcor))
biof70<- crop(biof70, extent(extcor))
names(biof50)<- names(biocpp)
names(biof70)<- names(biocpp)


ensf50<- ensemble(sdmodel, biof50, 
                  filename = "L.niloticus distribution2 at 50.img", 
                  setting=list(method="weighted", stat="AUC"))
ensf70<- ensemble(sdmodel, biof70, 
                  filename = "L.niloticus distributions at 70s.img", 
                  setting=list(method="weighted", stat="AUC"))

ensfrgrg<- ensemble(sdmodel, biof50, filename = "Distribution44.img",
                    setting=list(method="weighted", stat="AUC"))

ensf702<- ensemble(sdmodel, biof70, 
                   filename = "L.niloticus distribut at 70s.img", 
                   setting=list(method="weighted", stat="AUC"))

#==========================================================================
#Graphical output presentation
#-----------------------------

ensind<- stack(glme, rfe, game, svme, brte, ensp )

par(mfrow=c(2,3), mar=c(1, 1, 0.6, 0.6), family="Lucida Sans")

for (ens in c(glme, rfe, game, svme, brte, ensp)) {
  col<- colorRampPalette(c("grey", "coral4", "red", "purple"))
  plot(ens, cex=0.5, col=col(4))
}
mapview(glme)


par(mfrow=c(2,3), mar=c(1, 1, 0.6, 0.6), family="Lucida Sans")
for (ens in c(ensf702, ensfrgrg,ensf50, ensf70)) {
  col<- colorRampPalette(c("grey", "black", "red", "purple"))
  plot(ens, cex=0.5, col=col(4))
}
#=------------------------------------------------------------------
#Individual model algorithm ensembles at 50 and 70 time averages
#Ensemble modelling at 50s time averages
#=======================================================================
glm50<- getModelId(sdmodel, method = c("glm"))
GLM50<- ensemble(sdmodel, biof50, filename = "GLM50.img",
                 setting=list(method="weighted", stat="AUC", id=glm50))

rf50<- getModelId(sdmodel, method = c("rf"))
RF50<- ensemble(sdmodel, biof50, filename = "rf50.img",
                setting=list(method="weighted", stat="AUC", id=rf50))

gam50<- getModelId(sdmodel, method = c("gam"))
GAM50<- ensemble(sdmodel, biof50, filename = "gam50.img",
                 setting=list(method="weighted", stat="AUC", id=gam50))
svm50<- getModelId(sdmodel, method = c("svm"))
SVM50<- ensemble(sdmodel, biof50, filename = "svm50.img",
                 setting=list(method="weighted", stat="AUC", id=svm50))
brt50<- getModelId(sdmodel, method = c("brt"))
BRT50<- ensemble(sdmodel, biof50, filename = "brt50.img",
                 setting=list(method="weighted", stat="AUC", id=brt50))

#-----------------------------------------------------------------------
#Esemble for 70 time averages 
glm70<- getModelId(sdmodel, method = c("glm"))
GLM70<- ensemble(sdmodel, biof70, filename = "GLM70.img",
                 setting=list(method="weighted", stat="AUC", id=glm70))

rf70<- getModelId(sdmodel, method = c("rf"))
RF70<- ensemble(sdmodel, biof70, filename = "rf70.img",
                setting=list(method="weighted", stat="AUC", id=rf70))
gam70<- getModelId(sdmodel, method = c("gam"))
GAM70<- ensemble(sdmodel, biof70, filename = "gam70.img",
                 setting=list(method="weighted", stat="AUC", id=gam70))
svm70<- getModelId(sdmodel, method = c("svm"))
SVM70<- ensemble(sdmodel, biof70, filename = "svm70.img",
                 setting=list(method="weighted", stat="AUC", id=svm70))
brt70<- getModelId(sdmodel, method = c("brt"))
BRT70<- ensemble(sdmodel, biof70, filename = "brt70.img",
                 setting=list(method="weighted", stat="AUC", id=brt70))

#----------------------------------------------------------------------
#calcualting the differences
ensdiff1<- ensf50-ensp
ensdiff2<- ensf70-ensp
plot(stack(ensdiff1, ensdiff2))

#convert to presence absence
ev<- getEvaluation(sdmodel, stat = c("AUC", "TSS", "threshold"), opt = 2 )

meand<- mean(ev$TSS)
meand

ppan<- raster(ensp)

ppan[]<- ifelse(ensp[]>=0.65, 1, 0)

ppaf<- raster(ensf50)

ppaf[]<- ifelse(ensf50[]>=0.65, 1, 0)

ppa<- ppaf-ppan

plot(ppa)

colp<- colorRampPalette(c("black", "red", "blue"))
plot(ppa, col=colp(3))# Extinction and colonisation
save.image(file = "LatesAfrica.RData")

#====================================================================
#PART TWO
#========================================================================
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggsn)

#=========================================================================
coords<-   read.csv(file = "Perch.csv", header = TRUE, strip.white = T)%>%
  filter(!is.na(lat), continent%in%c("AFRICA"))


world <- ne_countries(scale = "medium", returnclass = "sf")
pointsw<- st_centroid(world)
world
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))


ggplot()+
  geom_sf(data = world, fill="white", color="black")+
  coord_sf(xlim = c(-19,48), ylim = c(-33,38))+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)
ggplot()+
  geom_sf(data = world, fill="white")+
  geom_sf(data = water, color="black")+
  geom_point(data = coords, aes(lon, lat), color="blue")+
  labs(x="Longitudes", y="Latitudes")+
  theme_bw()+
  theme(text = element_text(family = "Lucida Sans"))+
  coord_sf(xlim = c(-19,48), ylim = c(-33,38))
#geom_text(data = world, check_overlap = TRUE)

