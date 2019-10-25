
#calculate residuls for TA and DIC analysis
#libraries
library(seacarb)
library(lme4)
library(lmerTest)
library(cowplot)
library(tidyverse)
library(scales)
library(sf)
library(gstat)
library(stars)
library(lavaan)
library(piecewiseSEM)

#load data
Cdata<-read.csv('chemicaldata_maunalua.csv')
#remove rows with NAs
Cdata<-Cdata[complete.cases(Cdata),]

#calculate rest of carbonate params
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, T=Cdata$Temp_in, Patm=1, P=0, Pt=0, Sit=0,
          k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

Cdata[,c("CO2","HCO3","CO3","DIC","OmegaArag","OmegaCalcite","pCO2","fCO2")]<-
  CO2[,c("CO2","HCO3","CO3","DIC","OmegaAragonite","OmegaCalcite","pCO2","fCO2")]


#endmembers for Christina
#black point
BP.end.DIC<-3038
BP.end.TA<-2945
BP.end.Sal<-4.9
#wailupe
W.end.DIC<-1780
W.end.TA<-1743
W.end.Sal<-2

#hot end point from 05/24/2015
Hot.Sal<- 35.1961
Hot.pH<-8.069
Hot.TA<-2324
Hot.temp<-24.4490

offshore<-Cdata %>%
  group_by(Zone)%>%
  summarise_at(.vars = c("Salinity","Phosphate","Silicate","NN","Ammonia","pH","TA", "Temp_in"), .funs = mean)%>%
  filter(Zone == 'Offshore')

Hot.CO2<-carb(flag=8, Hot.pH, Hot.TA/1000000, S=Hot.Sal, T=Hot.temp, Patm=1, P=0, Pt=0, Sit=0,
              k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

Hot.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-Hot.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000


# offshore.CO2<-carb(flag=8, offshore$pH, offshore$TA/1000000, S=offshore$Salinity, T=offshore$Temp_in, Patm=1, P=offshore$Phosphate, Pt=0, Sit=0,
#               k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
# 
# offshore.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-offshore.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

# #With the offshore data
# BP.model<-lm(c(BP.end.TA,offshore.CO2$ALK)~c(BP.end.Sal,offshore.CO2$S))
# W.model<-lm(c(W.end.TA,offshore.CO2$ALK)~c(W.end.Sal,offshore.CO2$S))

#black point
BP.model<-lm(c(BP.end.TA,Hot.TA)~c(BP.end.Sal,Hot.Sal))
W.model<-lm(c(W.end.TA,Hot.TA)~c(W.end.Sal,Hot.Sal))


# #With the offshore data
# BP.model<-lm(c(BP.end.DIC,offshore.CO2$DIC)~c(BP.end.Sal,offshore.CO2$S))
# W.model<-lm(c(W.end.DIC,offshore.CO2$DIC)~c(W.end.Sal,offshore.CO2$S))


#DIC
BP.model.DIC<-lm(c(BP.end.DIC,Hot.CO2$DIC)~c(BP.end.Sal,Hot.Sal))
W.model.DIC<-lm(c(W.end.DIC,Hot.CO2$DIC)~c(W.end.Sal,Hot.Sal))


#predicted TA based on mixing line
BP.TA.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model$coefficients[2]+BP.model$coefficients[1]
W.TA.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model$coefficients[2]+W.model$coefficients[1]
Cdata$TA.pred<-c(W.TA.Pred, BP.TA.Pred)

#predicted DIC based on mixing line
BP.DIC.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model.DIC$coefficients[2]+BP.model.DIC$coefficients[1]
W.DIC.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model.DIC$coefficients[2]+W.model.DIC$coefficients[1]
Cdata$DIC.pred<-c(W.DIC.Pred, BP.DIC.Pred)

#calculate difference between mixing line and actual TA
BP.diff<- BP.TA.Pred-Cdata$TA[Cdata$Site=='BP'] #positive values are calcification and negative are dissolution
W.diff<- W.TA.Pred-Cdata$TA[Cdata$Site=='W'] #positive values are calcification and negative are dissolution
TA.diff<-c(W.diff,BP.diff)
Cdata$TA.diff<-TA.diff

#calculate difference between mixing line and actual DIC
BP.diff.DIC<- BP.DIC.Pred-Cdata$DIC[Cdata$Site=='BP'] #positive values are calcification and negative are dissolution
W.diff.DIC<- W.DIC.Pred-Cdata$DIC[Cdata$Site=='W'] #positive values are calcification and negative are dissolution
DIC.diff<-c(W.diff.DIC,BP.diff.DIC)
Cdata$DIC.diff<-DIC.diff


#### plot raw data and mixing line

ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Salinity, y = TA))+
  geom_line(aes(x=Salinity, y = TA.pred), col = 'blue')+
  ggtitle('Total Alkalinity')+
  facet_wrap(~Site)


ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Salinity, y = DIC))+
  geom_line(aes(x=Salinity, y = DIC.pred), col = 'blue')+
  ggtitle('Dissolved Inorganic Carbon')+
  facet_wrap(~Site)


##### Make a map of the deltas for each

# Cdata.BP<-Cdata %>%
#   dplyr::filter(Site =='W', Season == 'FALL', Tide == "L1")
# 
# coordinates(Cdata.BP) =~ Long+Lat
# 
# kriging_result = autoKrige(DIC.diff~1, Cdata.BP)
# plot(kriging_result, sp.layout = list(pts = list("sp.points", Cdata.BP)), justPosition = FALSE)

## pull out the results
# coords<-as.data.frame(kriging_result$krige_output@coords)
# coords$preds<-kriging_result$krige_output@data$var1.pred

# ggplot(coords, aes(x = x1, y = x2, col = preds))+
#   geom_tile()+
#  # scale_colour_gradient2(low = muted("red"), mid = "white",
#  #                        high = muted("blue"), midpoint = 0, space = "Lab",
# #                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
#   ggtitle('Low Tide Fall')


# Cdata.BP<-Cdata %>%
#   dplyr::filter(Site =='BP', Season == 'FALL', Tide == "H1")
# 
# coordinates(Cdata.BP) =~ Long+Lat
# 
# 
# kriging_result = autoKrige(TA.diff~1, Cdata.BP)
# plot(kriging_result, sp.layout = list(pts = list("sp.points", Cdata.BP)), justPosition = TRUE)
# 
# coords<-as.data.frame(kriging_result$krige_output@coords)
# coords$preds<-kriging_result$krige_output@data$var1.pred
# 
# ggplot(coords, aes(x = x1, y = x2, col = preds))+
#   geom_tile()+
#   scale_colour_gradient2(low = muted("red"), mid = "white",
#                          high = muted("blue"), midpoint = 0, space = "Lab",
#                          na.value = "grey50", guide = "colourbar", aesthetics = "colour")
# 
# #kriging_result$krige_output@data$var1.pred
# 
# #create a variogram
# # create coordinates
# coordinates(Cdata) = ~Long+Lat
# 
# v = variogram(TA.diff~1, data = Cdata)
# plot(v, plot.numbers = TRUE)
# 
# v.m = fit.variogram(v, vgm(1, "Exp", 2000,1))
# plot(v, v.m, plot.numbers = TRUE)


## plot the differences with points
ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Long, y = Lat, color = TA.diff, size = 1))+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = "blue", midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
  facet_wrap(Site~Tide+Season, scales = 'free')


# Just wailupe for DIC
Cdata %>%
  filter(Site=='W')%>%
ggplot()+
  geom_point(aes(x = Long, y = Lat, color = DIC.diff, size = 1))+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = "blue", midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
  facet_wrap(Season~Tide)+
  ggtitle('Wailupe')

# Just BP for DIC
Cdata %>%
  filter(Site=='BP')%>%
  ggplot()+
  geom_point(aes(x = Long, y = Lat, color = DIC.diff, size = 1))+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = "blue", midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
  facet_wrap(Season~Tide)+
  ggtitle('Black Point')


####################Anaysis########################### 

