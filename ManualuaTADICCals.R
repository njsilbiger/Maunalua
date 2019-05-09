
#calculate residuls for TA and DIC analysis
#libraries
library(seacarb)
library(lme4)
library(lmerTest)
library(effects)
library(cowplot)

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

Hot.CO2<-carb(flag=8, Hot.pH, Hot.TA/1000000, S=Hot.Sal, T=Hot.temp, Patm=1, P=0, Pt=0, Sit=0,
          k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

Hot.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-Hot.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#calculate a mixing line for TA
#Black point
par(mfrow=c(1,2))
plot(c(BP.end.Sal,Hot.Sal), c(BP.end.TA,Hot.TA), type='l', col = 'red', main = 'Black Point', xlab='Salinity', ylab='TA', xlim = c(0,36), ylim = c(2100,3000))
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='BP' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'])
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='BP' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'], pch=19)
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='BP' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'], pch=0)
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='BP' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'], pch=15)
#wailupe
plot(c(W.end.Sal,Hot.Sal), c(W.end.TA,Hot.TA), col='blue', main = 'Wailupe', type='l', xlab='Salinity', ylab='TA')
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='W' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'])
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='W' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'], pch=19)
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='W' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'], pch=0)
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'],Cdata$TA[Cdata$Site=='W' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'], pch=15)

#DIC
#Black point
par(mfrow=c(1,2))
plot(c(BP.end.Sal,Hot.Sal), c(BP.end.DIC, Hot.CO2$DIC), type='l', col = 'red', main = 'Black Point', xlab='Salinity', ylab='DIC', xlim = c(0,36), ylim = c(1800,3200))
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='BP' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'])
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='BP' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'], pch=19)
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='BP' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'], pch=0)
points(Cdata$Salinity[Cdata$Site=='BP' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='BP' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'], pch=15)
#wailupe
plot(c(W.end.Sal,Hot.Sal), c(W.end.DIC,Hot.CO2$DIC), col='blue', main = 'Wailupe', type='l', xlab='Salinity', ylab='DIC')
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='W' & Cdata$Tide=='L1' & Cdata$Season=='SPRING'])
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='W' & Cdata$Tide=='L2' & Cdata$Season=='SPRING'], pch=19)
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='W' & Cdata$Tide=='H1' & Cdata$Season=='SPRING'], pch=0)
points(Cdata$Salinity[Cdata$Site=='W' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'],Cdata$DIC[Cdata$Site=='W' & Cdata$Tide=='H2' & Cdata$Season=='SPRING'], pch=15)

#create the regression

#black point
BP.model<-lm(c(BP.end.TA,Hot.TA)~c(BP.end.Sal,Hot.Sal))
W.model<-lm(c(W.end.TA,Hot.TA)~c(W.end.Sal,Hot.Sal))
#DIC
BP.model.DIC<-lm(c(BP.end.DIC,Hot.CO2$DIC)~c(BP.end.Sal,Hot.Sal))
W.model.DIC<-lm(c(W.end.DIC,Hot.CO2$DIC)~c(W.end.Sal,Hot.Sal))

#predicted TA based on mixing line
BP.TA.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model$coefficients[2]+BP.model$coefficients[1]
W.TA.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model$coefficients[2]+W.model$coefficients[1]


#predicted DIC based on mixing line
BP.DIC.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model.DIC$coefficients[2]+BP.model.DIC$coefficients[1]
W.DIC.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model.DIC$coefficients[2]+W.model.DIC$coefficients[1]


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


####################Anaysis########################### 
# run anova to for Wailupe

#Make tide just high and low instead of H1, H2, L1, and L2
Cdata$Tide<-droplevels(Cdata$Tide) #this removes levels that don'e exist anymore (empty spaces for example)
levels(Cdata$Tide)<-c("H","H","L","L") # this makes H1 and H2 both H and same for L1 and L2

modW<-lmer(TA.diff~DIC.diff*Zone*Tide +(1|Season)+(1|Waypoint), data = Cdata[Cdata$Site=='W',])
anova(modW)


#make a prediction plot with tide on the sample plot
# easily to get a dataframe to do it yourself with ggplot2:
ne.effect <- data.frame(effect("DIC.diff*Zone*Tide", modW,
                               xlevels = list(DIC.diff = -150:350,
                                              Zone = unique(Cdata$Zone))))
#make the plot for wailupe
W.plot<-ggplot(ne.effect, aes(x = DIC.diff, y = fit, group = Tide, col = Tide)) +
  theme_bw()  +
  geom_line(aes(linetype = Tide), size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, linetype = Tide),
              color = "black", fill = "grey", alpha = 0.1) +
  geom_vline(xintercept = 0, lty=2)+ # add a vertical line at 0
  geom_hline(yintercept = 0, lty=2)+ # add a horizontal line at 0
  xlab('DIC residuals from mixing line')+
  ylab('TA residuals from mixing line')+
  ggtitle('Wailupe')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the gridlines
  facet_wrap(~ Zone)


#make a prediction plot with zone on the sample plot
# easily to get a dataframe to do it yourself with ggplot2:
ne.effect.zone <- data.frame(effect("DIC.diff*Zone*Tide", modW,
                               xlevels = list(DIC.diff = -150:350,
                                              Tide = unique(Cdata$Tide))))
#make the plot for wailupe
W.plot.zone<-ggplot(ne.effect, aes(x = DIC.diff, y = fit, group = Zone, col = Zone)) +
  theme_bw()  +
  geom_line(aes(linetype = Zone), size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, linetype = Zone),
              color = "black", fill = "grey", alpha = 0.1) +
  geom_vline(xintercept = 0, lty=2)+ # add a vertical line at 0
  geom_hline(yintercept = 0, lty=2)+ # add a horizontal line at 0
  xlab('DIC residuals from mixing line')+
  ylab('TA residuals from mixing line')+
  ggtitle('Wailupe')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the gridlines
  facet_wrap(~ Tide)

# run anova to for BP
modBP<-lmer(TA.diff~DIC.diff*Zone*Tide +(1|Season)+(1|Waypoint), data = Cdata[Cdata$Site=='BP',])
anova(modBP)

# easily to get a dataframe to do it yourself with ggplot2:
ne.effect.BP <- data.frame(effect("DIC.diff*Zone*Tide", modBP,
                               xlevels = list(DIC.diff = -150:350,
                                              Zone = unique(Cdata$Zone))))
#make the plot for wailupe with tide on the same plot
BP.plot<-ggplot(ne.effect.BP, aes(x = DIC.diff, y = fit, group = Tide, col = Tide)) +
  theme_bw()  +
  geom_line(aes(linetype = Tide), size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, linetype = Tide),
              color = "black", fill = "grey", alpha = 0.1) +
  geom_vline(xintercept = 0, lty=2)+ # add a vertical line at 0
  geom_hline(yintercept = 0, lty=2)+ # add a horizontal line at 0
  xlab('DIC residuals from mixing line')+
  ylab('TA residuals from mixing line')+
  ggtitle('Black Point')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the gridlines
  facet_wrap(~ Zone)


#make a prediction plot with zone on the sample plot
# easily to get a dataframe to do it yourself with ggplot2:
ne.effect.BP.zone <- data.frame(effect("DIC.diff*Zone*Tide", modBP,
                                    xlevels = list(DIC.diff = -150:350,
                                                   Tide = unique(Cdata$Tide))))
#make the plot for BP
BP.plot.zone<-ggplot(ne.effect.BP.zone, aes(x = DIC.diff, y = fit, group = Zone, col = Zone)) +
  theme_bw()  +
  geom_line(aes(linetype = Zone), size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, linetype = Zone),
              color = "black", fill = "grey", alpha = 0.1) +
  geom_vline(xintercept = 0, lty=2)+ # add a vertical line at 0
  geom_hline(yintercept = 0, lty=2)+ # add a horizontal line at 0
  xlab('DIC residuals from mixing line')+
  ylab('TA residuals from mixing line')+
  ggtitle('Black Point')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the gridlines
  facet_wrap(~ Tide)

##put both plots together using plot_grid
AllTide<-plot_grid(BP.plot, W.plot, labels=c('A', 'B'), nrow = 2)
#plot by zone
AllZone<-plot_grid(BP.plot.zone, W.plot.zone, labels=c('A', 'B'), nrow = 2)
#save the output in the output folder
ggsave(filename = 'Output/PlotsbyTide.pdf', plot = AllTide,device = 'pdf', width = 10, height = 8)
ggsave(filename = 'Output/PlotsbyZone.pdf', plot = AllZone,device = 'pdf', width = 8, height = 8)
