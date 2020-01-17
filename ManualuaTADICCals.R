
#calculate residuls for TA and DIC analysis
#libraries
library(seacarb)
library(lme4)
library(lmerTest)
library(effects)
library(cowplot)
library(lavaan)
library(piecewiseSEM)
library(semPlot)
library(DiagrammeR)
library(nlme)
library(brms)
library(tidybayes)
library(tidyverse)
library(ggthemes)

#load data
Cdata<-read.csv('chemicaldata_maunalua.csv')
#remove rows with NAs
Cdata<-Cdata[complete.cases(Cdata),]

#calculate rest of carbonate params
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, T=Cdata$Temp_in, Patm=1, P=Cdata$Phosphate/1000000, Pt=0, Sit=0,
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
BP.end.PO<-3.7
BP.end.NN<-163
BP.end.Si<-740

#wailupe
W.end.DIC<-1780
W.end.TA<-1743
W.end.Sal<-2
W.end.PO<-1.7
W.end.NN<-71
W.end.Si<-810

#hot end point from 05/24/2015
Hot.Sal<- 35.1961
Hot.pH<-8.069
Hot.DIC<-2006.3
Hot.TA<-2324
Hot.temp<-24.4490
Hot.PO<-0.11
Hot.NN<-0.01
Hot.Si<-1.17



#Hot.CO2<-carb(flag=8, Hot.pH, Hot.TA/1000000, S=Hot.Sal, T=Hot.temp, Patm=1, P=0, Pt=0, Sit=0,
#          k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

#Hot.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-Hot.CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#create the regression

#
#TA
BP.model<-lm(c(BP.end.TA,Hot.TA)~c(BP.end.Sal,Hot.Sal))
W.model<-lm(c(W.end.TA,Hot.TA)~c(W.end.Sal,Hot.Sal))
#DIC
BP.model.DIC<-lm(c(BP.end.DIC,Hot.DIC)~c(BP.end.Sal,Hot.Sal))
W.model.DIC<-lm(c(W.end.DIC,Hot.DIC)~c(W.end.Sal,Hot.Sal))
#NN
BP.model.NN<-lm(c(BP.end.NN,Hot.NN)~c(BP.end.Sal,Hot.Sal))
W.model.NN<-lm(c(W.end.NN,Hot.NN)~c(W.end.Sal,Hot.Sal))
#PO
BP.model.PO<-lm(c(BP.end.PO,Hot.PO)~c(BP.end.Sal,Hot.Sal))
W.model.PO<-lm(c(W.end.PO,Hot.PO)~c(W.end.Sal,Hot.Sal))


#predicted TA based on mixing line
BP.TA.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model$coefficients[2]+BP.model$coefficients[1]
W.TA.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model$coefficients[2]+W.model$coefficients[1]
Cdata$TA.pred<-c(W.TA.Pred, BP.TA.Pred)

#predicted DIC based on mixing line
BP.DIC.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model.DIC$coefficients[2]+BP.model.DIC$coefficients[1]
W.DIC.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model.DIC$coefficients[2]+W.model.DIC$coefficients[1]
Cdata$DIC.pred<-c(W.DIC.Pred, BP.DIC.Pred)

#predicted NN based on mixing line
BP.NN.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model.NN$coefficients[2]+BP.model.NN$coefficients[1]
W.NN.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model.NN$coefficients[2]+W.model.NN$coefficients[1]
Cdata$NN.pred<-c(W.NN.Pred, BP.NN.Pred)


#predicted NN based on mixing line
BP.PO.Pred<-Cdata$Salinity[Cdata$Site=='BP']*BP.model.PO$coefficients[2]+BP.model.PO$coefficients[1]
W.PO.Pred<-Cdata$Salinity[Cdata$Site=='W']*W.model.PO$coefficients[2]+W.model.PO$coefficients[1]
Cdata$PO.pred<-c(W.PO.Pred, BP.PO.Pred)

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

#calculate difference between mixing line and actual NN
BP.diff.NN<- BP.NN.Pred-Cdata$NN[Cdata$Site=='BP'] #positive values are calcification and negative are dissolution
W.diff.NN<- W.NN.Pred-Cdata$NN[Cdata$Site=='W'] #positive values are calcification and negative are dissolution
NN.diff<-c(W.diff.NN,BP.diff.NN)
Cdata$NN.diff<-NN.diff

#calculate difference between mixing line and actual PO
BP.diff.PO<- BP.PO.Pred-Cdata$Phosphate[Cdata$Site=='BP'] #positive values are calcification and negative are dissolution
W.diff.PO<- W.PO.Pred-Cdata$Phosphate[Cdata$Site=='W'] #positive values are calcification and negative are dissolution
PO.diff<-c(W.diff.PO,BP.diff.PO)
Cdata$PO.diff<-PO.diff

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

ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Salinity, y = NN))+
  geom_line(aes(x=Salinity, y = NN.pred), col = 'blue')+
  ggtitle('N+N')+
  facet_wrap(~Site)


ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Salinity, y = Phosphate))+
  geom_line(aes(x=Salinity, y = PO.pred), col = 'blue')+
  ggtitle('PO')+
  facet_wrap(~Site)

####################Anaysis########################### 
# run anova to for Wailupe

#Make tide just high and low instead of H1, H2, L1, and L2
Cdata$Tide<-droplevels(Cdata$Tide) #this removes levels that don'e exist anymore (empty spaces for example)
levels(Cdata$Tide)<-c("H","H","L","L") # this makes H1 and H2 both H and same for L1 and L2



# filter out the zones so that it is only diffures, ambient, and transition
Cdata<-Cdata %>% 
  dplyr::filter(Zone != 'Offshore')%>%
  droplevels()

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



## Look at relationship between feedbacks and salinity or silicate versus deltas
ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Silicate, y = pH))+
  facet_wrap(~Site)


Cdata2<-Cdata %>%
  dplyr::filter(TA.diff < 300) # remove possible outlier

## delta TA ~ pH
ggplot(Cdata2, aes(y = TA.diff, x = pH,group = Site))+
  geom_point()+
  geom_smooth(method = "lm", aes(y = TA.diff, x = pH, group = Site))+
  facet_wrap(~Site)

## delta TA ~ aragonite
ggplot(Cdata2, aes(y = TA.diff, x = Arag,group = Site))+
  geom_point()+
  geom_smooth(method = "lm", aes(y = TA.diff, x = Arag, group = Site))+
  facet_wrap(~Site)

Cdata2 %>%
dplyr::filter(Tide == "L")%>%  
ggplot(aes(x = DIC.diff, y = pH,group = Zone))+
  geom_point(aes( col = Zone))+
  geom_smooth(method = "lm", lwd = 2, se = FALSE, aes(x = DIC.diff, y = pH, group = Zone, col = Zone))

# delta DIC versus pH
Cdata2 %>%
 # dplyr::filter(Tide == "L")%>%  
  ggplot(aes(x = DIC.diff, y = pH))+
  geom_point()+
  geom_smooth(method = "lm", lwd = 2) +
  geom_vline(xintercept = 0, lty = 2, color = "grey")+
  xlab("Change in DIC from mixing line (Net Ecosystem Production)")+
  annotate("text", label = "Net photosynthesis", x = 20, y = 8.2, size = 6, hjust = 0)+
  annotate("text", label = "Net respiration", x = -100, y = 8.2, size = 6, hjust = 0)+
  ylab(expression(pH[t]))+
  theme(axis.title=element_text(size=16,face="bold"))+
 # facet_wrap(~Zone)+
  ggsave(filename = 'Output/DIC_pH.png')

# delta DIC versus delta TA
Cdata2 %>%
  ggplot(aes(x = DIC.diff, y = TA.diff, group = Tide))+
  geom_point(aes(size = percent_sgd))+
  geom_smooth(method = "lm", lwd = 2, aes(x = DIC.diff, y = TA.diff, group = Tide, col = Tide)) +
  geom_vline(xintercept = 0, lty = 2, color = "grey")+
  geom_hline(yintercept = 0, lty = 2, color = "grey")+
  xlab(expression("Change in DIC from mixing line \n(Net Ecosystem Production)"))+
 # annotate("text", label = "Net photosynthesis", x = 20, y = 150, size = 6, hjust = 0)+
#  annotate("text", label = "Net respiration", x = -100, y = 150, size = 6, hjust = 0)+
  ylab(expression('Change in TA from mixing line \n(Net Ecosystem Calcification)'))+
  labs(size = "% SGD")+
  theme(axis.title.x = element_text(margin=margin(30,0,0,0)),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"))+
  facet_wrap(~Zone)+
  ggsave(filename = 'Output/TA_DIC_LH.png')

a<-lm(pH~DIC.diff*percent_sgd, data = Cdata, subset = Cdata$Site=='W')
anova(a)
summary(a)

## Delta TA as a function of pH

# individual models by tide
TA.pH_low<-lm(TA.diff~pH, data = Cdata2, subset = Cdata2$Tide=='L')
r2TApH.low<-round(rsquared(TA.pH_low)$R.squared,2) # extract the r-squared
TA.pH_high<-lm(TA.diff~pH, data = Cdata2, subset = Cdata2$Tide=='H')
r2TApH.high<-round(rsquared(TA.pH_high)$R.squared,2) # extract the r-squared

# predictive power of pH on delta TA decreases (i.e. lower R2) when more SGD is present

# New facet label names for supp variable
tide.labs <- c(paste("High Tide: R2 =",r2TApH.high), paste("Low Tide: R2 =",r2TApH.low))
names(tide.labs) <- c("H", "L")

Cdata2 %>%
  ggplot(aes(x = pH, y = TA.diff))+
  geom_point(aes(size = percent_sgd, color  = Day_Night))+
  geom_smooth(method = "lm", lwd = 2, aes(x = pH, y = TA.diff)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey")+
  xlab(expression(pH[t]))+
  ylab(expression('Change in TA from mixing line \n(Net Ecosystem Calcification)'))+
  labs(size = "% SGD", color = "Day or Night")+
  theme(axis.title.x = element_text(margin=margin(30,0,0,0)),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"))+
  facet_wrap(~Tide, labeller = labeller(Tide = tide.labs))+
  ggsave(filename = 'Output/TA_pH.png')



# low salinities
Cdata2 %>%
  #  filter(Salinity<33)%>%
  ggplot(aes(x = pH, y = TA.diff/2, col = Salinity, group = Site))+
  geom_point(aes(shape = Tide, size = 1))+
  scale_colour_gradient2(low = "red", high = "blue3", space = "Lab",midpoint = 33,
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
  facet_wrap(~Site)+
  theme_bw()


#NN difference versus salinity
Cdata %>%
  #  filter(Salinity<33)%>%
  ggplot(aes(y = NN.diff, x = percent_sgd,  group = Site))+
  geom_point()+
  scale_colour_gradient2(low = "red", high = "blue3", space = "Lab",midpoint = 33,
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
  facet_wrap(~Site)+
  theme_bw()

#PO difference versus salinity
Cdata %>%
  #  filter(Salinity<33)%>%
  ggplot(aes(y = PO.diff, x = Salinity, col = Salinity, group = Site))+
  geom_point(aes(size = 1))+
  scale_colour_gradient2(low = "grey", high = "blue", mid = 'lightblue', midpoint = 31, space = "Lab", limits = c(25,35),
                         breaks = c(25,27,29,31,33,35))+
  facet_wrap(~Site)+
  theme_bw()

#SEM######################
# Test the effect of SGD (salinity) on pH which is mediated by N uptake and production rates. 
# Hypothesis: High SGD (low salinity/high silicate) increases N uptake of producers, which increases production (delta DIC), which increases pH

## standardize the data
Cdata$pH_std<-scale(Cdata$pH,center = TRUE, scale = TRUE)
Cdata$Si_std<-scale(Cdata$Silicate,center = TRUE, scale = TRUE)
Cdata$DIC.diff_std<-scale(Cdata$DIC.diff,center = TRUE, scale = TRUE)
Cdata$SGD_std<-scale(Cdata$percent_sgd,center = TRUE, scale = TRUE)
Cdata$NN_std<-scale(Cdata$NN,center = TRUE, scale = TRUE)
Cdata$NH4_std<-scale(Cdata$Ammonia,center = TRUE, scale = TRUE)
Cdata$PO_std<-scale(Cdata$Phosphate,center = TRUE, scale = TRUE)

SGD_prodFormula<-'
NH4_std ~ SGD_std
PO_std ~ SGD_std
NN_std ~ SGD_std
DIC.diff_std ~  NN_std + PO_std + NH4_std
pH_std ~ DIC.diff_std + SGD_std
'
#Day
SGD_prodModel_Day <- sem(SGD_prodFormula, data = Cdata[Cdata$Day_Night=='Day',], group = 'Tide')
summary(SGD_prodModel_Day, standardize = T)
semPaths(SGD_prodModel_Day, "std", edge.label.cex = 1.5,  
         intercepts =FALSE, layout = "tree", style = "lisrel",
         residuals = FALSE)

#Night
SGD_prodModel_Night <- sem(SGD_prodFormula, data = Cdata[Cdata$Day_Night=='Night',], group = 'Tide')
summary(SGD_prodModel_Night, standardize = T)
semPaths(SGD_prodModel_Night, "std", edge.label.cex = 1.5,  
         intercepts =FALSE, layout = "tree", style = "lisrel",
         residuals = FALSE)


## try SEM with random effect
## Filter our the Offshore data
Cdata3<-Cdata %>%
#  dplyr::filter(Zone!='Offshore')%>%
#  droplevels()%>%
  dplyr::mutate(Zone1 = as.numeric(Zone)) # apparently sem doesnt like characters

a<- lme(DIC.diff ~ Silicate*Zone1 , random = ~ 1 | Site , data = Cdata3)
b<- lme(pH ~ DIC.diff*Zone1 + Silicate*Zone1, random = ~ 1 | Site, data = Cdata3)       

SGD_mod<-psem(
    a,
    b
)

summary(SGD_mod)

### model a bayesian SEM
pH_mod <- bf(pH_std ~ Tide*(DIC.diff_std + SGD_std))
DIC_mod <- bf(DIC.diff_std ~ Tide*(NN_std + PO_std + NH4_std))
NN_mod<-bf(NN_std ~ SGD_std)
NH4_mod<-bf(NH4_std ~ SGD_std)
PO_mod<-bf(PO_std ~ SGD_std)

# full mediation
#pH_mod <- bf(pH ~ DIC.diff)
#DIC_mod <- bf(DIC.diff ~ NN )
#NN_mod<-bf(NN ~ percent_sgd)

k_fit_brms <- brm(pH_mod+
                    DIC_mod+ 
                    NN_mod +
                    NH4_mod+
                    PO_mod+
                    set_rescor(FALSE), 
                  data=Cdata,
                  cores=4, chains = 2)

# view the effect sized
fixef(k_fit_brms)

## all of the marginal effects plots
plot(marginal_effects(k_fit_brms), points = TRUE)

### let's predict pH with changes in SGD.  Note, non-terminal endogenous variables need to have an NA as their values.

newdata <- data.frame(percent_sgd = 50, DIC.diff=NA, NN = 20)

pH_pred <- fitted(k_fit_brms, newdata=newdata,
                     resp = "DICdiff", nsamples = 1000, 
                     summary = FALSE)

median(pH_pred)
posterior_interval(pH_pred)

## calculating the effect size from the interaction terms
## not sure if this is correct (Read section 7.2 https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/interactions.html)
post %>%
  transmute(gamma_Low    = b_pHstd_DIC.diff_std  + `b_pHstd_TideL:DIC.diff_std`,
            gamma_High = b_pHstd_DIC.diff_std ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(mean = mean(value))

post %>%
  transmute(gamma_Low    = b_pHstd_DIC.diff_std  + `b_pHstd_TideL:DIC.diff_std`,
            gamma_High = b_pHstd_DIC.diff_std ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  
  ggplot(aes(x = value, group = key, color = key, fill = key)) +
  geom_density(alpha = 1/4) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous('Standardized effect Size', expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL) +
  ggtitle("pH ~ DIC slopes",
          subtitle = "Blue = Low Tide, Green = High Tide") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none")
