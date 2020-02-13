## Run Bayesian SEM for Maunalua carbonate chemistry data
## By: Nyssa Silbiger
## Last updated: 2/10.2020
########################################################
#libraries
library(seacarb)
library(lme4)
library(lmerTest)
library(effects)
library(cowplot)
library(nlme)
library(brms)
library(lavaan)
library(piecewiseSEM)
library(tidybayes)
library(tidyverse)
library(ggthemes)
library(patchwork)
library(bayesplot)
library(dagitty)
library(ggdag)
library(ggtext)

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

Cdata <- Cdata %>%
  filter(TA.diff < 300) %>% # remove outlier
  mutate(log_NN = log(NN),
         log_PO = log(Phosphate),
         log_SGD = log(percent_sgd)) %>% # Need to log transform the NN, PO, and SGD data because it is highly left scewed
  mutate_at(.vars = c("pH", "Silicate","DIC.diff","log_SGD","log_NN","Ammonia","log_PO","TA.diff"), .funs = list(std = ~scale(.))) #standardize all the data


# each level is a model
## SGD drives changes in N and P directly; 
# N + P + NH4 directly drive changes in DIC diff (net production);
# DIC diff and SGD directly drive changes in pH;
# pH directly drives changes in TA diff (NEC)

### model a bayesian SEM (individual one for each site since the geo chemistry of the SGD is different)
## brms does some weird things with columns names that have non-alphanumerics. Here I am removing them all
colnames(Cdata)<-str_replace_all(colnames(Cdata), "[^[:alnum:]]", "")

# run one site at a time because they have different biogeochem in the SGD
NN_mod<-bf(logNNstd ~ logSGDstd) # NN ~ SGD which can change by site
PO_mod<-bf(logPOstd ~ logSGDstd) # PO ~ SGD which can change by site
DIC_mod <- bf(DICdiffstd ~ DayNight*(poly(logNNstd,2) + poly(logPOstd,2))) # DIC ~ nutrients, which can change by day/night (i.e. high nutrients could lead to high P during the day and high R at night what have opposite signs). It is also non-linear
pH_mod <- bf(pHstd ~ DICdiffstd + logSGDstd) # pH ~ NEP + SGD
TA_mod<-bf(TAdiffstd ~ pHstd) # NEC ~ pH, which can change by day/night

#NN_mod<-bf(NN_std ~ Site*SGD_std) # NN ~ SGD which can change by site
#NH4_mod<-bf(NH4_std ~ Site*SGD_std) # NH4 ~ SGD which can change by site
#PO_mod<-bf(PO_std ~ Site*SGD_std) # PO ~ SGD which can change by site

# Run the model first for Black Point
k_fit_brms <- brm(TA_mod+
                    pH_mod+
                    DIC_mod+ 
                    NN_mod +
                    PO_mod+
                     set_rescor(FALSE), 
                   data=Cdata[Cdata$Site=='BP',],
                   cores=4, chains = 3)

# view the effect sized
fixef(k_fit_brms)

#check it
#plot(k_fit_brms)

p1<-pp_check(k_fit_brms, resp="logNNstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("NN")

p2<-pp_check(k_fit_brms, resp="logPOstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("PO")

p3<-pp_check(k_fit_brms, resp="pHstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("pH")

p4<-pp_check(k_fit_brms, resp="DICdiffstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("DIC")

p5<-pp_check(k_fit_brms, resp="TAdiffstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("TA")

# view everything together using patchwork
p1+p2+p3+p4+p5+plot_layout(guides = "collect") +plot_annotation(title = 'Black Point Posterior Predictive Checkes', tag_levels = "A")+ggsave("Output/Posteriorchecks_BlackPoint.png")
## pp checks look good!

# plot some of the conditional effects
R1<-plot(conditional_effects(k_fit_brms, "logSGDstd", resp = "logNNstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal() 
R2<-plot(conditional_effects(k_fit_brms, "logSGDstd", resp = "logPOstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
R3<-plot(conditional_effects(k_fit_brms, "logSGDstd", resp = "pHstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
R4<-plot(conditional_effects(k_fit_brms, "logNNstd:DayNight", resp = "DICdiffstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
R5<-plot(conditional_effects(k_fit_brms, "logPOstd:DayNight", resp = "DICdiffstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
R6<-plot(conditional_effects(k_fit_brms, "DICdiffstd", resp = "pHstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
R7<-plot(conditional_effects(k_fit_brms, "pHstd", resp = "TAdiffstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()

R1+R2+R3+R4+R5+R6+R7+
  plot_annotation(title = 'Marginal Effects for all Black Point Models', tag_levels = "A")+
  plot_layout(guides = "collect")+
  ggsave("Output/marginaleffects_BlackPoint.png", width = 10, height = 7)

## get the posterior
post <- posterior_samples(k_fit_brms)

# plot the coefficients
Cof1<-post %>% 
  select(starts_with("b")) %>% 
  gather() %>% 
  ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, color = "firebrick4", alpha = 1/10) +
  stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                      size = 3/4, color = "firebrick4") +
  labs(title = "Coefficients",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank())

# plot the intercepts
Cof2<-post %>% 
  select(starts_with("Intercept")) %>% 
  gather() %>% 
  ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, color = "firebrick4", alpha = 1/10) +
  stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                      size = 3/4, color = "firebrick4") +
  labs(title = "Intercepts",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank())

# plot the variances
Cof3<-post %>% 
  select(starts_with("sigma")) %>% 
  gather() %>% 
  ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, color = "firebrick4", alpha = 1/10) +
  stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                      size = 3/4, color = "firebrick4") +
  labs(title = "Sigmas",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank())

Cof1+Cof2+Cof3+
  plot_annotation(title = 'Standardized coefficients for Black Point', tag_levels = "A")+
  ggsave("Output/coefficients_BlackPoint.png", width = 10, height = 7)

# pull out the estimates and set it up to join with the DAG
estimates<-data.frame(fixef(k_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name", "DayNight"), sep = "[^[:alnum:]]+") %>%
  select(name, to, 3:7)

### deal with interaction terms
interactions<-post %>%
  transmute(logNNstd_Night    = b_DICdiffstd_polylogNNstd21 + `b_DICdiffstd_DayNightNight:polylogNNstd21`,
            logNNstd_Day = b_DICdiffstd_polylogNNstd21,
            logNNstd2_Night    = b_DICdiffstd_polylogNNstd22 + `b_DICdiffstd_DayNightNight:polylogNNstd22`,
            logNNstd2_Day = b_DICdiffstd_polylogNNstd22,
            logPOstd_Night    = b_DICdiffstd_polylogPOstd21 + `b_DICdiffstd_DayNightNight:polylogPOstd21`,
            logPOstd_Day = b_DICdiffstd_polylogPOstd21,
            logPOstd2_Night    = b_DICdiffstd_polylogPOstd22 + `b_DICdiffstd_DayNightNight:polylogPOstd22`,
            logPOstd2_Day = b_DICdiffstd_polylogPOstd22
  ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value), 
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight")) %>%
  mutate(to = "DICdiffstd") %>%
  select(name,to, 2:5)

# remove all the DIC diff params from the model in names and replace with the calculated estimates for day and night
estimates<-estimates %>%
  filter(to != 'DICdiffstd') %>%
  filter(name != 'Intercept') %>%
  bind_rows(interactions) %>%
  mutate(DayNight = replace_na(DayNight, "Both"))


### Make a DAG ####

bigger_dag <- dagify(TAdiffstd ~ pHstd,
                     pHstd ~ DICdiffstd + logSGDstd,
                     DICdiffstd ~ logNNstd + logPOstd +logNNstd2 + logPOstd2,
                     logPOstd~logSGDstd,
                     logNNstd~logSGDstd,
                     exposure = "logSGDstd",
                     outcome = "TAdiffstd",
                     labels = c("TAdiffstd" = "NEC",
                                "pHstd" = "pH",
                                "logSGDstd" = "Log % SGD",
                                "DICdiffstd" ="NEP",
                                "logPOstd" = "log PO",
                                "logNNstd" = "log NN",
                                "logNNstd2" = "log NN^2",
                                "logPOstd2" = "log PO^2")) 



ggdag(bigger_dag, use_labels = "label", text = FALSE)+
  theme_dag()

DAGdata<-bigger_dag %>%
  dag_paths() %>% #### then left join these with the effect sizes
  as_tibble() %>%
  left_join(estimates) %>%
  mutate(DayNight = replace_na(DayNight, "Both"))


# DAG during the day for BP
Day_DAG<-DAGdata %>% 
  mutate(est_sign = as.character(sign(Estimate)))%>% # add a column for pos and negative
  mutate(edge_cols = case_when(est_sign == '-1' ~ "#B2182B", 
                               est_sign == '1' ~ "#2166AC")) %>%
   # filter(!is.na(Estimate))%>%
  filter(DayNight %in% c("Both", "Day"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point() +
  geom_dag_edges(aes(edge_width = abs(Estimate)/5, edge_colour = edge_cols)) +
  geom_dag_text_repel(aes(label = label))+
 # geom_dag_edges_diagonal() +
#  geom_dag_text(col = "white") +
  theme_dag() +
  ggtitle('Black Point Daytime')+
  theme(plot.title = element_text(hjust = 0.5))

# Black Point
Night_DAG<-DAGdata %>% 
  mutate(est_sign = as.character(sign(Estimate)))%>% # add a column for pos and negative
  mutate(edge_cols = case_when(est_sign == '-1' ~ "#B2182B", 
                               est_sign == '1' ~ "#2166AC")) %>%
  # filter(!is.na(Estimate))%>%
  filter(DayNight %in% c("Both", "Night"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point() +
  geom_dag_edges(aes(edge_width = abs(Estimate)/5, edge_colour = edge_cols)) +
  geom_dag_text_repel(aes(label = label))+
  # geom_dag_edges_diagonal() +
  #  geom_dag_text(col = "white") +
  theme_dag() +
  ggtitle('Black Point Nighttime')+
  theme(plot.title = element_text(hjust = 0.5))



DAGPlot_BP<-Day_DAG +Night_DAG+plot_annotation(title = 'Paths for Black Point', 
                                   tag_levels = "A",
                                   subtitle = "Line thickness is standardized effect size and color represents + (blue) and - (red) values",
                                  )
  ggsave("Output/DAGplots_BlackPoint.png", width = 12, height = 8)

###### Run Model for Wailupe ###############
  W_fit_brms <- brm(TA_mod+
                      pH_mod+
                      DIC_mod+ 
                      NN_mod +
                      PO_mod+
                      set_rescor(FALSE), 
                    data=Cdata[Cdata$Site=='W',],
                    cores=4, chains = 3)
  
  # view the effect sized
  fixef(W_fit_brms)

  #check it
  #plot(W_fit_brms)
  
  Wp1<-pp_check(W_fit_brms, resp="logNNstd") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("NN")
  
  Wp2<-pp_check(W_fit_brms, resp="logPOstd") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("PO")
  
  Wp3<-pp_check(W_fit_brms, resp="pHstd") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("pH")
  
  Wp4<-pp_check(W_fit_brms, resp="DICdiffstd") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("DIC")
  
  Wp5<-pp_check(W_fit_brms, resp="TAdiffstd") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("TA")
  
  # view everything together using patchwork
  Wp1+Wp2+Wp3+Wp4+Wp5+plot_layout(guides = "collect") +plot_annotation(title = 'Wailupe Posterior Predictive Checkes', tag_levels = "A")+ggsave("Output/Posteriorchecks_Wailupe.png")
  ## NN and PO need some work...
  
  # plot some of the conditional effects
  WR1<-plot(conditional_effects(W_fit_brms, "logSGDstd", resp = "logNNstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal() 
  WR2<-plot(conditional_effects(W_fit_brms, "logSGDstd", resp = "logPOstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
  WR3<-plot(conditional_effects(W_fit_brms, "logSGDstd", resp = "pHstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
  WR4<-plot(conditional_effects(W_fit_brms, "logNNstd:DayNight", resp = "DICdiffstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
  WR5<-plot(conditional_effects(W_fit_brms, "logPOstd:DayNight", resp = "DICdiffstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
  WR6<-plot(conditional_effects(W_fit_brms, "DICdiffstd", resp = "pHstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
  WR7<-plot(conditional_effects(W_fit_brms, "pHstd", resp = "TAdiffstd"), points = TRUE, plot = FALSE)[[1]] +theme_minimal()
  
  WR1+WR2+WR3+WR4+WR5+WR6+WR7+
    plot_annotation(title = 'Marginal Effects for all Wailupe Models', tag_levels = "A")+
    plot_layout(guides = "collect")+
    ggsave("Output/marginaleffects_Wailupe.png", width = 10, height = 7)
  
  ## get the posterior
  Wpost <- posterior_samples(W_fit_brms)
 
  
  # plot the coefficients
  WCof1<-Wpost %>% 
    select(starts_with("b")) %>% 
    gather() %>% 
    ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
    geom_vline(xintercept = 0, color = "firebrick4", alpha = 1/10) +
    stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                        size = 3/4, color = "firebrick4") +
    labs(title = "Coefficients",
         x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
          axis.text.y  = element_text(hjust = 0),
          axis.ticks.y = element_blank())
  
  # plot the intercepts
  WCof2<-Wpost %>% 
    select(starts_with("Intercept")) %>% 
    gather() %>% 
    ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
    geom_vline(xintercept = 0, color = "firebrick4", alpha = 1/10) +
    stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                        size = 3/4, color = "firebrick4") +
    labs(title = "Intercepts",
         x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
          axis.text.y  = element_text(hjust = 0),
          axis.ticks.y = element_blank())
  
  # plot the variances
  WCof3<-Wpost %>% 
    select(starts_with("sigma")) %>% 
    gather() %>% 
    ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
    geom_vline(xintercept = 0, color = "firebrick4", alpha = 1/10) +
    stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                        size = 3/4, color = "firebrick4") +
    labs(title = "Sigmas",
         x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
          axis.text.y  = element_text(hjust = 0),
          axis.ticks.y = element_blank())
  
 WCof1+WCof2+WCof3+
    plot_annotation(title = 'Standardized coefficients for Wailupe', tag_levels = "A")+
    ggsave("Output/coefficientsWailupe.png", width = 10, height = 7)
  
 # pull out the estimates and set it up to join with the DAG
 Westimates<-data.frame(fixef(W_fit_brms)) %>%
   rownames_to_column(var = "name")%>%
   separate(name, into = c("to","name", "DayNight"), sep = "[^[:alnum:]]+") %>%
   select(name, to, 3:7)
 
 ### deal with interaction terms
 Winteractions<-Wpost %>%
   transmute(logNNstd_Night    = b_DICdiffstd_polylogNNstd21 + `b_DICdiffstd_DayNightNight:polylogNNstd21`,
             logNNstd_Day = b_DICdiffstd_polylogNNstd21,
             logNNstd2_Night    = b_DICdiffstd_polylogNNstd22 + `b_DICdiffstd_DayNightNight:polylogNNstd22`,
             logNNstd2_Day = b_DICdiffstd_polylogNNstd22,
             logPOstd_Night    = b_DICdiffstd_polylogPOstd21 + `b_DICdiffstd_DayNightNight:polylogPOstd21`,
             logPOstd_Day = b_DICdiffstd_polylogPOstd21,
             logPOstd2_Night    = b_DICdiffstd_polylogPOstd22 + `b_DICdiffstd_DayNightNight:polylogPOstd22`,
             logPOstd2_Day = b_DICdiffstd_polylogPOstd22
   ) %>%
   gather(key, value) %>%
   group_by(key) %>%
   summarise(Estimate = mean(value), Est.Error = sd(value), 
             Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
             Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
   separate(key, into = c("name", "DayNight")) %>%
   mutate(to = "DICdiffstd") %>%
   select(name,to, 2:5)
 
 # remove all the DIC diff params from the model in names and replace with the calculated estimates for day and night
 Westimates<-Westimates %>%
   filter(to != 'DICdiffstd') %>%
   filter(name != 'Intercept') %>%
   bind_rows(Winteractions) %>%
   mutate(DayNight = replace_na(DayNight, "Both"))
 
 ## Make the Wailupe DAG.  The dag is the same, but the effect sizes are different
 WDAGdata<-bigger_dag %>%
   dag_paths() %>% #### then left join these with the effect sizes
   as_tibble() %>%
   left_join(Westimates) %>%
   mutate(DayNight = replace_na(DayNight, "Both"))
 
 
 # DAG during the day for Wailupe
 WDay_DAG<-WDAGdata %>% 
   mutate(est_sign = as.character(sign(Estimate)))%>% # add a column for pos and negative
   mutate(edge_cols = case_when(est_sign == '-1' ~ "#B2182B", 
                                est_sign == '1' ~ "#2166AC")) %>%
   # filter(!is.na(Estimate))%>%
   filter(DayNight %in% c("Both", "Day"))%>%
   ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
   geom_dag_point() +
   geom_dag_edges(aes(edge_width = abs(Estimate)/5, edge_colour = edge_cols)) +
   geom_dag_text_repel(aes(label = label))+
   # geom_dag_edges_diagonal() +
   #  geom_dag_text(col = "white") +
   theme_dag() +
   ggtitle('Wailupe Daytime')+
   theme(plot.title = element_text(hjust = 0.5))
 
 # Night wailpule
 WNight_DAG<-WDAGdata %>% 
   mutate(est_sign = as.character(sign(Estimate)))%>% # add a column for pos and negative
   mutate(edge_cols = case_when(est_sign == '-1' ~ "#B2182B", 
                                est_sign == '1' ~ "#2166AC")) %>%
   # filter(!is.na(Estimate))%>%
   filter(DayNight %in% c("Both", "Night"))%>%
   ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
   geom_dag_point() +
   geom_dag_edges(aes(edge_width = abs(Estimate)/5, edge_colour = edge_cols)) +
   geom_dag_text_repel(aes(label = label))+
   # geom_dag_edges_diagonal() +
   #  geom_dag_text(col = "white") +
   theme_dag() +
   ggtitle('Wailupe Nighttime')+
   theme(plot.title = element_text(hjust = 0.5))
 
DAGPlot_W<- WDay_DAG +WNight_DAG+plot_annotation(title = 'Paths for Wailupe', 
                                    tag_levels = "A",
                                  #  subtitle = "Line thickness is standardized effect size and color represents + (blue) and - (red) values",
 )
 ggsave("Output/DAGplotsWailupe.png", width = 12, height = 8)
 
 (DAGPlot_BP+plot_annotation(title = "Black Point"))/(DAGPlot_W+plot_annotation(title = "Wailupe")) +plot_annotation(tag_levels = "A")+
   ggsave("Output/DAGplotsBoth.png", width = 12, height = 13)
 