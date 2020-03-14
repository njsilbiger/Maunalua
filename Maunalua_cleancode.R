## Run Bayesian SEM for Maunalua carbonate chemistry data
## By: Nyssa Silbiger
## Last updated: 2/23/2020
########################################################
#libraries
library(seacarb)
library(lme4)
library(lmerTest)
library(effects)
library(cowplot)
library(nlme)
library(brms)
library(tidybayes)
library(tidyverse)
library(ggthemes)
library(patchwork)
library(bayesplot)
library(dagitty)
library(ggdag)
library(ggtext)
library(modelr)

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
#Spring
Hot.Sal<- 35.1961
Hot.pH<-8.069
Hot.DIC<-2006.3
Hot.TA<-2324
Hot.temp<-24.4490
Hot.PO<-0.11
Hot.NN<-0.01
Hot.Si<-1.17

#Fall (October 14th, 2015)
Hot.Sal.Fall<-35.2708
Hot.pH.Fall<-8.077
Hot.DIC.Fall<-1996.8
HOT.TA.Fall<-2320

#create the regression

#
#TA
# BP.modelspring<-lm(c(BP.end.TA,Hot.TA)~c(BP.end.Sal,Hot.Sal))
# W.modelspring<-lm(c(W.end.TA,Hot.TA)~c(W.end.Sal,Hot.Sal))
# BP.modelFall<-lm(c(BP.end.TA,HOT.TA.Fall)~c(BP.end.Sal,Hot.Sal.Fall))
# W.modelFall<-lm(c(W.end.TA,HOT.TA.Fall)~c(W.end.Sal,Hot.Sal.Fall))

# #DIC
# BP.model.DICspring<-lm(c(BP.end.DIC,Hot.DIC)~c(BP.end.Sal,Hot.Sal))
# W.model.DICspring<-lm(c(W.end.DIC,Hot.DIC)~c(W.end.Sal,Hot.Sal))
# BP.model.DICFall<-lm(c(BP.end.DIC,Hot.DIC.Fall)~c(BP.end.Sal,Hot.Sal.Fall))
# W.model.DICFall<-lm(c(W.end.DIC,Hot.DIC.Fall)~c(W.end.Sal,Hot.Sal.Fall))



# Cdata %>%
#   group_by(Site, Season)%>% # first calculate the mean values by site and season
#   summarise(Sal.mix = mean(Salinity, na.rm=TRUE),
#             TA.mix = mean(TA, na.rm=TRUE),
#             Si.mix = mean(Silicate, na.rm = TRUE)) %>%
#   left_join(Cdata) # join it with the original dataset

#predicted TA based on mixing line
## use cristina's methods  C1 = Cmix + (Cmix – Csgd)(((Smix – 35.2)/(Ssgd – Smix))  
Cdata<-Cdata %>% # calculate predicted data from mixing line based on salinity for each site and season
  #TA
  mutate(TA.pred = case_when(Site == 'BP'~ TA+(TA - BP.end.TA)*((Salinity - 35.2)/(BP.end.Sal - Salinity)),
                             Site == 'W' ~ TA+(TA - W.end.TA)*((Salinity - 35.2)/(W.end.Sal - Salinity))))%>%
  #DIC
  mutate(DIC.pred = case_when(Site == 'BP'~ DIC+(DIC - BP.end.DIC)*((Salinity - 35.2)/(BP.end.Sal - Salinity)),
                             Site == 'W' ~ DIC+(DIC - W.end.DIC)*((Salinity - 35.2)/(W.end.Sal - Salinity))))%>%
  #differences
  mutate(TA.diff = (Hot.TA-TA.pred)/2, #positive values are calcification and negative are dissolution
         DIC.diff = Hot.DIC - DIC) #positive values are net photosynthesis and negative are respiration

Cdata<-Cdata %>% # calculate predicted data from mixing line based on salinity for each site and season
  #TA
  mutate(TA.pred = case_when(Site == 'BP'~ TA+(TA - BP.end.TA)*((Silicate - Hot.Si)/(BP.end.Si - Silicate)),
                             Site == 'W' ~ TA+(TA - W.end.TA)*((Silicate - Hot.Si)/(W.end.Si - Silicate))))%>%
  #DIC
  mutate(DIC.pred = case_when(Site == 'BP'~ DIC+(DIC - BP.end.DIC)*((Silicate - Hot.Si)/(BP.end.Si - Silicate)),
                              Site == 'W' ~ DIC+(DIC - W.end.DIC)*((Silicate - Hot.Si)/(W.end.Si - Silicate))))%>%
  #differences
  mutate(TA.diff = (Hot.TA-TA.pred)/2, #positive values are calcification and negative are dissolution
         DIC.diff = Hot.DIC - DIC) #positive values are net photosynthesis and negative are respiration

#Spring
# Cdata<-Cdata %>% # calculate predicted data from mixing line based on salinity for each site and season
#  #TA
#   mutate(TA.pred = case_when(Site == 'BP' & Season =='SPRING' ~ Salinity*BP.modelspring$coefficients[2]+BP.modelspring$coefficients[1],
#                                 Site == 'BP' & Season =='FALL' ~ Salinity*BP.modelFall$coefficients[2]+BP.modelFall$coefficients[1],
#                                 Site == 'W' & Season =='SPRING' ~ Salinity*W.modelspring$coefficients[2]+W.modelspring$coefficients[1],
#                                 Site == 'W' & Season =='FALL' ~ Salinity*W.modelspring$coefficients[2]+W.modelspring$coefficients[1]))%>%
# #DIC
#   mutate(DIC.pred = case_when(Site == 'BP' & Season =='SPRING' ~ Salinity*BP.model.DICspring$coefficients[2]+BP.model.DICspring$coefficients[1],
#                              Site == 'BP' & Season =='FALL' ~ Salinity*BP.model.DICFall$coefficients[2]+BP.model.DICFall$coefficients[1],
#                              Site == 'W' & Season =='SPRING' ~ Salinity*W.model.DICspring$coefficients[2]+W.model.DICspring$coefficients[1],
#                              Site == 'W' & Season =='FALL' ~ Salinity*W.model.DICFall$coefficients[2]+W.model.DICFall$coefficients[1]))%>%
#   #differences
#   mutate(TA.diff = TA.pred - TA, #positive values are calcification and negative are dissolution
#          DIC.diff = DIC.pred - DIC) #positive values are net photosynthesis and negative are respiration
# 
# 


#### plot raw data and mixing line

ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Salinity, y = TA.pred, color = Day_Night))+ 
  coord_trans(x="log", y="log")+
  #geom_line(aes(x=Salinity, y = TA.pred), col = 'blue')+
  ggtitle('Total Alkalinity')+
  facet_wrap(~Site*Season)


ggplot(Cdata, aes(group = Site))+
  geom_point(aes(x = Salinity, y = DIC.pred, color = Day_Night))+
 # geom_line(aes(x=Salinity, y = DIC.pred), col = 'blue')+
  ggtitle('Dissolved Inorganic Carbon')+
  facet_wrap(~Site*Season)


####################Anaysis########################### 
# run anova to for Wailupe

#Make tide just high and low instead of H1, H2, L1, and L2
Cdata$Tide<-droplevels(Cdata$Tide) #this removes levels that don'e exist anymore (empty spaces for example)
levels(Cdata$Tide)<-c("H","H","L","L") # this makes H1 and H2 both H and same for L1 and L2

# filter out the zones so that it is only diffures, ambient, and transition
Cdata<-Cdata %>% 
  dplyr::filter(Zone != 'Offshore')%>%
  droplevels()


#SEM######################
# Test the effect of SGD (salinity) on pH which is mediated by N uptake and production rates. 
# Hypothesis: High SGD (low salinity/high silicate) increases N uptake of producers, which increases production (delta DIC), which increases pH

Cdata <- Cdata %>%
  filter(TA.diff < 150 & TA.diff > -10) %>% # remove outlier
  mutate(log_NN = log(NN),
         log_PO = log(Phosphate),
         log_SGD = log(percent_sgd),
         log_Salinity = log(Salinity)) %>% # Need to log transform the NN, PO, and SGD data because it is highly left scewed
  mutate_at(.vars = c("pH", "Silicate","DIC.diff","log_SGD","log_NN","Ammonia","log_PO","TA.diff", "Salinity", "log_Salinity", "Temp_in"), .funs = list(std = ~scale(.))) #standardize all the data

## change the factors for prettier names
Cdata<-Cdata %>%
  mutate(Tide = recode(Tide, L = "Low Tide",  
                       H = "High Tide"))  ## change the factors for prettier names
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
Temp_mod<-bf(Tempinstd ~ DayNight*logSGDstd*Season) ## SGD has cooler water and the intercept changes with season
#PO_mod<-bf(logPOstd ~ logSGDstd) # PO ~ SGD which can change by site
DIC_mod <- bf(DICdiffstd ~ DayNight*(poly(logNNstd,2)+poly(Tempinstd,2))) # DIC ~ nutrients and temperature, which can change by day/night (i.e. high nutrients could lead to high P during the day and high R at night what have opposite signs). It is also non-linear
#DIC_mod <- bf(DICdiffstd ~ logNNstd*Tempinstd) # DIC ~ nutrients and temperature, which can change by day/night (i.e. high nutrients could lead to high P during the day and high R at night what have opposite signs). It is also non-linear

pH_mod <- bf(pHstd ~ DICdiffstd + logSGDstd) # pH ~ NEP + SGD
TA_mod<-bf(TAdiffstd ~ pHstd+poly(Tempinstd,2)) # NEC ~ pH and temperature, which can change by Tide, because low has more nutrients it may distrupt this relationship (based on Silbiger et al. 2018)
#TA_mod<-bf(TAdiffstd ~ pHstd+Tempinstd +(1|Season)) # NEC ~ pH and temperature, which can change by Tide, because low has more nutrients it may distrupt this relationship (based on Silbiger et al. 2018)


# Run the model first for Black Point
k_fit_brms <- brm(TA_mod+
                    pH_mod+
                    DIC_mod+ 
                    Temp_mod+
                    NN_mod +
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

p2<-pp_check(k_fit_brms, resp="pHstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("pH")

p3<-pp_check(k_fit_brms, resp="DICdiffstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("DIC")

p4<-pp_check(k_fit_brms, resp="TAdiffstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("TA")

p5<-pp_check(k_fit_brms, resp="Tempinstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Temperature")

# view everything together using patchwork
p1+p2+p3+p4+p5+plot_layout(guides = "collect") +plot_annotation(title = 'Black Point Posterior Predictive Checks', tag_levels = "A")+ggsave("Output/Posteriorchecks_BlackPoint.png")
## pp checks look good!

# plot the conditional effects
conditions <- make_conditions(k_fit_brms, "Season") # for the three way interaction

R<-conditional_effects(k_fit_brms, "logSGDstd", resp = "logNNstd", method = "predict", resolution = 1000)
R1<-R$logNNstd.logNNstd_logSGDstd %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         lower = lower__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         upper = upper__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = exp(logSGD), y = exp(estimate)), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = exp(logSGD),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = percentsgd, y = NN), alpha = 0.1) +
  xlab("Percent SGD")+
  ylab(expression(paste("Nitrate + Nitrite (mmol L"^-1,")")))+
  coord_trans(x="log", y="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()

R<-conditional_effects(k_fit_brms, "logSGDstd:DayNight", resp = "Tempinstd",  conditions = conditions,method = "predict", resolution = 1000)
R2<-R$Tempinstd.Tempinstd_logSGDstd%>%
  mutate(estimate = estimate__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         lower = lower__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         upper = upper__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate, color = DayNight), lwd = 2)+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = percentsgd, y = Tempin, color = DayNight), alpha = 0.1) +
  xlab("Percent SGD")+
  ylab(expression(paste("Temperature (", degree, "C)")))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  theme_minimal()+
  facet_wrap(~Season)

R<-conditional_effects(k_fit_brms, "logSGDstd", resp = "pHstd", method = "predict", resolution = 1000)
R3<-R$pHstd.pHstd_logSGDstd %>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = percentsgd, y = pH), alpha = 0.1) +
  xlab("Percent SGD")+
  ylab(expression("pH"[t]))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  theme_minimal()

R<-conditional_effects(k_fit_brms, "logNNstd:DayNight", resp = "DICdiffstd", conditions = conditions, method = "predict", resolution = 1000)
R4<-R$`DICdiffstd.DICdiffstd_logNNstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         logNN = logNNstd*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logNN), y = estimate, group = DayNight, color = DayNight), lwd = 2)+
  geom_ribbon(aes(x = exp(logNN),ymin=lower, ymax=upper, group = DayNight, fill = DayNight), linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = NN, y = DICdiff, color = DayNight), alpha = 0.1) +
  xlab(expression(paste("Nitrate + Nitrite (mmol L"^-1,")")))+
  ylab(expression(paste("Net Ecosystem Production ( ", Delta, "DIC ", mu,"mol kg"^-1, ")")))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()

R<-conditional_effects(k_fit_brms, "Tempinstd:DayNight", resp = "DICdiffstd", conditions = conditions, method = "predict", resolution = 1000)
R5<-R$`DICdiffstd.DICdiffstd_Tempinstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate, group = DayNight, color = DayNight), lwd = 2)+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper, group = DayNight, fill = DayNight), linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = Tempin, y = DICdiff, color = DayNight), alpha = 0.1) +
  xlab(expression(paste("Temperature (", degree, "C)")))+
  ylab(expression(paste("Net Ecosystem Production ( ", Delta, "DIC ", mu,"mol kg"^-1, ")")))+
  theme_minimal()

R<-conditional_effects(k_fit_brms, "DICdiffstd", resp = "pHstd", method = "predict", resolution = 1000)
R6<-R$pHstd.pHstd_DICdiffstd%>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         DICdiff = DICdiffstd*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = DICdiff, y = estimate), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = DICdiff,ymin=lower, ymax=upper), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = DICdiff, y = pH), alpha = 0.1) +
  xlab(expression(paste("Net Ecosystem Production ( ", Delta, "DIC ", mu,"mol kg"^-1, ")")))+
  ylab(expression("pH"[t]))+
  theme_minimal()

R<-conditional_effects(k_fit_brms, "pHstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
R7<-R$TAdiffstd.TAdiffstd_pHstd%>%
  mutate(estimate = estimate__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         pH = pHstd*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = pH, y = estimate), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = pH,ymin=lower, ymax=upper), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = pH, y = TAdiff), alpha = 0.1) +
  xlab(expression("pH"[t]))+
  ylab(expression(paste("Net Ecosystem Calcification ( ", Delta, "TA ", mu,"mol kg"^-1, ")")))+
  theme_minimal()

R<-conditional_effects(k_fit_brms, "Tempinstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
R8<-R$TAdiffstd.TAdiffstd_Tempinstd%>%
  mutate(estimate = estimate__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate), lwd = 2, color = "blue")+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper),fill = "blue", linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = Tempin, y = TAdiff), alpha = 0.1) +
  xlab(expression(paste("Temperature (", degree, "C)")))+
  ylab(expression(paste("Net Ecosystem Calcification", " (",Delta, "TA ", mu,"mol kg"^-1, ")")))+
  theme_minimal()

R1+R2+R3+R4+R5+R6+R7+R8+
  plot_annotation(title = 'Marginal Effects for all Black Point Models', tag_levels = "A")+
  plot_layout(guides = "collect")+
  ggsave("Output/marginaleffects_BlackPoint.png", width = 12, height = 10)

## get the posterior
post <- posterior_samples(k_fit_brms)

# plot the coefficients
Cof1<-post %>% 
  select(starts_with("b"),-ends_with("Intercept")) %>%
  gather() %>% 
  group_by(key)%>%
  median_hdci()%>%
  mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
  separate(col = key,into = c("b", "dependent", "independent"),sep = "_")%>% #loose the b and bring the values back together
  mutate(key = paste(dependent, independent))%>%
  ggplot(aes(x = value, y = reorder(key, value), color = sig)) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point()+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  scale_color_manual(values = c("grey","firebrick4"))+
  # stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
  #                     size = 3/4, color = "firebrick4") +
  labs(title = "Black Point",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none")

# plot the intercepts
Cof2<-post %>% 
  select(starts_with("Intercept"),ends_with("Intercept")) %>% 
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
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(name != 'Intercept', to != 'DICdiffstd', to !="Tempinstd") # remove the intercepts and interaction terms (Calculated below) for the DAG
  
### deal with interaction terms. First for DIC models
InteractionsDIC<-post %>%
  transmute(logNNstd_Night    = `b_DICdiffstd_polylogNNstd21` + `b_DICdiffstd_DayNightNight:polylogNNstd21`,
            logNNstd_Day = `b_DICdiffstd_polylogNNstd21`,
            logNNstd2_Night    = `b_DICdiffstd_polylogNNstd22` + `b_DICdiffstd_DayNightNight:polylogNNstd22`,
            logNNstd2_Day = `b_DICdiffstd_polylogNNstd22`,
            Tempinstd_Night    = `b_DICdiffstd_polyTempinstd21` + `b_DICdiffstd_DayNightNight:polyTempinstd21`,
            Tempinstd_Day = `b_DICdiffstd_polyTempinstd21`,
            Tempinstd2_Night    = `b_DICdiffstd_polyTempinstd22` + `b_DICdiffstd_DayNightNight:polyTempinstd22`,
            Tempinstd2_Day = `b_DICdiffstd_polyTempinstd22`
          ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Season")) %>%
  mutate(to = "DICdiffstd") %>%
  select(name,to, 2:7)

# Interactions temp SGD
InteractionsSGD<-post %>%
  transmute(
            logSGDstd_Day_Spring = `b_Tempinstd_logSGDstd` +`b_Tempinstd_logSGDstd:SeasonSPRING`,
            logSGDstd_Night_Spring = `b_Tempinstd_logSGDstd` +`b_Tempinstd_DayNightNight:logSGDstd:SeasonSPRING` +`b_Tempinstd_DayNightNight`,
            logSGDstd_Day_Fall = `b_Tempinstd_logSGDstd`,
            logSGDstd_Night_Fall = `b_Tempinstd_logSGDstd`+`b_Tempinstd_DayNightNight`
  ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Season")) %>%
  mutate(to = "Tempinstd") %>%
  select(name,to, 2:7)

#Bind everything together
  estimates<-bind_rows(estimates, InteractionsDIC, InteractionsSGD) %>% # bind with the estimates above
  mutate(DayNight = replace_na(DayNight, "Both"),
         Season = replace_na(Season, "Both"))
  # change the "poly names'
  estimates$name[estimates$name=='polyTempinstd21'] = 'Tempinstd'
  estimates$name[estimates$name=='polyTempinstd22'] = 'Tempinstd2'


### Make a DAG ####
  bigger_dag <- dagify(TAdiffstd ~ pHstd+ Tempinstd+ Tempinstd2,
                       pHstd ~ DICdiffstd + logSGDstd,
                       Tempinstd ~ logSGDstd,
                       DICdiffstd ~ logNNstd + Tempinstd +logNNstd2 + Tempinstd2,
                       logNNstd ~ logSGDstd,
                       exposure = "logSGDstd",
                       outcome = "TAdiffstd",
                       labels = c("TAdiffstd" = "NEC",
                                  "pHstd" = "pH",
                                  "Tempinstd" = "Temperature",
                                  "Tempinstd2" = "Temperature^2",
                                  "logSGDstd" = "Log % SGD",
                                  "DICdiffstd" ="NEP",
                                  "logNNstd" = "log NN",
                                  "logNNstd2" = "log NN^2")) 
  
  
  #quick visual
  ggdag(bigger_dag, use_labels = "label", text = FALSE)+
    theme_dag()
  
  # join it with the estimates so that I can add colors and line thickness to related to effectsize
  DAGdata<-bigger_dag %>%
    dag_paths() %>% #### then left join these with the effect sizes
    as_tibble() %>%
    left_join(estimates) %>%
    mutate(DayNight = replace_na(DayNight, "Both"),
           Season = replace_na(Season, "Both"))%>%
    mutate(est_sign = as.character(sign(Estimate)))%>% # add a column for pos and negative
    mutate(edge_cols = case_when(est_sign == '-1' ~ "#d6604d",#"#fc8d59", 
                                 est_sign == '1' ~ "#4393c3")) %>%
    mutate(edge_lines = ifelse(sign(Q2.5)==sign(Q97.5),1,2),# add a dashed or solid line for significant effects (i.e. 95%CI does not overlap 0)
           edge_alpha = ifelse(sign(Q2.5)==sign(Q97.5),1,0.1)# make non-significant transparent
           ) 

  ## Basic DAG
  DAGdata %>% 
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Spring"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(
                       #edge_colour = 'grey'
                       #edge_linetype = edge_lines,
                       #edge_alpha = edge_alpha)
                       ) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggsave(filename = 'Output/BasicDAG.pdf', width = 6, height = 6)
    #ggtitle('Spring Daytime')+
    #theme(plot.title = element_text(hjust = 0.5))
  
  #Day Spring DAG
  DaySpring_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Spring"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_point() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Spring Daytime')+
    theme(plot.title = element_text(hjust = 0.5))
                                                                    
  
  #Night Spring DAG
  NightSpring_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Spring"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_point() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Spring Nighttime')+
    theme(plot.title = element_text(hjust = 0.5))
  
  #Day Fall DAG
  DayFall_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Fall"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_point() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Fall Daytime')+
    theme(plot.title = element_text(hjust = 0.5))

  #Night Fall DAG
  NightFall_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Fall"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
     geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_point() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Fall Nighttime')+
    theme(plot.title = element_text(hjust = 0.5)) 
  
  DAGPlot_BP<- (DaySpring_DAG +NightSpring_DAG)/(DayFall_DAG +NightFall_DAG)+
    plot_annotation(tag_levels = "A",title = "Black Point")+
    ggsave("Output/DAGplotsBP.pdf", width = 12, height = 13, useDingbats = FALSE)
  
  ###### Run Model for Wailupe (Need to add interaction terms with high and low tide in the model.. so many interactions) ###############

W_fit_brms <- brm(TA_mod+
                    pH_mod+
                    DIC_mod+ 
                    Temp_mod+
                    NN_mod +
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

Wp2<-pp_check(W_fit_brms, resp="pHstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("pH")

Wp3<-pp_check(W_fit_brms, resp="DICdiffstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("DIC")

Wp4<-pp_check(W_fit_brms, resp="TAdiffstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("TA")

Wp5<-pp_check(W_fit_brms, resp="Tempinstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Temperature")


# view everything together using patchwork
Wp1+Wp2+Wp3+Wp4+Wp5+plot_layout(guides = "collect") +plot_annotation(title = 'Wailupe Posterior Predictive Checks', tag_levels = "A")+ggsave("Output/Posteriorchecks_Wailupe.png")
## NN and PO need some work...

# plot some of the conditional effects
conditions <- make_conditions(W_fit_brms, "Season")

W<-conditional_effects(W_fit_brms, "logSGDstd", resp = "logNNstd", method = "predict", resolution = 1000)
WR1<-W$logNNstd.logNNstd_logSGDstd %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         lower = lower__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         upper = upper__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = exp(logSGD), y = exp(estimate)), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = exp(logSGD),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = percentsgd, y = NN), alpha = 0.1) +
  xlab("Percent SGD")+
  ylab(expression(paste("Nitrate + Nitrite (mmol L"^-1,")")))+
  coord_trans(x="log", y="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()


W<-conditional_effects(W_fit_brms, "logSGDstd:DayNight", resp = "Tempinstd",  conditions = conditions,method = "predict", resolution = 1000)
WR2<-W$Tempinstd.Tempinstd_logSGDstd%>%
  mutate(estimate = estimate__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         lower = lower__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         upper = upper__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate, color = DayNight), lwd = 2)+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = percentsgd, y = Tempin, color = DayNight), alpha = 0.1) +
  xlab("Percent SGD")+
  ylab(expression(paste("Temperature (", degree, "C)")))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  theme_minimal()+
  facet_wrap(~Season)


W<-conditional_effects(W_fit_brms, "logSGDstd", resp = "pHstd", method = "predict", resolution = 1000)
WR3<-W$pHstd.pHstd_logSGDstd %>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = percentsgd, y = pH), alpha = 0.1) +
  xlab("Percent SGD")+
  ylab(expression("pH"[t]))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  theme_minimal()

W<-conditional_effects(W_fit_brms, "logNNstd:DayNight", resp = "DICdiffstd", method = "predict", resolution = 1000)
WR4<-W$`DICdiffstd.DICdiffstd_logNNstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         logNN = logNNstd*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logNN), y = estimate, group = DayNight, color = DayNight), lwd = 2)+
  geom_ribbon(aes(x = exp(logNN),ymin=lower, ymax=upper, group = DayNight, fill = DayNight), linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = NN, y = DICdiff, color = DayNight), alpha = 0.1) +
  xlab(expression(paste("Nitrate + Nitrite (mmol L"^-1,")")))+
  ylab(expression(paste("Net Ecosystem Production ( ", Delta, "DIC ", mu,"mol kg"^-1, ")")))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()


W<-conditional_effects(W_fit_brms, "Tempinstd:DayNight", resp = "DICdiffstd", conditions = conditions, method = "predict", resolution = 1000)
WR5<-W$`DICdiffstd.DICdiffstd_Tempinstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate, group = DayNight, color = DayNight), lwd = 2)+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper, group = DayNight, fill = DayNight), linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = Tempin, y = DICdiff, color = DayNight), alpha = 0.1) +
  xlab(expression(paste("Temperature (", degree, "C)")))+
  ylab(expression(paste("Net Ecosystem Production ( ", Delta, "DIC ", mu,"mol kg"^-1, ")")))+
  theme_minimal()


W<-conditional_effects(W_fit_brms, "DICdiffstd", resp = "pHstd", method = "predict", resolution = 1000, conditions = conditions)
WR6<-Cdata %>%  ## conditional effects is cutting off some data do doing this the long way
  data_grid(DICdiffstd = seq_range(Cdata$DICdiffstd,2000), logSGDstd = median(Cdata$logSGDstd))%>%
  add_predicted_draws(W_fit_brms, resp = "pHstd", n = 2000, re_formula = NULL,
                      allow_new_levels = TRUE) %>%
  mutate(.prediction = .prediction*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         DICdiff = DICdiffstd*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center")
  )%>%
  median_qi(.width = c(.95)) %>%
  ggplot(aes(x = DICdiff, y = .prediction)) +
  geom_line(aes(y = .prediction), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = DICdiff,ymin=.prediction.lower, ymax=.prediction.upper), linetype=1.5, alpha=0.1, fill = "blue")+
  #stat_lineribbon(aes(y = .prediction), lwd = 2,.width = c(.95), alpha = 0.1, color = 'blue', fill = "blue") +
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = DICdiff, y = pH), alpha = 0.1) +
  xlab(expression(paste("Net Ecosystem Production ( ", Delta, "DIC ", mu,"mol kg"^-1, ")")))+
  ylab(expression("pH"[t]))+
  theme_minimal()

W<-conditional_effects(W_fit_brms, "pHstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
WR7<-Cdata %>% 
  data_grid(pHstd = seq_range(Cdata$pHstd,2000), Tempinstd = median(Cdata$Tempinstd))%>%
  add_predicted_draws(W_fit_brms, resp = "TAdiffstd", n = 2000, re_formula = NULL,
                      allow_new_levels = TRUE) %>%
  mutate(.prediction = .prediction*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         pH = pHstd*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center")
  )%>%
  median_qi(.width = c(.95)) %>%
  ggplot(aes(x = pH, y = .prediction)) +
  geom_line(aes(y = .prediction), lwd = 2, color = 'blue')+
  geom_ribbon(aes(x = pH,ymin=.prediction.lower, ymax=.prediction.upper), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = pH, y = TAdiff), alpha = 0.1) +
  ylab(expression(paste("Net Ecosystem Calcification ( ", Delta, "TA ", mu,"mol kg"^-1, ")")))+
  xlab(expression("pH"[t]))+
  theme_minimal()

W<-conditional_effects(W_fit_brms, "Tempinstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
WR8<-W$TAdiffstd.TAdiffstd_Tempinstd%>%
  mutate(estimate = estimate__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate), lwd = 2, color = "blue")+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper),fill = "blue", linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = Tempin, y = TAdiff), alpha = 0.1) +
  xlab(expression(paste("Temperature (", degree, "C)")))+
  ylab(expression(paste("Net Ecosystem Calcification", " (",Delta, "TA ", mu,"mol kg"^-1, ")")))+
  theme_minimal()

WR1+WR2+WR3+WR4+WR5+WR6+WR7+WR8+
  plot_annotation(title = 'Marginal Effects for all Wailupe Models', tag_levels = "A")+
  plot_layout(guides = "collect")+
  ggsave("Output/marginaleffects_Wailupe.png", width = 12, height = 10)

## get the posterior
Wpost <- posterior_samples(W_fit_brms)

# plot the coefficients
WCof1<-Wpost %>% 
  select(starts_with("b"),-ends_with("Intercept")) %>%
  gather() %>% 
  group_by(key)%>%
  median_hdci()%>%
  mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
  separate(col = key,into = c("b", "dependent", "independent"),sep = "_")%>% #loose the b and bring the values back together
  mutate(key = paste(dependent, independent))%>%
  ggplot(aes(x = value, y = reorder(key, value), color = sig)) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point()+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  scale_color_manual(values = c("grey","firebrick4"))+
  # stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
  #                     size = 3/4, color = "firebrick4") +
  labs(title = "Wailupe",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none")


# plot the intercepts
WCof2<-Wpost %>% 
  select(starts_with("Intercept"),ends_with("Intercept")) %>% 
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

## Plot of BP and Wailupe coefficients together
Cof1/WCof1+plot_annotation(tag_levels = "A")+
  ggsave("Output/coefficientsBoth.pdf", width = 10, height = 10)

# pull out the estimates and set it up to join with the DAG
# Westimates<-data.frame(fixef(W_fit_brms)) %>%
#   rownames_to_column(var = "name")%>%
#   separate(name, into = c("to","name", "DayNight"), sep = "[^[:alnum:]]+") %>%
#   select(name, to, 3:7)
# 
# ### deal with interaction terms
# Winteractions<-Wpost %>%
#   transmute(logNNstd_Night    = b_DICdiffstd_polylogNNstd21 + `b_DICdiffstd_DayNightNight:polylogNNstd21`,
#             logNNstd_Day = b_DICdiffstd_polylogNNstd21,
#             logNNstd2_Night    = b_DICdiffstd_polylogNNstd22 + `b_DICdiffstd_DayNightNight:polylogNNstd22`,
#             logNNstd2_Day = b_DICdiffstd_polylogNNstd22,
#             logPOstd_Night    = b_DICdiffstd_polylogPOstd21 + `b_DICdiffstd_DayNightNight:polylogPOstd21`,
#             logPOstd_Day = b_DICdiffstd_polylogPOstd21,
#             logPOstd2_Night    = b_DICdiffstd_polylogPOstd22 + `b_DICdiffstd_DayNightNight:polylogPOstd22`,
#             logPOstd2_Day = b_DICdiffstd_polylogPOstd22
#   ) %>%
#   gather(key, value) %>%
#   group_by(key) %>%
#   summarise(Estimate = mean(value), Est.Error = sd(value), 
#             Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
#             Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
#   separate(key, into = c("name", "DayNight")) %>%
#   mutate(to = "DICdiffstd") %>%
#   select(name,to, 2:5)
# 
# # remove all the DIC diff params from the model in names and replace with the calculated estimates for day and night
# Westimates<-Westimates %>%
#   filter(to != 'DICdiffstd') %>%
#   filter(name != 'Intercept') %>%
#   bind_rows(Winteractions) %>%
#   mutate(DayNight = replace_na(DayNight, "Both"))
