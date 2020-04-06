## Run Bayesian SEM for Maunalua carbonate chemistry data
## By: Nyssa Silbiger
## Last updated: 4/4/2020
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

# load bulk dissolution data
flow<-read.csv('Maunalua_data_flow.csv')

#join with Cdata
Cdata<-left_join(Cdata,flow)

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
                       H = "High Tide"),
         Season = recode(Season, SPRING = "Spring",
                         FALL = "Fall"))  ## change the factors for prettier names
# each level is a model
## SGD drives changes in N  directly; 
# N  directly drive changes in DIC diff (net production);
# DIC diff and SGD directly drive changes in pH;
# pH directly drives changes in TA diff (NEC)

### model a bayesian SEM (individual one for each site since the geo chemistry of the SGD is different)
## brms does some weird things with columns names that have non-alphanumerics. Here I am removing them all
colnames(Cdata)<-str_replace_all(colnames(Cdata), "[^[:alnum:]]", "")

# run one site at a time because they have different biogeochem in the SGD
NN_mod<-bf(logNNstd ~ logSGDstd) # NN ~ SGD which can change by site
Temp_mod<-bf(Tempinstd ~ DayNight*logSGDstd*Season) ## SGD has cooler water and the intercept changes with season
#PO_mod<-bf(logPOstd ~ logSGDstd) # PO ~ SGD which can change by site
DIC_mod <- bf(DICdiffstd ~ (Tide*DayNight*logNNstd*Season +Tempinstd)) # DIC ~ nutrients and temperature, which can change by day/night (i.e. high nutrients could lead to high P during the day and high R at night what have opposite signs). It is also non-linear
pH_mod <- bf(pHstd ~ DICdiffstd + logSGDstd) # pH ~ NEP + SGD
TA_mod<-bf(TAdiffstd ~ pHstd+Tempinstd) # NEC ~ pH and temperature, which can change by Tide, because low has more nutrients it may distrupt this relationship (based on Silbiger et al. 2018)

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
p1+p2+p3+p4+p5+plot_layout(guides = "collect") +
  plot_annotation(title = 'Black Point Posterior Predictive Checks', tag_levels = "A")+
  ggsave("Output/Posteriorchecks_BlackPoint.pdf", width = 5, height = 5)
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
  ylab(expression(atop("Nitrate + Nitrite", paste("(mmol L"^-1,")"))))+
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
  ylab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  theme_minimal()+
  facet_wrap(~Season)

#conditions <- make_conditions(k_fit_brms, "Tide") # for the three way interaction

R<-conditional_effects(k_fit_brms, "logSGDstd", resp = "pHstd",conditions = conditions, method = "predict", resolution = 1000)
R3<-R$pHstd.pHstd_logSGDstd %>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate), color = "blue", lwd = 2)+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper), fill = "blue", linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = percentsgd, y = pH), color = "blue", alpha = 0.1) +
  xlab("Percent SGD")+
  ylab(expression("pH"[t]))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  theme_minimal()

conditions <- make_conditions(k_fit_brms, c("Tide","Season")) # for the three way interaction

R<-conditional_effects(k_fit_brms, "logNNstd:DayNight",conditions = conditions,  resp = "DICdiffstd", method = "predict", resolution = 1000)
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
  xlab(expression(atop("Nitrate + Nitrite", paste("(mmol L"^-1,")"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,0.1,1,5,30))+
  theme_minimal()+
  facet_wrap(~Tide*Season)

conditions <- make_conditions(k_fit_brms, "Season") # for the three way interaction

R<-conditional_effects(k_fit_brms, "Tempinstd:Season", resp = "DICdiffstd",  method = "predict", resolution = 1000)
R5<-R$`DICdiffstd.DICdiffstd_Tempinstd:Season`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate, group = Season), lwd = 2, color = "blue")+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper, group = Season), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = Tempin, y = DICdiff), alpha = 0.1, color = "blue") +
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  facet_wrap(~Season)+
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
  xlab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
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
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
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
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
  theme_minimal()

R1+R2+R3+R4+R5+R6+R7+R8+
  plot_annotation(title = 'Marginal Effects for all Black Point Models', tag_levels = "A")+
  plot_layout(guides = "collect")+
  ggsave("Output/marginaleffects_BlackPoint.pdf", width = 12, height = 10)

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
estimates1<-data.frame(fixef(k_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(name != 'Intercept', to != 'DICdiffstd', to !="Tempinstd") # remove the intercepts and interaction terms (Calculated below) for the DAG

## pull out temp to DIC
estimates<-data.frame(fixef(k_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(to == 'DICdiffstd', name =="Tempinstd")%>% # remove the intercepts and interaction terms (Calculated below) for the DAG
  bind_rows(estimates1)
  
### deal with interaction terms. First for DIC models
InteractionsDIC<-post %>%
  transmute(logNNstd_Night_Fall_High    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd`,
            logNNstd_Day_Fall_High = `b_DICdiffstd_logNNstd`,
            
            logNNstd_Night_Fall_Low    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` +`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_TideLowTide:DayNightNight:logNNstd`,
            logNNstd_Day_Fall_Low = `b_DICdiffstd_logNNstd`+`b_DICdiffstd_TideLowTide:logNNstd`,
            
            logNNstd_Night_Spring_High    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` + `b_DICdiffstd_logNNstd:SeasonSpring` +`b_DICdiffstd_DayNightNight:logNNstd:SeasonSpring`,
            logNNstd_Day_Spring_High = `b_DICdiffstd_logNNstd`+ `b_DICdiffstd_logNNstd:SeasonSpring`,
            
            logNNstd_Night_Spring_Low    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` + `b_DICdiffstd_logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_TideLowTide:DayNightNight:logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:logNNstd:SeasonSpring`+`b_DICdiffstd_DayNightNight:logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:DayNightNight:logNNstd`,
            logNNstd_Day_Spring_Low = `b_DICdiffstd_logNNstd`+ `b_DICdiffstd_logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_TideLowTide:logNNstd:SeasonSpring`
            
               ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Season", "Tide")) %>%
  mutate(to = "DICdiffstd") %>%
  select(name,to, 2:8)

# Interactions temp SGD
InteractionsSGD<-post %>%
  transmute(
    logSGDstd_Day_Spring_Both = `b_Tempinstd_logSGDstd` +`b_Tempinstd_logSGDstd:SeasonSpring`,
    logSGDstd_Night_Spring_Both = `b_Tempinstd_logSGDstd` +`b_Tempinstd_logSGDstd:SeasonSpring` +`b_Tempinstd_DayNightNight:logSGDstd`+`b_Tempinstd_DayNightNight:logSGDstd:SeasonSpring`,
    logSGDstd_Day_Fall_Both = `b_Tempinstd_logSGDstd`,
    logSGDstd_Night_Fall_Both = `b_Tempinstd_logSGDstd`+`b_Tempinstd_DayNightNight:logSGDstd`
    
  ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Season", "Tide")) %>%
  mutate(to = "Tempinstd") %>%
  select(name,to, 2:8)

#Bind everything together
  estimates<-bind_rows(estimates, InteractionsDIC, InteractionsSGD) %>% # bind with the estimates above
  mutate(DayNight = replace_na(DayNight, "Both"),
         Season = replace_na(Season, "Both"),
         Tide = replace_na(Tide, "Both")  )

### Make a DAG ####
  bigger_dag <- dagify(TAdiffstd ~ pHstd+ Tempinstd,
                       pHstd ~ DICdiffstd + logSGDstd,
                       DICdiffstd ~ logNNstd + Tempinstd,
                       logNNstd ~ logSGDstd,
                       Tempinstd ~ logSGDstd,
                       exposure = "logSGDstd",
                       outcome = "TAdiffstd",
                       labels = c("TAdiffstd" = "NEC",
                                  "pHstd" = "pH",
                                  "Tempinstd" = "Temperature",
                                  "logSGDstd" = "Log % SGD",
                                  "DICdiffstd" ="NEP",
                                  "logNNstd" = "log NN")) %>%
    tidy_dagitty(layout = "tree", seed = 25)
  
  
  #quick visual
#  ggdag(bigger_dag, use_labels = "label", text = FALSE)+
#    theme_dag()

# join it with the estimates so that I can add colors and line thickness to related to effectsize
  DAGdata<-bigger_dag %>%
    dag_paths() %>% #### then left join these with the effect sizes
    as_tibble() %>%
    left_join(estimates) %>%
    mutate(DayNight = replace_na(DayNight, "Both"),
           Season = replace_na(Season, "Both"),
           Tide = replace_na(Tide, "Both"))%>%
    mutate(est_sign = as.character(sign(Estimate)))%>% # add a column for pos and negative
    mutate(edge_cols = case_when(est_sign == '-1' ~ "#d6604d",#"#fc8d59", 
                                 est_sign == '1' ~ "#4393c3")) %>%
    mutate(edge_lines = ifelse(sign(Q2.5)==sign(Q97.5),1,2),# add a dashed or solid line for significant effects (i.e. 95%CI does not overlap 0)
           edge_alpha = ifelse(sign(Q2.5)==sign(Q97.5),1,0.1)# make non-significant transparent
           ) 

  
  ## Basic DAG
  DAGdata %>%
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Spring"), Tide %in% c("Both","Low"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_node() +
    geom_dag_edges(
      #edge_colour = 'grey'
      #edge_linetype = edge_lines,
      #edge_alpha = edge_alpha)
    ) +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggsave(filename = 'Output/BasicDAG.pdf', width = 6, height = 6)

  #Day Spring High DAG
  DaySpringHigh_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Spring"), Tide %in% c("Both","High"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Spring Daytime High Tide')+
    theme(plot.title = element_text(hjust = 0.5))
                                                                    
  #Day Spring Low DAG
  DaySpringLow_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Spring"), Tide %in% c("Both","Low"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Spring Daytime Low Tide')+
    theme(plot.title = element_text(hjust = 0.5))
  
  #Night Spring High DAG
  NightSpringHigh_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Spring"), Tide %in% c("Both","High"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Spring Nighttime High Tide')+
    theme(plot.title = element_text(hjust = 0.5))
  
  #Night Spring Low DAG
  NightSpringLow_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Spring"), Tide %in% c("Both","Low"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Spring Nighttime Low Tide')+
    theme(plot.title = element_text(hjust = 0.5))
  
  #Day Fall High DAG
  DayFallHigh_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Fall"),Tide %in% c("Both","High"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Fall Daytime High Tide')+
    theme(plot.title = element_text(hjust = 0.5))

  #Day Fall Low DAG
  DayFallLow_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Fall"),Tide %in% c("Both","Low"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Fall Daytime Low Tide')+
    theme(plot.title = element_text(hjust = 0.5))
  
  #Night Fall High DAG
  NightFallHigh_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Fall"), Tide %in% c("Both","High"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
     geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Fall Nighttime High Tide')+
    theme(plot.title = element_text(hjust = 0.5)) 
  
  #Night Fall High DAG
  NightFallLow_DAG<-DAGdata %>% 
    filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Fall"), Tide %in% c("Both","Low"))%>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_width = abs(Estimate), 
                       edge_colour = edge_cols, 
                       #edge_linetype = edge_lines,
                       edge_alpha = edge_alpha)) +
    geom_dag_node() +
    geom_dag_text_repel(aes(label = label))+
    theme_dag() +
    ggtitle('Fall Nighttime Low Tide')+
    theme(plot.title = element_text(hjust = 0.5)) 
  
  DAGPlot_BP_High<- (DaySpringHigh_DAG +NightSpringHigh_DAG)/(DayFallHigh_DAG +NightFallHigh_DAG)+
    plot_annotation(tag_levels = "A",title = "Black Point")+
    ggsave("Output/DAGplotsBP_High.pdf", width = 12, height = 13, useDingbats = FALSE)
  
  DAGPlot_BP_Low<- (DaySpringLow_DAG +NightSpringLow_DAG)/(DayFallLow_DAG +NightFallLow_DAG)+
    plot_annotation(tag_levels = "A",title = "Black Point")+
    ggsave("Output/DAGplotsBP_Low.pdf", width = 12, height = 13, useDingbats = FALSE)
  
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
Wp1+Wp2+Wp3+Wp4+Wp5+plot_layout(guides = "collect") +
  plot_annotation(title = 'Wailupe Posterior Predictive Checks', tag_levels = "A")+
  ggsave("Output/Posteriorchecks_Wailupe.pdf", width = 5, height = 5)

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
  ylab(expression(atop("Nitrate + Nitrite", paste("(mmol L"^-1,")"))))+
  coord_trans(x="log", y="log")+
  scale_x_continuous(breaks = c(0,1,5,10,25))+
  scale_y_continuous(breaks = c(0,0.1,1,5,10))+
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
  ylab(expression(atop("Temperature",paste("(", degree, "C)"))))+
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

conditions <- make_conditions(W_fit_brms, c("Tide","Season")) # for the three way interaction

WR<-conditional_effects(W_fit_brms, "logNNstd:DayNight",conditions = conditions,  resp = "DICdiffstd", method = "predict", resolution = 1000)
WR4<-WR$`DICdiffstd.DICdiffstd_logNNstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         logNN = logNNstd*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logNN), y = estimate, group = DayNight, color = DayNight), lwd = 2)+
  geom_ribbon(aes(x = exp(logNN),ymin=lower, ymax=upper, group = DayNight, fill = DayNight), linetype=1.5, alpha=0.1)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = NN, y = DICdiff, color = DayNight), alpha = 0.1) +
  xlab(expression(atop("Nitrate + Nitrite", paste("(mmol L"^-1,")"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0,0.1,1,5,10))+
  theme_minimal()+
  facet_wrap(~Tide*Season)

conditions <- make_conditions(W_fit_brms, "Season") # for the three way interaction
WR<-conditional_effects(W_fit_brms, "Tempinstd:Season", resp = "DICdiffstd",  method = "predict", resolution = 1000)
WR5<-WR$`DICdiffstd.DICdiffstd_Tempinstd:Season`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate, group = Season), lwd = 2, color = "blue")+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper, group = Season), linetype=1.5, alpha=0.1, fill = "blue")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = Tempin, y = DICdiff), alpha = 0.1, color = "blue") +
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  facet_wrap(~Season)+
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
  xlab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
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
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
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
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
  theme_minimal()

WR1+WR2+WR3+WR4+WR5+WR6+WR7+WR8+
  plot_annotation(title = 'Marginal Effects for all Wailupe Models', tag_levels = "A")+
  plot_layout(guides = "collect")+
  ggsave("Output/marginaleffects_Wailupe.pdf", width = 12, height = 10)

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
  ggplot2::ggsave("Output/coefficientsBoth.pdf", width = 10, height = 10, useDingbats = FALSE)

Cof1+Cof2+Cof3+
  plot_annotation(title = 'Standardized coefficients for Black Point', tag_levels = "A")+
  ggsave("Output/coefficients_BlackPoint.png", width = 10, height = 7)

# pull out the estimates and set it up to join with the DAG
Westimates1<-data.frame(fixef(W_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(name != 'Intercept', to != 'DICdiffstd', to !="Tempinstd") # remove the intercepts and interaction terms (Calculated below) for the DAG

## pull out temp to DIC
Westimates<-data.frame(fixef(W_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(to == 'DICdiffstd', name =="Tempinstd")%>% # remove the intercepts and interaction terms (Calculated below) for the DAG
  bind_rows(Westimates1)


### deal with interaction terms. First for DIC models
WInteractionsDIC<-Wpost %>%
  transmute(logNNstd_Night_Fall_High    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd`,
            logNNstd_Day_Fall_High = `b_DICdiffstd_logNNstd`,
            
            logNNstd_Night_Fall_Low    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` +`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_TideLowTide:DayNightNight:logNNstd`,
            logNNstd_Day_Fall_Low = `b_DICdiffstd_logNNstd`+`b_DICdiffstd_TideLowTide:logNNstd`,
            
            logNNstd_Night_Spring_High    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` + `b_DICdiffstd_logNNstd:SeasonSpring` +`b_DICdiffstd_DayNightNight:logNNstd:SeasonSpring`,
            logNNstd_Day_Spring_High = `b_DICdiffstd_logNNstd`+ `b_DICdiffstd_logNNstd:SeasonSpring`,
            
            logNNstd_Night_Spring_Low    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` + `b_DICdiffstd_logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_TideLowTide:DayNightNight:logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:logNNstd:SeasonSpring`+`b_DICdiffstd_DayNightNight:logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:DayNightNight:logNNstd`,
            logNNstd_Day_Spring_Low = `b_DICdiffstd_logNNstd`+ `b_DICdiffstd_logNNstd:SeasonSpring`+`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_TideLowTide:logNNstd:SeasonSpring`
            
  ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Season", "Tide")) %>%
  mutate(to = "DICdiffstd") %>%
  select(name,to, 2:8)


# Interactions temp SGD
WInteractionsSGD<-Wpost %>%
  transmute(
    logSGDstd_Day_Spring_Both = `b_Tempinstd_logSGDstd` +`b_Tempinstd_logSGDstd:SeasonSpring`,
    logSGDstd_Night_Spring_Both = `b_Tempinstd_logSGDstd` +`b_Tempinstd_logSGDstd:SeasonSpring` +`b_Tempinstd_DayNightNight:logSGDstd`+`b_Tempinstd_DayNightNight:logSGDstd:SeasonSpring`,
    logSGDstd_Day_Fall_Both = `b_Tempinstd_logSGDstd`,
    logSGDstd_Night_Fall_Both = `b_Tempinstd_logSGDstd`+`b_Tempinstd_DayNightNight:logSGDstd`
    
  ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Season", "Tide")) %>%
  mutate(to = "Tempinstd") %>%
  select(name,to, 2:8)

#Bind everything together
Westimates<-bind_rows(Westimates, WInteractionsDIC, WInteractionsSGD) %>% # bind with the estimates above
  mutate(DayNight = replace_na(DayNight, "Both"),
         Season = replace_na(Season, "Both"),
         Tide = replace_na(Tide, "Both") )

# join it with the estimates so that I can add colors and line thickness to related to effectsize
WDAGdata<-bigger_dag %>%
  dag_paths() %>% #### then left join these with the effect sizes
  as_tibble() %>%
  left_join(Westimates) %>%
  mutate(DayNight = replace_na(DayNight, "Both"),
         Season = replace_na(Season, "Both"),
         Tide = replace_na(Tide, "Both")
  )%>%
  mutate(est_sign = as.character(sign(Estimate)))%>% # add a column for pos and negative
  mutate(edge_cols = case_when(est_sign == '-1' ~ "#d6604d",#"#fc8d59", 
                               est_sign == '1' ~ "#4393c3")) %>%
  mutate(edge_lines = ifelse(sign(Q2.5)==sign(Q97.5),1,2),# add a dashed or solid line for significant effects (i.e. 95%CI does not overlap 0)
         edge_alpha = ifelse(sign(Q2.5)==sign(Q97.5),1,0.1)# make non-significant transparent
  ) 

#Day Spring High DAG
WDaySpringHigh_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Spring"), Tide %in% c("Both", "High"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Spring Daytime High Tide')+
  theme(plot.title = element_text(hjust = 0.5))

#Day Spring Low DAG
WDaySpringLow_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Spring"), Tide %in% c("Both", "Low"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Spring Daytime Low Tide')+
  theme(plot.title = element_text(hjust = 0.5))

#Night Spring High DAG
WNightSpringHigh_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Spring"), Tide %in% c("Both", "High"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Spring Nighttime High Tide')+
  theme(plot.title = element_text(hjust = 0.5))

#Night Spring Low DAG
WNightSpringLow_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Spring"), Tide %in% c("Both", "Low"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Spring Nighttime Low Tide')+
  theme(plot.title = element_text(hjust = 0.5))

#Day Fall High DAG
WDayFallHigh_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Fall"), Tide %in% c("Both", "High"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Fall Daytime High Tide')+
  theme(plot.title = element_text(hjust = 0.5))

#Day Fall Low DAG
WDayFallLow_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Day"), Season %in% c("Both", "Fall"), Tide %in% c("Both", "Low"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Fall Daytime Low Tide')+
  theme(plot.title = element_text(hjust = 0.5))

#Night Fall High DAG
WNightFallHigh_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Fall"), Tide %in% c("Both", "High"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Fall Nighttime High Tide')+
  theme(plot.title = element_text(hjust = 0.5)) 

#Night Fall Low DAG
WNightFallLow_DAG<-WDAGdata %>% 
  filter(DayNight %in% c("Both", "Night"), Season %in% c("Both", "Fall"), Tide %in% c("Both", "Low"))%>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(aes(edge_width = abs(Estimate), 
                     edge_colour = edge_cols, 
                     #edge_linetype = edge_lines,
                     edge_alpha = edge_alpha)) +
  geom_dag_node() +
  geom_dag_text_repel(aes(label = label))+
  theme_dag() +
  ggtitle('Fall Nighttime Low Tide')+
  theme(plot.title = element_text(hjust = 0.5)) 

DAGPlot_WHigh<- (WDaySpringHigh_DAG +WNightSpringHigh_DAG)/(WDayFallHigh_DAG +WNightFallHigh_DAG)+
  plot_annotation(tag_levels = "A",title = "Wailupe")+
  ggsave("Output/DAGplotsWHigh.pdf", width = 12, height = 13, useDingbats = FALSE)

DAGPlot_WLow<- (WDaySpringLow_DAG +WNightSpringLow_DAG)/(WDayFallLow_DAG +WNightFallLow_DAG)+
  plot_annotation(tag_levels = "A",title = "Wailupe")+
  ggsave("Output/DAGplotsWLow.pdf", width = 12, height = 13, useDingbats = FALSE)

######################## Make some summary tables #############
Cdata %>%
  select(Site, Season, NN, pH, percentsgd,Salinity, Silicate,Tempin, DICdiff, TAdiff) %>%
  group_by(Site, Season)%>%
  summarise_all(.funs = list(~min(.), ~max(.), ~mean(.))) %>%
  pivot_longer(names_to = "Parameters", cols = "NN_min":"TAdiff_mean") %>%
  separate(Parameters, into = c("Parameters", "stat")) %>%
  arrange(desc(Parameters)) %>%
  pivot_wider(names_from = c(Site:Season, stat)) %>%
  write.csv(file = "SummaryTables/ranges.csv", row.names = FALSE)

# export the effect sizes
BPparams<-fixef(k_fit_brms) %>%
  as.data.frame() %>%
  rownames_to_column()%>%
  as_tibble() %>%
  rename_at(.vars = 2:5, paste0, "_BP") 

Wparams<-fixef(W_fit_brms) %>%
  as.data.frame() %>%
  rownames_to_column()%>%
  as_tibble() %>%
  rename_at(.vars = 2:5, paste0, "_W") 

params <-left_join(BPparams, Wparams)%>%
  write.csv(file = "SummaryTables/param_estimates.csv", row.names = FALSE)
