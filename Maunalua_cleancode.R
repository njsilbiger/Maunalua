## Run Bayesian SEM for Maunalua carbonate chemistry data
## By: Nyssa Silbiger
## Last updated: 10/19/2020
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
library(ggtext)
library(modelr)
library(ggfortify)

#load data
Cdata<-read.csv('chemicaldata_maunalua.csv', stringsAsFactors = TRUE)
#remove rows with NAs
Cdata<-Cdata[complete.cases(Cdata),]

#calculate rest of carbonate params
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, 
          T=Cdata$Temp_in, Patm=1, P=0, Pt=Cdata$Phosphate/1000000, Sit=Cdata$Silicate/1000000,
          k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propogation
er<-errors(flag=8, Cdata$pH, Cdata$TA/1000000, 
           S=Cdata$Salinity, T=Cdata$Temp_in, 
           Patm=1, P=0,Pt=Cdata$Phosphate/1000000,
           Sit=Cdata$Silicate/1000000,evar1 = 0.01, evar2 = 5e-6) 
      
#average error for DIC based on pH and TA
mean(er$DIC*1000000)
sd(er$DIC*1000000)/sqrt(nrow(er))

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

Cdata[,c("CO2","HCO3","CO3","DIC","OmegaArag","OmegaCalcite","pCO2","fCO2")]<-
  CO2[,c("CO2","HCO3","CO3","DIC","OmegaAragonite","OmegaCalcite","pCO2","fCO2")]

#endmembers from Christina
#black point
BP.end.DIC<-3038
BP.end.TA<-2946
BP.end.Sal<-4.9
BP.end.PO<-3.7
BP.end.NN<-163
BP.end.Si<-740

#wailupe
W.end.DIC<-1779
W.end.TA<-1754
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

#predicted TA based on mixing line
## use cristina's methods  C1 = Cmix + (Cmix – Csgd)(((Smix – 35.2)/(Ssgd – Smix))  

Cdata<-Cdata %>% # calculate predicted data from mixing line based on silicate for each site and season. I used Si because it is a better tracer for GW that Salinity
  #TA
  mutate(TA.pred = case_when(Site == 'BP'~ TA+(TA - BP.end.TA)*((Silicate - Hot.Si)/(BP.end.Si - Silicate)),
                             Site == 'W' ~ TA+(TA - W.end.TA)*((Silicate - Hot.Si)/(W.end.Si - Silicate))))%>%
  #DIC
  mutate(DIC.pred = case_when(Site == 'BP'~ DIC+(DIC - BP.end.DIC)*((Silicate - Hot.Si)/(BP.end.Si - Silicate)),
                              Site == 'W' ~ DIC+(DIC - W.end.DIC)*((Silicate - Hot.Si)/(W.end.Si - Silicate))))%>%
  #differences
  mutate(TA.diff = (Hot.TA-TA.pred)/2, #positive values are calcification and negative are dissolution
         DIC.diff = Hot.DIC - DIC) %>%#positive values are net photosynthesis and negative are respiration
  #Calculate percent SGD (Simix – SiSW)/ (SiGW – SiSW) with Christina's end members (L&O paper)
  mutate(percent_sgd = case_when(Site == 'BP'~0.1 +(100* (Silicate - 1.03)/(BP.end.Si - 1.03)),
                              Site == 'W' ~ 0.1 +(100* (Silicate - 1.03)/(W.end.Si - 1.03)))) #added 0.1 because the values will be log transformed
  

####################Anaysis########################### 
#Make tide just high and low instead of H1, H2, L1, and L2
Cdata$Tide<-droplevels(Cdata$Tide) #this removes levels that don'e exist anymore (empty spaces for example)
levels(Cdata$Tide)<-c("H","H","L","L") # this makes H1 and H2 both H and same for L1 and L2

# filter out the zones so that it is only diffuse, ambient, and transition. We did not get consistent data offshore because of Kayak disaster
Cdata<-Cdata %>% 
  dplyr::filter(Zone != 'Offshore')%>%
  droplevels()

#SEM######################
# Test the effect of SGD on ecosystem functioning mediated by changes in NN, Temperature, and pH. 
# Hypothesis: High SGD (low salinity/high silicate) increases N uptake of producers, which increases production (delta DIC), which increases pH, which increases NEC.
# But this relationship changes with Day/Night and Season (Temperature)

Cdata <- Cdata %>%
  filter(TA.diff < 150 & TA.diff > -10) %>% # remove 2 clear outliers
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
## SGD drives changes in N  and Temperature directly; 
# N and Temperature directly drive changes in DIC diff (net production);
# DIC diff and SGD directly drive changes in pH;
# pH directly drives changes in TA diff (NEC)

## N and P are highly colinear so I just used N to represent the nutrient gradient
# correlation coefficient for BP
BP.R <- round(cor(Cdata$NN[Cdata$Site =="BP"], Cdata$Phosphate[Cdata$Site =="BP"]),2)
## Cor coef for Wailupe
W.R <- round(cor(Cdata$NN[Cdata$Site =="W"], Cdata$Phosphate[Cdata$Site =="W"]),2)
#put it in a dataframe
cors<-data.frame(R = c(BP.R, W.R), Site = c("Black Point", "Wailupe"))
## get the slopes of the relationship between N and P to say something about N:P ratios at each site
BP.NP<-lm(NN~Phosphate, data = Cdata[Cdata$Site =="BP",])
W.NP<-lm(NN~Phosphate, data = Cdata[Cdata$Site =="W",])
summary(BP.NP)$coefficients[2]
summary(W.NP)$coefficients[2]

#make a plot
Cdata %>%
  mutate(Site = ifelse(Site == "BP", "Black Point", "Wailupe"))%>%
ggplot( aes(y = NN, x = Phosphate, group = Site, color = Site))+
  geom_point()+
  geom_smooth(method = "lm", formula = "y ~ x")+
  scale_color_manual(values = c("firebrick4", "orange"), name = " ")+
   #  geom_text(data = cors, aes(x = 7, y = 0.8, label = paste("Pearson's R =", R)))+
#  facet_wrap(~Site)+
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
  xlab(expression(atop("Phosphate", paste("(",mu, "mol L"^-1,")"))))+
  theme_bw()+
#  theme(strip.background = NULL)+
  ggsave("Output/NNvsPO.pdf", useDingbats = FALSE, width = 6, height = 4)


### model a bayesian SEM (individual one for each site since the geo chemistry of the SGD is different)
## brms does some weird things with columns names that have non-alphanumerics. Here I am removing them all
colnames(Cdata)<-str_replace_all(colnames(Cdata), "[^[:alnum:]]", "")

# These are all for a partially mediated model
# run one site at a time because they have different biogeochem in the SGD
NN_mod<-bf(logNNstd ~ logSGDstd) # NN ~ SGD which can change by site
Temp_mod<-bf(Tempinstd ~ DayNight*logSGDstd*Season) ## SGD has cooler water and relationship changes by DayNight and Seaso
#PO_mod<-bf(logPOstd ~ logSGDstd) # PO ~ SGD which can change by site
DIC_mod <- bf(DICdiffstd ~ DayNight*Tide*logNNstd +Season*Tempinstd) # DIC ~ nutrients and temperature, which can change by day/night (i.e. high nutrients could lead to high P during the day and high R at night what have opposite signs). Season also affects the relationship between NP and temp (light, flow, starting temp, etc)
pH_mod <- bf(pHstd ~ DICdiffstd + logSGDstd) # pH ~ NEP + SGD
TA_mod<-bf(TAdiffstd ~ pHstd+Tempinstd) # NEC ~ pH and temperature 

# compare to a fully mediated model which has direct pacths for SGD to NEP and NEC
TA_mod_full<-bf(TAdiffstd ~ pHstd+Tempinstd +logSGDstd) # NEC ~ pH and temperature 

# Run the model first for Black Point
# partial
k_fit_brms <- brm(TA_mod+
                    pH_mod+
                    DIC_mod+ 
                    Temp_mod+
                    NN_mod +
                    set_rescor(FALSE),
                  data=Cdata[Cdata$Site=='BP',],
                  cores=4, chains = 3)

# calculate LOO (leave one out) diagnostics
BP_loo<-loo(k_fit_brms, reloo = TRUE) # looks good!

# run the fully mediated model
k_fit_brms_full <- brm(TA_mod_full+
                         pH_mod+
                         DIC_mod+ 
                         Temp_mod+
                         NN_mod +
                         set_rescor(FALSE),
                       data=Cdata[Cdata$Site=='BP',],
                       cores=4, chains = 3)

# calculate LOO (leave one out) diagnostics
BP_looFull<-loo(k_fit_brms_full, reloo = TRUE) # looks good!

# difference  between partial and fill
loo(k_fit_brms, k_fit_brms_full, reloo = TRUE)
#there is almost no difference between the two models meaning that 
# no new information is added in the fully mediated model.  Will evaluate the partially mediated model

# view the effect sizes
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


# plot the conditional effects for all of the BP models

# Model 1
conditions <- make_conditions(k_fit_brms, "Season") # for the three way interaction

R<-conditional_effects(k_fit_brms, "logSGDstd", resp = "logNNstd", method = "predict", resolution = 1000)
R1<-R$logNNstd.logNNstd_logSGDstd %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         lower = lower__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         upper = upper__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = exp(logSGD), y = exp(estimate)), lwd = 1, color = 'grey')+
  geom_ribbon(aes(x = exp(logSGD),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = percentsgd, y = NN)) +
  xlab("% SGD")+
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
  ggtitle("Model 1")+
  coord_trans(x="log", y="log")+
  scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

#Model 2
R<-conditional_effects(k_fit_brms, "logSGDstd:DayNight", resp = "Tempinstd",  conditions = conditions,method = "predict", resolution = 1000)
R2<-R$`Tempinstd.Tempinstd_logSGDstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         lower = lower__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         upper = upper__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate, lty = DayNight, group = DayNight), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper, group = DayNight), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = percentsgd, y = Tempin, shape = DayNight) ) +
  xlab("% SGD")+
  ylab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  scale_shape_manual(values=c(0,15), name = "")+
  scale_linetype_manual(values = c(1,2),name = "")+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  ggtitle("Model 2")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')+
  facet_wrap(~Season)

#conditions <- make_conditions(k_fit_brms, "Tide") # for the three way interaction

#Model 3

R<-conditional_effects(k_fit_brms, "logSGDstd", resp = "pHstd",conditions = conditions, method = "predict", resolution = 1000)
R3<-R$pHstd.pHstd_logSGDstd %>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate), color = "grey", lwd = 1)+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper), fill = "grey", linetype=1.5, alpha=0.3)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = percentsgd, y = pH, color = DICdiff)) +
  xlab("% SGD")+
  ylab(expression("pH"[t]))+
  scale_y_continuous(limits = c(7.8, 8.4), breaks = seq(7.8, 8.4, by = 0.2))+
  labs(title = 'Model 3',
       color = "NEP")+
  coord_trans(x="log")+
  scale_color_gradient2(low = "#D8B365",
                         mid = "gray88",
                         high = "#5AB4AC",
                         midpoint = 0)+
  scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

R<-conditional_effects(k_fit_brms, "DICdiffstd", resp = "pHstd", method = "predict", resolution = 1000)
R6<-R$pHstd.pHstd_DICdiffstd%>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         DICdiff = DICdiffstd*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = DICdiff, y = estimate), lwd = 1, color = 'grey')+
  geom_ribbon(aes(x = DICdiff,ymin=lower, ymax=upper), linetype=1.5, fill = "grey", alpha = 0.3)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = DICdiff, y = pH, color = percentsgd)) +
  scale_color_gradient(name = "% SGD", trans = "log", breaks =c(0.2,1,5,10))+
  xlab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  ylab(expression("pH"[t]))+
  scale_y_continuous(limits = c(7.8, 8.4), breaks = seq(7.8, 8.4, by = 0.2))+
  labs(title = "Model 3")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

#Model 4

conditions <- make_conditions(k_fit_brms, c("Tide")) # for the three way interaction

R<-conditional_effects(k_fit_brms, "logNNstd:DayNight",conditions = conditions,  resp = "DICdiffstd", method = "predict", resolution = 1000)
R4<-R$`DICdiffstd.DICdiffstd_logNNstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         logNN = logNNstd*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logNN), y = estimate, group = DayNight, lty = DayNight), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = exp(logNN),ymin=lower, ymax=upper, group = DayNight), fill = "grey", linetype=1.5, alpha=0.3)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = NN, y = DICdiff, shape = DayNight, color = Tempin)) +
  xlab(expression(atop("Nitrate + Nitrite", paste("(mmol L"^-1,")"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  coord_trans(x="log")+
  scale_shape_manual(values=c(0,15), name = "")+
  scale_linetype_manual(values = c(1,2),name = "")+
  scale_x_continuous(breaks = c(0,0.1,1,5,30))+
  scale_y_continuous(limits = c(-200, 400), breaks = seq(-200, 500, by = 200))+
  scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
  labs(title = "Model 4")+
  theme_minimal()+
  facet_wrap(~Tide)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))


conditions <- make_conditions(k_fit_brms, "Season") # for the three way interaction

R<-conditional_effects(k_fit_brms, "Tempinstd:Season", resp = "DICdiffstd",  method = "predict", resolution = 1000)
R5<-R$`DICdiffstd.DICdiffstd_Tempinstd:Season`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>% ## only select the actual data to plot over so it does not over predict
  mutate(keep = case_when(Season == 'Spring' & Tempin <27 ~1,
                          Season == 'Fall' & Tempin >26 ~1))%>%
  filter(keep == 1)%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate, group = Season), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = Tempin, ymin=lower, ymax=upper, group = Season), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = Tempin, y = DICdiff, color = NN)) +
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  scale_color_gradient(name = "N+N", trans = "log", breaks =c(0.1,1,5,30), low = "lightgreen", high = "darkgreen")+
  labs(title = 'Model 4')+
  scale_y_continuous(limits = c(-200, 400), breaks = seq(-200, 500, by = 200))+
  theme_minimal()+
  facet_wrap(~Season)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

#Model 5
R<-conditional_effects(k_fit_brms, "pHstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
R7<-R$`TAdiffstd.TAdiffstd_pHstd`%>%
  mutate(estimate = estimate__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         pH = pHstd*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center")
  )%>%
  # mutate(keep = case_when(Season == 'Fall' & pH < 8.25 ~1,
  #                         Season == 'Spring' ~1))%>%
  # filter(keep == 1)%>%
  ggplot()+
  geom_line(aes(x = pH, y = estimate), lwd = 1, color = 'grey')+
  geom_ribbon(aes(x = pH,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = pH, y = TAdiff, color = Tempin)) +
  xlab(expression("pH"[t]))+
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
  scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
  labs(title = 'Model 5')+
  theme_minimal()+
#  facet_wrap(~Season)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))


R<-conditional_effects(k_fit_brms, "Tempinstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
R8<-R$`TAdiffstd.TAdiffstd_Tempinstd`%>%
  mutate(estimate = estimate__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
#  mutate(keep = case_when(Season == 'Spring' & Tempin <27 ~1,
#                          Season == 'Fall' & Tempin >26 ~1))%>%
#  filter(keep == 1)%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper),fill = "grey", linetype=1.5, alpha=0.3)+
  geom_point(data = Cdata[Cdata$Site=='BP',], aes(x = Tempin, y = TAdiff, color = pH)) +
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
  scale_color_gradient(low = "peachpuff", high = "sienna", name = expression("pH"[t]))+
  labs(title = "Model 5")+
  theme_minimal()+
#  facet_wrap(~Season)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

## bring them all together in patchwork
(R1|R2)/(R3|R6)/(R4|R5)/(R7|R8)+plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "A")+
  ggsave("Output/marginaleffects_BlackPoint.pdf", width = 10, height = 15, useDingbats = FALSE)

## get the posterior
post <- posterior_samples(k_fit_brms)

# plot the coefficients
BPCof<-post %>% 
  select(starts_with("b"),-ends_with("Intercept")) %>%
  gather() %>% 
  group_by(key)%>%
  median_hdci()%>%
  mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
  separate(col = key,into = c("b", "dependent", "independent"),sep = "_")%>% #loose the b and bring the values back together
 # mutate(key = paste(dependent, independent))%>%
  mutate(dependent = factor(dependent, levels = c("logNNstd","Tempinstd","pHstd","DICdiffstd","TAdiffstd")))%>%
  mutate(dependent  = recode(dependent, logNNstd = "NN", Tempinstd = "Temperature", pHstd = "pH", DICdiffstd = "NEP", TAdiffstd  ="NEC"),
         Site = "Black Point") 

Cof1<- BPCof %>%
 ggplot(aes(x = value, y = reorder(independent, value), alpha = sig, color = "firebrick4")) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point()+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  scale_alpha_manual(values = c(0.2,1))+
  scale_color_manual(values = "firebrick4")+
  # stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
  #                     size = 3/4, color = "firebrick4") +
  labs(title = "Black Point",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "none")+
  facet_grid(~dependent, scales = "free_y", space='free')

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

# pull out the estimates and set it up to join with the DAG
estimates1<-data.frame(fixef(k_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(name != 'Intercept', to != 'DICdiffstd', to !="Tempinstd") # remove the intercepts and interaction terms (Calculated below) for the DAG

## pull out main effects of season, day/night and tide
estimates<-data.frame(fixef(k_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(name %in% c("TideLowTide","SeasonSpring","DayNightNight"))%>% # remove the intercepts and interaction terms (Calculated below) for the DAG
  bind_rows(estimates1)

### deal with interaction terms. First for DIC models
InteractionsDIC<-post %>%
  transmute(logNNstd_Night_High_Both    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd`,
            logNNstd_Day_High_Both = `b_DICdiffstd_logNNstd`,
            
            logNNstd_Night_Low_Both    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` +`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_DayNightNight:TideLowTide:logNNstd`,
            logNNstd_Day_Low_Both = `b_DICdiffstd_logNNstd`+`b_DICdiffstd_TideLowTide:logNNstd`,
            
            Temperature_Both_Both_Spring    = `b_DICdiffstd_Tempinstd` + `b_DICdiffstd_SeasonSpring:Tempinstd`,
            Temperature_Both_Both_Fall = `b_DICdiffstd_Tempinstd`,
            
  ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Tide", "Season")) %>%
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

write.csv(estimates, "Output/BlackPointEstimates.csv")

  ###### Run Model for Wailupe ###############
# partially mediated model
W_fit_brms <- brm(TA_mod+
                    pH_mod+
                    DIC_mod+ 
                    Temp_mod+
                    NN_mod +
                    set_rescor(FALSE), 
                  data=Cdata[Cdata$Site=='W',],
                  cores=4, chains = 3)
# calculate LOO (leave one out) diagnostics
W_loo<-loo(W_fit_brms, reloo = TRUE) # looks good!


#fully mediated model
W_fit_brms_full <- brm(TA_mod_full+
                    pH_mod+
                    DIC_mod+ 
                    Temp_mod+
                    NN_mod +
                    set_rescor(FALSE), 
                  data=Cdata[Cdata$Site=='W',],
                  cores=4, chains = 3)

# calculate LOO (leave one out) diagnostics
W_looFull<-loo(W_fit_brms_full, reloo = TRUE) # looks good!

# difference  between partial and fill
loo(W_fit_brms, W_fit_brms_full, reloo = TRUE)
## difference is 5.3 +/- 3.5 suggesting some new information is added at Wailupe

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

#Model 1
conditions <- make_conditions(W_fit_brms, "Season") # for the three way interaction

W<-conditional_effects(W_fit_brms, "logSGDstd", resp = "logNNstd", method = "predict", resolution = 1000)
WR1<-W$logNNstd.logNNstd_logSGDstd %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         lower = lower__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         upper = upper__*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = exp(logSGD), y = exp(estimate)), lwd = 1, color = 'grey')+
  geom_ribbon(aes(x = exp(logSGD),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = percentsgd, y = NN)) +
  xlab("% SGD")+
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
  ggtitle("Model 1")+
  coord_trans(x="log", y="log")+
  scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

#Model 2
W<-conditional_effects(W_fit_brms, "logSGDstd:DayNight", resp = "Tempinstd",  conditions = conditions,method = "predict", resolution = 1000)
WR2<-W$`Tempinstd.Tempinstd_logSGDstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         lower = lower__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         upper = upper__*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate, lty = DayNight, group = DayNight), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper, group = DayNight), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = percentsgd, y = Tempin, shape = DayNight) ) +
  xlab("% SGD")+
  ylab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  scale_shape_manual(values=c(0,15), name = "")+
  scale_linetype_manual(values = c(1,2),name = "")+
  coord_trans(x="log")+
  scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  ggtitle("Model 2")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')+
  facet_wrap(~Season)

#Model 3
W<-conditional_effects(W_fit_brms, "logSGDstd", resp = "pHstd",conditions = conditions, method = "predict", resolution = 1000)
WR3<-W$pHstd.pHstd_logSGDstd %>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         logSGD = logSGDstd*attr(Cdata$logSGDstd,"scaled:scale")+attr(Cdata$logSGDstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = exp(logSGD), y = estimate), color = "grey", lwd = 1)+
  geom_ribbon(aes(x = exp(logSGD),ymin=lower, ymax=upper), fill = "grey", linetype=1.5, alpha=0.3)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = percentsgd, y = pH, color = DICdiff)) +
  xlab("% SGD")+
  ylab(expression("pH"[t]))+
  scale_y_continuous(limits = c(7.7, 8.4), breaks = seq(7.8, 8.4, by = 0.2))+
  labs(title = 'Model 3',
       color = "NEP")+
  coord_trans(x="log")+
  scale_color_gradient2(low = "#D8B365",
                        mid = "gray88",
                        high = "#5AB4AC",
                        midpoint = 0)+
  scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))


W<-conditional_effects(W_fit_brms, "DICdiffstd", resp = "pHstd", method = "predict", resolution = 1000)
WR6<-W$pHstd.pHstd_DICdiffstd%>%
  mutate(estimate = estimate__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         lower = lower__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         upper = upper__*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
         DICdiff = DICdiffstd*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center")
  )%>%
  ggplot()+
  geom_line(aes(x = DICdiff, y = estimate), lwd = 1, color = 'grey')+
  geom_ribbon(aes(x = DICdiff,ymin=lower, ymax=upper), linetype=1.5, fill = "grey", alpha = 0.3)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = DICdiff, y = pH, color = percentsgd)) +
  scale_color_gradient(name = "% SGD", trans = "log", breaks =c(0.2,1,5,10))+
  xlab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  ylab(expression("pH"[t]))+
  scale_y_continuous(limits = c(7.7, 8.4), breaks = seq(7.8, 8.4, by = 0.2))+
  labs(title = "Model 3")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

conditions <- make_conditions(W_fit_brms, c("Tide")) # for the three way interaction

#Model 4
W<-conditional_effects(W_fit_brms, "logNNstd:DayNight",conditions = conditions,  resp = "DICdiffstd", method = "predict", resolution = 1000)
WR4<-W$`DICdiffstd.DICdiffstd_logNNstd:DayNight`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         logNN = logNNstd*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center")
  )%>%
  mutate(keep = case_when(Tide == 'High Tide' & logNN <log(2) ~1,
                          Tide == 'Low Tide' ~1))%>%
  filter(keep == 1)%>%
  ggplot()+
  geom_line(aes(x = exp(logNN), y = estimate, group = DayNight, lty = DayNight), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = exp(logNN),ymin=lower, ymax=upper, group = DayNight), fill = "grey", linetype=1.5, alpha=0.3)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = NN, y = DICdiff, shape = DayNight, color = Tempin)) +
  xlab(expression(atop("Nitrate + Nitrite", paste("(mmol L"^-1,")"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  coord_trans(x="log")+
  scale_shape_manual(values=c(0,15), name = "")+
  scale_linetype_manual(values = c(1,2),name = "")+
  scale_x_continuous(breaks = c(0,0.1,1,5,30))+
  scale_y_continuous(limits = c(-200, 400), breaks = seq(-200, 500, by = 200))+
  scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
  labs(title = "Model 4")+
  theme_minimal()+
  facet_wrap(~Tide)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))


conditions <- make_conditions(W_fit_brms, "Season") # for the three way interaction

W<-conditional_effects(W_fit_brms, "Tempinstd:Season", resp = "DICdiffstd",  method = "predict", resolution = 1000)
WR5<-W$`DICdiffstd.DICdiffstd_Tempinstd:Season`%>%
  mutate(estimate = estimate__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>% ## only select the actual data to plot over so it does not over predict
  mutate(keep = case_when(Season == 'Spring' & Tempin <26.5 ~1,
                          Season == 'Fall' & Tempin >24 ~1))%>%
  filter(keep == 1)%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate, group = Season), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = Tempin, ymin=lower, ymax=upper, group = Season), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = Tempin, y = DICdiff, color = NN)) +
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Production", paste("(", Delta, "DIC ", mu,"mol kg"^-1, ")"))))+
  scale_color_gradient(name = "N+N", trans = "log", breaks =c(0.1,1,5,30), low = "lightgreen", high = "darkgreen")+
  labs(title = 'Model 4')+
  scale_y_continuous(limits = c(-200, 400), breaks = seq(-200, 500, by = 200))+
  theme_minimal()+
  facet_wrap(~Season)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

#Model 5
W<-conditional_effects(W_fit_brms, "pHstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
WR7<-W$`TAdiffstd.TAdiffstd_pHstd`%>%
  mutate(estimate = estimate__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         pH = pHstd*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center")
  )%>%
  # mutate(keep = case_when(Season == 'Fall' & pH < 8.25 ~1,
  #                         Season == 'Spring' ~1))%>%
  # filter(keep == 1)%>%
  ggplot()+
  geom_line(aes(x = pH, y = estimate), lwd = 1, color = 'grey')+
  geom_ribbon(aes(x = pH,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = pH, y = TAdiff, color = Tempin)) +
  xlab(expression("pH"[t]))+
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
  scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
  labs(title = 'Model 5')+
  theme_minimal()+
  #  facet_wrap(~Season)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

W<-conditional_effects(W_fit_brms, "Tempinstd", resp = "TAdiffstd", method = "predict", resolution = 1000)
WR8<-W$`TAdiffstd.TAdiffstd_Tempinstd`%>%
  mutate(estimate = estimate__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         lower = lower__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         upper = upper__*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center"),
         Tempin = Tempinstd*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center")
  )%>%
  #  mutate(keep = case_when(Season == 'Spring' & Tempin <27 ~1,
  #                          Season == 'Fall' & Tempin >26 ~1))%>%
  #  filter(keep == 1)%>%
  ggplot()+
  geom_line(aes(x = Tempin, y = estimate), lwd = 1, color = "grey")+
  geom_ribbon(aes(x = Tempin,ymin=lower, ymax=upper),fill = "grey", linetype=1.5, alpha=0.3)+
  geom_point(data = Cdata[Cdata$Site=='W',], aes(x = Tempin, y = TAdiff, color = pH)) +
  xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  ylab(expression(atop("Net Ecosystem Calcification",paste("(",Delta, "TA/2 ", mu,"mol kg"^-1, ")"))))+
  scale_color_gradient(low = "peachpuff", high = "sienna", name = expression("pH"[t]))+
  labs(title = "Model 5")+
  theme_minimal()+
  #  facet_wrap(~Season)+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
        legend.box = "horizontal")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

# Bring it together in patchwork
(WR1|WR2)/(WR3|WR6)/(WR4|WR5)/(WR7|WR8)+plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "A")+
  ggsave("Output/marginaleffects_Wailupe.pdf", width = 10, height = 15, useDingbats = FALSE)


## get the posterior
Wpost <- posterior_samples(W_fit_brms)

# Get and plot the coefficients for both Wailupe and Black Point
WCof<-Wpost %>% 
  select(starts_with("b"),-ends_with("Intercept")) %>%
  gather() %>% 
  group_by(key)%>%
  median_hdci()%>%
  mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
  separate(col = key,into = c("b", "dependent", "independent"),sep = "_")%>% #loose the b and bring the values back together
  # mutate(key = paste(dependent, independent))%>%
  mutate(dependent = factor(dependent, levels = c("logNNstd","Tempinstd","pHstd","DICdiffstd","TAdiffstd")))%>%
  mutate(dependent  = recode(dependent, logNNstd = "NN", Tempinstd = "Temperature", pHstd = "pH", DICdiffstd = "NEP", TAdiffstd  ="NEC"),
         Site = "Wailupe") %>%
  ##Bind with the Black Point coefficients
  bind_rows(BPCof)

#Make the plot
CoefPlot<-WCof%>%
  ggplot(aes(x = value, y = reorder(independent, value), alpha = sig, color = Site)) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point(size = 2)+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  scale_alpha_manual(values = c(0.2,1))+
  scale_color_manual(values = c("firebrick4", "orange"), name = " ")+
  # stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
  #                     size = 3/4, color = "firebrick4") +
  labs(#title = "Black Point",
       x = NULL, y = NULL) +
  theme_bw() +
  guides(alpha = FALSE)+
  theme(panel.grid   = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
            )+
  facet_grid(~dependent, scales = "free_y", space='free')+
  ggplot2::ggsave("Output/coefficientsBoth.pdf", width = 10, height = 5, useDingbats = FALSE)


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


# pull out the estimates and set it up to join with the DAG
estimates1<-data.frame(fixef(W_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(name != 'Intercept', to != 'DICdiffstd', to !="Tempinstd") # remove the intercepts and interaction terms (Calculated below) for the DAG

## pull out main effects of season, day/night and tide
estimates<-data.frame(fixef(W_fit_brms)) %>%
  rownames_to_column(var = "name")%>%
  separate(name, into = c("to","name"), sep = "_") %>%
  select(name, to, 3:6) %>%
  filter(name %in% c("TideLowTide","SeasonSpring","DayNightNight"))%>% # remove the intercepts and interaction terms (Calculated below) for the DAG
  bind_rows(estimates1)

### deal with interaction terms. First for DIC models
InteractionsDIC<-Wpost %>%
  transmute(logNNstd_Night_High_Both    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd`,
            logNNstd_Day_High_Both = `b_DICdiffstd_logNNstd`,
            
            logNNstd_Night_Low_Both    = `b_DICdiffstd_logNNstd` + `b_DICdiffstd_DayNightNight:logNNstd` +`b_DICdiffstd_TideLowTide:logNNstd` +`b_DICdiffstd_DayNightNight:TideLowTide:logNNstd`,
            logNNstd_Day_Low_Both = `b_DICdiffstd_logNNstd`+`b_DICdiffstd_TideLowTide:logNNstd`,
            
            Temperature_Both_Both_Spring    = `b_DICdiffstd_Tempinstd` + `b_DICdiffstd_SeasonSpring:Tempinstd`,
            Temperature_Both_Both_Fall = `b_DICdiffstd_Tempinstd`,
            
  ) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(Estimate = mean(value), Est.Error = sd(value),
            Q2.5 = mean(value) - qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n()),
            Q97.5 = mean(value) + qt(1- 0.05/2, (n() - 1))*sd(value)/sqrt(n())) %>%
  separate(key, into = c("name", "DayNight", "Tide", "Season")) %>%
  mutate(to = "DICdiffstd") %>%
  select(name,to, 2:8)


# Interactions temp SGD
InteractionsSGD<-Wpost %>%
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

write.csv(estimates, "Output/WailupeEstimates.csv")

### Make a DAG (see DAG.R) ####

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

## look at some values also by Day NIght and Zones as well

Cdata %>%
  select(Site, Zone, DayNight, Season, NN, pH, percentsgd,Salinity, Silicate,Tempin, DICdiff, TAdiff) %>%
  group_by(Site, Zone, DayNight, Season)%>%
  summarise_all(.funs = list(~mean(.))) %>%
  pivot_longer(names_to = "Parameters", cols = "NN":"TAdiff") %>%
  separate(Parameters, into = c("Parameters", "stat")) %>%
  arrange(desc(Parameters)) %>%
  pivot_wider(names_from = c(Site:Season, stat)) %>%
  View()

Cdata %>%
  select(Site, DayNight, Season, Tide, NN) %>%
  group_by(Site, DayNight, Tide, Season)%>%
  summarise_all(.funs = list(~max(.)))

## Make a PCA of all the environmental data. Facet it by side and show differences in day/night, tide, and season
df_BP <- Cdata %>%
  filter(Site =="BP")%>%
  mutate(season_day = paste(Season, DayNight)) # makes it easier to code


pca_res <- prcomp(df_BP%>% # makes it easier to code
                    select(logSGDstd, logNNstd, pHstd, Tempinstd, DICdiffstd, TAdiffstd))

# make plot for Black Point
PCA_BP<-autoplot(pca_res, data = df_BP, colour = 'Tide', shape = "season_day",
         loadings.label.label = c("% SGD", "N+N","pH","Temperature","NEP","NEC"), 
         loadings.label.repel=T,
         loadings = TRUE, loadings.colour = 'lightblue',loadings.label.colour = 'dodgerblue',
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.vjust = 1.2)+
  scale_shape_manual(values = c(1,16,2,17))+
  scale_color_manual(values = c("black","grey"))+
  theme_few()+
  theme(legend.title=element_blank(),
        legend.position = "none")

#Wailupe
df_W <- Cdata %>%
  filter(Site =="W")%>%
  mutate(season_day = paste(Season, DayNight)) # makes it easier to code


pca_res_W <- prcomp(df_W%>% # makes it easier to code
                    select(logSGDstd, logNNstd, pHstd, Tempinstd, DICdiffstd, TAdiffstd))

# make plot for Black Point
PCA_W<-autoplot(pca_res_W, data = df_W, colour = 'Tide', shape = "season_day",
                 loadings.label.label = c("% SGD", "N+N","pH","Temperature","NEP","NEC"), 
                 loadings.label.repel=T,
                 loadings = TRUE, loadings.colour = 'lightblue',loadings.label.colour = 'dodgerblue',
                 loadings.label = TRUE, loadings.label.size = 4, loadings.label.vjust = 1.2)+
  scale_shape_manual(values = c(1,16,2,17))+
  scale_color_manual(values = c("black","grey"))+
  theme_few()+
  theme(legend.title=element_blank())

# bring together in patchwork
PCA_BP+PCA_W + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A")+
  ggsave("Output/PCAplots.pdf", height = 5, width = 10)
