## used to make predictions for 50% change (increase and decrease) in SGD

# load additional libraries
library(matrixStats)
library(facetscales)

# Need to source the Maunalua Clean code before running this.

## Get predictions for a 50% (increase and decrease) change in SGD
set.seed(131)

# Add a 50% increase in SGD to every value, log it and then standardize it
Cdata_BP<-Cdata%>%
  select(Waypoint, Site,Tide, Season, DayNight, percentsgd, logNN, Tempin, pH, DICdiff, TAdiff,logSGD, logSGDstd, logNNstd, Tempinstd, pHstd, DICdiffstd, TAdiffstd)%>%
  filter(Site == "BP")%>%
  mutate(Site = "Kupikipiki'o",
         Obs_pred = "Observed",
     #    SGD10 = percentsgd+percentsgd, # add a 50% increase in SGD, log it, and standardize it
    #     logSGD10 = log(SGD10),
         # put it on the same scale as the original dataset
         #(logSGD10 - attr(Cdata$logSGDstd,"scaled:center")/attr(Cdata$logSGDstd,"scaled:scale"))
         logSGD50std = log(exp(logSGDstd)+0.5*exp(logSGDstd)),  # add a 50% increase in SGD, log it, and standardize it
         logSGDneg50std = log(exp(logSGDstd)-0.5*exp(logSGDstd)) # 50% decrease in SGD
    )
         #logSGD10std = scale(logSGD10, scale = TRUE, center = TRUE))

# wailupe
Cdata_W<-Cdata%>%
  select(Waypoint, Site,Tide, Season, DayNight, percentsgd, logNN, Tempin, pH, DICdiff, TAdiff,logSGD, logSGDstd, logNNstd, Tempinstd, pHstd, DICdiffstd, TAdiffstd)%>%
  filter(Site == "W")%>%
  mutate(Site = "Wailupe",
         Obs_pred = "Observed",
   #      SGD10 = percentsgd+percentsgd, # add a 50% increase in SGD, log it, and standardize it
    #     logSGD10 = log(SGD10),
         # put it on the same scale as the original dataset
         logSGD50std = log(exp(logSGDstd)+0.5*exp(logSGDstd)),  # add a 50% increase in SGD, log it, and standardize it
         logSGDneg50std = log(exp(logSGDstd)-0.5*exp(logSGDstd)) # 50% decrease in SGD
         ) 


## THis only works for predictions that are directly related to the changed variable. It does not calculate downstream changes
Predictions_nochange<-predicted_draws(k_fit_brms, newdata=Cdata_BP) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median) %>% # get the median of each draw
  mutate(Site = "Kupikipiki'o")

# get predictions for NN and Temp for increase in 50%
Predictions_logNN_Temp<-predicted_draws(k_fit_brms, newdata=Cdata_BP %>% select(!logSGDstd) %>% rename(logSGDstd = logSGD50std))%>%
  filter(.category %in% c("Tempinstd","logNNstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

#predictions for -50%
Predictions_logNN_Temp_neg50<-predicted_draws(k_fit_brms, newdata=Cdata_BP %>% select(!logSGDstd) %>% rename(logSGDstd = logSGDneg50std))%>%
  filter(.category %in% c("Tempinstd","logNNstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for NEP +50%
newdata_DIC<-Predictions_logNN_Temp %>%
  select(-c(logNNstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  select(Season, Tide, DayNight, logNNstd, Tempinstd, TAdiffstd, pHstd, DICdiffstd, logSGDstd)

Predictions_DIC<-predicted_draws(k_fit_brms, newdata=newdata_DIC)%>%
  filter(.category %in% c("DICdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)


## get predictions for NEP -50%
newdata_DICneg50<-Predictions_logNN_Temp_neg50 %>%
  select(-c(logNNstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  select(Season, Tide, DayNight, logNNstd, Tempinstd, TAdiffstd, pHstd, DICdiffstd, logSGDstd)

Predictions_DICneg50<-predicted_draws(k_fit_brms, newdata=newdata_DICneg50)%>%
  filter(.category %in% c("DICdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for pH +50%
newdata_pH<-newdata_DIC %>%
  select(-DICdiffstd) %>%
  bind_cols(Predictions_DIC %>%
               select(.prediction, .category) %>%
               pivot_wider(names_from = .category, values_from = .prediction) %>%
               select(DICdiffstd))%>%
  ungroup()%>%
  select(-c(.row, DayNight1, Season1))

Predictions_pH<-predicted_draws(k_fit_brms, newdata=newdata_pH)%>%
  filter(.category %in% c("pHstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for pH -50%
newdata_pHneg50<-newdata_DICneg50 %>%
  select(-DICdiffstd) %>%
  bind_cols(Predictions_DICneg50 %>%
              select(.prediction, .category) %>%
              pivot_wider(names_from = .category, values_from = .prediction) %>%
              select(DICdiffstd))%>%
  ungroup()%>%
  select(-c(.row, DayNight1, Season1))

Predictions_pHneg50<-predicted_draws(k_fit_brms, newdata=newdata_pHneg50)%>%
  filter(.category %in% c("pHstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for NEC +50
newdata_NEC<-Predictions_pH %>%
  select(-c(pHstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  bind_cols(newdata_DIC %>%
              select(Tempinstd)) %>%
  select(-c(.row, .draw, .iteration))

Predictions_NEC<-predicted_draws(k_fit_brms, newdata=newdata_NEC)%>%
  filter(.category %in% c("TAdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for NEC -50
newdata_NECneg50<-Predictions_pHneg50 %>%
  select(-c(pHstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  bind_cols(newdata_DICneg50 %>%
              select(Tempinstd)) %>%
  select(-c(.row, .draw, .iteration))

Predictions_NECneg50<-predicted_draws(k_fit_brms, newdata=newdata_NECneg50)%>%
  filter(.category %in% c("TAdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

# backcalculate the predictions to the original scale
Predictions_nochange<-Predictions_nochange %>%
  mutate(predict_backcalc  = case_when(.category  == 'logNNstd' ~ .prediction*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
                                       .category  == 'pHstd' ~ .prediction*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
                                       .category  == 'Tempinstd' ~ .prediction*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
                                       .category  == 'DICdiffstd' ~ .prediction*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
                                       .category  == 'TAdiffstd' ~ .prediction*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center")
                                     ),
         Obs_pred = "Original SGD")  %>%
  select(.category,.prediction,predict_backcalc, Obs_pred, Season, DayNight, Tide)
                                       
                                       
                                       
## bring the predictions at higher SGD together for +50%
PredictionsChange<-Predictions_logNN_Temp %>%
  select(.category,.prediction, Season, DayNight, Tide) %>%
  bind_rows(Predictions_DIC%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_pH%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_NEC%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  mutate(predict_backcalc  = case_when(.category  == 'logNNstd' ~ .prediction*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
                                       .category  == 'pHstd' ~ .prediction*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
                                       .category  == 'Tempinstd' ~ .prediction*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
                                       .category  == 'DICdiffstd' ~ .prediction*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
                                       .category  == 'TAdiffstd' ~ .prediction*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center")
                                      ),
         Obs_pred = "Increased SGD") 
  

## bring the predictions at higher SGD together for -50%
PredictionsChangeneg50<-Predictions_logNN_Temp_neg50 %>%
  select(.category,.prediction, Season, DayNight, Tide) %>%
  bind_rows(Predictions_DICneg50%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_pHneg50%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_NECneg50%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  mutate(predict_backcalc  = case_when(.category  == 'logNNstd' ~ .prediction*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
                                       .category  == 'pHstd' ~ .prediction*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
                                       .category  == 'Tempinstd' ~ .prediction*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
                                       .category  == 'DICdiffstd' ~ .prediction*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
                                       .category  == 'TAdiffstd' ~ .prediction*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center")
  ),
  Obs_pred = "Decreased SGD") 

## bring both together
All_Predictions<-bind_rows(Predictions_nochange, PredictionsChange, PredictionsChangeneg50) %>%
  mutate(Site = "Kupikipiki'o") 


## Now for Wailupe ##########
Predictions_nochangeW<-predicted_draws(W_fit_brms, newdata=Cdata_W) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median) %>% # get the median of each draw
  mutate(Site = "Wailupe")%>%
  mutate(predict_backcalc  = case_when(.category  == 'logNNstd' ~ .prediction*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
                                       .category  == 'pHstd' ~ .prediction*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
                                       .category  == 'Tempinstd' ~ .prediction*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
                                       .category  == 'DICdiffstd' ~ .prediction*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
                                       .category  == 'TAdiffstd' ~ .prediction*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center")
  ),
  Obs_pred = "Original SGD")  %>%
  select(.category,.prediction,predict_backcalc, Obs_pred, Site, Tide)

# get predictions for NN and Temp +50%
Predictions_logNN_TempW<-predicted_draws(W_fit_brms, newdata=Cdata_W %>% select(!logSGDstd) %>% rename(logSGDstd = logSGD50std))%>%
  filter(.category %in% c("Tempinstd","logNNstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

# get predictions for NN and Temp -50%
Predictions_logNN_TempWneg50<-predicted_draws(W_fit_brms, newdata=Cdata_W %>% select(!logSGDstd) %>% rename(logSGDstd = logSGDneg50std))%>%
  filter(.category %in% c("Tempinstd","logNNstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for NEP +50%
newdata_DICW<-Predictions_logNN_TempW %>%
  select(-c(logNNstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  select(Season, Tide, DayNight, logNNstd, Tempinstd, TAdiffstd, pHstd, DICdiffstd, logSGDstd)

Predictions_DICW<-predicted_draws(W_fit_brms, newdata=newdata_DICW)%>%
  filter(.category %in% c("DICdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for NEP -50%
newdata_DICWneg50<-Predictions_logNN_TempWneg50 %>%
  select(-c(logNNstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  select(Season, Tide, DayNight, logNNstd, Tempinstd, TAdiffstd, pHstd, DICdiffstd, logSGDstd)

Predictions_DICWneg50<-predicted_draws(W_fit_brms, newdata=newdata_DICWneg50)%>%
  filter(.category %in% c("DICdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for pH +50%
newdata_pHW<-newdata_DICW %>%
  select(-DICdiffstd) %>%
  bind_cols(Predictions_DICW %>%
              select(.prediction, .category) %>%
              pivot_wider(names_from = .category, values_from = .prediction) %>%
              select(DICdiffstd))%>%
  ungroup()%>%
  select(-c(.row, DayNight1, Season1))

Predictions_pHW<-predicted_draws(W_fit_brms, newdata=newdata_pHW)%>%
  filter(.category %in% c("pHstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for pH -50%
newdata_pHWneg50<-newdata_DICWneg50 %>%
  select(-DICdiffstd) %>%
  bind_cols(Predictions_DICWneg50 %>%
              select(.prediction, .category) %>%
              pivot_wider(names_from = .category, values_from = .prediction) %>%
              select(DICdiffstd))%>%
  ungroup()%>%
  select(-c(.row, DayNight1, Season1))

Predictions_pHWneg50<-predicted_draws(W_fit_brms, newdata=newdata_pHWneg50)%>%
  filter(.category %in% c("pHstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for NEC +50%
newdata_NECW<-Predictions_pHW %>%
  select(-c(pHstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  bind_cols(newdata_DICW %>%
              select(Tempinstd)) %>%
  select(-c(.row, .draw, .iteration))

Predictions_NECW<-predicted_draws(W_fit_brms, newdata=newdata_NECW)%>%
  filter(.category %in% c("TAdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

## get predictions for NEC -50%
newdata_NECWneg50<-Predictions_pHWneg50 %>%
  select(-c(pHstd,Tempinstd))%>%
  pivot_wider(names_from = .category, values_from = .prediction) %>%
  ungroup()%>%
  bind_cols(newdata_DICWneg50 %>%
              select(Tempinstd)) %>%
  select(-c(.row, .draw, .iteration))

Predictions_NECWneg50<-predicted_draws(W_fit_brms, newdata=newdata_NECWneg50)%>%
  filter(.category %in% c("TAdiffstd")) %>%
  group_by(.row,.category, Season, DayNight, Tide)%>%
  summarise_if(is.numeric,median)

# backcalculate the predictions to the original scale

## bring the predictions at higher SGD together +50%
PredictionsChangeW<-Predictions_logNN_TempW %>%
  select(.category,.prediction, Season, DayNight, Tide) %>%
  bind_rows(Predictions_DICW%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_pHW%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_NECW%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  mutate(predict_backcalc  = case_when(.category  == 'logNNstd' ~ .prediction*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
                                       .category  == 'pHstd' ~ .prediction*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
                                       .category  == 'Tempinstd' ~ .prediction*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
                                       .category  == 'DICdiffstd' ~ .prediction*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
                                       .category  == 'TAdiffstd' ~ .prediction*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center")
  ),
  Obs_pred = "Increased SGD",
  Site = "Wailupe") 

## bring the predictions at higher SGD together -50%
PredictionsChangeWneg50<-Predictions_logNN_TempWneg50 %>%
  select(.category,.prediction, Season, DayNight, Tide) %>%
  bind_rows(Predictions_DICWneg50%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_pHWneg50%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  bind_rows(Predictions_NECWneg50%>%
              select(.category,.prediction, Season, DayNight, Tide))%>%
  mutate(predict_backcalc  = case_when(.category  == 'logNNstd' ~ .prediction*attr(Cdata$logNNstd,"scaled:scale")+attr(Cdata$logNNstd,"scaled:center"),
                                       .category  == 'pHstd' ~ .prediction*attr(Cdata$pHstd,"scaled:scale")+attr(Cdata$pHstd,"scaled:center"),
                                       .category  == 'Tempinstd' ~ .prediction*attr(Cdata$Tempinstd,"scaled:scale")+attr(Cdata$Tempinstd,"scaled:center"),
                                       .category  == 'DICdiffstd' ~ .prediction*attr(Cdata$DICdiffstd,"scaled:scale")+attr(Cdata$DICdiffstd,"scaled:center"),
                                       .category  == 'TAdiffstd' ~ .prediction*attr(Cdata$TAdiffstd,"scaled:scale")+attr(Cdata$TAdiffstd,"scaled:center")
  ),
  Obs_pred = "Decreased SGD",
  Site = "Wailupe") 
## bring both together
All_PredictionsW<-bind_rows(Predictions_nochangeW, PredictionsChangeW,PredictionsChangeWneg50)

## Bind Both BP and Wailupe
Final_Predictions<-bind_rows(All_Predictions, All_PredictionsW) %>%
  mutate(Obs_pred = fct_relevel(Obs_pred, levels =c("Original SGD","Increased SGD", "Decreased SGD")))

Final_Predictions %>% 
  ggplot(aes(x = predict_backcalc, fill = Obs_pred), alpha = 0.2)+
  geom_density(alpha = 0.2)+
  facet_wrap(~.category*Site, scales = "free", ncol = 2)+
  theme_bw()

# make a list so that I can selectively scale the y axes
scales_y <- list(
TAdiffstd = scale_y_continuous(),
pHstd = scale_y_continuous(),
DICdiffstd = scale_y_continuous(),
Tempinstd = scale_y_continuous(),
logNNstd = scale_y_log10()
)

cat.labs<-c("NEC umol kg-1", "pH","NEP umol kg-1","Temperature deg C","log N+N umol L-1")
names(cat.labs)<-c("TAdiffstd","pHstd","DICdiffstd","Tempinstd","logNNstd")

Final_Predictions %>% 
  ggplot(aes(x = Site, y = predict_backcalc, fill = Obs_pred), alpha = 0.2)+
  geom_violin(alpha = 0.2)+
  stat_summary(fun = "median", geom = "point",
               colour = "black", size = 4, position = position_dodge(width = 0.9), show.legend = FALSE) +
  scale_fill_manual(values = c("lightcoral","lightblue3", "lightgreen"))+
 # facet_wrap(~.category, scales = "free", ncol = 2)+
  theme_bw()+
  ylab("Predicted value")+
  xlab("")+
  theme(legend.title = element_blank())+
  facet_grid_sc(rows = vars(.category), 
                scales = list(y = scales_y),
                labeller = labeller(.category = cat.labs)) # this changes the axis to what I want



# calculate medians
Meds<-Final_Predictions %>%
  ungroup()%>%
  group_by(Site,.category, Obs_pred)%>%
  summarise(predict_backcalc = median(predict_backcalc))

Final_Predictions %>% 
 # filter(Site =="Kupikipiki'o")%>%
  ggplot(aes(x = predict_backcalc, fill = Obs_pred), alpha = 0.2)+
  geom_density(alpha = 0.2)+
#  stat_summary(fun = "median", geom = "point",
 #              colour = "black", size = 4, position = position_dodge(width = 0.9), show.legend = FALSE) +
  scale_fill_manual(values = c("lightcoral","lightblue3", "lightgreen"))+
  # facet_wrap(~.category, scales = "free", ncol = 2)+
  theme_bw()+
  xlab("Predicted value")+
  ylab("")+
  geom_vline(data = Meds, aes(xintercept =predict_backcalc, color = Obs_pred, group = .category))+ # add a line at the medians
 # theme(legend.title = element_blank())+
  facet_wrap(.category~Site,scales = "free", ncol = 2,
#  facet_grid_sc(cols = vars(.category), 
#                #rows = vars(Site),
 #               scales = list(x = scales_x, y = "free"),# this changes the axis to what I want
                labeller = labeller(.category = cat.labs), strip.position = "bottom") +
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        legend.title = element_blank()) +
  ggsave("Output/PredictionsDouble.pdf", width = 6, height = 8)

# calculate % change for each variable
Meds %>%
  group_by(Site, .category)%>%
  mutate(per.change = 100*(predict_backcalc[Obs_pred=="Increased SGD"] -
                         predict_backcalc[Obs_pred=="Original SGD"])/predict_backcalc[Obs_pred=="Original SGD"]) %>%
  mutate(per.changeneg50 = 100*(predict_backcalc[Obs_pred=="Decreased SGD"] -
                             predict_backcalc[Obs_pred=="Original SGD"])/predict_backcalc[Obs_pred=="Original SGD"]) %>%
  filter(Obs_pred =="Original SGD")
#  select(-c(Obs_pred, predict_backcalc)) %>%
 # pivot_wider(names_from = c(Site), values_from = c(per.change, per.changeneg50))


# calculate change in current and 50% increase predictions
predictionchange<-Final_Predictions %>%
  ungroup()%>%
  group_by(Site,.category, .row)%>%
  mutate(difference = predict_backcalc[Obs_pred =="Increased SGD"]-predict_backcalc[Obs_pred =="Original SGD"],
         differencepercent = 100*(difference/abs(predict_backcalc[Obs_pred =="Original SGD"]))) %>% 
  mutate(differenceneg50 = predict_backcalc[Obs_pred =="Decreased SGD"]-predict_backcalc[Obs_pred =="Original SGD"],
         differencepercentneg50 = 100*(differenceneg50/abs(predict_backcalc[Obs_pred =="Original SGD"]))) %>% 
    filter(Obs_pred =="Original SGD") # these are now differences so only pull out one set so it doesnt repeat twice
  
predictionchange %>% 
  ggplot(aes(x = differencepercent, fill = Site), alpha = 0.2)+
  geom_density(alpha = 0.2)+
  geom_vline(aes(xintercept = 0), lty =2)+
  scale_fill_manual(values = c("firebrick4","gold"))+
  facet_wrap(~.category, scales = "free", nrow = 1,
             labeller = labeller(.category = cat.labs), strip.position = "bottom") +
    xlab("Difference in predictions between current and 50% increase in SGD")+
  ylab("Density")+
  theme_few()+
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        legend.title = element_blank()) +
  ggsave("Output/predictiondifference50.pdf", width = 13, height = 4, useDingbats = FALSE)
  
  # count the # of times values increased or decreases or stayed the same 
# predictionchange %>%
#   group_by(.category, Site,) %>%
#   summarise(n_samples = n(),
#             n_lessthan = sum(difference < 0),
#             p_decreased = round(100*(n_lessthan / n_samples)),2 )%>%
#   select(.category, Site, p_decreased) %>%
#   pivot_wider(names_from = .category, values_from = p_decreased) %>%

  # join with the lat and long data
.row <- c(seq(1,nrow(Cdata_BP),1),seq(1,nrow(Cdata_W),1))
mappredictions<-Cdata %>%
  select(Lat, Long, Site) %>%
  mutate(.row = .row) %>%
  mutate(Site = case_when(Site=="BP"~"Kupikipiki'o",
                          Site=="W"~"Wailupe")) %>%
  left_join(predictionchange) %>% # get average change by location
  group_by(Site, Lat, Long, .category)%>%
  summarise(mean.val = median(difference),
            mean.val_per = median(differencepercent),
            mean.valneg50 = median(differenceneg50),
            mean.val_perneg50 = median(differencepercentneg50),)
  
mappredictions %>%
  group_by(.category, Site,) %>%
  summarise(n_samples = n(),
            n_lessthan = sum(mean.val < 0),
            p_decreased = round(100*(n_lessthan / n_samples))) %>%
  select(.category, Site, p_decreased) %>%
  pivot_wider(names_from = Site, values_from = p_decreased) %>%
  na.omit()

mappredictions %>%
  group_by(.category, Site) %>%
  summarise(min50 = min(mean.val_per,na.rm=TRUE), max50 = max(mean.val_per,na.rm=TRUE),
            minneg50 = min(mean.val_perneg50,na.rm=TRUE), maxneg50 = max(mean.val_perneg50,na.rm=TRUE))
  