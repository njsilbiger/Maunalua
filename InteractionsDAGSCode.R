### pull out all the code I wrote for estimates, interactions, and DAGS. Not needed anymore... right now

## Black Point #############
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



######## Wailupe ################
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

