# Use data from Nelson et al. 2015 to assess the relationship between radon and silicate

#load libraries
library(tidyverse)

# read in the data
SI<-c(2.2,
      3.9,
      7.1,
      13,
      23,
      42,
      75,
      134,
      241,
      432,
      774)

Radon <- c(1.1,
           2.1,
           4.3,
           8.7,
           18,
           36,
           53,
           80,
           120,
           180,
           270)

Data<-data.frame(SI,Radon)

ggplot(Data, aes(x = Radon, y = SI))+
  geom_point()+
 # coord_trans(x="log", y="log")+
  scale_x_continuous(trans='log') +
  scale_y_continuous(trans='log')+
  geom_smooth(method = "lm")+
  xlab("log(Radon)")+
  ylab("log(Silicate)")+
  theme_bw()

mod<-lm(log(SI)~log(Radon), data = Data)
anova(mod)
summary(mod)
