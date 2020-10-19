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
#  coord_trans(x="log10", y="log10")+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')+
  geom_smooth(method = "lm")+
  xlab(expression("log(Radon) DPM m"^{-3}))+
  ylab(expression(paste("log(Silicate) ", mu,"mol L"^{-1})))+
  theme_bw()+
  ggsave("Output/radon_silicate.png")

mod<-lm(log(SI)~log(Radon), data = Data)
anova(mod)
summary(mod)
