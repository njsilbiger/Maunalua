
#calculate residuls for TA and DIC analysis
#libraries
library(seacarb)

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


#plot relationship between sgd and delta TA
model.SGDTA<-lm(Cdata$TA.diff~Cdata$percent_sgd)

par(mfrow=c(1,1))
plot(Cdata$percent_sgd[Cdata$Site=='BP'],Cdata$TA.diff[Cdata$Site=='BP'], xlab='Percent SGD', ylab='delta TA')
points(Cdata$percent_sgd[Cdata$Site=='W'],Cdata$TA.diff[Cdata$Site=='W'], col='red')
lines(seq(0,100,0.01),seq(0,100,0.01)*model.SGDTA$coefficients[2]+model.SGDTA$coefficients[1] )
legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)

par(mfrow=c(1,1))
plot(Cdata$percent_sgd[Cdata$Site=='BP' & Cdata$Tide=='L2'],Cdata$TA.diff[Cdata$Site=='BP'& Cdata$Tide=='L2'])
points(Cdata$percent_sgd[Cdata$Site=='W'& Cdata$Tide=='L2'],Cdata$TA.diff[Cdata$Site=='W'& Cdata$Tide=='L2'], col='red')
legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)


#DIC
model.SGDDIC<-lm(Cdata$DIC.diff~Cdata$percent_sgd)

par(mfrow=c(1,1))
plot(Cdata$percent_sgd[Cdata$Site=='BP'],Cdata$DIC.diff[Cdata$Site=='BP'], xlab='Percent SGD', ylab='delta DIC')
points(Cdata$percent_sgd[Cdata$Site=='W'],Cdata$DIC.diff[Cdata$Site=='W'], col='red')
#lines(seq(0,100,0.01),seq(0,100,0.01)*model.SGDTA$coefficients[2]+model.SGDTA$coefficients[1] )
legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)


#regression between pH and delta pH
pH2<-Cdata$pH^2
model.pHTA<-lm(Cdata$TA.diff~Cdata$pH)

par(mfrow=c(1,1))
plot(Cdata$pH[Cdata$Site=='BP'],Cdata$TA.diff[Cdata$Site=='BP'], xlab = 'pH', ylab = 'Delta TA')
points(Cdata$pH[Cdata$Site=='W'],Cdata$TA.diff[Cdata$Site=='W'], col='red')
lines(seq(7.8,8.8,0.01),seq(7.8,8.8,0.01)*model.pHTA$coefficients[2]+model.pHTA$coefficients[1] )
legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)

#DIC
#model.pHDIC<-lm(Cdata$DIC.diff~Cdata$pH)

#par(mfrow=c(1,1))
#plot(Cdata$pH[Cdata$Site=='BP'],Cdata$DIC.diff[Cdata$Site=='BP'], xlab = 'pH', ylab = 'Delta DIC')
#points(Cdata$pH[Cdata$Site=='W'],Cdata$DIC.diff[Cdata$Site=='W'], col='red')
#lines(seq(7.8,8.8,0.01),seq(7.8,8.8,0.01)*model.pHDIC$coefficients[2]+model.pHDIC$coefficients[1] )
#legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)

#pH as a function of DIC
model.DICpH<-lm(Cdata$pH~Cdata$DIC.diff)

par(mfrow=c(1,1))
plot(Cdata$DIC.diff[Cdata$Site=='BP'],Cdata$pH[Cdata$Site=='BP'], ylab = 'pH', xlab = 'Delta DIC')
points(Cdata$DIC.diff[Cdata$Site=='W'],Cdata$pH[Cdata$Site=='W'], col='red')
lines(seq(-100,300,0.1),seq(-100,300,0.1)*model.DICpH$coefficients[2]+model.DICpH$coefficients[1] )
legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)


model.pHSGD<-lm(Cdata$pH~Cdata$percent_sgd)

# nutrients
par(mfrow=c(1,1))
plot(Cdata$NN[Cdata$Site=='BP'],Cdata$TA.diff[Cdata$Site=='BP'], xlab = 'NN', ylab = 'Delta TA')
points(Cdata$NN[Cdata$Site=='W'],Cdata$TA.diff[Cdata$Site=='W'], col='red')
legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)

#TA vs DIC plots
plot(Cdata$DIC,Cdata$TA, xlab='DIC', ylab='TA')
plot(Cdata$DIC.diff,Cdata$TA.diff, xlab='delta DIC', ylab='delta TA')


#just low tide
#pH as a function of DIC
model.DICpH<-lm(Cdata$pH~Cdata$DIC.diff)

par(mfrow=c(1,1))
plot(Cdata$DIC.diff[Cdata$Site=='BP' & Cdata$Tide=='L1'],Cdata$pH[Cdata$Site=='BP'& Cdata$Tide=='L1'], ylab = 'pH', xlab = 'Delta DIC')
points(Cdata$DIC.diff[Cdata$Site=='W'& Cdata$Tide=='L1'],Cdata$pH[Cdata$Site=='W'& Cdata$Tide=='L1'], col='red')

points(Cdata$DIC.diff[Cdata$Site=='BP' & Cdata$Tide=='L2'],Cdata$pH[Cdata$Site=='BP'& Cdata$Tide=='L2'], ylab = 'pH', xlab = 'Delta DIC', pch=19)
points(Cdata$DIC.diff[Cdata$Site=='W'& Cdata$Tide=='L2'],Cdata$pH[Cdata$Site=='W'& Cdata$Tide=='L2'], col='red', pch=19)

lines(seq(-100,300,0.1),seq(-100,300,0.1)*model.DICpH$coefficients[2]+model.DICpH$coefficients[1] )
legend('topright', legend = c('Black Point', 'Wailupe'), col=c('black','red'), pch=21)
