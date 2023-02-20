library(tidyverse)
library(effects)
require(pscl)
library(MASS)

#Data contains variables year, latitude and longitude of station, Depth category
#eDNA concentration (ng/microlitre) and acoustic registration (NASC - nautical area scattering coefficient)
#Presence absence is assigned based on presence of acoustic registration
dat <- read_csv("capelin_eDNA_acoustic.csv") %>% 
          mutate(pa = ifelse(SA > 0, 1, 0),
                 Depth = factor(Depth))
          
view(dat)

#------------- Generalized Linear Model ----------------------------

#------------- Modelling Presence Absence -------------------
#Research question
#Is there presence of fish when eDNA is present (binomial)?

#Binomial model
glm.pa <- glm(factor(pa) ~  sqrt(eDNA), 
              data = dat ,
              family=binomial)
summary(glm.pa)
par(mfrow=c(2,2))
plot(glm.pa)
deviance(glm.pa)


all.effects <- allEffects(mod = glm.pa)

png(filename = "plot_binomial.png", width = 10, height = 6, units = 'in', res = 300)
plot(all.effects, type = "response", ylim = c(0, 1), ylab = "Presence Absence of capelin", xlab = "eDNA concentration")
dev.off()



#--------------- Modelling concentration ----------------------------
#Is SA higher when eDNA concentration is higher?

#------------- Plotting data distribution ------------------

p1 <- ggplot(dat %>% 
               filter(SA < 20000)) +
  geom_histogram(aes(SA), binwidth = 1000) 

p2 <- ggplot(dat %>% 
               filter(SA < 20000)) +
  geom_histogram(aes(log(SA+1))) 

p3 <- ggplot(dat %>% 
               filter(SA < 20000)) +
  geom_histogram(aes(SA^2)) 


gridExtra::grid.arrange(p1,p2,p3)

p4 <- ggplot(dat)  +
  geom_histogram(aes(eDNA)) 

p5 <- ggplot(dat) +
  geom_histogram(aes(log(eDNA+1))) 

p6 <- ggplot(dat) +
  geom_histogram(aes(eDNA^2)) 


gridExtra::grid.arrange(p4,p5,p6)

#Normal distribution not achieved because data has lots of zeros.
#Other transformations such as log+1 and box-cox doesn't normalize the data distribution either
#However the Gaussian link assumes normal distribution of errors from the model.
#There are no assumptions necessary about the distribution of original data.
#Model adequacy was verified by examination of residuals

#------------------ Gaussain model ---------------------------------

glm.conc <- glm(sqrt(SA) ~  sqrt(eDNA), 
           data = dat ,
           family=gaussian)

summary(glm.conc)
par(mfrow=c(2,2))
plot(glm.conc)

p1.conc <- plot(predictorEffect("eDNA", glm.conc), xlab ="eDNA concentration",
           axes=list(grid=TRUE,
                     x=list(rug=FALSE),
                     y=list(lab="Strength of acoustic registration")))

png(filename = "glm_concentration.png", width = 10, height = 6, units = 'in', res = 300)
p1.conc
dev.off()

#Three ways to check normality of residual
#Residual versus fitted
#If the red trend line is approximately flat and close to zero one can assume the residual are normally distributed.

png(file="residual_fitted.png", width = 150 ,height = 100, units = 'mm', res = 300)
plot(glm.conc, which = 1)
dev.off()


#Histogram of residuals
resid_hist <- data.frame(glm.conc$residuals) %>% 
  ggplot(aes(glm.conc.residuals)) + 
  #geom_histogram(aes(x=glm$residuals), binwidth = 20)
  geom_histogram(aes(y=..density..,), position="identity", alpha=0.5, col='grey50', fill='grey50')+
  geom_density(aes(glm.conc.residuals),adjust=2.5,color='blue') +
  ylab("Density") + xlab("Model residuals")

png(file="residual_histogram.png", width = 150 ,height = 100, units = 'mm', res = 300)
resid_hist
dev.off()

#Boxplot - if median is around 0 normality can be assumed
png(file="residual_boxplot.png", width = 150 ,height = 100, units = 'mm', res = 300)
boxplot(glm.conc$residuals, main="Residual Box Plot")
dev.off()

#Another test for normality is to plot residuals against the dependent variable which should show 
#a linear trend where capelin is present (zeros will group in one corner).
data <- data.frame(x = dat %>% dplyr::select(SA), residuals = glm.conc$residuals)

#Plot y against residuals
ggplot(data, 
       aes(x = SA, y = residuals)) + geom_point()
#This is highlighted by removing two high outliers in the acoustic registration
#Plot y against residuals
ggplot(data %>% filter(SA < 20000), 
       aes(x = SA, y = residuals)) + geom_point() 
  

#------------------ Negative binomial model--------------------------------------

#Exploring GLM with negative binomial distribution because of zeros in the data
#Really high acoustic registrations are excluded to enable model convergence
glm.conc.nb <- glm.nb(SA ~  eDNA, 
                   data = dat %>% filter(SA < 20000), link = "sqrt",
                   control = glm.control(maxit = 500), init.theta = 1)


summary(glm.conc.nb)
par(mfrow=c(2,2))
plot(glm.conc.nb)


#all.effects.nb <- allEffects(mod = glm.conc.nb)
#plot(all.effects.nb, type = "link")

pl.conc.nb <- plot(predictorEffect("eDNA", glm.conc.nb),xlab="eDNA concentration",
                axes=list(grid=TRUE,
                          x=list(rug=TRUE),
                          y=list(lab="Strength of acoustic registration")))

png(file="glm_negbin.png", width = 150 ,height = 100, units = 'mm', res = 300)
pl.conc.nb
dev.off()


