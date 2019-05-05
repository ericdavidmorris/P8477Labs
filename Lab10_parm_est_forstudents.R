## Lab 10: parameter estimation
## serial interval and R0
## install the packages first
require(flexsurv) # for survivial analysis
require(survival) # for survivial analysis
require(R0); # for estimating R0

####################################################################################
## Part 1. Estimate the serial inteval for SARS using data from Lipsitch et al. 2003
####################################################################################
## read in data
#getwd()  ## use "getwd()" to see your current directory and the path
#setwd('YOUR DIRECTORY SAVING THE DATA')
da.sars=read.csv('SARS_serial_intervals.csv')
## each row is 1 case
## 1st column 'tm1' is the lower bound of the interval
## 2nd column 'tm2' is the upper bound of the interval
## they are the same for this dataset, as we only have a single value for each case
## 3rd column 'event' is the The status indicator, normally 
## 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). 
## here we have event=1 (showing symptoms)

## before doing the analysis, plot the dataset to see how it looks
par(mar=c(3,3,1,1),cex=1.2,mgp=c(1.5,.5,0))
hist(da.sars[,1],breaks=20,col='grey',ylim=c(0,25), main='',
     xlab='Serial Interval (days)',ylab='Number of cases')

## to use the survival function in the package (either flexsurv or survival)
## we have to first covert the data to a survival object
## to do so, run the command:
xsurv=Surv(da.sars$tm1,da.sars$tm2,da.sars$event,type='interval')

####################################################
## Fit to the Weibull distribution
####################################################
## First try the survival function, with a Weibull distribution:
surv1=flexsurvreg(xsurv~1,dist='weibull')

surv1; # check the model output
surv1$res;  # model parm estimates are save in 'res'

## Calculate the mean based on the estimates
## for Weibull distribution:
surv1.mean=surv1$res['scale','est']*gamma(1+1/surv1$res['shape','est']); # the mean for weibull
surv1.sd=sqrt(surv1$res['scale','est']^2*(gamma(1+2/surv1$res['shape','est'])-(gamma(1+1/surv1$res['shape','est']))^2))
print(paste(round(surv1.mean,1),'+/-',round(surv1.sd,1)))

## compute the model fit:
tm=seq(min(da.sars[,'tm1']),max(da.sars[,'tm2']),by=.2) # more time points than the observed, so we can fill the gap
# compute the number of cases per the survival model
fit1=nrow(da.sars)*dweibull(tm,scale=surv1$res['scale','est'],
                            shape=surv1$res['shape','est']);

# Super-impose the model fit on the data for comparison
par(mar=c(3,3,1,1),cex=1.2,mgp=c(1.5,.5,0))
hist(da.sars[,1],breaks=20,col='grey',ylim=c(0,25), main='',
     xlab='Serial Interval (days)',ylab='Number of cases')
lines(tm,fit1,col='red',lwd=2)
legend('topright',cex=.9,seg.len = .8,
       legend=c('Observed','Fitted (Weibull)'),
       lty=c(0,1),pch=c(22,NA),lwd=c(NA,2),pt.bg = c('grey',NA),
       col=c('grey','red'),bty='n')


####################################################
## Fit to the Exponential distribution
####################################################
surv2=flexsurvreg(xsurv~1,dist='exponential')

surv2; # check the model output
surv2$res; # model parm estimates are save in 'res'

# compute the number of cases per the model
fit2=nrow(da.sars)*dexp(tm,rate=surv2$res['rate','est'])


####################################################
## Fit to the log-normal distribution
####################################################
surv3=flexsurvreg(xsurv~1,dist='lognormal')
surv3;  # check the model output
surv3$res;  # model parm estimates are save in 'res'
# compute the number of cases per the model
fit3=nrow(da.sars)*dlnorm(tm,meanlog=surv3$res['meanlog','est'],sdlog = surv3$res['sdlog','est'])



####################################################
# Plot results all together:
par(mar=c(3,3,1,1),cex=1.2,mgp=c(1.5,.5,0))
hist(da.sars[,1],breaks=20,col='grey',ylim=c(0,25), main='',xlim=c(0,24),
     xlab='Serial Interval (days)',ylab='Number of cases')
lines(tm,fit1,col='red',lwd=2)
lines(tm,fit2,col='blue',lwd=2)
lines(tm,fit3,col='orange',lwd=2)
legend('topright',cex=.9,seg.len = .8,
       legend=c('Observed','Fitted (Weibull)','Fitted (Exponential)','Fitted (Log-normal)'),
       lty=c(0,1,1,1),pch=c(22,NA,NA,NA),lwd=c(NA,2,2,2),pt.bg = c('grey',NA,NA,NA),
       col=c('grey','red','blue','orange'),bty='n')




############################################################
## LQ2: compute the mean, sd, and AIC and compare
############################################################
## Weibull
surv1.mean=surv1$res['scale','est']*gamma(1+1/surv1$res['shape','est']); # the mean for weibull
surv1.sd=sqrt(surv1$res['scale','est']^2*(gamma(1+2/surv1$res['shape','est'])-(gamma(1+1/surv1$res['shape','est']))^2))
surv1.AIC=surv1$AIC

## Exponential:
surv2.mean=1/surv2$res['rate','est'];
surv2.sd=1/surv2$res['rate','est'];
surv2.AIC=surv2$AIC

## Log-normal:
## NOTE: IT IS ON LOG SCALE
surv3.mean=exp(surv3$res['meanlog','est'])
surv3.sd=exp(surv3$res['sdlog','est'])
surv3.AIC=surv3$AIC



####################################################################################
## Part 2: Estimating R0 from the exponential growth phase of the epidemic
####################################################################################
## data: daily incidence during 1918 influenza pandemic in Germany (from 'R0' library)
da.flu=read.csv('data_1918pandemic_Germany.csv')


## Always plot and check the data first
par(mar=c(3,3,1,1),cex=1.2,mgp=c(1.5,.5,0))
plot(da.flu[,1],da.flu[,2],cex=2,xlab='',ylab='Cases')
## cumpute the cumulative incidence using the cumsum function
cumI=cumsum(da.flu[,2]); 

# plot and see:

plot(da.flu[,1],log(cumI),cex=2,xlab='',ylab='Log(Cumulative Incidence)')



# Fit the first 7, 14, 21 days:
D=3; # set the generation time to 3 days
Ndays=21;  # ADJUST THE NUMBER OF DAYS INCLUDED IN THE FIT HERE
tm1=1:Ndays;
fit1=lm(log(cumI[1:Ndays])~tm1)
summary(fit1)

# compute R0 based on the model-fit
R=1+fit1$coefficients[2]*D   
#slope:fit1$coefficients[2]



####################################################################################
## Part 3: R0 - Maximum Likelihood Estimation (MLE)
####################################################################################
data("Germany.1918") # Request the data (it's from the package)
Germany.1918 # print the data to see the structure

# First we need the distribution of generation time (i.e. serial interval)
mGT<-generation.time("gamma", c(2.45, 1.38))

# Maximum Likelihood Estimation using the est.R0.ML function
est.R0.ML(Germany.1918, # the data
          mGT, # the distribution of serial interval
          begin=1, # the start of the data
          end=14, # ADJUST THE NUMBER OF DAYS TO INCLUDE IN THE MODEL HERE
          range=c(0.01,50) # the range of possible values to test
          )



