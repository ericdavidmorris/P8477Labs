# Lab week 8: vector-borne infections

library(deSolve)

# function for a simple mosquito borne disease
SIRMosVec = function(time, state, parms) {
  with(as.list(c(state, parms)), {
    # infection in humans
    dSH = vH - SH * r * (THM * IM) - muH * SH
    dIH = SH * r * (THM * IM) - gamma * IH - muH * IH
    # dRH = gamma * IH - muH * RH
    
    # infection in mosquitoes
    dSM = vM - SM * r * (TMH * IH) - muM * SM
    dIM = SM * r * (TMH * IH) - muM * IM
    list(c(dSH, dIH, dSM, dIM))
  })
}

############################################################
## PART 1 Model Dynamics
############################################################
# inital states and parameters
NH=1e7; IH = 1; SH=NH-IH;
muH=1/50/365; # human life span: 50 yr
vH=NH*muH;
TMH = 0.8; # prob infection from human to mosquito;
THM = 0.5; # prob infection from mosquito to human;
gamma=1/7; # infectious period: 7 days
NM=1e8; IM=1; SM=NM-IM;
muM = 1/7; # 1 week life span for mosquito
vM=NM*muM; # births in mosquito population
b=.5; # number of bite per mosquito per day;
r = b / NH; # bite rate per human per mosquito

parameters = c(muH = muH, muM = muM,
               vH = vH, vM = vM, 
               THM = THM, TMH = TMH, 
               gamma = gamma, r = b / NH)

state = c(SH = SH, IH = 1,SM = SM,  IM = 1)

times=1:(365*100);
## solve the odes using R ode sovler:
sim=ode(y=state,times=times,func=SIRMosVec,parms = parameters)

# plot results, e.g. for Humans:
par(mfrow=c(2,1),mar=c(3,3,1,1), cex=1, mgp=c(1.8,.5,0))
matplot(sim[,'time'],sim[,c('SH','IH')],type='l', 
        log='y', # NOTE: THE Y-AXIS IS ON LOG SCALE
        lwd=1,col=c('blue','red'),lty=1, main='Humans', cex.main=1,
        ylab='Numbers (on log scale)',xlab='Time (days)')
legend('bottomright',c('SH','IH'),col=c('blue','red'),
       lty=1, cex=1, lwd=1, bty='n')

# now plot results for Mosquitoes:

par(mfrow=c(2,1),mar=c(3,3,1,1), cex=1, mgp=c(1.8,.5,0))
matplot(sim[,'time'],sim[,c('SM','IM')],type='l', 
        log='y', # NOTE: THE Y-AXIS IS ON LOG SCALE
        lwd=1,col=c('blue','red'),lty=1, main='Mosquitoes', cex.main=1,
        ylab='Numbers (on log scale)',xlab='Time (days)')
legend('bottomright',c('SM','IM'),col=c('blue','red'),
       lty=1, cex=1, lwd=1, bty='n')

summary(sim[,'IH'])

summary(sim[,'IM'])

# plot % immune vs time
par(mfrow=c(1,1),mar=c(3,3,.5,.5), cex=1.2, mgp=c(1.8,.5,0))
plot(sim[,'time'],(NH-sim[,'SH']-sim[,'IH'])/NH*100,type='l', lwd=2,
     lty=1, main = 'Humans', cex.main = 0.5, ylab='% Immune',xlab='Time (days)',ylim=c(0,100))

############################################################
## PART 2 Factors shaping the model Dynamics
############################################################
## Calculating R0
## check the equatioin for R0 in the lecture slides

r0 = ((b^2 * THM * TMH) / (muM * (gamma + muH))) * (NM/NH)

#####################
## Test NM/NH vs R0 vs epidemic dynamics
## use the R Shiny App: ShinyApp_VectorBorneDis.R

 