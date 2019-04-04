## Lab week 6: ODE models with diff risk groups

library(deSolve)

########################################################################
## PART1: dynamics for a 2 risk group model
########################################################################
## ODE model for 2 risk groups
SIS2riskGrs <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSH = - betaHH * SH * IH - betaHL * SH * IL + gamma * IH
    dIH = betaHH * SH * IH + betaHL * SH * IL - gamma * IH
    dSL = - betaLL * SL * IL - betaLH * SL * IH + gamma * IL
    dIL = betaLL * SL * IL + betaLH * SL * IH - gamma * IL
    
    list(c(dSH, dIH, dSL, dIL))
  })
}
# Parameters and initial conditions.
parameters <- c(betaHH = 10.0, betaHL = 0.1, betaLH = 0.1,
                betaLL = 1, gamma = 1)
NH = 0.2; NL=1-NH;
IH=1e-8; SH=NH-IH;
IL=1e-5; SL=NL-IL;
state = c(SH=SH, IH = IH, SL=SL, IL = IL)

# Run the model, plot and show: 


times = seq(0, 30, by = 0.2)

sim_SIS2riskGrs = ode(y = state, times = times, func = SIS2riskGrs, parms = parameters);

fracs = sim_SIS2riskGrs[, 2:5]

matplot(fracs, type = 'l', lty = 1, col = c('red', 'red', 'blue', 'blue'))

# 1) the prevalence in each group over time

par(mfrow = c(1,1), cex = 1.2, mgp = c(2,.5,0), mar = c(3,3,1,1))
matplot(fracs[,c(2,4)], log = 'y', type = 'l', lty = 1, lwd = 2, col = c('red','blue'), 
        ylab = '% Infectious individuals', xlab = "Time", ylim = c(1e-8,1e-1),
        xlim=c(0,30*5),xaxt = 'n') # time step is .2, so need to rescale
axis(1,at = 0:30*5, lab = 0:30)
legend('bottomright',c('high risk','low risk'), col = c('red','blue'), lty=1,bty='n')


# 2) the %S in each group over time
par(mfrow=c(1,1),cex=1.2,mgp=c(2,.5,0),mar=c(3,3,1,1))
matplot(fracs[,c(1,3)],log='y',type='l',lty=1,lwd=2,col=c('red','blue'), 
        ylab='% S',xlab="Time", 
        xlim=c(0,30*5),xaxt='n') # time step is .2, so need to rescale
axis(1,at=0:30*5,lab=0:30)
legend('right',c('high risk','low risk'),col=c('red','blue'),lty=1,bty='n')

########################################################################
## PART2: Calculating R0 for the risk structured model 
## and test how R0 change with degree of assortative mixing
########################################################################
## 1. calculate R0 for the two beta matrix

beta=matrix(c(1, 1.5, 1.5, 0.5),2,2) #matrix(c(10, 0.1, 0.1, 1),2,2)   # FILL IN YOUR BETA MATRIX HERE
n=c(NH,NL)      # n is the vector storing th proportion in each group
n.matrix=diag(n,2,2)  # matrix related to the population size in each group
# to see it:
View(n.matrix)

gamma=1;
R.matrix=n.matrix %*% beta / gamma
# to see the output of the eigen function:
eigen(R.matrix)
## To find R0
R0=eigen(R.matrix)$values[1]

## or directly:
R0=eigen(n.matrix %*% beta)$values[1]/gamma

## 2. proportion to both the contact numbers of the source group and the sink group
beta=matrix(0,2,2)
Ncontact=c(5,1); # number of contact in each group
Nk=c(.2,.8); # proportion of populaiton in each group
M=sum(Ncontact*Nk); # mean number of contact 
## to consider assortative mixing, use a to change the level of assortativeness
## a is the proportion of within group mixing
a=1.0; # change a to 0, .2, .4, .6, .8, 1
for(i in 1:2){
  for (j in 1:2){
    if(i==j) {
      beta[i,j]=a*Ncontact[i]/Nk[i]+(1-a)*Ncontact[i]^2/M
    } else {
      beta[i,j]=(1-a)*Ncontact[i]*Ncontact[j]/M
    }
  }
}

R0=eigen(diag(Nk, 2, 2) %*% beta)$values[1]/gamma
R0

a_values = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
r0_values = c(3.222222, 3.491356, 3.808143, 4.169353, 4.568974, 5)

matplot(a_values, r0_values, type='b',lty=1,lwd=2, pch = 20, col=c('red'), 
        ylab='R0 Values',xlab="a values") 

########################################################################
## Part 3: Simulating HIV early epidemic
########################################################################
# the HIV model with 5 risk groups
HIV5riskGrs <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dS1 = nu1-sum(BETA[1,] * c(I1,I2,I3,I4,I5))* S1 - mu*S1
    dS2 = nu2-sum(BETA[2,] * c(I1,I2,I3,I4,I5))* S2 - mu*S2
    dS3 = nu3-sum(BETA[3,] * c(I1,I2,I3,I4,I5))* S3 - mu*S3
    dS4 = nu4-sum(BETA[4,] * c(I1,I2,I3,I4,I5))* S4 - mu*S4
    dS5 = nu5-sum(BETA[5,] * c(I1,I2,I3,I4,I5))* S5 - mu*S5
    
    dI1 = sum(BETA[1,] * c(I1,I2,I3,I4,I5))* S1 - mu*I1 - gamma*I1
    dI2 = sum(BETA[2,] * c(I1,I2,I3,I4,I5))* S2 - mu*I2 - gamma*I2
    dI3 = sum(BETA[3,] * c(I1,I2,I3,I4,I5))* S3 - mu*I3 - gamma*I3
    dI4 = sum(BETA[4,] * c(I1,I2,I3,I4,I5))* S4 - mu*I4 - gamma*I4
    dI5 = sum(BETA[5,] * c(I1,I2,I3,I4,I5))* S5 - mu*I5 - gamma*I5
    
    dA1 = d*gamma*I1 - mu*A1 - m*A1
    dA2 = d*gamma*I2 - mu*A2 - m*A2
    dA3 = d*gamma*I3 - mu*A3 - m*A3
    dA4 = d*gamma*I4 - mu*A4 - m*A4
    dA5 = d*gamma*I5 - mu*A5 - m*A5
    
    # cumulative incidence
    dcumA1 = d*gamma*I1
    dcumA2 = d*gamma*I2
    dcumA3 = d*gamma*I3
    dcumA4 = d*gamma*I4
    dcumA5 = d*gamma*I5
    
    list(c(dS1,dS2,dS3,dS4,dS5, dI1,dI2,dI3,dI4,dI5, 
           dA1,dA2,dA3,dA4,dA5, dcumA1,dcumA2,dcumA3,dcumA4,dcumA5))
  })
}
mu=.0312;  # in years
N=c(.06,.31,.52,.08,.03); # population in each group
NU=mu*N; # recruits to each group
gamma=.2; m=1; d=.3; # in years
beta=.0217; # scaling factor
BETA=matrix(c(rep(0,5),
              0,.65,2.15,12.9,21.5,
              0,2.15,7.17,43.1,71.8,
              0,12.9,43.1,258,431,
              0,21.5,71.8,431,718),5,5,byrow=T)*beta
I0=c(0,0,0,0,1e-5);
S0=N-I0;
A0=cumA0=rep(0,5);
parameters=c(BETA=BETA,mu=mu,gamma=gamma,m=m,d=d,
             nu1=NU[1],nu2=NU[2],nu3=NU[3],nu4=NU[4],nu5=NU[5])

state=c(S1=S0[1],S2=S0[2],S3=S0[3],S4=S0[4],S5=S0[5],
        I1=I0[1],I2=I0[2],I3=I0[3],I4=I0[4],I5=I0[5],
        A1=A0[1],A2=A0[2],A3=A0[3],A4=A0[4],A5=A0[5],
        cumA1=cumA0[1],cumA2=cumA0[2],cumA3=cumA0[3],cumA4=cumA0[4],cumA5=cumA0[5])
times=0:100;

simHIV=ode(y=state,times=times,func=HIV5riskGrs,parms=parameters);
# incidence per 1000 for each group
Incidence=(simHIV[-1,c('cumA1','cumA2','cumA3','cumA4','cumA5')]-
             simHIV[-length(times),c('cumA1','cumA2','cumA3','cumA4','cumA5')])/
          matrix(N,100,5,byrow=T)*1000
# total incidence per 1000 population
totInci=rowSums((simHIV[-1,c('cumA1','cumA2','cumA3','cumA4','cumA5')]-
                   simHIV[-length(times),c('cumA1','cumA2','cumA3','cumA4','cumA5')]))*1e3

S=simHIV[,c('S1','S2','S3','S4','S5')]/matrix(N,101,5,byrow=T)*100
totS=rowSums(simHIV[,c('S1','S2','S3','S4','S5')])*100

## plot results
par(mfrow=c(1,1),cex=1.2,mar=c(3,3,1,1),mgp=c(1.8,.5,0))
matplot(Incidence[,2:5],type='l',lty=2:5,col='black','red','blue','green','orange',lwd=1.5,ylim=c(0,35),
        ylab='Incidence AIDS per year per 1000',xlab='Time (years)')
lines(totInci,lwd=2)
legend('topright',leg=c('1-5','6-50','51-100','100+','total'),lty=c(2:5,1),
       col='black',lwd=c(rep(1.5,4),2),cex=.9,bty='n')

#plot sus
par(mfrow=c(1,1),cex=1.2,mar=c(3,3,1,1),mgp=c(1.8,.5,0))
matplot(S[,2:5],type='l',lty=2:5,col='black','red','blue','green','orange',lwd=1.5,ylim=c(0,100),
        ylab='% Susceptible AIDS per year per 1000',xlab='Time (years)')
lines(totS,lwd=2)
legend('right',leg=c('1-5','6-50','51-100','100+','total'),lty=c(2:5,1),
       col='black',lwd=c(rep(1.5,4),2),cex=.9,bty='n')


# [Q6] literature review on HIV/ADIS dynamics and (if you'd like, redo the simulation w. updated parms)