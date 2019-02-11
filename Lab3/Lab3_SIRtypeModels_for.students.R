## Lab 3: SIR type models

library(deSolve)

#######################################################
## PART 1 SIR model with demography
#######################################################
# HERE IS THE LOVELY SIR MODEL WE USED LAST WEEK
SIR = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # rate of change
    dS = -beta*S*I/N;
    dI = beta*S*I/N - gamma*I;
    
    # NOTE: EQUATION FOR CUMULATIVE INCIDENCE IS SHOWN BELOW
    # dcumInic=beta*S*I/N;
    
    # return the rate of change
    list(c(dS,dI))  
  }) # end with(as.list...)
}

## SIR model with demography
SIRdem = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # rate of change
    dS = mu*N - beta*S*I/N - mu*S;
    dI = beta*S*I/N - gamma*I - mu*I;
    
    # return the rate of change
    list(c(dS,dI))
  }) # end with(as.list...)
}

## READ THE CODE AND COMPARE THE DIFFERENCES
## RUN BOTH WITH THE FOLLOWING 
## INITIAL CONDITIONS/PARAMETERS
N = 1e5; I0 = 20; S0 = 0.2*N;
state = c(S = S0, I = I0);
# parameters
beta = 520/365; gamma = 1/7; mu = 1/(70*365);
parameters.SIR = c(beta = beta, gamma = gamma); # ASSEMBLE YOUR PARAMETER SET FOR EACH MODEL
parameters.SIRdem = c(beta = beta, gamma = gamma, mu = mu); # ASSEMBLE YOUR PARAMETER SET FOR EACH MODEL

times = seq(0,365*100, by = 1);

sim.SIR = ode(y = state, times = times, func = SIR, parms = parameters.SIR);
s.SIR = sim.SIR[,'S']/N  # %S
i.SIR = sim.SIR[,'I']/N  # %I

# for the SIRdem model
sim.SIRdem = ode(y = state,times = times, func = SIRdem, parms = parameters.SIRdem);
s.SIRdem = sim.SIRdem[,'S']/N
i.SIRdem = sim.SIRdem[,'I']/N

# plot results: S & I vs time
par(mfrow = c(2,1),mar = c(3,3,1,1), mgp = c(1.8,.5,0))
plot(sim.SIRdem[,'time'],s.SIRdem,ylab = '% S',xlab = 'Time (day)', type = 'l',col = 'blue',lwd = 1)
plot(sim.SIRdem[,'time'],i.SIRdem,ylab = '% I',xlab = 'Time (day)', type = 'l',col = 'red',lwd = 1)


par(mfrow = c(2,1),mar = c(3,3,1,1), mgp = c(1.8,.5,0))
plot(sim.SIR[,'time'],s.SIR,ylab = '% S',xlab = 'Time (day)', type = 'l',col = 'blue',lwd = 1)
plot(sim.SIR[,'time'],i.SIR,ylab = '% I',xlab = 'Time (day)', type = 'l',col = 'red',lwd = 1)


## to see just the *first* few years: 
# set the limit for x-axis using xlim=c(xmin,xmax)
plot(sim.SIRdem[,'time'],s.SIRdem,ylab = '% S',xlab = 'Time',xlim = c(1,365*10), 
     type = 'l',col = 'blue')
plot(sim.SIRdem[,'time'],i.SIRdem,ylab = '% I',xlab = 'Time',xlim = c(1,365*10),
     type = 'l',col = 'red')

plot(sim.SIR[,'time'],s.SIR,ylab = '% S',xlab = 'Time',xlim = c(1,365*10), 
     type = 'l',col = 'blue')
plot(sim.SIR[,'time'],i.SIR,ylab = '% I',xlab = 'Time',xlim = c(1,365*10),
     type = 'l',col = 'red')

## to see just the *last* few years: 
# set the limit for x-axis using xlim=c(xmin,xmax)
plot(sim.SIRdem[,'time'],s.SIRdem,ylab = '% S',xlab = 'Time',xlim = c(365*90+1,365*100), 
     type = 'l',col = 'blue')
plot(sim.SIRdem[,'time'],i.SIRdem,ylab = '% I',xlab = 'Time',xlim=c(365*90+1,365*100),
     type = 'l',col = 'red')

#######################################################
## PART 2: SIS model
#######################################################
SIS = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    
    # EQUATIONS for the SIS:
    dS = -beta * S * I/N + gamma * I;
    dI = beta * S * I/N - gamma * I;  # FILL IN THE EQUATION HERE
    
    # return the rate of change
    list(c(dS,dI))
  }) # end with(as.list...)
}

## AFTER FINISHING YOUR SIS MODEL, RUN WITH THE FOLLOWING 
## INITIAL CONDITIONS/PARAMETERS
N = 1e5; I0 = 10; S0 = N - I0;
state = c(S = S0, I = I0);
parameters.SIS = c(beta = .5,gamma = .3);

times = seq(0,100, by = 1);
sim.SIS = ode(y = state, times = times, func = SIS, parms = parameters.SIS);
s.SIS = sim.SIS[,'S']
i.SIS = sim.SIS[,'I']

# plot results: S & I vs time
par(mfrow = c(2,1),mar = c(3,3,1,1), mgp = c(1.8,.5,0))
plot(sim.SIS[,'time'],s.SIS,ylab = '% S',xlab = 'Time (day)', type = 'l',col = 'blue',lwd = 1)
plot(sim.SIS[,'time'],i.SIS,ylab = '% I',xlab = 'Time (day)', type = 'l',col = 'red',lwd = 1)

#######################################################
## PART 3: SEIR MODEL (VS SIR)
#######################################################
## AGAIN, THE SIR MODEL FOR YOUR REFERENCE
## NOTE THE EXTRA EQN FOR CUMMULATIVE INCIDENCE

SIR = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dS = -beta*S*I/N;
    dI = beta*S*I/N-gamma*I;
    
    # cumulative incidence
    dcumInci = beta*S*I/N;
    
    # return the rate of change
    list(c(dS,dI,dcumInci))
  }) # end with(as.list...)
}


SEIR=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dS = -beta*S*I/N;
    dE = beta * S * I / N - alpha*E;
    dI = alpha * E - gamma *I;
    
    # FILL IN YOUR SEIR MODEL HERE

    #cumulative incidence
    dcumInci = beta * S * I /N;
    # return the rate of change
    # RETURN ALL STATE VARIABLES: 
    list(c(dS,dE,dI,dcumInci))
  }) # end with(as.list...)
}

## NOW TEST YOUR MODELS USING THE FOLLOWING PARMS/INI CONDS:
times=0:100
N=1e5; E0=10; I0=0; S0=N-E0-I0;
stateSEIR=c(S=S0,E=E0,I=I0, cumInci=0);
paramSEIR=c(beta=.6,alpha=1/1.5,gamma=1/3); 
simSEIR=ode(y=stateSEIR,times=times,func=SEIR,parms=paramSEIR) # run SEIR

stateSIR=c(S=S0,I=E0, cumInci=E0);
paramSIR=c(beta=.4,gamma=1/4.5);  # to keep the same serial interval and R0
simSIR=ode(y=stateSIR,times=times,func=SIR,parms=paramSIR) # run SEIR

s_SEIR=simSEIR[,'S']/N; i_SEIR=simSEIR[,'I']/N; 
s_SIR=simSIR[,'S']/N; i_SIR=simSIR[,'I']/N;

# TO CALCULATE THE INCIDENCE USING THE CUMMULATIVE INCIDENCE
# RECALL: X[-1] MEANS DELETE THE FIRST ELEMENT IN X
newi_SEIR=simSEIR[-1,'cumInci']/N-simSEIR[-length(times),'cumInci']/N; # new cases
newi_SIR=simSIR[-1,'cumInci']/N-simSIR[-length(times),'cumInci']/N; # new cases


# TO UNDERSTAND THE LINE 'PAR', READ ABOUT WHAT EACH ARGUEMENT IS DOING
# FROM, E.G. http://www.statmethods.net/advgraphs/parameters.html
# GOOGLE IF YOU'D LIKE TO LEARN MORE
par(mfrow=c(3,1),cex=.8,mgp=c(1.8,.5,0),mar=c(3,3,1,1))
plot(times,s_SEIR,ylim=c(0,1),
     ylab='% S',xlab='Time', type='l',col='blue')
lines(times,s_SIR,col='blue',lty=2)
legend('topright',cex = .8,
       legend = c('SEIR Model','SIR Model'),
       col = c('blue','blue'), lty = c(1,2),
       bty = 'n')
# REMEMBER TO ADD A LEGEND INDICATING WHAT YOU'RE PLOTTING
ymax=max(i_SEIR,i_SIR)*1.05
plot(times,i_SEIR,ylim=c(0,ymax),col='red',ylab='% I',xlab='Time', type='l')
lines(times,i_SIR,col='red',lty=2)
legend('topright',cex = .8,
       legend = c('SEIR Model','SIR Model'),
       col = c('red','red'), lty = c(1,2),
       bty = 'n')
# REMEMBER TO ADD A LEGEND INDICATING WHAT YOU'RE PLOTTING
ymax=max(newi_SEIR,newi_SIR)*1.05
plot(newi_SEIR,ylim=c(0,ymax),col='green',ylab='Incidence',xlab='Time', type='l')
lines(newi_SIR,col='green',lty=2)
legend('topright',cex = .8,
       legend = c('SEIR Model','SIR Model'),
       col = c('green','green'), lty = c(1,2),
       bty = 'n')
# REMEMBER TO ADD A LEGEND INDICATING WHAT YOU'RE PLOTTING
