## Lab 3: SIR type models

library(deSolve)


#######################################################
## PART 3: SEIR MODEL (VS SIR)
#######################################################
## AGAIN, THE SIR MODEL FOR YOUR REFERENCE
## NOTE THE EXTRA EQN FOR CUMMULATIVE INCIDENCE


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
times=0:200
N=90000; E0=500; I0=0; S0=N-E0-I0;
stateSEIR=c(S=S0,E=E0,I=I0, cumInci=0);
paramSEIR1=c(beta=1.4,alpha=1/4.3,gamma=1/2.5); 
paramSEIR2=c(beta=0.692,alpha=1/4.3,gamma=1/2.5); 
paramSEIR3=c(beta=0.52,alpha=1/4.3,gamma=1/2.5); 
paramSEIR4=c(beta=0.4114,alpha=1/4.3,gamma=1/3.5); 

simSEIR1=ode(y=stateSEIR,times=times,func=SEIR,parms=paramSEIR1) # run SEIR
simSEIR2=ode(y=stateSEIR,times=times,func=SEIR,parms=paramSEIR2) # run SEIR
simSEIR3=ode(y=stateSEIR,times=times,func=SEIR,parms=paramSEIR3) # run SEIR
simSEIR4=ode(y=stateSEIR,times=times,func=SEIR,parms=paramSEIR4) # run SEIR

s_SEIR1=simSEIR1[,'S']/N; i_SEIR1=simSEIR1[,'I']/N; 
s_SEIR2=simSEIR2[,'S']/N; i_SEIR2=simSEIR2[,'I']/N; 
s_SEIR3=simSEIR3[,'S']/N; i_SEIR3=simSEIR3[,'I']/N; 
s_SEIR4=simSEIR4[,'S']/N; i_SEIR4=simSEIR4[,'I']/N; 


# TO CALCULATE THE INCIDENCE USING THE CUMMULATIVE INCIDENCE
# RECALL: X[-1] MEANS DELETE THE FIRST ELEMENT IN X
newi_SEIR1=simSEIR1[-1,'cumInci']/N-simSEIR1[-length(times),'cumInci']/N; # new cases
newi_SEIR2=simSEIR2[-1,'cumInci']/N-simSEIR2[-length(times),'cumInci']/N; # new cases
newi_SEIR3=simSEIR3[-1,'cumInci']/N-simSEIR3[-length(times),'cumInci']/N; # new cases
newi_SEIR4=simSEIR4[-1,'cumInci']/N-simSEIR4[-length(times),'cumInci']/N; # new cases


# TO UNDERSTAND THE LINE 'PAR', READ ABOUT WHAT EACH ARGUEMENT IS DOING
# FROM, E.G. http://www.statmethods.net/advgraphs/parameters.html
# GOOGLE IF YOU'D LIKE TO LEARN MORE
par(mfrow=c(3,1),cex=.8,mgp=c(1.8,.5,0),mar=c(3,3,1,1))
plot(times,s_SEIR1,ylim=c(0,1),
     ylab='% S',xlab='Time', type='l',col='blue')
lines(times,s_SEIR2,col='green')
lines(times,s_SEIR3,col='red')
lines(times,s_SEIR4,col='orange')
legend('bottomleft',cex = .8,
       legend = c('R0=3.5','R0=1.73', 'R0=1.3', 'R0=1.44'),
       col = c('blue','green', 'red', 'orange'), lty = c(1,1,1,1),
       bty = 'n')
# REMEMBER TO ADD A LEGEND INDICATING WHAT YOU'RE PLOTTING
ymax=max(i_SEIR,i_SIR)*1.05
plot(times,i_SEIR1,ylim=c(0,ymax),col='blue',ylab='% I',xlab='Time', type='l')
lines(times,i_SEIR2,col='green')
lines(times,i_SEIR3,col='red')
lines(times,i_SEIR4,col='orange')
legend('topright',cex = .8,
       legend = c('R0=3.5','R0=1.73', 'R0=1.3', 'R0=1.44'),
       col = c('blue','green', 'red', 'orange'), lty = c(1,1,1,1),
       bty = 'n')
# REMEMBER TO ADD A LEGEND INDICATING WHAT YOU'RE PLOTTING
ymax=max(newi_SEIR,newi_SIR)*1.05
plot(newi_SEIR1,ylim=c(0,ymax),col='blue',ylab='Incidence',xlab='Time', type='l')
lines(newi_SEIR2,col='green')
lines(newi_SEIR3,col='red')
lines(newi_SEIR4,col='orange')
legend('topright',cex = .8,
       legend = c('R0=3.5','R0=1.73', 'R0=1.3', 'R0=1.44'),
       col = c('blue','green', 'red', 'orange'), lty = c(1,1,1,1),
       bty = 'n')
# REMEMBER TO ADD A LEGEND INDICATING WHAT YOU'RE PLOTTING