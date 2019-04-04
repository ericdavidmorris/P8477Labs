## Lab 7 Stochastic models
## Ebola: Legrand et al. 2007
## A transition rate, depending only on the present state of the population, 
## is allocated to each transition labmda_i (see Table 2). 
## At each iteration of the algorithm, a time tau_i is drawn from an exponential distribution 
## with parameter labmda_i for each transition. 
## The next transition m is the transition that has the minimum time to
# occurence (tau_m). Counts in each compartment are updated accordingly.
# rate of change
# rate.se=S/N*(beta.I*I+beta.H*H+beta.F*FF);  # S->E
# rate.ei=alpha*E;  # E->I
# rate.ih=gamma.h*theta1*I; # I->H
# rate.hf=gamma.dh*delta2*H;  # H->F
# rate.fr=gamma.f*FF;  # F->R
# rate.ir=gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
# rate.if=delta1*(1-theta1)*gamma.d*I; # I->F
# rate.hr=gamma.ih*(1-delta2)*H; # H->R


################################################################################
# PART 1: STOCHASTIC APPROACH: TIME TO NEXT EVENT
################################################################################
# RUN THE MODEL WITHOUT INTERVENTION
# parameters (no control): time in week
alpha=1; # incubation period: 7 days = 1 week
gamma.h=7/5; # from onset to hopspitalization: 5 days
gamma.d=7/9.6; # from onset to death: 9.6 days
gamma.i=7/10; # from onset to end of infectiousness for survivors: 10 days
gamma.f=7/2; # from death to traditional burial 2 days
gamma.ih=7/(10-5); # from hospitalization to end of infectiousness for survivors
gamma.dh=7/(9.6-5); # from hospitalization to death
theta1=.67; # proportion infectious in the hospital
delta1=.8; # CFR for unhospitalized
delta2=.8; # CFR for hospitalize
beta.I=.588; # transmission rate in the community
beta.H=.794; # transmission rate in the hospital
beta.F=7.653; # transmission rate at funerals; 

# INTIAL CONDITIONS
N=2e5; E0=H0=FF0=R0=0; I0=3; S0=N-I0; 
tm=0; # in week
S=S0; E=E0; I=I0; H=H0; FF=FF0; R=R0; cumI=I0;
res=c(tm,NA,S,E,I,H,FF,R,cumI); 
num_wk=10; # NUMBER OF WEEK, CHANGE ACCORDINGLY
## RUN THE MODEL FOR num_wk WEEKS WITHOUT CONTROL MEASURES
## DEPENDING ON WHETHER THE EPID TAKES OFF AFTER THE FIRST FEW CASES,
## THE TIME TO COMPLETE THE RUN VARIES
## IT WILL TAKE MUCH LONGER TO RUN IF YOU INCREASE THE NUM_WK TO LARGER NUMBER
while((I>0 | E>0 | H>0 | FF>0) & S>0 & tm<num_wk){ # simulate epid for the first 20 week without ctrl
  # step 1: compute the transition rates
  rate.se=S/N*(beta.I*I+beta.H*H+beta.F*FF);  # S->E
  rate.ei=alpha*E;  # E->I
  rate.ih=gamma.h*theta1*I; # I->H
  rate.hf=gamma.dh*delta2*H;  # H->F
  rate.fr=gamma.f*FF;  # F->R
  rate.ir=gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
  rate.if=delta1*(1-theta1)*gamma.d*I; # I->F
  rate.hr=gamma.ih*(1-delta2)*H; # H->R
  rate.tot=rate.se+rate.ei+rate.ih+rate.hf+rate.fr+rate.ir+rate.if+rate.hr
  
  # step 2: draw a random number, u1, betw 0 and 1 and 
  # calcualte the time after which the next transition occurs
  u1=runif(1); # draw a random number 
  tau_i=-log(u1)/rate.tot

  # step 3: compute the probability that each type of transition will occur based on the rates
  # use this to calculate the range in which a number drawn at random 
  # must lie for a given transition to occur
  p.se=rate.se/rate.tot;  # prob the next event is S->E
  p.ei=rate.ei/rate.tot;  # prob the next event is E->I
  p.ih=rate.ih/rate.tot;  # prob the next event is I->H
  p.hf=rate.hf/rate.tot;  # prob the next event is H->F
  p.fr=rate.fr/rate.tot;  # prob the next event is F->R
  p.ir=rate.ir/rate.tot;  # prob the next event is I->R
  p.if=rate.if/rate.tot;  # prob the next event is I->F
  p.hr=rate.hr/rate.tot;  # prob the next event is H->R
  
  # step 4: draw a random number, u2, to determine the transition event which occurs next.
  u2=runif(1); # draw a random number 
  
  # step 5: use the result from step 4 to update the number of state varibles
  if(u2<p.se) { # u2 lies in [0,p.se), the next event is S->E
    S=S-1; E=E+1;
  } else if (u2<p.se+p.ei) { # u2 lies in [p.se,p.se+p.ei), the next event is E->I
    E=E-1; I=I+1;
    cumI=cumI+1; # record cumulative incidence
  } else if (u2<p.se+p.ei+p.ih) { # u2 lies in [p.se+p.ei,p.se+p.ei+p.ih), the next event is I->H
    I=I-1; H=H+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf) { # u2 lies in [p.se+p.ei+p.ih,p.se+p.ei+p.ih+p.hf), the next event is H->F
    H=H-1; FF=FF+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr){ # FF->R
    FF=FF-1; R=R+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir){ # I->R
    I=I-1; R=R+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir+p.if) { # I->FF
    I=I-1; FF=FF+1;
  } else { # H->R
    H=H-1; R=R+1;
  }
  tm=tm+tau_i; # update time
  # save the results
  res=rbind(res,c(tm,tau_i,S,E,I,H,FF,R,cumI))
} 
colnames(res)=c('time','tm_step','S','E','I','H','FF','R','cumI')


## [LQ1]: run 3 times and plot results for each run
# run and save 3 times

plot(res[, 'time'],res[, 'cumI'],type='l', lwd=2, col=c('red'),
     lty=1, main = 'Third Run - Cumulative Incidence vs. Time', cex.main = 1, ylab='Cumulative Incidence',xlab='Time (weeks)')


## [LQ2]: time step and # event per week
## 1. run it for 10 weeks
## 2. compile the results to weekly interval 
## and compute # events for each week
wkly.res=matrix(NA,num_wk,8); colnames(wkly.res)=c('time','S','E','I','H','FF','R','cumI')
num_events=NULL;
for(wk in 1:num_wk){
  idx=tail(which(res[,'time']<wk),1);
  # idx also records the number of events happened from time 0 to the current week
  num_events=c(num_events,ifelse(wk==1,idx,idx-sum(num_events)));
  wkly.res[wk,]=res[idx,c('time','S','E','I','H','FF','R','cumI')];
}
tmsteps=res[-1,'tm_step']


plot(wkly.res[, 'time'],num_events,type='l', lwd=2, col=c('red'),
     lty=1, main = 'Q2 - # of Events over 10 weeks', cex.main = 1, ylab='Number of Events',xlab='Time (weeks)')


plot(tmsteps,type='l', lwd=2, col=c('red'),
     lty=1, main = 'Q2 - Time Steps', cex.main = 1, ylab='Time Steps',xlab='Time Step #')



## [LQ3] RUN TIME VS. SIMULATION TIME
## Run the model for 9, 11, 13, 15, and 17 weeks, respectively, 
## and record the computing time for each simulation. 
## You can use the App (“RShinyApp_StochModelRunTime.R”) to do so; 
## computing time is shown on the top of the plot. 
## How dose the computing time change with number of weeks simulated? 
## Does it change linearly (which is the case for simulation with fixed time step)? (1pt)



plot(c(9, 11, 13, 15, 17), c(0.057, 0.069, 2, 8.7, 100),type='l', lwd=2, col=c('red'),
     lty=1, main = 'Time to run simulation' , cex.main = 1, ylab='Seconds',xlab='Time (weeks) of epidemic', ylim = c(0,100))




################################################################################
## PART 2. UNCERTAINTY TEST
################################################################################
## RUN THE MODEL WITH CONTROL MEASURES AT WEEK 9
## DEPENDING ON WHETHER THE EPID TAKES OFF AFTER THE FIRST FEW CASES,
## THE TIME TO COMPLETE THE RUN VARIES
## RUN THE MODEL ONCE TO SEE HOW IT WORKS
z=1-.12; # intervention in the community
z.H=1; # intervention in the hospital
z.F=1; # intervention in safe burial
N=2e5; E0=H0=FF0=R0=0; I0=3; S0=N-I0; 
tm=0; # in week
S=S0; E=E0; I=I0; H=H0; FF=FF0; R=R0; cumI=I0;
res=c(tm,S,E,I,H,FF,R,cumI); 
while((I>0 | E>0 | H>0 | FF>0) & S>0){ # simulate epid for the first 9 weeks without ctrl
  # step 1: compute the transition rates
  if (tm<9){ # no intervention
    beta.I=.588; # transmission rate in the community
    beta.H=.794; # transmission rate in the hospital
    beta.F=7.653; # transmission rate at funerals; 
  } else { # with intervention on the three transmission routes
    beta.I=.588 * (1-z); # transmission rate in the community
    beta.H=.794 * (1-z.H); # transmission rate in the hospital
    beta.F=7.653 * (1-z.F); # transmission rate at funerals
  }
  rate.se=S/N*(beta.I*I+beta.H*H+beta.F*FF);  # S->E
  rate.ei=alpha*E;  # E->I
  rate.ih=gamma.h*theta1*I; # I->H
  rate.hf=gamma.dh*delta2*H;  # H->F
  rate.fr=gamma.f*FF;  # F->R
  rate.ir=gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
  rate.if=delta1*(1-theta1)*gamma.d*I; # I->F
  rate.hr=gamma.ih*(1-delta2)*H; # H->R
  rate.tot=rate.se+rate.ei+rate.ih+rate.hf+rate.fr+rate.ir+rate.if+rate.hr
  
  # step 2: draw a random number, u1, betw 0 and 1 and 
  # calcualte the time after which the next transition occurs
  u1=runif(1); # draw a random number 
  tau_i=-log(u1)/rate.tot
  
  # step 3: compute the probability that each type of transition will occur based on the rates
  # use this to calculate the range in which a number drawn at random 
  # must lie for a given transition to occur
  p.se=rate.se/rate.tot;  # prob the next event is S->E
  p.ei=rate.ei/rate.tot;
  p.ih=rate.ih/rate.tot;
  p.hf=rate.hf/rate.tot;
  p.fr=rate.fr/rate.tot;
  p.ir=rate.ir/rate.tot;
  p.if=rate.if/rate.tot;
  p.hr=rate.hr/rate.tot;
  
  # step 4: draw a random number, u2, to determine the transition event which occurs next.
  u2=runif(1); # draw a random number 
  
  # step 5: use the result from step 4 to update the number of state varibles
  if(u2<p.se) { # u2 lies in [0,p.se), the next event is S->E
    S=S-1; E=E+1;
  } else if (u2<p.se+p.ei) { # u2 lies in [p.se,p.se+p.ei), the next event is E->I
    E=E-1; I=I+1;
    cumI=cumI+1; # record cumulative incidence
  } else if (u2<p.se+p.ei+p.ih) { # u2 lies in [p.se+p.ei,p.se+p.ei+p.ih), the next event is I->H
    I=I-1; H=H+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf) { # u2 lies in [p.se+p.ei+p.ih,p.se+p.ei+p.ih+p.hf), the next event is H->F
    H=H-1; FF=FF+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr){ # FF->R
    FF=FF-1; R=R+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir){ # I->R
    I=I-1; R=R+1;
  } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir+p.if) { # I->FF
    I=I-1; FF=FF+1;
  } else { # H->R
    H=H-1; R=R+1;
  }
  tm=tm+tau_i; # update time
  # save the results
  res=rbind(res,c(tm,S,E,I,H,FF,R,cumI))
} 
colnames(res)=c('time','S','E','I','H','FF','R','cumI')
## compile the results to weekly interval
num_wk=ceiling(res[nrow(res),'time'])
wkly.res=matrix(NA,num_wk,8); colnames(wkly.res)=c('time','S','E','I','H','FF','R','cumI')
num_events=NULL;
for(wk in 1:num_wk){
  idx=tail(which(res[,'time']<wk),1);
  # idx also records the number of events happened from time 0 to the current week
  num_events=c(num_events,ifelse(wk==1,idx,idx-sum(num_events)));
  wkly.res[wk,]=res[idx,];
}






## RUN THE MODEL FOR 1000 TIMES
## NOTE: THIS WILL TAKE A WHILE (WAIT... WAIT... WAIT...)
## PATIENT...PATIENT...PATIENT...
num_run=1000;
out=matrix(0,num_run,4); colnames(out)=c('num_run','duration','peak','final.size');
out[,1]=1:num_run;
for(ir in 1:num_run){
  N=2e5; E0=H0=FF0=R0=0; I0=3; S0=N-I0; 
  tm=0; # in week
  S=S0; E=E0; I=I0; H=H0; FF=FF0; R=R0; cumI=I0;
  res=c(tm,S,E,I,H,FF,R,cumI); 
  while((I>0 | E>0 | H>0 | FF>0) & S>0){ # simulate epid for the first 20 week without ctrl
    # step 1: compute the transition rates
    if (tm<9){ # no intervention
      beta.I=.588; # transmission rate in the community
      beta.H=.794; # transmission rate in the hospital
      beta.F=7.653; # transmission rate at funerals; 
    } else { # with intervention on the three transmission routes
      beta.I=.588 * (1-z); # transmission rate in the community
      beta.H=.794 * (1-z.H); # transmission rate in the hospital
      beta.F=7.653 * (1-z.F); # transmission rate at funerals
    }
    rate.se=S/N*(beta.I*I+beta.H*H+beta.F*FF);  # S->E
    rate.ei=alpha*E;  # E->I
    rate.ih=gamma.h*theta1*I; # I->H
    rate.hf=gamma.dh*delta2*H;  # H->F
    rate.fr=gamma.f*FF;  # F->R
    rate.ir=gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
    rate.if=delta1*(1-theta1)*gamma.d*I; # I->F
    rate.hr=gamma.ih*(1-delta2)*H; # H->R
    rate.tot=rate.se+rate.ei+rate.ih+rate.hf+rate.fr+rate.ir+rate.if+rate.hr
    
    # step 2: draw a random number, u1, betw 0 and 1 and 
    # calcualte the time after which the next transition occurs
    u1=runif(1); # draw a random number 
    tau_i=-log(u1)/rate.tot
    
    # step 3: compute the probability that each type of transition will occur based on the rates
    # use this to calculate the range in which a number drawn at random 
    # must lie for a given transition to occur
    p.se=rate.se/rate.tot;  # prob the next event is S->E
    p.ei=rate.ei/rate.tot;
    p.ih=rate.ih/rate.tot;
    p.hf=rate.hf/rate.tot;
    p.fr=rate.fr/rate.tot;
    p.ir=rate.ir/rate.tot;
    p.if=rate.if/rate.tot;
    p.hr=rate.hr/rate.tot;
    
    # step 4: draw a random number, u2, to determine the transition event which occurs next.
    u2=runif(1); # draw a random number 
    
    # step 5: use the result from step 4 to update the number of state varibles
    if(u2<p.se) { # u2 lies in [0,p.se), the next event is S->E
      S=S-1; E=E+1;
    } else if (u2<p.se+p.ei) { # u2 lies in [p.se,p.se+p.ei), the next event is E->I
      E=E-1; I=I+1;
      cumI=cumI+1; # record cumulative incidence
    } else if (u2<p.se+p.ei+p.ih) { # u2 lies in [p.se+p.ei,p.se+p.ei+p.ih), the next event is I->H
      I=I-1; H=H+1;
    } else if (u2<p.se+p.ei+p.ih+p.hf) { # u2 lies in [p.se+p.ei+p.ih,p.se+p.ei+p.ih+p.hf), the next event is H->F
      H=H-1; FF=FF+1;
    } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr){ # FF->R
      FF=FF-1; R=R+1;
    } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir){ # I->R
      I=I-1; R=R+1;
    } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir+p.if) { # I->FF
      I=I-1; FF=FF+1;
    } else { # H->R
      H=H-1; R=R+1;
    }
    tm=tm+tau_i; # update time
    # save the results
    res=rbind(res,c(tm,S,E,I,H,FF,R,cumI))
  } 
  colnames(res)=c('time','S','E','I','H','FF','R','cumI')
  ## compile the results to weekly interval
  num_wk=ceiling(res[nrow(res),'time']) # that is the epidemic duration
  wkly.res=matrix(NA,num_wk,8); colnames(wkly.res)=c('time','S','E','I','H','FF','R','cumI')
  num_events=NULL;
  for(wk in 1:num_wk){
    idx=tail(which(res[,'time']<wk),1);
    # idx also records the number of events happened from time 0 to the current week
    num_events=c(num_events,ifelse(wk==1,idx,idx-sum(num_events)));
    wkly.res[wk,]=res[idx,];
  }
  wklyInci=wkly.res[-1,'cumI']-wkly.res[-nrow(wkly.res),'cumI']
  
  ## save result for this run
  out[ir,'duration']=num_wk;
  out[ir,'peak']=max(wklyInci);
  out[ir,'final.size']=wkly.res[nrow(wkly.res),'cumI'];
  print(c(paste('run #',ir),out[ir,2:4]));
}

hist(out[,'final.size'], col=c('red'), main = 'Final Epidemic Size', cex.main = 1, ylab = 'Frequency of Final Size', xlab = 'Final Epidemic Size', ylim = c(0, 350), xlim = c(0, 2500))

hist(out[,'peak'], col=c('red'), main = 'Epidemic Peak', cex.main = 1, ylab = 'Frequency of Epidemic Peak', xlab = 'Epidemic Peak', ylim = c(0, 350), xlim = c(0, 600))

## NOTE: the matrix "out" records the results 
## (epidemic duration, peak magnitude, and final epidemic size)
## from each run in each row of the matrix

## PLOT THE DISTRIBUTION
## HINT: USE THE FUNCITON: hist() to plot histogram

