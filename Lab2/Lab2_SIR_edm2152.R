## Lab 2: SIR MODEL
## SCRIPT TO TEST THE SIR MODEL

## if you haven't had it installed yet, type in the following first:
## install.packages('deSolve')
library(deSolve)

## PRE-LAB: TRY THE ODE SOLVER
## SIMPLE ODE: dy/dx=x^2
myfunc = function(x,y,parms){
  dy = x^2
  list(dy)
}
xs = seq(0,10,by = .2)
state = c(y = 0); 
out = ode(y = state, times = xs, func = myfunc, parms = NULL)

y.true = 1/3*xs^3

## check if the ode is working
par(mar = c(3,3,1,1), mgp = c(1.8,.5,0), cex = 1.2)
plot(out[,1], out[,2], xlab = 'x',ylab = 'y', type = 'l')
points(xs,y.true,col = 'red', pch = 20)
legend('topleft', legend = c('Analytic solution: y.ture','Numerical solution: ODE'),
       col = c('red','black'), lty = c(NA,1), pch = c(20,NA), bty = 'n')

## try the SIR
SIR = function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dS = -beta*S*I/N;
    dI = beta*S*I/N - gamma*I;
    dR = gamma * I;
    # return the rate of change
    list(c(dS,dI, dR))
  }) # end with(as.list...)
}


# specify initial conditions/parameters
N = 1e5; # population size
I0 = 10; # initial No. of Infectious people
S0 = N - I0;  # initial No. of Susceptible people
state = c(S = S0, I = I0);  # store the inital conditions I & S together 
parameters = c(beta = .5, gamma = .3);  # store the model parameters

times = seq(0,100, by = 1);
sim = ode(y = state, times = times, func = SIR, parms = parameters);

View(sim)  # to see what's in the model output 'sim'
s = sim[,'S']/N
i = sim[,'I']/N


################################################################################
## YOU CAN USE THE FOLLOWING SECTION TO SEE THE CHANGES IN S AND I STEP BY STEP
################################################################################
## CHANGE IN NUMBER OF INFECTIOUS
numWk = 15; tmstep = 7; # weekly increase
wklyS = wklyI = rep(0,numWk + 1);
wklyS[1] = S0; wklyI[1] = I0;
wk = 1; # first week
times = seq(tmstep*(wk - 1),tmstep*wk, by = 1);
state = c(S = wklyS[wk], I = wklyI[wk])
sim = ode(y = state, times = times, func = SIR, parms = parameters);
wklyS[wk + 1] = tail(sim[,'S'],1); # save new state
wklyI[wk + 1] = tail(sim[,'I'],1);
plot(sim[,'time'], sim[,'I'], ylim = c(0, N*.12), xlim = c(0,tmstep*numWk),
     ylab = 'I', xlab = 'Time', pch = 20)
for (wk in 2:numWk) {
  times = seq(tmstep*(wk - 1), tmstep*wk, by = 1);
  state = c(S = wklyS[wk], I = wklyI[wk])
  sim = ode(y = state, times = times, func = SIR, parms = parameters);
  wklyS[wk + 1] = tail(sim[,'S'], 1); # save new state
  wklyI[wk + 1] = tail(sim[,'I'], 1);
  points(tail(sim[,'time'], tmstep), tail(sim[,'I'],tmstep), pch = 20)
  readline("Press <return to continue") ## PRESS ENTER TO PROCEED
}

## CHANGE IN THE NUMBER OF SUSCEPTIBLES
wklyS = wklyI = rep(0, numWk + 1);
wklyS[1] = S0; wklyI[1] = I0;
wk = 1; # first week
times = seq(tmstep*(wk - 1), tmstep*wk,by = 1);
state = c(S = wklyS[wk], I = wklyI[wk])
sim = ode(y = state, times = times, func = SIR, parms = parameters);
wklyS[wk + 1] = tail(sim[,'S'], 1); # save new state
wklyI[wk + 1] = tail(sim[,'I'], 1);
plot(sim[,'time'], sim[,'S']/N, ylim = c(0,1), xlim = c(0,tmstep*numWk),
     ylab = '% Susceptible', xlab = 'Time', pch = 20)
for (wk in 2:numWk) {
  times = seq(tmstep*(wk - 1), tmstep*wk, by = 1);
  state = c(S = wklyS[wk], I = wklyI[wk])
  sim = ode(y = state, times = times, func = SIR, parms = parameters);
  wklyS[wk + 1] = tail(sim[,'S'], 1); # save new state
  wklyI[wk + 1] = tail(sim[,'I'], 1);
  points(tail(sim[,'time'], tmstep), tail(sim[,'S']/N,tmstep), pch = 20)
  readline("Press <return to continue") ## PRESS ENTER TO PROCEED
}


## LAB REPORT EXERCISES:


################################################################################
## Q1: What do we find R (# recovered)?

# Compute and plot R over time for the following epidemic setting. What is the final epidemic size (i.e. R at the end of the epidemic)?
################################################################################

# Specifying the initial conditions/parameters:

N = 1e5; # N, the population size
I0 = 10; # I, the initial # of Infectious people
S0 = N - I0;  # S, the initial # of Susceptible people
R0 = N - I0 - S0; # Calculating R0 with the R = N - I - S equation
state = c(S = S0, I = I0, R = R0);  # Storing the inital conditions I & S & R together 
parameters = c(beta = .5, gamma = .3);  # Storing the model parameters beta and gamma

times = seq(0,365, by = 1); # Running for 365 days (a year)
sim = ode(y = state, times = times, func = SIR, parms = parameters);

View(sim)
s = sim[, 'S']/N
i = sim[, 'I']/N
r = sim[, 'R']/N

#plot R over time 

plot(times, sim[, 'R'],
     xlab = 'Day',
     ylab = 'Recovered (R)',
     type = 'l', lty = 1,
     col = 'red', lwd = 2,
     main = 'Plotting Recovered over the course of Epidemic Year')

sim[366, 'R']

# R = N-sim[,"I"]-sim[,"S"] if you don't edit the SIR

################################################################################
## Q2: TEST HOW BETA AND GAMMA (R0) AFFECT THE EPIDEMIC 
################################################################################

## CHANGE THE PARAMETERS SPECIFIC IN THE LAB ASSIGNMENT
#(1) beta=0.5; gamma=0.3 (2) beta=0.5; gamma=0.4 (3) beta=0.5; gamma=0.5 (4) beta=0.5; gamma=0.6

## ERIC NOTES: (beta/gamma = r_0, so as gamma increases, r_0 decreases)

## RUN THE MODEL AND SEE THE DIFFERENCE
## REPORT YOUR FINDINGS


N = 1e5; # N, the population size
I0 = 10; # I, the initial # of Infectious people
S0 = N - I0;  # S, the initial # of Susceptible people
R0 = N - I0 - S0; # Calculating R0 with the R = N - I - S equation
state = c(S = S0, I = I0, R = R0);  # Storing the inital conditions I & S & R together 
parameters = c(beta = .5, gamma = .6);  # Storing the model parameters beta and gamma

times = seq(0,365, by = 1); # Running for 365 days (a year)
sim = ode(y = state, times = times, func = SIR, parms = parameters);

plot(times, sim[, 'I'],
     xlab = 'Day',
     ylab = 'I',
     type = 'l', lty = 1,
     col = 'red', lwd = 2,
     main = 'Plotting I over the course of the Year, beta = .5, gamma = .6')


s = sim[, 'S']/N
i = sim[, 'I']/N
r = sim[, 'R']/N


sim[366, 'R']



################################################################################
## Qs3&4: TEST EXPONENTIAL PERIOD: the first few weeks
################################################################################


N = 1e5; # N, the population size
I0 = 10; # I, the initial # of Infectious people
S0 = N - I0;  # S, the initial # of Susceptible people
R0 = N - I0 - S0; # Calculating R0 with the R = N - I - S equation
state = c(S = S0, I = I0, R = R0);  # Storing the inital conditions I & S & R together 
parameters = c(beta = .5, gamma = .3);  # Storing the model parameters beta and gamma

times = seq(0,100, by = 1); # Running for 365 days (a year)
sim = ode(y = state, times = times, func = SIR, parms = parameters);

s = sim[,'S']/N
i = sim[,'I']/N


WkExp = 5; # CHANGE THE NUMBER OF WEEK HERE
Iexp = sim[seq(1, length = WkExp, by = 7), 'I']; 
# NOTE: THE TIME STEP IN THE SIMULATION IS DAY, BUT WE ARE LOOKING AT WEEK HERE
# SO WE USE 'seq(1,length=WkExp,by=7)' TO EXTRACT THE CORRESPONDING DATE FOR EACH WEEK
lnI = log(Iexp); # TAKE THE LOG
tt = seq(1, length = WkExp, by = 7)
fit = lm(lnI~tt) # LINEAR REGRESSION
slope = fit$coeff[2] ## extract the slope for the linear regression

par(mfrow = c(1,1), cex = 1.5)
plot(tt, lnI, xlab = 'Time', main = 'Number of Weeks = 8')
lines(tt, fit$coeff[1] + fit$coeff[2]*tt, col = 'red', lwd = 2)

## NOW TRY HOW LONG THE EXPONENTIAL PERIOD CAN LAST BY CHANING WkExp by yourself

################################################################################
## Q5: TEST EPIDEMIC BURNOUT
################################################################################

N = 1e5; # N, the population size
I0 = 10; # I, the initial # of Infectious people
S0 = N - I0;  # S, the initial # of Susceptible people
R0 = N - I0 - S0; # Calculating R0 with the R = N - I - S equation
state = c(S = S0, I = I0, R = R0);  # Storing the inital conditions I & S & R together 
parameters = c(beta = .5, gamma = .3);  # Storing the model parameters beta and gamma

times = seq(0,600, by = 1); # Running for 365 days (a year)
sim = ode(y = state, times = times, func = SIR, parms = parameters);

s = sim[,'S']/N
i = sim[,'I']/N

sim[601, 'S']

sim[601, 'I']

tail(sim,1)