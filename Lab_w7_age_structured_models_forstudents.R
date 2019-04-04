## Lab week 7: age-structured model

library(deSolve)

SIR2ageGrs <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSC = nu - SC * (betaCC * IC + betaCA * IA) - muC * SC - lC * SC
    dIC = SC * (betaCC * IC + betaCA * IA) - gamma * IC - muC * IC - lC * IC
    dSA = lC * SC - SA * (betaAC * IC + betaAA * IA) - muA * SA
    dIA = lC * IC + SA * (betaAC * IC + betaAA * IA) - gamma * IA - muA * IA
    
    list(c(dSC, dIC, dSA, dIA))
  })
}

########################################################################################
## PART 1: Two AGE GROUP MODEL
########################################################################################
## parameters/inital conditions
betaCC = 100 * 11; betaCA = betaAC = 10 * 11; betaAA = 20 * 11;
gamma = 365 / 14; # 2 weeks  
lC = 0.066667; # 1/15; 
muC = 0; muA = 0.0166667; #1/(75-15); # mean life expectancy=75yr; subtract the 15 yrs in youth with 0 death rate
nC = muA/(lC + muA); nA = 1 - nC;
nu = (lC + muA) * nC;

IC0 = IA0 = .0001; SC0 = SA0 = .1; 

## try the two age group model
parameters=c(betaCC=betaCC,betaCA=betaCA, betaAC=betaAC, betaAA=betaAA,
             gamma=gamma,lC=lC,muC=muC,muA=muA)
state=c(SC=SC0,IC=IC0,SA=SA0,IA=IA0);
times=seq(1,100,by=1/365) # TIME IS IN YEAR
sim=ode(times = times, func = SIR2ageGrs,y = state, parms = parameters)

## proportion susceptible in the two groups
fSC=sim[,'SC']/nC; # note: it is normalized by population size in each group
fSA=sim[,'SA']/nA;
## fS=sim[,c('SC','SA')]/matrix(c(nC,nA),nrow(sim),2,byrow=T) # do it altogether
## proportion infectious in the two groups
fIC=sim[,'IC']/nC;
fIA=sim[,'IA']/nA;

fS=sim[,c('SC','SA')]/matrix(c(nC,nA),nrow(sim),2,byrow=T)

## plot # people infected
totI=rowSums(sim[,c('IC','IA')]);  # total infectious
totS=rowSums(sim[,c('SC','SA')]);  # total susecptible


#par(mfrow = c(1,1), cex = 1.2, mgp = c(2,.5,0), mar = c(3,3,1,1))
matplot(fS, type = 'l', lty = 1, lwd = 2, col = c('red','blue'), 
        ylab = 'Fraction of Susceptibles (Si/ni)', xlab = "Time in days (100 years)",
        xlim=c(0,365*100)) 
legend('topright',c('Children','Adults'), col = c('red','blue'), lty=1,bty='n')

# plot Si/ni


#R0 calculation

beta=matrix(c(betaCC, betaCA, betaAC, betaAA),2,2) #matrix(c(10, 0.1, 0.1, 1),2,2)   # FILL IN YOUR BETA MATRIX HERE
n=c(nC,nA)      # n is the vector storing th proportion in each group
n.matrix=diag(n,2,2)  # matrix related to the population size in each group
# to see it:
View(n.matrix)

gamma=365 / 14;
R.matrix=n.matrix %*% beta / gamma
# to see the output of the eigen function:
eigen(R.matrix)
## To find R0
R0=eigen(R.matrix)$values[1]

## or directly:
R0=eigen(n.matrix %*% beta)$values[1]/gamma


########################################################################################
## PART 2: VACCINATION
########################################################################################
## FIRST: MODIFY THE SIR2ageGrs MODEL TO INCLUDE VACCINATION
## USE THE EQUATIONS PROVIDED IN THE SLIDES
SIR2ageGrsVac <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    
    dSC = nu*(1 - p) - SC * (betaCC * IC + betaCA * IA) - muC * SC - lC * SC
    dIC = SC * (betaCC * IC + betaCA * IA) - gamma * IC - muC * IC - lC * IC
    dSA = lC * SC - SA * (betaAC * IC + betaAA * IA) - muA * SA
    dIA = lC * IC + SA * (betaAC * IC + betaAA * IA) - gamma * IA - muA * IA
    
    list(c(dSC, dIC, dSA, dIA))
  })
}

state = c(SC = SC0,IC = IC0,SA = SA0,IA = IA0)
parametersVac = c(betaCC = betaCC,betaCA = betaCA, betaAC = betaAC, betaAA = betaAA,
                gamma = gamma,lC = lC,muC = muC,muA = muA,p = 0.5)
times = seq(1,100, by = 1/365)
sim = ode(times = times, func = SIR2ageGrsVac,y = state, parms = parametersVac)

#finding Is
fIC = sim[,'IC']/nC;
fIA = sim[,'IA']/nA;

fI = sim[,c('IC','IA')]/matrix(c(nC,nA),nrow(sim),2,byrow = T)

#plot
matplot(fI, type = 'l', lty = 1, lwd = 2, col = c('red','blue'), 
        ylab = 'Fraction of Infectious (Ii/ni)', xlab = "Time in days (100 years)",
        xlim = c(0,365*100)) 
legend('topright',c('Children','Adults'), col = c('red','blue'), lty = 1,bty = 'n')



####################
## TEST DIFFERENT VACCINATION RATES
## use the R Shiny App: ShinyApp_Vac.R
