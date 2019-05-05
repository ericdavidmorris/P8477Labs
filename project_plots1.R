library(deSolve)

plague = function(t,state,parameters) {
  with(as.list(c(state,parameters)), {
    
    # rate of change
    dS =  -S/N * (beta * I);
    dE = S/N * (beta * I) - alpha * E
    dI = alpha * E - I * (gamma.int * theta + gamma.d * (1 - theta));
    dINT = I * (gamma.int * theta) - INT * (gamma.i.r * (1 - delta.int) + gamma.int.d * delta.int) 
    dD = INT * (gamma.int.d * delta.int) + I * (gamma.d * (1 - theta))
    dR =  INT * (gamma.i.r * (1 - delta.int))
    
    dcumInci = alpha * E # cumulative incidence
    
    # return the rate of change
    list(c(dS,  dE,  dI,  dINT, dD, dR, dcumInci))
  }) # end with(as.list...)
}

####################################################################################
## 
####################################################################################
# Initial conditions:
N = 90000; E0 = INT0 = D0 = R0 = 0; I0 = 500; cumInci0 = 3; S0 = N - I0;
state = c(S = S0, E = E0, I = I0, INT = INT0, R = R0, D = D0, cumInci = cumInci0);

alpha = 1/4.3; # incubation period: 4.3 days
gamma.int = 1; # from I to intervention: 1 days
gamma.d = 1/1.9; # from I to death: 1.9 days
gamma.i.r = 1/7; # from intervention to recovery: 7 days
gamma.int.d = 1/3; # from intervention to death: 3 days
theta1 = 0 # proportion survey said would participate in intervention: 0, 0.25, 0.50, 0.75, 1 blue
theta2 = 0.25 #red
theta3 = 0.50 #green
theta4 = 0.75 #orange
theta5 = 1 #purple
delta.int_1 = 0.055 # CFR for intervention
delta.int_2 = 0.71; 
delta.noint = 1; # CFR for no intervention
#beta1 = 0.41;# transmission rate in the community, 0.41, 0.52, 1.4, 3.0
beta2 = 1.5;


parms_beta2_theta1_delta1 = c(alpha = alpha, # incubation period: 4.3 days
                gamma.int = gamma.int, # from I to intervention: 1 days
                gamma.d = gamma.d, # from I to death: 1.9 days
                gamma.i.r = gamma.i.r, # from intervention to recovery: 7 days
                gamma.int.d = gamma.int.d, # from intervention to death: 3 days
                theta = theta1, # proportion survey said would participate in intervention
                delta.int = delta.int_1, # CFR for intervention
                delta.noint = delta.noint, # CFR for no intervention
                beta = beta2 # transmission rate in the community
);

parms_beta2_theta2_delta1 = c(alpha = alpha, # incubation period: 4.3 days
                              gamma.int = gamma.int, # from I to intervention: 1 days
                              gamma.d = gamma.d, # from I to death: 1.9 days
                              gamma.i.r = gamma.i.r, # from intervention to recovery: 7 days
                              gamma.int.d = gamma.int.d, # from intervention to death: 3 days
                              theta = theta2, # proportion survey said would participate in intervention
                              delta.int = delta.int_1, # CFR for intervention
                              delta.noint = delta.noint, # CFR for no intervention
                              beta = beta2 # transmission rate in the community
);

parms_beta2_theta3_delta1 = c(alpha = alpha, # incubation period: 4.3 days
                              gamma.int = gamma.int, # from I to intervention: 1 days
                              gamma.d = gamma.d, # from I to death: 1.9 days
                              gamma.i.r = gamma.i.r, # from intervention to recovery: 7 days
                              gamma.int.d = gamma.int.d, # from intervention to death: 3 days
                              theta = theta3, # proportion survey said would participate in intervention
                              delta.int = delta.int_1, # CFR for intervention
                              delta.noint = delta.noint, # CFR for no intervention
                              beta = beta2 # transmission rate in the community
);

parms_beta2_theta4_delta1 = c(alpha = alpha, # incubation period: 4.3 days
                              gamma.int = gamma.int, # from I to intervention: 1 days
                              gamma.d = gamma.d, # from I to death: 1.9 days
                              gamma.i.r = gamma.i.r, # from intervention to recovery: 7 days
                              gamma.int.d = gamma.int.d, # from intervention to death: 3 days
                              theta = theta4, # proportion survey said would participate in intervention
                              delta.int = delta.int_1, # CFR for intervention
                              delta.noint = delta.noint, # CFR for no intervention
                              beta = beta2 # transmission rate in the community
);

parms_beta2_theta5_delta1 = c(alpha = alpha, # incubation period: 4.3 days
                              gamma.int = gamma.int, # from I to intervention: 1 days
                              gamma.d = gamma.d, # from I to death: 1.9 days
                              gamma.i.r = gamma.i.r, # from intervention to recovery: 7 days
                              gamma.int.d = gamma.int.d, # from intervention to death: 3 days
                              theta = theta5, # proportion survey said would participate in intervention
                              delta.int = delta.int_1, # CFR for intervention
                              delta.noint = delta.noint, # CFR for no intervention
                              beta = beta2 # transmission rate in the community
);



####################################
## RUN FOR 100 DAYS
####################################
times = 1:100 # run for 100 days
sim_theta1 = ode(y = state, times = times, func = plague, parms = parms_beta2_theta1_delta1 )
sim_theta2 = ode(y = state, times = times, func = plague, parms = parms_beta2_theta2_delta1 )
sim_theta3 = ode(y = state, times = times, func = plague, parms = parms_beta2_theta3_delta1 )
sim_theta4 = ode(y = state, times = times, func = plague, parms = parms_beta2_theta4_delta1 )
sim_theta5 = ode(y = state, times = times, func = plague, parms = parms_beta2_theta5_delta1 )

inci_theta1 = sim_theta1[seq(7, nrow(sim_theta1), by = 7),'cumInci'] - c(0, sim_theta1[seq(7, nrow(sim_theta1) - 7, by = 7),'cumInci']) # get weekly incidence
inci_theta2 = sim_theta2[seq(7, nrow(sim_theta2), by = 7),'cumInci'] - c(0, sim_theta2[seq(7, nrow(sim_theta2) - 7, by = 7),'cumInci']) # get weekly incidence
inci_theta3 = sim_theta3[seq(7, nrow(sim_theta3), by = 7),'cumInci'] - c(0, sim_theta3[seq(7, nrow(sim_theta3) - 7, by = 7),'cumInci'])
inci_theta4 = sim_theta4[seq(7, nrow(sim_theta4), by = 7),'cumInci'] - c(0, sim_theta4[seq(7, nrow(sim_theta4) - 7, by = 7),'cumInci'])
inci_theta5 = sim_theta5[seq(7, nrow(sim_theta5), by = 7),'cumInci'] - c(0, sim_theta5[seq(7, nrow(sim_theta5) - 7, by = 7),'cumInci'])


par(mfrow = c(2,1), mar = c(3,3,1,1), mgp = c(1.8,.5,0), cex = 1)
plot(inci_theta1, ylab = 'Weekly incidence', xlab = 'Weeks', ylim = c(0,30000), type = 'l', lwd = 2, col = "blue", main = 'Transmission rate of 1.5 and levels of intervention participation')
lines(inci_theta2, col='red')
lines(inci_theta3, col='green')
lines(inci_theta4, col='yellow')
lines(inci_theta5, col='purple')
legend('topright', c('0%','25%', '50%', '75%', '100%'),
       lty = c(1, 1, 1), pch = c(NA, NA, NA), col = c('blue','red', 'green', 'yellow', 'purple'), cex = .6, bty = 'n')

plot(sim_theta1[,'cumInci'], ylab = 'Cumulative incidence', xlab = 'Days', type = 'l', lwd = 2, col = "blue")
lines(sim_theta2[,'cumInci'], col='red')
lines(sim_theta3[,'cumInci'], col='green')
lines(sim_theta4[,'cumInci'], col='yellow')
lines(sim_theta5[,'cumInci'], col='purple')
legend('bottomright', c('0%','25%', '50%', '75%', '100%'),
       lty = c(1, 1, 1), pch = c(NA, NA, NA), col = c('blue','red', 'green', 'yellow', 'purple'), cex = .6, bty = 'n')
