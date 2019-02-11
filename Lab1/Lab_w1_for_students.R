# Script for R tutorial

# numbers
a = 1; 
b = 2;

# strings:
aa = "I'm a string"
bb = "Columbia"

# logical/Boolean:
flag1 = TRUE; # all in cap
flag2 = T; # or just T
flag3 = FALSE; # all in cap
flag4 = F # or just T

# vectors
x = c(1,2,3,4,5,6,7,8,9,10);
y = c('SIR','SEIR','SIRS')


# sequence
z1 = 1:100; 
z2 = seq(1,100, by = 2);


# matrix
X = matrix(1:12, 3, 4)
Y = matrix(1:12, 3, 4, byrow = T)


# array: when your dataset has >2 dimensions
arr = array(0,c(2,3,4))


# List
LL = list(model = 'SIR',
        Susceptible = seq(1e5,.9e5, by = -500),
        Infected = seq(1,1000, by = 500))

## Operations
x*5
10 + x*5
z = 100
x*z # not xz

## Inputs and outputs
# da=read.csv('~/Documents/Teaching/IDmodeling2017/Lab_w1_R101/Lab1data.csv')
da = read.csv(file.choose())


head(da) # by default: the first 6 rows
head(da,2); # the first 2 rows
tail(da); # by default: the last 6 rows
tail(da,1); # the last row

##########################################
## PLOT
##########################################
## scatter plot
# the basic
plot(x = da[,'time'], y = da[,'I'])
plot(x = da[,'time'], y = da[,'I'], pch = 20)
# label things
plot(da[,'time'],da[,'I'],
     xlab = 'Time',
     ylab = 'Number Infected',
     type = 'l', lty = 1,
     col = 'red', lwd = 2)

# make it prettier
par(cex = 1.2, mar = c(3,3,1,1), mgp = c(2,.5,0)) # set the parameters
plot(da[,'time'],da[,'I'],
     xlab = 'Time',
     ylab = 'Number Infected',
     type = 'l', lty = 1, lwd = 2)

par(cex = 1.2, mar = c(3,3,1,1), mgp = c(2,.5,0)) # set the parameters
plot(da[,'time'],da[,'S'],
     xlab = 'Time',
     ylab = 'Number Susceptible',
     type = 'l', lty = 1, lwd = 2,
     col = 'red')

# plot multiple lines together
par(cex = 1, mar = c(3,3,1,1), mgp = c(2,.5,0)) # set the parameters
matplot(da[,'time'],da[,c('S','I')],
        xlab = 'Time',ylab = '# people',
        type = 'l', lty = 1,
        col = c('blue','red'))
legend('topright',cex = .8,
       legend = c('Susceptible','Infected'),
       col = c('blue','red'),lty = c(1,1),
       bty = 'n')

par(cex = 1, mar = c(3,3,1,1), mgp = c(2,.5,0)) # set the parameters
matplot(da[,'time'],da[,c('S','I')],
        xlab = 'Time', ylab = '# people',
        type = 'l',col = 'black',lty = c(1,2))
legend('topright',cex = .8,
       legend = c('Susceptible','Infected'),
       lty = c(1,2),bty = 'n')


# histogram
par(cex = 1,mar = c(3,3,3,3),mgp = c(2,.5,0)) # set the parameters
hist(da[,'I'],xlab = 'Infected',
     main = 'Histogram of the Infected')

# boxplot
par(cex = 1.2,mar = c(3,3,1,1),mgp = c(2,.5,0)) # set the parameters
boxplot(X,xlab = 'Groups',ylab = 'Value')

# barplot: let's use a built-in dataset called VADeaths
# too see the dataset, use the command: View(VADeaths)
barplot(VADeaths,beside = T,legend = rownames(VADeaths))

# plot multiple plot together
par(mfrow = c(2,2),cex = .7,mar = c(3,3,2,1),mgp = c(2,.5,0)) # set the parameters
plot(da[,'time'],da[,'I'],main = 'Infected', xlab = 'Time',
     ylab = '#Infected',type = 'l',lty = 1,lwd = 2)
plot(da[,'time'],da[,'S'],main = 'Susceptible', 
     xlab = 'Time', ylab = '#Susceptible',
     type = 'l',lty = 1, lwd = 2,col = 'blue')
hist(da[,'I'],xlab = 'Infected',
     main = 'Hist of the Infected')
boxplot(X,main = 'boxplot',xlab = 'Groups',ylab = 'Value')


