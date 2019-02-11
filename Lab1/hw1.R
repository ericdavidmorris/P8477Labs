## Importing Data

hw1 = read.csv(file.choose())

## Question 1 

nrow(hw1)

ncol(hw1)

hw1[2, 6]  

hw1[30, 5:14]

hw1[,5]

hw1[36, 'New.York.City']

summary(hw1$New.York.City)

## Plot for Question 2

plot(hw1$New.York.City, type = 'l', lty = 1, main = 'Influenza-like illness in NYC', 
     xaxt = 'n',  # this is telling R not to plot x-axis
     xlab = 'Date', ylab = 'ILI per 100,000')
axis(1, at = 1:nrow(hw1), labels = hw1$Date)  # add an axis


## Question 3 plot

par(cex = 1, mar = c(3,3,1,1), mgp = c(2,.5,0)) # set the parameters
matplot(hw1[,c('New.York.City','New.York')],
        xaxt = 'n',
        xlab = 'Date', ylab = 'ILI per 100,000',
        type = 'l', lty = 1,
        col = c('blue','red'))
axis(1, at = 1:nrow(hw1), labels = hw1$'Date') 
legend('topright',cex = .8,
       legend = c('New York City','New York State'),
       col = c('blue','red'),lty = c(1,1),
       bty = 'n')

## Question 4 boxplot

boxplot(t(hw1[1:52, 2:55]), 
        xlab = 'Week', ylab = "ILI of all locations")

## Question 5 

a = 2
b = 10
c = 3

x_1 = (-b + sqrt((b^2) - (4*a*c))) / (2*a)

x_2 = (-b - sqrt((b^2) - (4*a*c))) / (2*a)  

d = 5 
e = -6 
f = 1

x_3 = (-e + sqrt((e^2) - (4*d*f))) / (2*d)

x_4 = (-e - sqrt((e^2) - (4*d*f))) / (2*d)

# Question 6 Quad Function

Fn_sol_quadratic = function(a, b, c){
  
  pos_root = (-b + sqrt((b^2) - (4*a*c))) / (2*a)
  
  neg_root = (-b - sqrt((b^2) - (4*a*c))) / (2*a)
  
  print(pos_root)
  
  print(neg_root)
}

Fn_sol_quadratic(2, 10, 3)