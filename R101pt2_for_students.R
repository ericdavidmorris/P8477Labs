# Script for R tutorial (cont'd)
# Control Statements

##########################################
## If statements
##########################################
# if (something is TURE){
#   do plan-A
# } else {
#   do plan-B
# }
input=readline('please enter a number: ')
if(input>100){
  print('the number you enter is greater than 100')
} else {
  print('the number you enter is less or equal to 100')
}

p=runif(1) # randomly pick a number from 0 to 1
if (p<.33) {
  print(paste('p=',p),quote = F)
  print('small chance')
} else if (p<.66) {
  print(paste('p=',p),quote = F)
  print('medium chance')
} else {
  print(paste('p=',p),quote = F)
  print('great chance')
}

##########################################
## for loops
##########################################
# EX1. Use the looping variable as a counter but the task may be linked to the looping variable 
for(i in 1:10){
  print(i)
}

# EX2. Variables inside the loop depend on the looping variable and hence the result
X=1:10; Y=10:1 
for(i in 1:10){
  x=X[i]; # note: each time we grab a different x from the vector indexed by the looping variable
  y=Y[i]; # note: each time we grab a different y from the vector indexed by the looping variable
  z=x+y;  # each time z would based on the new x and y
  print(paste0('step', i,': x=',x,'; y=',y,'; x+y=',z))
}

# EX3. Use the looping variable as a counter for repeated task
# use the if-else code above to pick 10 random number (p) 
# and print the probabilities
for(i in 1:10){
  p=runif(1) # randomly pick a number from 0 to 1
  if (p<.33) {
    print(paste('p=',p,'small chance'),quote = F)
  } else if (p<.66) {
    print(paste('p=',p, 'medium chance'),quote = F)
  } else {
    print(paste('p=',p, 'great chance'),quote = F)
  }
}

## nested for loop: for-loop in a for-loop
for (i in 1:5){ 
  for (j in letters[1:3]){ 
    print(c(i,j))
  }
}

M=matrix(0,5,5) # initialize a matrix with 0
for (i in 1:5){ # go through rows
  for (j in 1:5){ # go through columns
    # fill the upper triangle with 1
    if (i<=j) M[i,j]=1
  }
}
M

##########################################
## while loop
##########################################
total=0
while (total<10){
  r=runif(1,0,3) # draw a random number, r, from [0,3]
  total=total+r # add r to the total
  print(total)
}




##########################################
## combination of diff control statements
##########################################
# Example:
# ask user to enter their age, and check if the entry is valid
# if <0, print error message 3 times (just to be annoying, haha)
# if >160, not likely, ask to enter again
# if between 0-160: move on

age = -1
while (age < 0 | age >160){
  age = as.numeric(readline('Please enter your age: '))
  if (age < 0) {
    for(i in 1:3){
      print ('Error!!! Your age should be >=0!!!')
    }
  } else if (age > 160){
    print ('Are you sure you are that old?')
    print ('Think again...')
  } else {
    print (paste('Thanks for entering. Your age is',age))
  }
}



