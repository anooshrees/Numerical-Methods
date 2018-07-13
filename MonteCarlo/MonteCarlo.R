##############################################
# Author: Anooshree Sengupta
# Created on: 9/13/17
# Description: Code using Monte Carlo method
#              to predict the value of pi
##############################################

##### RUN ME #####
# Constants used in Monte Carlo approximation
radius <- 1
n_runs <- 10^4
print_freq <- 10^2

# First, check out limits of runif() -- R's random number generator
random_num_hist(10^5)

# Get vector of difference between estimates and pi given constants
mc_error <- montecarlo(radius, n_runs, print_freq)

# Plot of error values; should converge to 0
plot(1:n_runs, mc_error,
     ylim=c(min(mc_error),max(mc_error)), xlim=c(0,n_runs), col="blue", 
     ylab = "Difference Error", xlab="Run Number", 
     main=paste("Error During", n_runs, "Runs"), cex = 0.5) 
# Line at zero to aid user in seeing convergence
abline(h=0)

#### FUNCTIONS ####

##############################################
# Function that uses Monte Carlo approximation
# to approximate the value of pi.
# INPUTS: radius, the size of the radius of the 
# 			     inscribed circle within a square that 
#			     is used for this method
#		 num_runs, the number of times a “dart” is thrown	
#		 print_out, the interval at which the 
#               estimate will be printed
# OUTPUT: vector containing difference between estimate
#         and pi for each run
##############################################
montecarlo <- function(radius, num_runs, print_out){
  count = 1
  y_vals <- rep(0, num_runs)
  n_circle <- 0
  
  while(count <= num_runs){
    x <- runif(1, 0, radius)
    y <- runif(1, 0, radius)
    if((x^2 + y^2) <= radius^2){
      n_circle = n_circle+1
    }
    estimate <- acos(-1) - (n_circle/count)*4
    y_vals[count] <- estimate
    if(count%%print_out == 0){
      print(paste("On run", count, "difference between pi and estimate is:", estimate))
      print(paste("n_circle is", n_circle))
    }
    count = count+1
  }
  return(y_vals)
}

##############################################
# Function that plots a histogram with the
# outputs of the random number generator runif()
# INPUTS: num_runs, the number of times a random
#                   number between 0 and 1 is 
#                   generated
# OUTPUT: a histogram of random values and frequency
##############################################
random_num_hist <- function(num_runs){
  vec <- rep(0, num_runs)
  for(i in 1:num_runs){
    vec[i] <- runif(1, 0, 1)
  }
  hist(vec, plot = TRUE, main = "Random Number Output", xlab = "Value")
}
