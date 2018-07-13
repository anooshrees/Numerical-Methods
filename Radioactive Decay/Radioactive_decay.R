##############################################
# Author: Anooshree Sengupta
# Created on: 11/13/17
# Description: Code that provides numerical 
#              prediction of the solution for
#              the half life differential equation
#              N'(t) = -N(t)/t using the approximation
#              N(t + delta_t) = N(t)*[1-(delta_t/tau)]
#              over a user-provided domain and with
#              a user-designated delta_t for a given
#              material with a user-provided tau and
#              initial amount of user-provided N
##############################################

##### RUN ME ####

# Provide the constants for N(t)=N*e^(-t/tau) and step size
n_initial = 100
delta_t = 3
tau = 29 # tau for strontium 90, in years

# Provide end of domain -- it always starts at t=0
end = 500

# And run the function! It will automatically graph approximation
# and return the (x,y) points for future reference
points <- half_life(n_initial, delta_t, tau, end)

# Running RMS error function will return RMS error and add
# a red plot with close form values to the graph
rms_error(points, n_initial, tau)

#### FUNCTIONS ####

##############################################
# Function that provides a numerical approximation
# of the function values of N(t) = N*e^(-t/tau)
# using the estimate:
# N(t + delta_t) = N(t)*[1-(delta_t/tau)]
# INPUTS: n_initial, the value of N in the function
#		      delta_t, the step size between function approximations
#         tau, the tau constant in the function
#         start_time, the beginning of the domain for the function
#         end_time, the end of the domain
# OUTPUT: data frame containing (x,y) points with the 
#         approximation for the function with the 
#         user-provided constants. Points are in domain
#         and x values are delta_t apart.
#         Also plots the approximation.
# PRECONDITION: function exists on provided domain
#               delta_t is not larger than domain
##############################################

half_life <- function(n_initial, delta_t, tau, end_time){
  time_vector <- seq(0, end_time, by = delta_t)
  n_vector <- rep(0, length(time_vector))
  previous_n <- n_initial
  
  for(i in 1:length(n_vector)){
   n_vector[i] <- previous_n*(1-(delta_t/tau))
   previous_n <- n_vector[i]
  }
  
  plot(time_vector, n_vector,
       ylim=c(min(n_vector),max(n_vector)), xlim=c(0,end_time), col="blue", 
       ylab = "Material Remaining", xlab="Time", 
       main="Radioactive Decay for Given Material", cex = 0.5) 
  
  return(data.frame(time_vector, n_vector))
}

#######################################
# Returns the RMS error between numerical
# and approximated values of radioactive decay
# function with form N(t) = N*e^(-t/tau)
# INPUT: estimate, a data frame containing (x,y) points 
#                  for the estimated values of the function
#        n_initial, the value for the N constant in the 
#                   function
#        tau, the tau constant in the function
# OUTPUT: the RMS error between the values of estimate 
#         and close form function values
#         A plot of the close form values on the domain 
#         provided by estimate's x values 
# PRECONDITION: estimate is the estimate for the function
#               with the constants n_initial and tau
#               both estimate and close-form calculation exist
#               on domain provided by x values in first column
#               of estimate.
#######################################

rms_error <- function(estimate, n_initial, tau){
  time <- estimate[,1]
  actual <- rep(0,nrow(estimate))
  error_sum <- 0
  for(i in 1:nrow(estimate)){
    actual[i] <- (n_initial*exp(-time[i]/tau))
    error_sum = ((actual[i] - estimate[i,2])^2) + error_sum
  }
  
  # Add actual points of graph to existant graph
  points(time, actual, pch=5, col="red", cex = 0.5)
  
  print("The RMS error is:")
  return(sqrt(error_sum/nrow(estimate)))
}
