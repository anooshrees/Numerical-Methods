##############################################
# Author: Anooshree Sengupta
# Created on: 11/17/17
# Description: Code that provides numerical 
#              prediction of the angular displacement
#              theta for a simple pendulum using
#              the following approximations:
#              θ(t+Δt)=θ(t)+ω(t+Δt)*Δt
#              and ω(t+Δt) = ω(t)-(g/l)*sin(θ(t))Δt
#              over a user-provided domain and delta_t.
#
#              Code will also require initial values of 
#              the angular displacement θ and the angular 
#              velocity ω and the value of the gravitational 
#              constant g, and pendulum length l
#
#              For exported graphs, approx. is the blue
#              curve, and the close-form is red
##############################################

##### RUN ME ####

# Provide the constants for formula listed in header
theta = pi/5 # will work better for smaller angles; use radians!
omega = 0 
g = 9.8 # gravitational constant is 9.8 on earth
l = 3
delta_t = 0.01

# Provide end of domain -- it always starts at t=0
end = 12

# And run the function! It will automatically graph approximation
# and return the (x,y) points for future reference
approximation <- simple_pendulum(theta, omega, g, l, delta_t, end)

# Running RMS error function will return RMS error and add
# a red plot with close form values to the graph
rms_error(approximation, theta, g, l)

#### FUNCTIONS ####

##############################################
# Function that provides a numerical approximation
# of the function values of theta in θ'' + (g/l)sin(θ)=0,
# as the function predicts the angle of a simple pendulum
# using the estimate:
# θ(t+Δt)=θ(t)+ω(t+Δt)*Δt
# and ω(t+Δt) = ω(t)-(g/l)*sin(θ(t))Δt
# INPUTS: theta_initial, the starting angle of the pendulum
#         omega_initial, the starting angular velocity of the pendulum
#         gravity, the gravitational constant of the planet
#         length, the length of the pendulum's string
#		  delta_t, the step size between function approximations
#         end_time, the end of the domain, assuming start at 0
# OUTPUT: data frame containing (x,y) points with the 
#         approximation for theta with the 
#         user-provided constants. Points are in domain
#         and time values are delta_t apart.
#         Also plots the approximation.
# PRECONDITION: function exists on provided domain
#               delta_t is not larger than domain
##############################################

simple_pendulum <- function(theta_initial, omega_initial, gravity, length, delta_t, end_time){
  time_vector <- seq(0, end_time, by = delta_t)
  
  theta_vector <- rep(0, length(time_vector))
  theta_vector[1] <- theta_initial
  
  omega_vector <- rep(0, length)
  omega_vector[1] <- omega_initial
  
  for(i in 2:length(time_vector)){
    omega_vector[i] <- omega_vector[i-1]-(gravity/length)*sin(theta_vector[i-1])*delta_t
    theta_vector[i] <- theta_vector[i-1]+omega_vector[i]*delta_t
  }
  
  # pdf(file = "~/Desktop/pendulum.pdf")
  plot(time_vector, theta_vector,
       ylim=c(min(theta_vector),max(theta_vector)), xlim=c(0,end_time), col="blue", 
       ylab = "Angle (Radians)", xlab="Time", 
       main="Angle Over Time of Simple Pendulum", cex = 0.25) 
  
  return(data.frame(time_vector, theta_vector))
}

#######################################
# Returns the RMS error between numerical
# and approximated values of theta derived from
# function θ(t+Δt)=θ(t)+ω(t+Δt)*Δt
# INPUT: estimate, a data frame containing (x,y) points 
#                  for the estimated values of the function
#        theta_initial, the initial angle of the pendulum
#        g, the gravitational constant
#        l, the length of the pendulum
# OUTPUT: the RMS error between the values of estimate 
#         and close form function values
#         A plot of the close form values on the domain 
#         provided by estimate's x values 
# PRECONDITION: estimate is the estimate for the function
#               with the constants theta_inital, g, and l
#               both estimate and close-form calculation exist
#               on domain provided by time values in first column
#               of estimate.
#######################################

rms_error <- function(estimate, theta_initial, g, l){
  time <- estimate[,1]
  actual <- rep(0,nrow(estimate))
  error_sum <- 0
  for(i in 1:nrow(estimate)){
    actual[i] <- (theta_initial*cos(sqrt(g/l)*time[i]))
    error_sum = ((actual[i] - estimate[i,2])^2) + error_sum
  }
  
  # Add actual points of graph to existant graph
  points(time, actual, pch=5, col="red", cex = 0.25)
  
  print("The RMS error is:")
  return(sqrt(error_sum/nrow(estimate)))
}
