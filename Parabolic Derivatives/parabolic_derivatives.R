##############################################
# Author: Anooshree Sengupta
# Created on: 9/13/17
# Description: Code to find values from the function
#              e^(-x^2), then to compare the results
#              from the parabolic and numerical
#              derivatives
##############################################

####### RUN ME #######

# Random values for min, max, and step, set
# to values that follow preconditions and input
# definitions of create_table
# Feel free to define your own!
min = -5
max = 5
step = .5

#gauss <- gaussian_fun(min, max, step) # gaussian curve
gauss <- sinc_function(min, max, step)
para <- para_derivative(gauss)  # parabolic derivative
# dgauss <- gaussian_derivative(min(para$x_vals), max(para$x_vals), step) # numerical derivative
 dgauss <- sinc_derivative(min(para$x_vals), max(para$x_vals), step)

plot(gauss$x_vals, gauss$y_vals,
     ylim=c(-1,1), xlim=c(min,max), pch=17, col="blue", 
     ylab = "y", xlab="x") # gaussian function will be blue pts
points(para$x_vals, para$y_vals, 
       pch=16, col="green") # parametric derivative will be green pts
points(dgauss$x_vals, dgauss$y_vals, 
       pch=1, col="red") # numeric derivative will be red pts

# RMS error
print(paste("The RMS error is", signif(rms_error(dgauss, para))))

############ FUNCTIONS ##############

#######################################
# Creates a table of (x,y) pairs for the 
# gaussian function e^(-x^2)
# INPUT: min, the minimum x value
#        max, the maximum x value
#        step, the difference between consecutive
#              x values
# OUTPUT: a dataframe with x and y values
# PRECONDITIONS: function exists over all x values
# POSTCONDITIONS: data frame returned is sorted
#######################################

gaussian_fun <- function(min, max, step){
  e <- exp(1)
  x_vals <- seq(min, max, by = step)
  y_vals <- rep(0, length(seq(min, max, by = step)))
  for(i in 1:length(seq(min, max, by = step))){
    y_vals[i] <- signif(e^(-(x_vals[i]^2)))
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Creates a table of (x,y) pairs for the 
# function -2x*e^(-x^2)
# INPUT: min, the minimum x value
#        max, the maximum x value
#        step, the difference between consecutive
#              x values
# OUTPUT: a dataframe with x and y values,
#         with x values being the midpoints
#         between values in the provided dataframe
# PRECONDITIONS: function exists over all x values
#                provided dataframe is sorted
# POSTCONDITIONS: data frame returned is sorted
#######################################

gaussian_derivative <- function(min, max, step){
  e <- exp(1)
  x_vals <- seq(min, max, by = step)
  y_vals <- rep(0, length(seq(min, max, by = step)))
  for(i in 1:length(seq(min, max, by = step))){
    y_vals[i] <- signif(-2*x_vals[i]*e^(-(x_vals[i]^2)))
  }
  return(data.frame(x_vals, y_vals))
}


#######################################
# Creates a table of (x,y) pairs for the 
# approximated derivative of function represented in x
# using a parametric fit
# INPUT: x, a data frame containing sorted (x,y) points 
#           for the function e^(-x^2), or another function
# OUTPUT: a dataframe with x and y values from the three pt
#         derivative that follows the following formula:
#         y = ax+b, given the parametric fit of 
#         y = ax^2 +bx + c to the original function
# PRECONDITION: data in x is sorted
#               function contained in x is well-behaved
#               x has at least 3 (x,y) pairs
# POSTCONDITION: data frame returned is sorted
#######################################
para_derivative <- function(x){
  x_vals <- rep(0, nrow(x)-2)
  y_vals <- rep(0, nrow(x)-2)
  for(k in 1:(nrow(x)-2)){
    # Define the three points for 
    # a parametric fit
    x1 <- x[k,1]
    x2 <- x[k+1,1]
    x3 <- x[k+2,1]
    
    y1 <- x[k,2]
    y2 <- x[k+1,2]
    y3 <- x[k+2,2] 
    
    # Formulas for a and b found using a closed
    # form calculation in Mathematica
    a <- (x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1))/((x1-x2)*(x1-x3)*(x2-x3))
    b <- ((x3^2)*(y1-y2)+(x1^2)*(y2-y3)+(x2^2)*(y3-y1))/((x1-x2)*(x1-x3)*(x2-x3))
    
    x_vals[k] <- x2
    y_vals[k] <- 2*a*x2 + b
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Returns the RMS error between numerical
# and approximated derivatives
# INPUT: expected, a data frame containing sorted (x,y) points 
#                  for the numerical derivative of a function
#        observed, a data frame containing sorted (x,y) points
#                  for the approximated derivative of a function
# OUTPUT: the RMS error between the y values of expected 
#         and observed
# PRECONDITION: data in expected and observed is sorted
#               expected and observed have the same x_vals
#               expected and observed have the same number
#               of (x,y) points
#######################################
rms_error <- function(expected, observed){
  error_sum <- 0
  for(i in 1:nrow(expected)){
    error_sum = ((expected[i,2] - observed[i,2])^2) + error_sum
  }
  return(sqrt(error_sum/nrow(expected)))
}