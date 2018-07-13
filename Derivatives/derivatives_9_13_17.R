##############################################
# Author: Anooshree Sengupta
# Created on: 9/4/17
# Description: Code to find values from gaussian and sinc
#              functions, then to compare the results
#              from the approximated and numerical
#              derivatives
##############################################

####### RUN ME ######

# Random values for min, max, and step, set
# to values that follow preconditions and input
# definitions of create_table
min = -5
max = 5
step = .5

# table of values for e^(-x^2)
 function_vals <- gaussian_fun(min, max, step)

# table of values for sinc
# function_vals <- sinc_function(min, max, step)

# table of values for approximated derivative 
# read function comments for description of each
# CHOOSE ONE; comment out the rest

# approximate_derivative <- derivative_table(function_vals)
# approximate_derivative <- three_pt_derivative(function_vals)
 approximate_derivative <- five_point_stencil(function_vals)

# # table of values for -2x*e^(-x^2); if used, 
# # comment out above numerical_derivative below
numerical_derivative <- gaussian_derivative(
  min(approximate_derivative$x_vals),
  max(approximate_derivative$x_vals),
  step)

# table of values for sin(x)/x; if used,
# comment out above numerical_derivative above
# numerical_derivative <- sinc_derivative(
#   min(approximate_derivative$x_vals),
#   max(approximate_derivative$x_vals),
#   step)
#   
# table of absolute differences between y values of derivatives
comparison <- compare_methods(approximate_derivative, numerical_derivative)

# Plot 
plot(function_vals$x_vals, function_vals$y_vals,
     ylim=c(-1,1), xlim=c(min,max), pch=1, col="blue", 
     ylab = "y", xlab="x")
points(approximate_derivative$x_vals, approximate_derivative$y_vals, 
       pch=2, col="green") # approximated derivative will be green pts
points(numerical_derivative$x_vals, numerical_derivative$y_vals,
       pch=3, col="purple")

# Standard deviation
print(paste("The standard deviation is", signif(sd(comparison$abs_difference))))

# RMS error
print(paste("The RMS error is", signif(rms_error(numerical_derivative, approximate_derivative))))

########### FUNCTIONS ################

#######################################
# Creates a table of (x,y) pairs for the 
# function e^(-x^2)
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
# function sin(x)/x
# INPUT: min, the minimum x value
#        max, the maximum x value
#        step, the difference between consecutive
#              x values
# OUTPUT: a dataframe with x and y values
# PRECONDITIONS: function exists over all x values
# POSTCONDITIONS: data frame returned is sorted
#######################################

sinc_function <- function(min, max, step){
  x_vals <- seq(min, max, by = step)
  y_vals <- rep(0, length(seq(min, max, by = step)))
  for(i in 1:length(seq(min, max, by = step))){
    if(x_vals[i] == 0){
      y_vals[i] <- 1
    }
    else{
      y_vals[i] <- sin(x_vals[i])/x_vals[i]
    }
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Creates a table of (x,y) pairs for the 
# function [cos(x)/x] - [sin(x)/x^2]
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

sinc_derivative <- function(min, max, step){
  x_vals <- seq(min, max, by = step)
  y_vals <- rep(0, length(seq(min, max, by = step)))
  for(i in 1:length(seq(min, max, by = step))){
    if(x_vals[i] == 0){
      y_vals[i] <- 0
    }
    else{
      y_vals[i] <- (cos(x_vals[i])/x_vals[i]) - (sin(x_vals[i])/(x_vals[i]^2))
    }
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Creates a table of (x,y) pairs for the 
# approximated derivative of e^(-x^2)
# INPUT: x, a data frame containing sorted (x,y) points 
#           for the function e^(-x^2)
# OUTPUT: a dataframe with x and y values,
#         with x values being the midpoints
#         between values in the provided dataframe,
#         and y values being the slope between
#         two values in the x column of the provided
#         dataframe
# PRECONDITION: data in x is sorted
#               function contained in x is well-behaved
# POSTCONDITION: data frame returned is sorted
#######################################

derivative_table <- function(x){
  x_vals <- rep(0, nrow(x)-1)
  y_vals <- rep(0, nrow(x)-1)
  for(i in 1:nrow(x)-1){
    x_vals[i] <- (x[i,1]+ x[i+1,1])/2
    y_vals[i] <- signif((x[i+1,2] - x[i,2])/(x[i+1,1] - x[i,1]))
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Creates a table of (x,y) pairs for the 
# approximated derivative of function represented in x
# INPUT: x, a data frame containing sorted (x,y) points 
#           for the function e^(-x^2), or another function
# OUTPUT: a dataframe with x and y values from the five,
#         point stencil that follows the following formula:
#         f'(x) = (-f(x+2h)+8f(x+h)-8f(x-h)+f(x-2h))/12h
# PRECONDITION: data in x is sorted
#               function contained in x is well-behaved
#               x has at least 5 (x,y) pairs
# POSTCONDITION: data frame returned is sorted
#######################################
five_point_stencil <- function(x){
  x_vals <- rep(0, nrow(x)-4)
  y_vals <- rep(0, nrow(x)-4)
  for(k in 1:(nrow(x)-4)){
    i <- k+2
    x_vals[k] <- x[i,1]
    y_vals[k] <- signif((-x[i+2, 2] + 8*x[i+1,2] - 8*x[i-1,2] + x[i-2,2])/(12*(x[i+1,1]-x[i,1])))
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Creates a table of (x,y) pairs for the 
# approximated derivative of function represented in x
# INPUT: x, a data frame containing sorted (x,y) points 
#           for the function e^(-x^2), or another function
# OUTPUT: a dataframe with x and y values from the three pt
#         derivative that follows the following formula:
#         f'(x) = (-f(x+2h)+8f(x+h)-8f(x-h)+f(x-2h))/12h
# PRECONDITION: data in x is sorted
#               function contained in x is well-behaved
#               x has at least 3 (x,y) pairs
# POSTCONDITION: data frame returned is sorted
#######################################
three_pt_derivative <- function(x){
  x_vals <- rep(0, nrow(x)-2)
  y_vals <- rep(0, nrow(x)-2)
  for(k in 1:(nrow(x)-2)){
    i <- k+1
    x_vals[k] <- x[i,1]
    y_vals[k] <- signif((x[i+1,2] - x[i-1,2])/(x[i+1,1] - x[i-1,1]))
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Compares the results from numerical and
# approximated derivatives
# INPUT: d1, a data frame containing sorted (x,y) points 
#             for the function -2x*e^(-x^2)
#        d2, a data frame containing sorted (x,y) points
#               for the approximated derivative of e^(-x^2)
# NOTE: positions of inputs can switch
# OUTPUT: a dataframe with x values and the abs. differences
#         between the y values of d1 and d2
# PRECONDITION: data in d1 and d2 is sorted
#               d1 and d2 have the same x_vals
#######################################
compare_methods <- function(d1, d2){
  results <- data.frame(d1$x_vals, abs(d1$y_vals-d2$y_vals))
  colnames(results) <- c("x_vals", "abs_difference")
  return(results)
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