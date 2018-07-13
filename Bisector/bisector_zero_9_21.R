##############################################
# Author: Anooshree Sengupta
# Created on: 9/19/17
# Description: Code to find zeroes of any 
#              given function using the 
#              bisector method
##############################################

############ RUN ME ################

# Random values for min, max, step, intersection
# set to values that follow preconditions and input
# definitions of functions below 

# Feel free to define your own!
min = -10
max = 10
step = .1
intersection = 4

line_table <- line_w_zero(min, max, step, intersection)
gauss <- gaussian_fun(min, max, step) # gaussian curve

bisection(min, max, gauss, eps)

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
# Creates a table of (x,y) pairs for a line 
# with slope 1 and x-intercept of intersection
# formula: y = x-intersection
# INPUT: min, the minimum x value
#        max, the maximum x value
#        step, the difference between consecutive
#              x values
#        intersection, where the function intersects
#                      the x-axis 
# OUTPUT: a dataframe with x and y values
# PRECONDITIONS: function exists over all x values
# POSTCONDITIONS: data frame returned is sorted
#######################################

line_w_zero <- function(min, max, step, intersection){
  x_vals <- seq(min, max, by = step)
  y_vals <- rep(0, length(seq(min, max, by = step)))
  for(i in 1:length(seq(min, max, by = step))){
    y_vals[i] <- x_vals[i]-intersection
  }
  return(data.frame(x_vals, y_vals))
}

#######################################
# Find the root of a function from a table of (x,y) 
# pairs using the bisection method
# Code adopted from the following image:
# https://en.wikipedia.org/wiki/Bisection_method#/media/File:Bisection_method.svg
# INPUT: start, the beginning of the domain 
#               containing a zero
#        end, the end of the dorm containing
#             a zero
#        x, the dataframe containing (x,y) pairs
#           for a given function on [a,b]
# OUTPUT: the x value of a zero in the given domain
# PRECONDITIONS: function exists over all x values
#######################################

bisection <- function(start, end, x, episilon){
  a <- start
  b <- end
  a_val <- x[which(x$x_vals == a), 2]
  b_val <- x[which(x$x_vals == b), 2]
  
  if (a_val == 0){
    print(paste("The value of the root is:", a))
    return(a)
  }
  if (b_val == 0){
    print(paste("The value of the root is:", b))
    return(b)
  }
  if (a_val * b_val >= 0){
    print("You provided a function that doesn't CROSS 0 in the given domain, returning NULL")
    return(NULL)
  }

  c <- a
  while ((b-a) >= episilon){
    if (a_val == 0){
      print("The value of the root is:")
      return(a)
    }
    if (b_val == 0){
      print("The value of the root is:")
      return(b)
    }
    c <- (a+b)/2;
    if(c %in% x[,1]){
      c_val <- x[which(x$x_vals == c), 2]
    }
    else{
      c_val <- (a_val+b_val)/2
    }
    
    if (c_val == 0){
      print("The value of the root is:")
      return(c)
    }
    else{
      if(c_val*a_val < 0){
        b <- c
        b_val <- c_val
      }
      else{
        a <- c
        a_val <- c_val
      }
    }
  }
  
  print("Stopped due to step size smaller than episilon, stopped at:")
  return(c)
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
# Find the root of a function from a table of (x,y) 
# pairs using the Newton-Raphson method
# INPUT: start, the beginning of the domain 
#               containing a zero
#        end, the end of the dorm containing
#             a zero
#        x, the dataframe containing (x,y) pairs
#           for a given function on [a,b]
# OUTPUT: the x value of a zero in the given domain
# PRECONDITIONS: function exists over all x values
#######################################
newton_raphson <- function(x){
  x1 <- x[1,1] #set to the very first x value
  num_rows <- nrow(x)
  
  while(x1 != 0){
    x1 <- x1 - 
  }
}


