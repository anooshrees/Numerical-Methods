##############################################
# Author: Anooshree Sengupta
# Created on: 9/19/17
# Description: Code to find zeroes of any 
#              given function using the 
#              bisector method
##############################################

############ RUN ME ################

bisection(-10, 10, eps)
newton_raphson(-10, 10, eps)

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
# Creates a table of (x,y) pairs for x^3
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
cubic <- function(min, max, step){
  x_vals <- seq(min, max, by = step)
  y_vals <- rep(0, length(seq(min, max, by = step)))
  
  for(i in 1:length(seq(min, max, by = step))){
    y_vals[i] <- x_vals[i]^3
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

bisection <- function(start, end, episilon){
  a <- start
  b <- end
  a_val <- test_fun(a)
  b_val <- test_fun(b)
  
  if (a_val == 0){
    print(paste("The value of the root is:", a))
    return(a)
  }
  if (b_val == 0){
    print(paste("The value of the root is:", b))
    return(b)
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
    c_val <- test_fun(c)
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
  
  print(paste(
    "Stopped due to step size smaller than episilon, stopped at the point listed below.",
    "If the y value is less than 2.220446e-16 (machine episilon), assume 0.",
    "Otherwise, please doublecheck that your function crosses 0 in the given interval by graphing."
  ))
  return(paste(c,",", c_val))
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
# pairs using the Newton-Raphson method, which
# follows the formula: x1 = x0 + (f(x0)/f'(x0))
# INPUT: start, the beginning of the domain 
#               containing a zero
#        end, the end of the dorm containing
#             a zero
#        eps, the machine constant episilon
#        iterations, the number of times the program
#                    will use the newton-raphson method
# OUTPUT: the x value of a zero in the given domain
# PRECONDITIONS: function exists over all x values
#######################################
newton_raphson <- function(start, end, eps, iterations){
  x <- start #set to the very first x value
  count <- 0
  # Use newton raphson method iterations times
  while(count<iterations){
    x <- (x+eps) - test_fun(x+eps)/test_der(x+eps)
    count <- count + 1
  }
  # Check the received point -- does it work?
  if(abs(test_fun(x)) <= eps && start<=x && x<=end){
    print("Found a zero (value smaller than episilon) at the following point:")
    return(paste(x, ", ", test_fun(x), sep = ""))
  }
  else{
    print("Could not find a root in given interval after the given number of iterations")
    return()
  }
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

bisection_min <- function(start, end, episilon, step){
  a <- start
  b <- end
  a_val <- three_pt_derivative(a, step)
  b_val <- three_pt_derivative(b, step)

  if (abs(a_val) <= eps && three_pt_second(a, step)>episilon){
    print(paste("A minimum is:", a))
    return(a)
  }
  if (abs(b_val) <= eps && three_pt_second(b, step)>episilon){
    print(paste("A minimum is:", b))
    return(b)
  }
  c <- a
  while ((b-a) >= episilon){
    if (abs(a_val) <= eps && three_pt_second(a, step)>episilon){
      print("Found a minimum at:")
      return(paste(a, ", ", test_fun(a), sep = ""))
    }
    if (abs(b_val) <= eps && three_pt_second(b, step)>episilon){
      print("Found a minimum at:")
      return(paste(b, ", ", test_fun(b), sep = ""))
    }
    c <- (a+b)/2;
    c_val <- three_pt_derivative(c, step)
    
    if (abs(c_val) <= eps && three_pt_second(c, step)>episilon){
      print("Found a minimum at:")
      return(paste(c, ", ", test_fun(c), sep = ""))
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
  return("Could not find a local minimum")
}


###############################
# Test function of cos(5x) + 2x
###############################
test_fun <- function(x){
  return(x*e^(x))
}

###############################
# Derivative of test function,
# 2 - 5sin(5x)
###############################
test_der <- function(x){
  return(2*x)
}

##########################
# 2nd Derivative of test
##########################
test_2_der <- function(x){
  return(2)
}

#######################################
# Creates a table of (x,y) pairs for the 
# approximated derivative of function represented in x
# INPUT: a point and step size 
# OUTPUT: a y value for the three point derivative approximation
#         based on the provided point, + step, and - step
# PRECONDITION: function exists at points and those
#               plus and minus step
#######################################
three_pt_derivative <- function(x_val, step){
  y_val <- signif((test_fun(x_val+step)-test_fun(x_val-step))/(2*step))
  return(y_val)
}

#######################################
# Creates a table of (x,y) pairs for the 
# approximated second derivative of function represented in x
# INPUT: a point and step size 
# OUTPUT: a y value for the three point derivative approximation
#         based on the provided point, + step, and - step
# PRECONDITION: function exists and is diff. at points and those
#               plus and minus step
#######################################
three_pt_second <- function(x_val, step){
  y_val <- signif((three_pt_derivative(x_val+step, step)-three_pt_derivative(x_val-step, step))/(2*step))
  return(y_val)
}
