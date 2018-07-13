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
bisection(min, max, line_table, eps)

############ FUNCTIONS ##############

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
# Code adopted from pseudocode at 
# described on Wikipedia page describing method
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
    print("You provided a function that doesn't CROSS 0 in the given domain")
    return(NULL)
  }

  c <- a
  while ((b-a) >= episilon){
    if (a_val == 0){
      print(paste("The value of the root is:", a))
      return(a)
    }
    if (b_val == 0){
      print(paste("The value of the root is:", b))
      return(b)
    }
    c <- (a+b)/2;
    c_val <- (a_val+b_val)/2
    if (c_val == 0){
      print(paste("The value of the root is:", c))
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
  
  print(paste("The value of the root is:", c))
  return(c)
}
