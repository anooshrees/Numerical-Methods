##############################################
# Author: Anooshree Sengupta
# Created on: 10/18/17
# Description: Functional fit for linear,
#              parabolic, and gaussian
##############################################

######## RUN ME #########

# Linear table test
line_table <- line_w_zero(-10, 10, 0.1, 3)
write.csv(line_table, file = "~/Desktop/School 2017-2018/Numerical Methods/line_table.csv")

# Gaussian table test
gaussian_table <- gaussian_fun(-10, 10, 0.1)
write.csv(gaussian_table, file = "~/Desktop/School 2017-2018/Numerical Methods/gaussian_table.csv")

# Actual test of function
start = 0 # beginning of domain
end = 41 # end of domain
step = 5 # step used for two-point derivative
q_vector = c(275000,12,6,5) # beginning guess for constants
file_name = "GeigerHisto.txt"

functional_fit(start, end, step, q_vector, file_name)

############# FUNCTIONS #######################

##############################################
# Function that computes the constants for a 
# functional fit of a given set of data points
# contained in a file by reducing the error 
# function across each of the constants. The 
# method can fit points to either a line, 
# two-degree polynomial, or a guassian function.
# It guesses the form based on the number of 
# constants provided
# INPUT: start, the beginning of the domain of points
#        end, the end of the domain of points
#        file_name, the file which contains the data 
#        step, the step size that will be used when finding
#              the two-point derivative when reducing error
#        q_vector, the vector of constants for the functional 
#                  form. Refer to functional_form for details
# OUTPUT: a vector with the constants for a given function form
#         that provides fit for set of data points
# PRECONDITION: points can be fit to one of three forms
#               step is not larger than domain
##############################################
functional_fit <- function(start, end, step, q_vector, file_name){
  # Read in the table containing the data values, assumes tsv
  df <- read.csv(file = paste("~/Desktop/School 2017-2018/Numerical Methods/", file_name, sep=""), sep="\t")
  df <- df[intersect(which(df[,1]<=end), which(df[,1]>=start)),]

  # Initial values for each of the vectors and variables used in error calculation
  sum_vector = rep(0, length(q_vector))
  error_vector = rep(0, length(q_vector))
  lambda = 1
  past_error = rep(1, length(q_vector))
  past_sum = rep(1, length(q_vector))
  count = 0

  # While loop continuing for as many counts as possible (before RStudio crashes)
  while(count <= 10000){
    # Calculate vectors with error and partial derivatives of each q
    for(q in 1:length(q_vector)){
      error = 0
      sum = 0
      for(i in 1:nrow(df)){
        sum = sum + (df[i,2]-functional_form(df[i,1], q_vector))*partial_derivative(q_vector, df[i,1], step)[q]
        error = error + (df[i,2]-functional_form(df[i,1], q_vector))^2
      }
      # Divide error by 2 to avoid the factor of 2 in derivative
      error = error/(2)
      error_vector[q] <- error
      sum_vector[q] <- sum
    }
    
    # Normalize the partial derivative vector
    sum_vector <- sum_vector/sqrt(sum(sum_vector^2)) 

    # Now go through both vectors and change the value of q
    curr_lambda = lambda # Hold previous lambda constant throughout loop
    for(q in 1:length(q_vector)){
      if(error_vector[q]>past_error[q]){
        lambda = curr_lambda/2
      }
      if(error_vector[q]<past_error[q]){
        lambda = curr_lambda*1.5
      }
      # Change the value in q_vector
      delta_q <- lambda*sum_vector[q]
      q_vector[q] <- q_vector[q] + delta_q
    }
    # Update vectors and count
    past_error = error_vector
    past_sum = sum_vector
    count = count+1
    
    # Print to keep user updated
    if(count%%250 == 0){
      print(paste("Currently on run", count))
      print(paste("The constants are", q_vector))
      print(paste("The error is", error_vector))
      print(paste("The learning factor is", lambda))
    }
  }
  
  # print out function approximation and q_vector
  if(length(q_vector) == 2){
    q_vector <- signif(q_vector)
    print(paste("Our fit is y = ", signif(q_vector[1]), "*x + (", signif(q_vector[2]), ")", sep = ""))
    print(paste("error is:", past_error))
    
    # graph the function and the points
    curve(functional_form(x, q_vector), from=start, to=end, , xlab="x", ylab="y")
    points(df[,1], df[,2], 
           pch=2, col="green")
    return(q_vector)
  }
  if(length(q_vector) == 3){
    q_vector <- signif(q_vector)
    print(paste("Our fit is y = ", signif(q_vector[1]), "*(x^2) + ", signif(q_vector[2]), "*x + (", signif(q_vector[3]), ")", sep = ""))
    print(paste("error is:", past_error))
    
    # graph the function and the points
    curve(functional_form(x, q_vector), from=start, to=end, , xlab="x", ylab="y")
    points(df[,1], df[,2], 
           pch=2, col="green")
    return(q_vector)
  }
  if(length(q_vector == 4)){
    q_vector <- signif(q_vector)
    print(paste("Our fit is y = ", q_vector[1], "*e^-(", -q_vector[2], "+x)^2)/(", q_vector[3], "^2) + ", q_vector[4], sep=""))
    print(paste("error is:", past_error))
    
    # graph the function and the points
    curve(functional_form(x, q_vector), from=start, to=end, , xlab="x", ylab="y")
    points(df[,1], df[,2], 
           pch=2, col="green")
    return(q_vector)
  }
  if(length(q_vector == 7)){
    q_vector <- signif(q_vector)
    print(paste("Our fit is y = ", q_vector[1], "*e^-(", q_vector[2], "-x)^2)/(", q_vector[3], "^2) + ", q_vector[4], "*e^-(", q_vector[5], "-x)^2)/(", q_vector[6], "^2) + ", q_vector[7], sep=""))
    print(paste("error is:", past_error))
    
    # graph the function and the points
    curve(functional_form(x, q_vector), from=start, to=end, , xlab="x", ylab="y")
    points(df[,1], df[,2], 
           pch=2, col="green")
    return(q_vector)
  }
}


##############################################
# Takes the partial derivative relative to each
# of the constants defined by the functional forms
# in functional_form, based on the length of the 
# q_vector provided by the user
# INPUT: point, the point at which derivative is calculated
#        step, the step used between points when calculating
#              the partial derivative
#        q_vector, the vector of constants for the functional 
#                  form. Refer to functional_form for details
# OUTPUT: a vector with the partial derivatives for each 
#         constant in q vector at the point, point
# PRECONDITION: points can be fit to one of three forms
#               step is not larger than domain
##############################################
partial_derivative <- function(q_vector, point, step){
  partial_vector <- rep(0, length(q_vector))
  for(q in 1:length(q_vector)){
    greater <- q_vector
    greater[q] <- greater[q]+step
    partial_vector[q] <- (functional_form(point, greater) - functional_form(point, q_vector))/step
  }
  return(partial_vector)
}
##############################################
# Function that provides the form of the 
# function the points are being fit to.
# The form is either a line, two-degree polynomial,
# or a gaussian function, depending on the number
# of constants provided:
# 2: line, 3: polynomial, 4: gaussian
# INPUT: q_vector, the vector of constants for the functional 
#                  form. 
#        x, the value at which the function is being calculated
# OUTPUT: a vector with the constants for a given function form
#         that provides fit for set of data points
# PRECONDITION: points can be fit to one of three forms
#               step is not larger than domain
##############################################
functional_form <- function(x, q_vector){
  # linear form
  if(length(q_vector) == 2){
    return(q_vector[1]*x + q_vector[2])
  }
  # parabolic form
  else if(length(q_vector) == 3){
    return(q_vector[1]*x^2 + q_vector[2]*x + q_vector[3])
  }
  # gaussian form
  else if(length(q_vector) == 4){
    return(q_vector[1]*e^(-(((x-q_vector[2])^2)/(q_vector[3]^2)))+q_vector[4])
  }
  # double gaussian form
  else{
    return(q_vector[1]*exp(-((x-q_vector[2])^2)/(q_vector[3])) + q_vector[4]^2*exp(-((x-q_vector[5])^2)/(q_vector[6]^2)) + q_vector[7])
  }
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