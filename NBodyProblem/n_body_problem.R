##############################################
# Author: Anooshree Sengupta
# Created on: 1/2/17
# Description:
##############################################

n_body_movement <-  function(num_objects, initial_x, initial_y, initial_z,
                             initial_vx, initial_vy, initial_vz, delta_t,
                             end_time, masses){
  
  time_vector <- seq(0, end, by = delta_t)
  
  #############################
  # Create velocity data frames
  #############################
  vx <- data.frame(matrix(NA, nrow = length(time_vector), ncol = num_objects))
  rownames(velocities) <- time_vector
  colnames <- masses
  vx[1,] <- initial_vx
  
  vy <- data.frame(matrix(NA, nrow = length(time_vector), ncol = num_objects))
  rownames(velocities) <- time_vector
  colnames <- masses
  vy[1,] <- initial_vy
  
  vz <- data.frame(matrix(NA, nrow = length(time_vector), ncol = num_objects))
  rownames(velocities) <- time_vector
  colnames <- masses
  vz[1,] <- initial_vz
  
  ############################
  # Create position data frame
  ############################
  x <- data.frame(matrix(NA, nrow = length(time_vector), ncol = num_objects))
  rownames(velocities) <- time_vector
  colnames <- masses
  x[1,] <- initial_x
  
  y <- data.frame(matrix(NA, nrow = length(time_vector), ncol = num_objects))
  rownames(velocities) <- time_vector
  colnames <- masses
  y[1,] <- initial_y
  
  z <- data.frame(matrix(NA, nrow = length(time_vector), ncol = num_objects))
  rownames(velocities) <- time_vector
  colnames <- masses
  z[1,] <- initial_z
  
  for(i in 1:(length(time_vector)-1)){
    for(j in 1:(num_objects)){
      sum_vx = 0 
      sum_vy = 0
      sum_vz = 0
      for(k in 1:num_objects){
        sum_vx = sum_vx + difference_eqtn()
        sum_vy = sum_vy + difference_eqtn()
        sum_vz = sum_vz + difference_eqtn()
      }
      vx[i+1, j] = vx[i,j]-sum_vx
      vy[i+1, j] = vy[i,j]-sum_vy
      vz[i+1, j] = vz[i,j]-sum_vz
      
      x[i+1, j] = x[i,j]+(vx[i+1,j]*delta_t)
      y[i+1, j] = y[i,j]+(vy[i+1,j]*delta_t)
      z[i+1, j] = z[i,j]+(vz[i+1,j]*delta_t)
    }
  }
    
  

}