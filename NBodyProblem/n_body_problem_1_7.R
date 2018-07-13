##############################################
# Author: Anooshree Sengupta
# Created on: 1/7/17
# Description: Models the movement of objects over
#              a user-specified amount of time based
#              on Newton's Universal Law of Gravition,
#              as represented by the following differential equations
#              taken from the course handout on Athena:
#              vj(t+Δt) = vj(t)-∑((mi*G)/|rj(t)-ri(t))∣^3)*(rj(t)-ri(t))Δt
#              rj(t+Δt) = rj(t)+vj(t+Δt)Δt
#
#              These vector and position values are 3D vectors, so the 
#              code calculates the x, y, and z position of each object
#              and simultaneously graphs the changing positions while
#              reporting their values to the user.
#
#              This code can be used to understand interactions between 
#              objects, but has greater applications, such as modeling
#              the movement of stars in a galaxy. If objects collide,
#              the code merges them into one object with their combined mass, and
#              the velocity required by the conversation of momentum.
#
#              The RUN ME sections offers suggested simulations from the handout,
#              but the user can change those values and must provide the time
#              domain and intervals, initials velocities, and initial positions.
#
#              For specifics on graphing and exporting graphs to pdf files,
#              the user should go directly to the bottom of the nbody_movement
#              function for specifics on handling each simulation.
##############################################

###### RUN ME ######
library(scatterplot3d) # used for plotting objects

# define the universal gravitational constant G
G =  6.754*10^(-11)

################################################################
# Run first test case: two 100 kg objects placed 10 meters apart
# Note: numerator values are in meters; operations following them
#       convert them to solar masses, parsecs, and km/s.
#       This is to make it easier for the user to change the values
#       by meters for this small-scale interaction
################################################################
num_objects <- 2

# 10 meters apart
initial_x <- c(0,0) # initial x positions (parsecs)
initial_y <- c(0,10/(3.086*10^16)) # initial y positions (parsecs)
initial_z <- c(0,0) # initial z positions (parsecs)

# Begin at rest
initial_vx <- c(0,0) # initial x velocities (km/s)
initial_vy <- c(0,0) # initial y velocities (km/s)
initial_vz <- c(0,0) # initial z velocities (km/s)

# delta_t can be easily changed by user
delta_t <- 100 # 100 seconds
end_time <- 345600 # about four days in seconds, collides at ~3.4 days

# both objects have a mass of 100 kg
masses <- c(100, 100)/(1.989*10^30) # masses of the objects (solar masses)

# both objects have a radius of 1 m
radii <- c(1,1)/(3.086*10^16) # radii of the objects (parsecs)

n_body_movement(num_objects, initial_x, initial_y, initial_z,
                initial_vx, initial_vy, initial_vz, delta_t,
                end_time, masses, radii)

################################################################
# Run second test case: the earth rotating around the sun
################################################################
num_objects <- 2

#  0.000004848 parsecs apart (mean distance), earth is first object and sun is the second
initial_x <- c(0.000004848, 0) # initial x positions 
initial_y <- c(0,0) # initial y positions 
initial_z <- c(0,0) # initial z positions 

# Begin with earth moving at 30 km/s and sun at rest
# NOTE: Realistically, the sun moves at 200 km/s around the galaxy.
#       However, since the earth is orbiting the sun and
#       thus also takes on those 200 km/s, we can eliminate it
#       because relative to the earth, the sun is fairly still
initial_vx <- c(0,0) # initial x velocities 
initial_vy <- c(30,0) # initial y velocities 
initial_vz <- c(0,0) # initial z velocities 

# delta_t can be easily changed by user
delta_t <- 86400 # 1 day
end_time <- 31536000 # approximately one year in seconds

# mass of earth and sun found through Google (same with velocities and masses)
masses <- c(0.000003003, 1) # masses of the objects (solar masses)

# radius of earth and sun found via Google
radii <- c(2.065*10^-10, 2.25461*10^-8) # radii of the objects (parsecs)

n_body_movement(num_objects, initial_x, initial_y, initial_z,
                initial_vx, initial_vy, initial_vz, delta_t,
                end_time, masses, radii)

################################################################
# Run third test case: build a galaxy
################################################################

# maximum number of galaxies this code can handle (closest mult. of 10)
# For pure observation and numerical output, system handles up to 2000
# For graphing, use 100 for a clear, understandable graph
num_objects <- 100

initial_x <- rep(0, num_objects)
initial_y <- rep(0, num_objects)
initial_z <- rep(0, num_objects)

initial_vx <- rep(0, num_objects)
initial_vy <- rep(0, num_objects)
initial_vz <- rep(0, num_objects)

masses <- rep(0, num_objects)
radii <- rep(0, num_objects)

# set very first object at the origin; makes random generation
# of positions of following objects easier

# at origin
initial_x[1] <- 0
initial_y[1] <- 0
initial_z[1] <- 0 

# a multiplier that controls how large
# the velocity and position of an object
# can become. 0.5773503 is approximately
# the cube root of 1/3, which means the 
# magnitude of velocity won't exceed 200 km/s
# and position won't exceed 30K parsecs away
# from the next object.
mult = 0.5773503

# give a random velocity
initial_vx[1] <- runif(1, -mult, mult)*200
initial_vy[1] <- runif(1, -mult, mult)*200
initial_vz[1] <- runif(1, -mult, mult)*200

masses[1] <- runif(1, 0.025, 10)

for(particle in 2:num_objects){
  # position the stars between  0-30K parsecs away from each other 
  # (30K parsecs is the size of milky way)
  initial_x[particle] <- initial_x[particle-1]+runif(1, -mult, mult)*30000
  initial_y[particle] <- initial_y[particle-1]+runif(1, -mult, mult)*30000
  initial_z[particle] <- initial_z[particle-1]+runif(1, -mult, mult)*30000
  
  # give a random velocity (based on 200 km/s because that is the vel of the sun)
  initial_vx[particle] <- runif(1, -mult, mult)*200
  initial_vy[particle] <- runif(1, -mult, mult)*200
  initial_vz[particle] <- runif(1, -mult, mult)*200
  
  # assign a random mass: from a brown dwarf or to red giant
  masses[particle] <- runif(1, 0.025, 10)
  # radius creates stars of sizes within a magnitude of ten from
  # the sun's radius. As always, this random number generation
  # can be changed; this is simply for demonstration purposes
  radii[particle] <- runif(1, 1, 4)*(10^runif(1, -9, -7))
}

# delta_t can be easily changed by user
delta_t <- 630720000000000 # 20 million years in seconds
end_time <- 6.3072*10^15 # 200 million years in seconds

n_body_movement(num_objects, initial_x, initial_y, initial_z,
                initial_vx, initial_vy, initial_vz, delta_t,
                end_time, masses, radii)

########################################
# User Input Option
# Values can be inputted manually or 
# read from a file; the choice is yours!
########################################
num_objects <- 

initial_x <- c() # initial x positions (parsecs)
initial_y <- c() # initial y positions (parsecs)
initial_z <- c() # initial z positions (parsecs)

# Begin at rest
initial_vx <- c() # initial x velocities (km/s)
initial_vy <- c() # initial y velocities (km/s)
initial_vz <- c() # initial z velocities (km/s)

# delta_t can be easily changed by user
delta_t <-  # seconds
end_time <-  # seconds

# both objects have a mass of 100 kg
masses <- c() # masses of the objects (solar masses)

# both objects have a radius of 1 m
radii <- c() # radii of the objects (parsecs)

n_body_movement(num_objects, initial_x, initial_y, initial_z,
                initial_vx, initial_vy, initial_vz, delta_t,
                end_time, masses, radii)

###### FUNCTIONS ######

################################################################
# Calculates numerical prediction of the velocity and 
# position for an object using the following differential equations 
# derived from Newton's Universal Law of Gravitation over a user-
# designated period of time and graphs the objects in 3D
# (taken from the "n-body problem" document on Schoology):
# 1. vj(t+Δt) = vj(t)-∑((mi*G)/|rj(t)-ri(t))∣^3)*(rj(t)-ri(t))Δt
# 2. rj(t+Δt) = rj(t)+vj(t+Δt)Δt
#
# INPUTS: num_objects, the number of objects
#         initial_x, a vector containing the inital x positions of the objects
#         initial_y, a vector containing the inital y positions of the objects
#         initial_z, a vector containing the inital z positions of the objects
#
#         initial_vx, a vector containing the inital x velocities of the objects
#         initial_vy, a vector containing the inital y velocities of the objects
#         initial_vz, a vector containing the inital z velocities of the objects
#         
#         delta_t, the user-provided time interval between calculations
#         end_time, the end time of the simulation
#         masses, a vector containing the masses of the objects
#         radii, a vector containing the radii of the objects
#
# OUTPUT: the final positions and velocities of the objects
#         prints out the positions throughout the code
#         continually graphs the motion of the objects in 3D
#
# PRECONDITION: function exists on provided domain
#               delta_t is not larger than domain
################################################################
n_body_movement <-  function(num_objects, initial_x, initial_y, initial_z,
                             initial_vx, initial_vy, initial_vz, delta_t,
                             end_time, masses, radii){
  
  time_vector <- seq(0, end_time, by = delta_t)
  
  #####################################################
  # Convert to MKS from solar masses, parsecs, and km/s
  #####################################################
  # solar masses to kg
  masses <- masses*1.989*10^30
  # parsecs to m
  initial_x <- initial_x*3.086*10^16
  initial_y <- initial_y*3.086*10^16
  initial_z <- initial_z*3.086*10^16
  radii <- radii*3.086*10^16
  # km/s to m/s
  initial_vx <- initial_vx*1000
  initial_vy <- initial_vy*1000
  initial_vz <- initial_vz*1000
  
  #############################
  # Create velocity vectors
  #############################
  vx <- initial_vx
  vy <- initial_vy
  vz <- initial_vz
  
  ############################
  # Create position vectors
  ############################
  x <- initial_x
  y <- initial_y
  z <- initial_z
  
  ##################################
  # Perform calculations and graph
  ##################################
  for(i in 1:length(time_vector)){
    
    #################################
    # Check for any collisions first 
    #################################
    if(num_objects > 1){
      for(obj in 1:(num_objects-1)){
        for(other_obj in (obj+1):(num_objects)){
          # find distance at which objects would collide
          dist = radii[obj]+radii[other_obj]
          # then check whether or not they are 
          # touching -- if so, merge into one object
          if(abs(x[obj]-x[other_obj]) <= dist+eps && 
            abs(y[obj]-y[other_obj]) <= dist+eps &&
            abs(z[obj]-z[other_obj]) <= dist+eps){
            # save masses
            m1 <- masses[obj]
            m2 <- masses[other_obj]
            # alert the user
            print("*****************************************")
            print(paste("Objects of mass ", m1, " and mass ", m2, " have collided"))
            print("*****************************************")
            # change velocities
            vx[obj] <- (m1*vx[obj] + m2*vx[other_obj])/(m1+m2)
            vx[other_obj] <- 0
            vy[obj] <- (m1*vy[obj] + m2*vy[other_obj])/(m1+m2)
            vy[other_obj] <- 0
            vz[obj] <- (m1*vz[obj] + m2*vz[other_obj])/(m1+m2)
            vz[other_obj] <- 0
            # change masses
            masses[obj] <- m1+m2
            masses <- masses[-other_obj]
            # change position and update num_objects
            x <- x[-other_obj]
            y <- y[-other_obj]
            z <- z[-other_obj]
            num_objects <- num_objects - 1
          }
        }
      }
    }
    
    # store current values in vectors to prevent the code 
    # from using updated positions in the velocity calculation
    previous_x <- x
    previous_y <- y
    previous_z <- z
    
    ##################################################
    # Iterate through objects and perform calculations
    ##################################################
    for(j in 1:(num_objects)){
      # These sums will be added to previous velocities in 
      # order to complete difference equations
      sum_vx = 0 
      sum_vy = 0
      sum_vz = 0
      
      # another for loop in order to complete summation in 
      # diff. eq.
      for(k in 1:num_objects){
        if(k != j){
          # denominator in summation portion is
          # the magnitude of the difference between
          # locations, cubed
          constant =
            sqrt(
              (previous_x[j]-previous_x[k])^2 +
                (previous_y[j]-previous_y[k])^2 +
                (previous_z[j]-previous_z[k])^2)^3
          
          sum_vx = sum_vx + difference_eqtn(j, k, previous_x, masses, delta_t, constant)
          sum_vy = sum_vy + difference_eqtn(j, k, previous_y, masses, delta_t, constant)
          sum_vz = sum_vz + difference_eqtn(j, k, previous_z, masses, delta_t, constant)
        }
       }
      
      # update the velocities of the given object by subtracting
      # the sum calculated earlier
      vx[j] = vx[j]-sum_vx
      vy[j] = vy[j]-sum_vy
      vz[j] = vz[j]-sum_vz
      
      # update the positions using the velocities from above
      x[j] = x[j]+(vx[j]*delta_t)
      y[j] = y[j]+(vy[j]*delta_t)
      z[j] = z[j]+(vz[j]*delta_t)
      
    }
    # continually print out the positions of each of the object 
    # Recommend commenting this out for galaxy simulation
    print(paste("At time:", time_vector[i]))
    for(object in 1:(num_objects)){
      print(paste("Object of mass ", masses[object], " is at (", x[object], ", ", y[object], ", ", z[object], ")", sep = ""))
    }
    
    # plot the positions of each of the objects in order to update user
    # every number or iterations (too often and graphs don't show up on
    # the user viewing panel)
    
    # The x, y, and z limits for the plots are "magic numbers" because
    # R makes it difficult to pass parameters to those values, and giving
    # those values to the user would not be effective as they cannot predict
    # gravitational movement for large objects to that large of a degree.
    # That being said, the recommended plots for the three simulations
    # are included and labelled below; comment out the ones you don't need.
    
    # two 100 kg objects
    # pdf and dev.off() line should be commented out if user just
    # wants graph in the RStudio viewer
    if(i%%1000 == 0){
      pdf(file = paste("~/Desktop/two100kg_time_", time_vector[i], "seconds.pdf", sep = ""),width=7,height=5)
      scatterplot3d(x, y, z,
                    xlim = c(-4,4), ylim = c(-1,11), zlim = c(-1,1))
      dev.off()
    }

    # earth and sun, creates a pdf file with position once every week
    # pdf and dev.off() line should be commented out if user just
    # wants graph in the RStudio viewer
    if(i%%(7) == 0){
      pdf(file = paste("~/Desktop/earthandsun_week_", i/7, ".pdf", sep = ""),width=7,height=5)
      scatterplot3d(x, y, z,
                    xlab = "X Position (m)", ylab = "Y Position (m)", zlab = "Z Position (m)",
                    xlim = c(-2.5*10^11, 2.5*10^11), ylim = c(-2.5*10^11, 2.5*10^11), zlim = c(-1,1))
      dev.off()
    }
    
    # galaxy exports to pdf for better understanding, but pdf lines can
    # be commented out for just scatterplot in RStudio
    pdf(file = paste("~/Desktop/galaxy100", i, ".pdf", sep = ""),width=7,height=5)
      scatterplot3d(x, y, z,
                    xlab = "X Position (m)", ylab = "Y Position (m)", zlab = "Z Position (m)",
                    xlim = c(min(x), max(x)), ylim = c(min(y), max(y)), zlim = c(min(z), max(z)))
    dev.off()
  }
}


################################################################
# Calculates one iteration of the sum that is subtracted from the   
# previous velocity in the difference equation to predict the future 
# velocity of the object
# INPUTS: current, the index of the object whose velocity is being calculated
#         other, the index of another object in the simulation
#         dimension, a vector of positions for one dimension (x,y, or z) 
#         masses, a vector containing hte masses of each object
#         delta_t, a user-designated time-step
#         constant, a calculated value that serves as the denominator
# OUTPUT: One calculation for the summation portion of the diff. eqt.
# PRECONDITION: The "current" and "other" objects are not the same
################################################################
difference_eqtn <- function(current, other, dimension, masses, delta_t, constant){
  # return part of summation if not same object
  return((masses[other]*G)*(dimension[current]-dimension[other])*delta_t/constant)
}
