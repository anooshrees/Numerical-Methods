##############################################
# Author: Anooshree Sengupta
# Created on: 11/30/17
# Description: Code that provides numerical 
#              prediction of the velocity and 
#              position for a projectile spherical object using
#              the following approximations using
#              the Euler-Cromer approach (derived in document):
#              vx(t+Δt)=vx(t)−(1/m)*2π(R^2)ρv(t)vx(t)Δt,
#              x(t+Δt)=x(t)+vx(t+Δt)Δt,
#              vy(t+Δt)=vy(t)−(1/m)*(2π(R^2)ρv(t)vy(t)+g)Δt,
#              y(t+Δt)=y(t)+vy(t+Δt)Δt 
#              over a user-provided domain and delta_t.
#
#              Code will also require initial values of 
#              the mass m, radius R of the spherical object,
#              initial velocity v, and initial angle θ
##############################################

###### RUN ME ######
# define constants used in all models
theta = pi/4 # must be in radians
radius = 0.155 # 15.5 centimeters, from handout
mass = 43.2 # 43.2 kgs, from handout
initial_v = 827 # 827 m/s, from handout

# decide on interval size and domain (assumed start at 0)
delta_t = 0.0005
end = 50

# set the initial position and density
initial_rho = 1.225 # 1.225 kg/m^3, from handout
initial_x = 0 # m, 0 is recommended
initial_y = 0 # m

###########################
# Constant air density model
###########################
# uses all of the predefined constants, designates density type as
# "Constant," and sets all of the Adiabatic-specific constants to 0
constant_model <- projectile_motion(theta, radius, mass, delta_t,
                                    end, initial_rho, initial_y, 
                                    initial_x, initial_v, "Constant",
                                    0, 0, 0)
# Create a new plot with the position of the constant model
plot(constant_model$x, constant_model$y,
     ylim=c(min(constant_model$y),max(constant_model$y)), xlim=c(0,max(constant_model$x)), col="blue", 
     ylab = "Y Position", xlab="X Position", 
     main="Position of a Projectile Object", cex = 0.25) 
# OR (must choose one) add points to existing plot (there must 
# be an existing plot)
points(constant_model$x, constant_model$y, pch=5, col="blue", cex = 0.25)

#####################
# Zero air density model
#####################
# Same as "Constant" model, but sets initial_rho to 0 (assumes no
# air density)
zero_model <- projectile_motion(theta, radius, mass, delta_t,
                                    end, 0, initial_y, 
                                    initial_x, initial_v, "Constant",
                                    0, 0, 0)
# Create a new plot
plot(zero_model$x, zero_model$y,
     ylim=c(min(zero_model$y),max(zero_model$y)), xlim=c(0,max(zero_model$x)), col="red", 
     ylab = "Y Position", xlab="X Position", 
     main="Position of a Projectile Object", cex = 0.25) 
# OR (must choose one) add points to existing plot (there must 
# be an existing plot)
points(zero_model$x, zero_model$y, pch=5, col="red", cex = 0.25)

###########################
# Adiabatic Density Profile
###########################
# using calculation: ρ=ρ0(1−(ay/T0))^α, 
# so need to set constants
a = 6.5*(10^-3) # K/m, from handout
alpha = 2.5 # from handout
T0 = 291.15 # K, from handout

# uses all of the predefined constants, designates density type as
# "Adiabatic," and uses all of the Adiabatic-specific constants 
adiabatic_model <- projectile_motion(theta, radius, mass, delta_t,
                                    end, initial_rho, initial_y, 
                                    initial_x, initial_v, "Adiabatic",
                                    a, alpha, T0)

plot(adiabatic_model$x, adiabatic_model$y,
     ylim=c(min(adiabatic_model$y),max(adiabatic_model$y)), xlim=c(0,max(adiabatic_model$x)), col="purple", 
     ylab = "Y Position", xlab="X Position", 
     main="Position of a Projectile Object", cex = 0.25) 
# OR (must choose one) add points to existing plot (there must 
# be an existing plot)
points(adiabatic_model$x, adiabatic_model$y, pch=5, col="purple", cex = 0.25)

############################
# Isothermal Density Profile
############################
# using calculation: ρ=ρ0*e^(−y/y0)
initial_y = 10^4 # m, from handout

# uses constants except initial_y, designates density type as
# "Isothermal," and sets all of the Adiabatic-specific constants to 0
isothermal_model <- projectile_motion(theta, radius, mass, delta_t,
                                    end, initial_rho, initial_y, 
                                    initial_x, initial_v, "Isothermal",
                                    0, 0, 0)

plot(isothermal_model$x, isothermal_model$y,
     ylim=c(min(isothermal_model$y),max(isothermal_model$y)), xlim=c(0,max(isothermal_model$x)), col="green", 
     ylab = "Y Position", xlab="X Position", 
     main="Position of a Projectile Object", cex = 0.25) 
# OR (must choose one) add points to existing plot (there must 
# be an existing plot)
points(isothermal_model$x, isothermal_model$y, pch=5, col="green", cex = 0.25)

###### FUNCTIONS ######

################################################################
# Calculates numerical prediction of the velocity and 
# position for a projectile spherical object using
# the following approximations using the Euler-Cromer approach 
# (derived in the "projectile motion" document on Schoology):
# 1. vx(t+Δt)=vx(t)−(1/m)*2π(R^2)ρv(t)vx(t)Δt,
# 2. x(t+Δt)=x(t)+vx(t+Δt)Δt,
# 3. vy(t+Δt)=vy(t)−(1/m)*(2π(R^2)ρv(t)vy(t)+g)Δt,
# 4. y(t+Δt)=y(t)+vy(t+Δt)Δt 
#
# INPUTS: theta, the firing angle of spherical object
#         radius, the radius of the spherical object
#         mass, the mass of the object
#         delta_t, the step size between function approx.
#         end, the end time of the motion (assumed start at 0)
#         initial_rho, the initial air density
#         initial_y, the initial y position
#         initial_x, the initial x position
#         type_density, the type of density model used
#         a, a constant used in the Adiabatic density model
#         alpha, another constant in the Adiabatic model
#         initial_temp, initial temperature, used in Adiabatic
#
# OUTPUT: a dataframe structure containing the x position, y position,
#         x velocity, and y velocity for each time (designated by delta_t
#         and end)
#         dataframe can be used in RUN ME section for graphing motion
#         and position values
#
# PRECONDITION: function exists on provided domain
#               delta_t is not larger than domain
################################################################

projectile_motion <-  function(theta, radius, mass, delta_t, end, 
                                 initial_rho, initial_y, initial_x, 
                                 initial_v, type_density, a, alpha,
                                 initial_temp){
  time_vector <- seq(0, end, by = delta_t)
  
  # Velocity vectors
  vx <- rep(0, length(time_vector))
  vx[1] <- initial_v*cos(theta)
  vy <- rep(0, length(time_vector))
  vy[1] <- initial_v*sin(theta)
  
  # Position vectors
  x <- rep(0, length(time_vector))
  x[1] <- initial_x
  y <- rep(0, length(time_vector))
  y[1] <- initial_y
  
  for(i in 1:(length(time_vector)-1)){
    past_v <- sqrt(vx[i]^2 + vy[i]^2)
    
    # For Constant and Zero options:
    rho <- initial_rho
    
    # Change rho if a function of y position
    if(type_density == "Adiabatic"){
      rho <- initial_rho*(1-(a*y[i]/initial_temp))^alpha
    }
    if(type_density == "Isothermal"){
      rho <- initial_rho*exp(-y[i]/y[1])
    }
    
    # approximate future velocity
    vx[i+1] <- vx[i]-(1/mass)*2*pi*(radius^2)*rho*past_v*vx[i]*delta_t
    vy[i+1] <- vy[i]-((1/mass)*2*pi*(radius^2)*rho*past_v*vy[i]+g)*delta_t
    
    # approximate future position
    x[i+1] <- x[i]+vx[i+1]*delta_t
    y[i+1] <- y[i]+vy[i+1]*delta_t
  }
  
  # returns a dataframe of the results
  return(data.frame(vx, vy, x, y))
}