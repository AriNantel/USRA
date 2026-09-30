library(readxl)
library(deSolve)
library(sindyr)
library(pracma)
library(palinsol)

datasets <- clean_data()

ice_volume_clean <- datasets[[1]]
co2_clean <- datasets[[2]]
ocean_temp_clean <- datasets[[3]]
SST_clean <- datasets[[5]]
GAST_clean <- datasets[[4]]
dust_clean <- datasets[[6]]
#air_temp_clean <- datasets[[5]]

dt <- 0.01

# Time span

# Here we define our dimensionless time for generated data define by t = t_star * 10 from the barry saltzman book
t_star <- seq(0, 50, by = 0.1)
times <- t_star * 10

# Here we defien a common_timescale for the datasets
common_timescale <- seq(400, 11, by = -dt)
#common_timescale <- seq(11, 400, dt)
model_time <- seq(from = 0, by = dt, length.out = length(common_timescale))

# Create data on the same timescale for Ice extent, Co2 concenntration and deep ocean temp

smooth_ice_fit <- smooth.spline(ice_volume_clean$Age,  ice_volume_clean$Ice_Volume)
smooth_ice <- predict(smooth_ice_fit, common_timescale)$y

smooth_co2_fit <- smooth.spline(co2_clean$Age,  co2_clean$CO2)
smooth_co2 <- predict(smooth_co2_fit, common_timescale)$y

smooth_ocean_temp_fit <- smooth.spline(ocean_temp_clean$Age,  ocean_temp_clean$Ocean_Temp)
smooth_ocean_temp <- predict(smooth_ocean_temp_fit, common_timescale)$y

smooth_SST_fit <- smooth.spline(SST_clean$Age, SST_clean$SST_Temp)
smooth_SST <- predict(smooth_SST_fit, common_timescale)$y

smooth_GAST_fit <- smooth.spline(GAST_clean$Age, GAST_clean$GAST)
smooth_GAST <- predict(smooth_GAST_fit, common_timescale)$y

smooth_Dust_fit <- smooth.spline(dust_clean$Age, dust_clean$Dust)
smooth_Dust <- predict(smooth_Dust_fit, common_timescale)$y

#smooth_air_temp_fit <- smooth.spline(air_temp_clean$Age, air_temp_clean$air_temp)
#smooth_air_temp <- predict(smooth_air_temp_fit, common_timescale)$y

normalized_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

normalized <- function(x) {
  (x - mean(x)) / sd(x)
}

solar_radiation <- function(common_timescale) {

  # Palinsol requires ages in years before present but our common_timescale is in ky BP so we multiply it by 1000. 
  # They also consider present as 0 and negatives as going back in time so we multiply our time axis to reflect the paleoclimate data
  orbital_time <- -common_timescale * 1000

  #print(head(orbital_time))

  # We want the oldest time to be our first timestep so we must reverse the order to go from 500 kt BP to 0
  #orbital_time <- rev(orbital_time)

  latitude <- 65 * pi/180 # (radians)

  insolation <- function(times, astrosol=ber78,...)
  sapply(times, function(tt) Insol(orbit=astrosol(tt), lon = pi/2, lat = latitude))

  # Daily mean incoming solar radiation at TOA (W/m2)
  isl <- insolation(orbital_time, ber78)
  #print(head(isl))
  #print(head(orbital_time))
  isl_df <- data.frame(age=orbital_time, isl=isl)
  return(isl_df)
}

# Used for generated data
isl_df = solar_radiation(common_timescale = times)

# Used for dataset data
#isl_df = solar_radiation(common_timescale = common_timescale)

isl <- as.matrix(isl_df$isl)
isl_normalized <- (isl - mean(isl)) / sd(isl)
# Here t_star refers to the non dimentional time, when generating synthetic data
R_interp <- approxfun(t_star, isl_normalized, rule = 2)

# Here we interpolate the insolation data based on the common_timescale for the datasets
# R_interp <- approxfun(model_time, isl_normalized, rule = 2)



xs <- data.frame(
  x = smooth_ice,
  y = smooth_co2,
  z = smooth_ocean_temp,
  w = smooth_GAST,
  v = smooth_SST,
  s = smooth_Dust
  #v = smooth_air_temp
)

# Get rid of NA values
xs <- na.omit(xs)

xs_normalized <- as.data.frame(lapply(xs, normalized))

# print(nrow(xs_normalized))
# print(head(xs_normalized))
# print(tail(xs))

# Define parameters (change these to explore different behavior)
p <- 1
q <- 2.5
r <- 1.3
s <- 0.6
v <- 0.2
u <- 0.5

# Define the system of ODEs for generated data with external forcing, if we do not want external forcing change to u <- 0
generate_system_Milank <- function(t, state, parameters) {
  x <- state[1]
  y <- state[2]
  z <- state[3]

  Rt <- R_interp(t)

  dx <- -x - y - v * z - u * Rt
  dy <- -p * z + r * y - s * y^2 - y^3
  dz <- -q * (x + z)

  list(c(dx, dy, dz))
}

# Initial conditions for generated data
initial_conditions <- c(0.6, -0.8, -0.2)

# Inital conditions for dataset data, we use the first row of the datasets
#initial_conditions <- as.numeric(xs_normalized[1,])
#initial_conditions_2D <- initial_conditions[1:2]
#initial_conditions_3D <- initial_conditions[1:3]
#initial_conditions_4D <- initial_conditions[1:4]

# Solve ODEs for each initial condition for the generated data
solutions_generated_data <- ode(
  y = initial_conditions,
  times = t_star,
  func = generate_system_Milank,
  parms = NULL,
  method = "lsoda"
)

df <- data.frame(solutions_generated_data)
colnames(df) <- c("time", "ice", "co2", "temp")

xs_generated <- data.frame(df$ice, df$co2, df$temp)
colnames(xs_generated) <- c("x", "y", "z")

xs_gen_normalized <- as.data.frame(lapply(xs_generated, normalized))

# lambda_values = list(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07)

# for (lam in lambda_values) {
#   # Using SINDyr library
#   sindy.obj = sindyc(xs = xs_normalized, u = isl_normalized, dt = dt, lambda = lam)
#   cat("\n========================\n")
#   cat("Lambda =", lam, "\n")
#   cat("========================\n")
#   print(get_equations(sindy.obj$B))
# }

#sindy.obj = sindyc(xs = xs_normalized, u = isl_normalized, dt = dt, lambda = 0.05)
#print(sindy.obj$B)
#print(get_equations(sindy.obj$B))

#sindy.obj = sindy(xs = xs_normalized, dt = dt, lambda = 0.015)
#print(sindy.obj$B)
#print(get_equations(sindy.obj$B))


#test = sindyc(xs = xs_normalized, dt = dt, lambda = 0.0)
#print(test$B)

#result <- sindyc(xs_gen_normalized, u = isl_normalized, dt = 0.1, lambda = 0.1)
#print(result$B)