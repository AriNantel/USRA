library(readxl)
library(deSolve)
library(sindyr)
library(pracma)
library(palinsol)

dt <- 0.01

# Time span

# Here we define our dimensionless time for generated data define by t = t_star * 10 from the barry saltzman book
t_star <- seq(0, 50, by = 0.1)
times <- 500 - t_star * 10

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

isl <- as.matrix(isl_df$isl)
isl_normalized <- (isl - mean(isl)) / sd(isl)
# Here t_star refers to the non dimentional time, when generating synthetic data
R_interp <- approxfun(t_star, isl_normalized, rule = 2)

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

added_nonsense <- sin(2 * pi * t_star / 5)

xs_gen_with_sin <- xs_gen_normalized
xs_gen_with_sin$n <- normalized(added_nonsense)

#result <- sindyc(xs_gen_normalized, u = isl_normalized, dt = 0.1, lambda = 0.1)
result <- sindyc(xs_gen_with_sin, u = isl_normalized, dt = 0.1, lambda = 0.1)

print(result$B)
print(get_equations(result$B))