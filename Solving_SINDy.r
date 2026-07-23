library(readxl)
library(deSolve)
library(sindyr)
library(pracma)
library(palinsol)

clean_data <- function() {
  ocean_10_original <- read_excel(
    "sosdian2009.xls",
    sheet = "23-24 MgCa BWT",
    skip = 5
  )

  ocean_150_original <- read_excel(
    "sosdian2009.xls",
    sheet = "607 MgCa BWT",
    skip = 5
  )

  ocean_10_original <- as.data.frame(ocean_10_original)
  ocean_150_original <- as.data.frame(ocean_150_original)

  ocean_temp_original <- rbind(ocean_10_original, ocean_150_original)

  ocean_temp_clean <- ocean_temp_original[, c("Age (ka)", "BWT (°C)")]
  names(ocean_temp_clean) <- c("Age", "Ocean_Temp")
  #deep_ocean_temp_clean$Age <- deep_ocean_temp_clean$Age * 1000
  ocean_temp_clean <- ocean_temp_clean[ocean_temp_clean$Age >= 10.37 & ocean_temp_clean$Age <= 400, ]
  ocean_temp_clean <- ocean_temp_clean[order(ocean_temp_clean$Age, decreasing = TRUE),]

  ice_volume_original <- read.csv(
    "rohling2021-wh-main.txt",
    header = TRUE,
    sep = "\t",
    comment.char = "#"
  )

  ice_volume_clean <- ice_volume_original[, c("t.ka.", "V_NH.SH")]
  names(ice_volume_clean) <- c("Age", "Ice_Volume")
  ice_volume_clean$Age <- ice_volume_clean$Age * -1
  ice_volume_clean <- ice_volume_clean[order(ice_volume_clean$Age), ]
  ice_volume_clean <- ice_volume_clean[ice_volume_clean$Age >= 0 & ice_volume_clean$Age <= 400, ]
  ice_volume_clean <- ice_volume_clean[order(ice_volume_clean$Age, decreasing = TRUE),]

  co2_original <- read_excel(
    "antarctica2015co2.xls",
    sheet = "CO2 Composite",
    skip = 14
  )

  co2_original <- as.data.frame(co2_original)

  co2_clean <- co2_original[, c("Gasage (yr BP)", "CO2 (ppmv)")]
  names(co2_clean) <- c("Age", "CO2")
  co2_clean$Age <- co2_clean$Age / 1000
  co2_clean <- co2_clean[co2_clean$Age >= 0 & co2_clean$Age <= 400, ]
  co2_clean <- co2_clean[order(co2_clean$Age, decreasing = TRUE),]

  return(list(ice_volume_clean, co2_clean, ocean_temp_clean))
}

datasets <- clean_data()

ice_volume_clean <- datasets[[1]]
co2_clean <- datasets[[2]]
ocean_temp_clean <- datasets[[3]]

dt <- 0.01

# Time span

# Here we define our dimensionless time for generated data define by t = t_star * 10 from the barry saltzman book
t_star <- seq(0, 50, by = 0.1)
times <- t_star * 10

# Here we defien a common_timescale for the datasets
common_timescale <- seq(400, 11, by = -dt)
model_time <- seq(from = 0, by = dt, length.out = nrow(xs_normalized))

# Create data on the same timescale for Ice extent, Co2 concenntration and deep ocean temp

smooth_ice_fit <- smooth.spline(ice_volume_clean$Age,  ice_volume_clean$Ice_Volume)
smooth_ice <- predict(smooth_ice_fit, common_timescale)$y

smooth_co2_fit <- smooth.spline(co2_clean$Age,  co2_clean$CO2)
smooth_co2 <- predict(smooth_co2_fit, common_timescale)$y

smooth_ocean_temp_fit <- smooth.spline(ocean_temp_clean$Age,  ocean_temp_clean$Ocean_Temp)
smooth_ocean_temp <- predict(smooth_ocean_temp_fit, common_timescale)$y

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
  print(head(orbital_time))
  isl_df <- data.frame(age=orbital_time, isl=isl)
  return(isl_df)
}

# Used for generated data
# isl_df = solar_radiation(common_timescale = times)

# Used for dataset data
isl_df = solar_radiation(common_timescale = common_timescale)

isl <- as.matrix(isl_df$isl)
isl_normalized <- (isl - mean(isl)) / sd(isl)
# Here t_star refers to the non dimentional time, when generating synthetic data
# R_interp <- approxfun(t_star, isl_normalized, rule = 2)

# Here we interpolate the insolation data based on the common_timescale for the datasets
R_interp <- approxfun(model_time, isl_normalized, rule = 2)

xs <- data.frame(
  x = smooth_ice,
  y = smooth_co2,
  z = smooth_ocean_temp
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
# initial_conditions <- c(0.6, -0.8, -0.2)

# Inital conditions for dataset data, we use the first row of the datasets
initial_conditions <- as.numeric(xs_normalized[1,])

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

xs_plus_forcing <- cbind(xs_normalized,isl_normalized)
names(xs_plus_forcing) <- c("x", "y", "z", "u")

#Theta = features(xs_plus_forcing, polyorder = 2)

# Using SINDyr library
sindy.obj = sindyc(xs = xs_normalized, u = isl_normalized, dt = dt, lambda = 0.015)
print(sindy.obj$B)

#test = sindyc(xs = xs_normalized, dt = dt, lambda = 0.0)
#print(test$B)

#result <- sindyc(xs_gen_normalized, u = isl_normalized, dt = 0.1, lambda = 0.1)
#print(result$B)