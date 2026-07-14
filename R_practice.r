library(readxl)
library(deSolve)
library(sindyr)
library(pracma)
library(palinsol)

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

dt <- 0.01

# Creating a common timescale by jumps of 5000 as that is largest gap
common_timsecale <- seq(11, 400, by = dt)

# Create data on the same timescale for Ice extent, Co2 concenntration and deep ocean temp

smooth_ice_fit <- smooth.spline(ice_volume_clean$Age,  ice_volume_clean$Ice_Volume)
smooth_ice <- predict(smooth_ice_fit, common_timsecale)$y

smooth_co2_fit <- smooth.spline(co2_clean$Age,  co2_clean$CO2)
smooth_co2 <- predict(smooth_co2_fit, common_timsecale)$y

smooth_ocean_temp_fit <- smooth.spline(ocean_temp_clean$Age,  ocean_temp_clean$Ocean_Temp)
smooth_ocean_temp <- predict(smooth_ocean_temp_fit, common_timsecale)$y

solar_radiation <- function(common_timsecale){

  orbital_time = -common_timsecale * 1000

  latitude <- 65 * pi/180 # (radians)

  insolation <- function(times, astrosol=ber78,...)
  sapply(times, function(tt) Insol(orbit=astrosol(tt), lon = pi/2, lat = latitude))

  # Daily mean incoming solar radiation at TOA (W/m2) 
  isl <- insolation(orbital_time, ber78)
  return(isl)
}

isl = solar_radiation(common_timsecale = common_timsecale)

xs <- data.frame(
  x = smooth_ice,
  y = smooth_co2,
  z = smooth_ocean_temp,
  isl = isl
)

# Get rid of NA values
xs <- na.omit(xs)

# Scale my matrix of data on common timescale
# xs_scaled <- as.data.frame(scale(xs))

normalized <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

xs_normalized <- as.data.frame(lapply(xs, normalized))

# print(nrow(xs_normalized))
# print(head(xs_normalized))
# print(tail(xs))

# Building Theta matrix
# where x is the ice extent, y is the atmospheric CO2 and z deep ocean temp

x <- xs_normalized$ice
y <- xs_normalized$co2
z <- xs_normalized$temp

# Define parameters (change these to explore different behavior)
#p <- 1
#q <- 2.5
#r <- 1.3
#s <- 0.6
#v <- 0.2

# Define the system of ODEs for generated data
#generate_system <- function(t, state, parameters) {
#  x <- state[1]
#  y <- state[2]
#  z <- state[3]
  
#  dx <- -x - y - v * z
#  dy <- -p * z + r * y - s * y^2 - y^3
#  dz <- -q * (x + z)
  
#  list(c(dx, dy, dz))
#}

# Initial conditions
# initial_conditions <- c(0.1, 0.8, -0.2)

# Time span
#dt <- 0.01
#times <- seq(0, 100, by = dt)

# Solve ODEs for each initial condition for the generated data
#solutions_generated_data <- ode(
#  y = initial_conditions,
#  times = times,
#  func = generate_system,
#  parms = NULL,
#  method = "lsoda"
#)

#df <- data.frame(solutions_generated_data)
#colnames(df) <- c("time", "ice", "co2", "temp")

#xs_generated <- data.frame(df$ice, df$co2, df$temp)
#colnames(xs_generated) <- c("ice", "co2", "temp")

#set.seed(20)
# whitenoise_x = ts(rnorm(79, mean = 0.18))
# whitenoise_y = ts(rnorm(79, mean = 0.43))
# whitenoise_z = ts(rnorm(79, mean = 0.03))

# whitenoise_xs <- data.frame(x = whitenoise_x, y = whitenoise_y, z = whitenoise_z)

# Using SINDyr library
sindy.obj = sindy(xs = xs_normalized, dt = 0.01, lambda = 0.25)
print(sindy.obj$B)