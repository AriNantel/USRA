library(readxl)
library(deSolve)
library(sindyr)
library(pracma)

deep_ocean_temp_original <- read.csv(
  "miller2024-sealevel.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

deep_ocean_temp_clean <- deep_ocean_temp_original[, c("Age", "pT_0")]
names(deep_ocean_temp_clean) <- c("Age", "Deep_Ocean_Temp")
#deep_ocean_temp_clean$Age <- deep_ocean_temp_clean$Age * 1000
deep_ocean_temp_clean <- deep_ocean_temp_clean[deep_ocean_temp_clean$Age >= 0 & deep_ocean_temp_clean$Age <= 400, ]

# Sea level and oxygen isotopes can inidcate more or less ice sheets but i not a direct proxy for sea ice extent/ ice volume
#ice_extent_original <- read.csv(
#  "lisiecki2005-d18o-stack-noaa.txt",
#  header = TRUE,
#  sep = "\t",
#  comment.char = "#"
#)

#ice_extent_clean <- ice_extent_original[, c("age_calkaBP", "d18O_benthic")]
#names(ice_extent_clean) <- c("Age", "Ice_Extent")
#ice_extent_clean$Age <- ice_extent_clean$Age * 1000
#ice_extent_clean <- ice_extent_clean[ice_extent_clean$Age >= 0 & ice_extent_clean$Age <= 400000, ]

#sea_level_original <- read_excel(
#  "1-s2.0-S0012821X15003404-mmc2.xlsx",
#  sheet = "Stacks"
#)

#sea_level_original <- as.data.frame(sea_level_original)

#sea_level_clean <- sea_level_original[, c("Age", "Detrended sea-level equivalent")]
#names(sea_level_clean) <- c("Age", "Sea_Level")
#sea_level_clean$Age <- sea_level_clean$Age * 1000
#sea_level_clean <- sea_level_clean[sea_level_clean$Age >= 0 & sea_level_clean$Age <= 400000, ]

ice_volume_original <- read.csv(
  "bintanja2008-noaa.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

ice_volume_clean <- ice_volume_original[, c("age_calkaBP", "Ice_nam")]
names(ice_volume_clean) <- c("Age", "Ice_Volume")
#ice_volume_clean$Age <- ice_volume_clean$Age * 1000
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

# print(head(deep_ocean_temp_clean))
# print(head(ice_extent_clean))
# print(head(co2_clean))

dt <- 5

# Creating a common timescale by jumps of 5000 as that is largest gap
common_timsecale <- seq(0, 400, by = dt)

# Create data on the same timescale for Ice extent, Co2 concenntration and deep ocean temp


#interpolated_ice_extent <- approx(
#  x = ice_extent_clean$Age,
#  y = ice_extent_clean$Ice_Extent,
#  xout = common_timsecale
#)

#interpolated_sea_level <- approx(
#  x = sea_level_clean$Age,
#  y = sea_level_clean$Sea_Level,
#  xout = common_timsecale
#)

interpolated_Ice_Volume <- approx(
  x = ice_volume_clean$Age,
  y = ice_volume_clean$Ice_Volume,
  xout = common_timsecale
)

interpolated_co2 <- approx(
  x = co2_clean$Age,
  y = co2_clean$CO2,
  xout = common_timsecale
)

interpolated_deep_ocean_temp <- approx(
  x = deep_ocean_temp_clean$Age,
  y = deep_ocean_temp_clean$Deep_Ocean_Temp,
  xout = common_timsecale
)

xs <- data.frame(
  x = interpolated_Ice_Volume$y,
  y = interpolated_co2$y,
  z = interpolated_deep_ocean_temp$y
)

xs <- na.omit(xs)

# Scale my matrix of data on common timescale
xs_scaled <- as.data.frame(scale(xs))

normalized <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

xs_normalized <- as.data.frame(lapply(xs, normalized))

print(nrow(xs_normalized))
# print(head(xs_normalized))
# print(tail(xs))

# Building Theta matrix
# where x is the ice extent, y is the atmospheric CO2 and z deep ocean temp

x <- xs_normalized$x
y <- xs_normalized$y
z <- xs_normalized$z

theta <- cbind(
  "(Intercept)" = 1,
  x = x,
  y = y,
  z = z,
  "x:x" = x^2,
  "x:y" = x * y,
  "x:z" = x * z,
  "y:y" = y^2,
  "y:z" = y * z,
  "z:z" = z^2,
  "x:x:x" = x^3,
  "x:x:y" = x^2 * y,
  "x:x:z" = x^2 * z,
  "x:y:y" = y^2 * x,
  "x:y:z" = x * y * z,
  "x:z:z" = z^2 * x,
  "y:y:y" = y^3,
  "y:y:z" = y^2 * z,
  "y:z:z" = z^2 * y,
  "z:z:z" = z^3
)

# Create integrated differential matrix for soft SINDy
dx <- sweep(
  xs_normalized,
  2,
  as.numeric(xs_normalized[1, ]),
  FUN = "-"
)
dx_matrix <- as.matrix(dx)

theta_int <- matrix(0, nrow(theta), ncol(theta))

theta_int[1, ] <- dt * theta[1, ]

for (n in 2:nrow(theta)) {
  theta_int[n, ] <- theta_int[n - 1, ] + dt * theta[n, ]
}

colnames(theta_int) <- colnames(theta)

# Define parameters (change these to explore different behavior)
p <- 1
q <- 2.5
r <- 1.3
s <- 0.6
v <- 0.2

# Define the system of ODEs
system <- function(t, state, parameters) {
  x <- state[1]
  y <- state[2]
  z <- state[3]
  
  dx <- -x - y - v * z
  dy <- -p * z + r * y - s * y^2 - y^3
  dz <- -q * (x + z)
  
  list(c(dx, dy, dz))
}

# Initial conditions
initial_conditions <- c(0.1, 0.8, -0.2)

# Time span

dt <- 0.01
times <- seq(0, 100, by = dt)

# Solve ODEs for each initial condition
solutions <- ode(
  y = initial_conditions,
  times = times,
  func = system,
  parms = NULL,
  method = "lsoda"
)

df <- data.frame(solutions)
colnames(df) <- c("time", "x", "y", "z")

xs <- data.frame(df$x, df$y, df$z)
colnames(xs) <- c("x", "y", "z")

set.seed(20)
whitenoise_x = ts(rnorm(79, mean = 0.18))
whitenoise_y = ts(rnorm(79, mean = 0.43))
whitenoise_z = ts(rnorm(79, mean = 0.03))

whitenoise_xs <- data.frame(x = whitenoise_x, y = whitenoise_y, z = whitenoise_z)

# Using SINDyr library
sindy.obj = sindy(xs = xs_normalized, dt = 5, lambda = 0.01)
print(sindy.obj$B)
