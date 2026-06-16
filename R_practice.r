library(readxl)
library(deSolve)
library(sindyr)

deep_ocean_temp_original <- read.csv(
  "miller2024-sealevel.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

deep_ocean_temp_clean <- deep_ocean_temp_original[, c("Age", "pT_0")]
names(deep_ocean_temp_clean) <- c("Age", "Deep_Ocean_Temp")
deep_ocean_temp_clean$Age <- deep_ocean_temp_clean$Age * 1000
deep_ocean_temp_clean <- deep_ocean_temp_clean[deep_ocean_temp_clean$Age >= 0 & deep_ocean_temp_clean$Age <= 400000, ]

ice_extent_original <- read.csv(
  "lisiecki2005-d18o-stack-noaa.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

ice_extent_clean <- ice_extent_original[, c("age_calkaBP", "d18O_benthic")]
names(ice_extent_clean) <- c("Age", "Ice_Extent")
ice_extent_clean$Age <- ice_extent_clean$Age * 1000
ice_extent_clean <- ice_extent_clean[ice_extent_clean$Age >= 0 & ice_extent_clean$Age <= 400000, ]

sea_level_original <- read_excel(
  "1-s2.0-S0012821X15003404-mmc2.xlsx",
  sheet = "Stacks"
)

sea_level_original <- as.data.frame(sea_level_original)

sea_level_clean <- sea_level_original[, c("Age", "Detrended sea-level equivalent")]
names(sea_level_clean) <- c("Age", "Sea_Level")
sea_level_clean$Age <- sea_level_clean$Age * 1000
sea_level_clean <- sea_level_clean[sea_level_clean$Age >= 0 & sea_level_clean$Age <= 400000, ]

co2_original <- read_excel(
  "antarctica2015co2.xls",
  sheet = "CO2 Composite",
  skip = 14
)

co2_original <- as.data.frame(co2_original)

co2_clean <- co2_original[, c("Gasage (yr BP)", "CO2 (ppmv)")]
names(co2_clean) <- c("Age", "CO2")
co2_clean <- co2_clean[co2_clean$Age >= 0 & co2_clean$Age <= 400000, ]

# print(head(deep_ocean_temp_clean))
# print(head(ice_extent_clean))
# print(head(co2_clean))

# Creating a common timescale by jumps of 5000 as that is largest gap
common_timsecale <- seq(0, 400000, by = 5000)

# Create data on the same timescale for Ice extent, Co2 concenntration and deep ocean temp
interpolated_ice_extent <- approx(
  x = ice_extent_clean$Age,
  y = ice_extent_clean$Ice_Extent,
  xout = common_timsecale
)

interpolated_sea_level <- approx(
  x = sea_level_clean$Age,
  y = sea_level_clean$Sea_Level,
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
  x = interpolated_sea_level$y,
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

# print(head(xs_normalized))
# print(tail(xs))

# Building Theta matrix
# where x is the ice extent, y is the atmospheric CO2 and z deep ocean temp

x <- xs_scaled$x
y <- xs_scaled$y
z <- xs_scaled$z

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
  xs_scaled,
  2,
  as.numeric(xs_scaled[1, ]),
  FUN = "-"
)
dx_matrix <- as.matrix(dx)

theta_int <- matrix(0, nrow(theta), ncol(theta))

theta_int[1, ] <- 5000 * theta[1, ]

for (n in 2:nrow(theta)) {
  theta_int[n, ] <- theta_int[n - 1, ] + 5000 * theta[n, ]
}

colnames(theta_int) <- colnames(theta)


# Using SINDyr library
sindy.obj = sindy(xs = xs_scaled, dt = 5000, lambda = 1e-5)
print(sindy.obj$B)

# Implementing least square regression on my own