# Load required package
library(deSolve)
library(sindyr)

source("R_practice.R")

# Define parameters (change these to explore different behavior)
p <- 1
q <- 2.5
r <- 1.3
s <- 0.6
v <- 0.2

# Define the system of ODEs according to barry saltzman equations
saltzman_system <- function(t, state, parameters) {
  x <- state[1]
  y <- state[2]
  z <- state[3]

  dx <- -x - y - v * z
  dy <- -p * z + r * y - s * y^2 - y^3
  dz <- -q * (x + z)

  list(c(dx, dy, dz))
}

# Initial conditions
initial_conditions <- as.numeric(xs_normalized[1, ])

# Time span

dt <- 0.01
#times <- seq(0, 30, by = dt)
times <- common_timsecale


# Solve ODEs for each initial condition
solutions <- ode(
  y = initial_conditions,
  times = times,
  func = saltzman_system,
  parms = NULL,
  method = "lsoda"
)

# Plot
df.a <- data.frame(solutions)
colnames(df.a) <- c("time", "x", "y", "z")


# Define the system of ODEs
system <- function(t, state, parameters) {
  x <- state[1]
  y <- state[2]
  z <- state[3]

  #dx <- -0.968 * x + -0.993 * y + -0.172 * z
  #dy <-  0.384 * x +  1.342 * y + -0.683 * z + -0.555 * y^2 + -0.916 * y^3 +  0.012 * z^3
  #dz <-  0.515 * y + -0.423 * z +  0.143 * y^2 +  0.240 * y^3

  #dx <- 0.303 - 1.171 * x - 1.122 * y + 0.234 * z + 1.687 * x^2 + 1.973 * x*y - 0.666 *x*z + 1.824 * y^2 - 0.584 *y*z - 0.853 *x*x*x - 0.946 *x*x*y + 0.345 *x*x*z - 1.628 *x*y*y + 1.200 *x*y*z - 0.824 *y*y*y + 0.101 *y*y*z
  #dy <- -0.371 + 1.156 * x + 1.028 * y - 0.989 * x^2 - 1.812 *x*y - 0.560 * y^2 - 1.062 *y*z + 0.475 * z^2 + 0.712 *x*x*y + 0.485 *x*x*z + 0.424 *x*y*y + 0.851 *x*y*z - 0.900 *x*z*z + 0.605 *y*y*z
  #dz <- -0.583 + 1.711 * x - 0.535 * y + 1.740 * z - 1.552 * x^2 + 0.831 *x*y - 3.175 *x*z + 1.735 * y^2 - 1.324 *y*z - 0.982 *z*z + 0.409 *x*x*x - 0.443 *x*x*y + 1.390 *x*x*z - 1.307 *x*y*y + 1.095 *x*y*z + 0.926 *x*z*z - 0.750 *y*y*y - 0.145 *y*y*z + 0.715 *y*z*z + 0.106 *z*z*z
  
  # Normal sindy with lambda = 0.25
  dx <- -0.3408959*x - 0.6553861*y + 0.5914701*z + 0.9908404*x^2 + 1.2947956*x*y - 1.4918978*x*z + 1.6883207*y^2 - 1.1672142*y*z - 0.6813752*x^3 - 0.7988121*x^2*y + 0.8126251*x^2*z - 1.6405513*x*y^2 + 1.8816830*x*y*z - 0.8566249*y^3 + 0.3341154*y^2*z
  dy <- -0.3710603 + 1.1563477*x + 1.0282703*y - 0.9894505*x^2 - 1.8122942*x*y - 0.5600973*y^2 - 1.0615566*y*z + 0.4749107*z^2 + 0.7124899*x^2*y + 0.4846215*x^2*z + 0.4238264*x*y^2 + 0.8508179*x*y*z - 0.8996630*x*z^2 + 0.6048002*y^2*z
  dz <- -0.5859008 + 1.6880933*x - 0.5573774*y + 1.7447657*z - 1.4640817*x^2 + 0.9334177*x*y - 3.3018644*x*z + 1.8605015*y^2 - 1.5307828*y*z - 0.8136639*z^2 + 0.3576680*x^3 - 0.5424893*x^2*y + 1.4410937*x^2*z - 1.4551431*x*y^2 + 1.3245459*x*y*z + 0.9100749*x*z^2 - 0.8440160*y^3 + 0.6907300*y*z^2

  list(c(dx, dy, dz))
}

# Solve ODEs for each initial condition
solutions <- ode(
  y = initial_conditions,
  times = times,
  func = system,
  parms = NULL,
  method = "lsoda"
)

df.b <- data.frame(solutions)
colnames(df.b) <- c("time", "x", "y", "z")

plot_result <- function(common_timsecale, xs_normalized, recovered_equations, title = "SINDy ODE extraction from datasets"){
  plot(
    x = common_timsecale,
    y = xs_normalized$x,
    type = "l",
    ylim = c(-2, 2),
    ylab = "(-)",
    xlab = "Time (10 ky)",
    xlim = c(11, 400),
    main = title
  )

  lines(x = common_timsecale,
        recovered_equations$x,
        col = "orange",
        lty = 3, lwd = 4)
  lines(x = common_timsecale,
        recovered_equations$y,
        col = "orange",
        lty = 3, lwd = 4)
  lines(x = common_timsecale,
        recovered_equations$z,
        col = "orange",
        lty = 3, lwd = 4)

  #lines(x = df.a$time, df.a$x, col = "black")
  #lines(x = df.a$time, df.a$y, col = "red")
  #lines(x = df.a$time, df.a$z, col = "blue")

  lines(x = common_timsecale, xs_normalized$x, col = "black")
  lines(x = common_timsecale, xs_normalized$y, col = "red")
  lines(x = common_timsecale, xs_normalized$z, col = "blue")


  abline(h = 0, lty = 2)
  legend(
    "top",
    c("x = Ice", "y = CO2", "z = T_Ocean"),
    col = c("black", "red", "blue"),
    pch = 16,
    bty = "n",
    horiz = TRUE
  )
  legend(
    "bottom",
    c("x, y, z from Weak SINDy"),
    col = c("orange"),
    pch = 16,
    bty = "n",
    horiz = TRUE
  )
}

plot_scatter <- function(raw_x, raw_y, smooth_y, timescale, y_label, title = NULL) {
  plot(
    x = raw_x,
    y = raw_y,
    type = "p",
    col = "red",
    ylab = y_label,
    xlab = "Time (10 ky)",
    xlim = c(0, 400),
    main = title
  )

  lines(x = timescale, y = smooth_y, col = "black")

  legend(
    "top",
    c("x = raw data", "y = smooth data"),
    col = c("red", "black"),
    pch = c(16, NA),
    lty = c(NA, 1),
    bty = "n",
    horiz = TRUE
  )

  invisible(NULL)
}

plot_scatter(ice_volume_clean$Age, ice_volume_clean$Ice_Volume, smooth_ice, common_timsecale, "ice_volume", "Raw Ice data vs smoothed curve")
plot_scatter(co2_clean$Age, co2_clean$CO2, smooth_co2, common_timsecale, "co2", "Raw co2 data vs smoothed curve")
plot_scatter(ocean_temp_clean$Age, ocean_temp_clean$Ocean_Temp, smooth_ocean_temp, common_timsecale, "Ocean Temp", "Raw ocean temp data vs smoothed curve")

plot_result(common_timsecale = common_timsecale, xs_normalized = xs_normalized, recovered_equations = df.b)

print("i finishied running")