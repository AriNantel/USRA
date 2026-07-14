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
  #isl <- state[4]

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

  # Normal sindy with lambda = 0.25 + isl
  #dx <- 0.3057160 - 0.9142313*x - 0.5402391*y + 1.7062704*x^2 - 0.9999887*x*z + 0.4260319*y^2 + 0.4916019*z^2 - 0.8476447*isl - 1.2256619*x^3 + 1.1392757*x^2*z + 2.0181136*x*y*z - 1.1069817*x*z^2 + 0.8936664*x*isl + 1.7407349*y*isl + 0.4430091*z^3 - 1.3411695*y*z^2 - 0.6169967*x*y*isl + 0.2512526*y*z*isl - 0.4487506*y^2*isl + 0.6696400*isl^2 - 0.6448090*x*isl^2 - 1.2963118*y*isl^2
  #dy <- -0.7054636 + 1.7712878*x + 1.8607147*y - 1.2172802*x^2 - 3.5326022*x*y + 1.2067898*x*z - 1.3664124*y^2 - 0.9770882*y*z + 0.5260653*z^2 + 0.6448168*isl + 1.8102388*x^2*y - 1.2568418*x^2*z + 1.5721412*x*y^2 - 0.7621718*x*y*z - 2.1769685*x*isl + 0.9497213*y^2*z + 0.6958018*y*z^2 - 0.8460432*z^3 - 1.8335723*z*isl + 1.5140968*x^2*isl + 0.5867626*x*y*isl + 1.1925783*x*z*isl - 0.9914675*y*z*isl + 1.8608833*z^2*isl + 1.1573759*isl^2 - 0.4578916*x*isl^2 - 0.5689652*z*isl^2 - 0.4281389*isl^3
  #dz <- -0.7045642 + 3.4702170*x + 0.9342548*y + 1.2276720*z - 5.0678030*x^2 - 4.5672862*x*y - 1.8498788*x*z - 1.2487683*y*z - 1.1396437*z^2 - 1.6373735*isl + 2.3561814*x^3 + 3.7464498*x^2*y + 0.4580969*x^2*z + 1.5851805*x*y^2 + 0.5804718*x*y*z + 1.2705964*x*z^2 + 3.2873763*x*isl - 0.2664695*y^3 + 1.2120315*y*z^2 + 2.3606284*y*isl - 1.8242324*x^2*isl - 2.6197066*x*y*isl - 1.0730341*x*z*isl - 0.8593215*y^2*isl - 1.1667054*y*z*isl + 0.4434772*x*isl^2 + 0.4963072*y*isl^2 - 0.2801780*isl^3
  #disl <- 0

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

plot_combined_result <- function(common_timsecale, xs_normalized, recovered_equations, title = "SINDy ODE extraction from datasets"){
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
    c("x, y, z from SINDy"),
    col = c("orange"),
    pch = 16,
    bty = "n",
    horiz = TRUE
  )
}

plot_stacked_result <- function(common_timsecale, xs_normalized, recovered_equations, title = NULL){

  par(mfrow = c(2, 2))

  plot(
    x = common_timsecale,
    y = xs_normalized$x,
    type = "l",
    ylim = c(-2, 2),
    ylab = "(-)",
    xlab = "Time (10 ky)",
    xlim = c(11, 400),
    main = "SINDy ODE extraction from datasets for ice extent"
  )

  lines(x = common_timsecale,
        recovered_equations$x,
        col = "orange",
        lty = 3, lwd = 4)

  legend(
    "topright",
    c("Ice Data", "Recovered"),
    col = c("black", "orange"),
    lty = c(1, 3),
    lwd = c(1, 4),
    bty = "n"
  )

  plot(
    x = common_timsecale,
    y = xs_normalized$y,
    col = "red",
    type = "l",
    ylim = c(-2, 2),
    ylab = "(-)",
    xlab = "Time (10 ky)",
    xlim = c(11, 400),
    main = "SINDy ODE extraction from datasets for CO2 concentration"
  )

  lines(x = common_timsecale,
        recovered_equations$y,
        col = "orange",
        lty = 3, lwd = 4)

  legend(
    "topright",
    c("CO2 Dataset", "Recovered"),
    col = c("red", "orange"),
    lty = c(1, 3),
    lwd = c(1, 4),
    bty = "n"
  )

  plot(
    x = common_timsecale,
    y = xs_normalized$z,
    col = "blue",
    type = "l",
    ylim = c(-2, 2),
    ylab = "(-)",
    xlab = "Time (10 ky)",
    xlim = c(11, 400),
    main = "SINDy ODE extraction from datasets for Ocean temperature"
  )
  lines(x = common_timsecale,
        recovered_equations$z,
        col = "orange",
        lty = 3, lwd = 4)

  legend(
    "topright",
    c("Dataset", "Recovered"),
    col = c("blue", "orange"),
    lty = c(1, 3),
    lwd = c(1, 4),
    bty = "n"
  )

  abline(h = 0, lty = 2)

  par(mfrow = c(1, 1))
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



#plot_scatter(ice_volume_clean$Age, ice_volume_clean$Ice_Volume, smooth_ice, common_timsecale, "ice_volume", "Raw Ice data vs smoothed curve")
#plot_scatter(co2_clean$Age, co2_clean$CO2, smooth_co2, common_timsecale, "co2", "Raw co2 data vs smoothed curve")
#plot_scatter(ocean_temp_clean$Age, ocean_temp_clean$Ocean_Temp, smooth_ocean_temp, common_timsecale, "Ocean Temp", "Raw ocean temp data vs smoothed curve")

plot_stacked_result(common_timsecale = common_timsecale, xs_normalized = xs_normalized, recovered_equations = df.b)

print("i finishied running")