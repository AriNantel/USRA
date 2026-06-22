# Load required package
library(deSolve)
library(sindyr)

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
times <- seq(0, 30, by = dt)


# Solve ODEs for each initial condition
solutions <- ode(
  y = initial_conditions,
  times = times,
  func = system,
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
  
dx <- -0.968 * x + -0.993 * y + -0.172 * z

dy <-  0.384 * x +  1.342 * y + -0.683 * z + -0.555 * y^2 + -0.916 * y^3 +  0.012 * z^3

dz <-  0.515 * y + -0.423 * z +  0.143 * y^2 +  0.240 * y^3
  
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


plot(
  x = df.a$time,
  y = df.a$x,
  type = "l",
  ylim = c(-2, 2),
  ylab = "(-)",
  xlab = "Time (10 ky)",
  xlim = c(0, 30),
  main = "SINDy ODE extraction for simulated data with 5% noise"
)

lines(x = df.b$time,
      df.b$x,
      col = "orange",
      lty = 3, lwd = 4)
lines(x = df.b$time,
      df.b$y,
      col = "orange",
      lty = 3, lwd = 4)
lines(x = df.b$time,
      df.b$z,
      col = "orange",
      lty = 3, lwd = 4)

lines(x = df.a$time, df.a$x, col = "black")
lines(x = df.a$time, df.a$y, col = "red")
lines(x = df.a$time, df.a$z, col = "blue")


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

print("i finishied running")