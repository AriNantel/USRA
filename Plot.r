# Load required package
library(deSolve)
library(sindyr)

source("Solving_SINDy.R")

# Define parameters (change these to explore different behavior)
p <- 1
q <- 2.5
r <- 1.3
s <- 0.6
v <- 0.2
u <- 0.5

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

# Solve ODEs for each initial condition
saltzman_solutions <- ode(
  y = initial_conditions,
  times = times,
  func = saltzman_system,
  parms = NULL,
  method = "lsoda"
)

# Plot
df.a <- data.frame(saltzman_solutions)
colnames(df.a) <- c("time", "x", "y", "z")


# Define the system of ODEs
system <- function(t, state, parameters) {
  x <- state[1]
  y <- state[2]
  z <- state[3]
  
  u <- parameters$R_interp(t)

  #dx <- -0.968 * x + -0.993 * y + -0.172 * z
  #dy <-  0.384 * x +  1.342 * y + -0.683 * z + -0.555 * y^2 + -0.916 * y^3 +  0.012 * z^3
  #dz <-  0.515 * y + -0.423 * z +  0.143 * y^2 +  0.240 * y^3

  #dx <- 0.303 - 1.171 * x - 1.122 * y + 0.234 * z + 1.687 * x^2 + 1.973 * x*y - 0.666 *x*z + 1.824 * y^2 - 0.584 *y*z - 0.853 *x*x*x - 0.946 *x*x*y + 0.345 *x*x*z - 1.628 *x*y*y + 1.200 *x*y*z - 0.824 *y*y*y + 0.101 *y*y*z
  #dy <- -0.371 + 1.156 * x + 1.028 * y - 0.989 * x^2 - 1.812 *x*y - 0.560 * y^2 - 1.062 *y*z + 0.475 * z^2 + 0.712 *x*x*y + 0.485 *x*x*z + 0.424 *x*y*y + 0.851 *x*y*z - 0.900 *x*z*z + 0.605 *y*y*z
  #z <- -0.583 + 1.711 * x - 0.535 * y + 1.740 * z - 1.552 * x^2 + 0.831 *x*y - 3.175 *x*z + 1.735 * y^2 - 1.324 *y*z - 0.982 *z*z + 0.409 *x*x*x - 0.443 *x*x*y + 1.390 *x*x*z - 1.307 *x*y*y + 1.095 *x*y*z + 0.926 *x*z*z - 0.750 *y*y*y - 0.145 *y*y*z + 0.715 *y*z*z + 0.106 *z*z*z
  
  # Normal sindy with lambda = 0.25
  #dx <- -0.3408959*x - 0.6553861*y + 0.5914701*z + 0.9908404*x^2 + 1.2947956*x*y - 1.4918978*x*z + 1.6883207*y^2 - 1.1672142*y*z - 0.6813752*x^3 - 0.7988121*x^2*y + 0.8126251*x^2*z - 1.6405513*x*y^2 + 1.8816830*x*y*z - 0.8566249*y^3 + 0.3341154*y^2*z
  #dy <- -0.3710603 + 1.1563477*x + 1.0282703*y - 0.9894505*x^2 - 1.8122942*x*y - 0.5600973*y^2 - 1.0615566*y*z + 0.4749107*z^2 + 0.7124899*x^2*y + 0.4846215*x^2*z + 0.4238264*x*y^2 + 0.8508179*x*y*z - 0.8996630*x*z^2 + 0.6048002*y^2*z
  #dz <- -0.5859008 + 1.6880933*x - 0.5573774*y + 1.7447657*z - 1.4640817*x^2 + 0.9334177*x*y - 3.3018644*x*z + 1.8605015*y^2 - 1.5307828*y*z - 0.8136639*z^2 + 0.3576680*x^3 - 0.5424893*x^2*y + 1.4410937*x^2*z - 1.4551431*x*y^2 + 1.3245459*x*y*z + 0.9100749*x*z^2 - 0.8440160*y^3 + 0.6907300*y*z^2

  # Normal sindy with lambda = 0.25 + isl
  #dx <- 0.3057160 - 0.9142313*x - 0.5402391*y + 1.7062704*x^2 - 0.9999887*x*z + 0.4260319*y^2 + 0.4916019*z^2 - 0.8476447*isl - 1.2256619*x^3 + 1.1392757*x^2*z + 2.0181136*x*y*z - 1.1069817*x*z^2 + 0.8936664*x*isl + 1.7407349*y*isl + 0.4430091*z^3 - 1.3411695*y*z^2 - 0.6169967*x*y*isl + 0.2512526*y*z*isl - 0.4487506*y^2*isl + 0.6696400*isl^2 - 0.6448090*x*isl^2 - 1.2963118*y*isl^2
  #dy <- -0.7054636 + 1.7712878*x + 1.8607147*y - 1.2172802*x^2 - 3.5326022*x*y + 1.2067898*x*z - 1.3664124*y^2 - 0.9770882*y*z + 0.5260653*z^2 + 0.6448168*isl + 1.8102388*x^2*y - 1.2568418*x^2*z + 1.5721412*x*y^2 - 0.7621718*x*y*z - 2.1769685*x*isl + 0.9497213*y^2*z + 0.6958018*y*z^2 - 0.8460432*z^3 - 1.8335723*z*isl + 1.5140968*x^2*isl + 0.5867626*x*y*isl + 1.1925783*x*z*isl - 0.9914675*y*z*isl + 1.8608833*z^2*isl + 1.1573759*isl^2 - 0.4578916*x*isl^2 - 0.5689652*z*isl^2 - 0.4281389*isl^3
  #dz <- -0.7045642 + 3.4702170*x + 0.9342548*y + 1.2276720*z - 5.0678030*x^2 - 4.5672862*x*y - 1.8498788*x*z - 1.2487683*y*z - 1.1396437*z^2 - 1.6373735*isl + 2.3561814*x^3 + 3.7464498*x^2*y + 0.4580969*x^2*z + 1.5851805*x*y^2 + 0.5804718*x*y*z + 1.2705964*x*z^2 + 3.2873763*x*isl - 0.2664695*y^3 + 1.2120315*y*z^2 + 2.3606284*y*isl - 1.8242324*x^2*isl - 2.6197066*x*y*isl - 1.0730341*x*z*isl - 0.8593215*y^2*isl - 1.1667054*y*z*isl + 0.4434772*x*isl^2 + 0.4963072*y*isl^2 - 0.2801780*isl^3
  #disl <- 0

  # Generated data with milankovich cycles age not handeeled properly
  #dx <- -0.9743054 * x - 0.9905567 * y - 0.1695750 * z - 0.6530760 * u
  #dy <-  1.3918707 * y - 0.9470937 * z - 0.5413422 * y^3
  #dz <- -2.626116 * x - 2.489054 * z

  #dx <- -0.1050817 + 0.9452510*x - 1.2341000*y + 0.8407446*z + 0.2892553*u - 0.1836959*x^2 - 0.3210831*x*z + 0.1994932*x*u + 0.2323529*y*u - 0.8115401*x^3 - 1.6594451*x^2*y - 0.5806966*x^2*z - 0.3151262*x*y^2 - 0.9783775*x*y*z + 0.1663541*x*y*u + 0.8925698*x*z^2 - 0.2343162*x*z*u + 0.4264578*y^3 + 0.3605472*y^2*z + 0.2126087*y*z^2 + 1.0944890*z^3 - 0.2747702*z^2*u
  #dy <- 1.3918752*y - 0.9470934*z - 0.5413448*y^3
  #dz <- -2.626120*x - 2.489058*z

  #dx <- 0.2002533 + 1.2093752*x - 0.5578123*y + 0.3589730*z + 0.3801035*x^2 + 0.4604642*x*y + 0.5045046*x*z + 0.4966401*x*u + 0.5561062*y*z + 0.4612660*z*u + 2.3730036*x^3 + 0.1823070*x^2*y + 8.6887071*x^2*z - 0.5427830*x*y^2 + 1.9026007*x*y*z - 0.6833693*x*y*u + 10.3215821*x*z^2 + 0.9345179*x*z*u + 0.4061236*y^2*z + 1.0780984*y*z^2 - 0.4142158*y*z*u + 4.7997624*z^3 + 0.6979408*z^2*u
  #dy <- 1.3175486*y - 0.9415787*z - 0.4901178*y^3
  #dz <- -2.606445*x - 2.487482*z

  # From generated normalized data with milank forcing, times = t_star * 10 and rev(orbital_forcing)
  #dx <- -0.7521633*x - 0.9463480*y - 0.6415690*u - 0.2516033*x^3 - 0.6365085*x^2*z - 0.6159260*x*z^2 - 0.2146560*z^3
  #dy <- 1.1438281*y - 0.9578595*z + 0.3521423*x^3 + 0.3638186*x^2*y + 1.3898780*x^2*z + 1.0405844*x*y*z + 1.5228534*x*z^2 - 0.4361185*y^3 + 0.7635863*y*z^2 + 0.4940146*z^3
  #dz <- -2.606444*x - 2.487481*z

  # From generated normalized data with Milank forcing, times = t_star * 10, rev(orbital_forcing) and reveresed xs_gen_normalized
  #dx <- -0.2002370 - 1.2094880*x + 0.5577688*y - 0.3590866*z - 0.3800624*x^2 - 0.4604083*x*y - 0.5044614*x*z - 0.4966272*x*u - 0.5560662*y*z - 0.4612516*z*u - 2.3729734*x^3 - 0.1825837*x^2*y - 8.6882477*x^2*z + 0.5426606*x*y^2 - 1.9028461*x*y*z + 0.6833494*x*y*u - 10.3207904*x*z^2 - 0.9344692*x*z*u - 0.4061986*y^2*z - 1.0780707*y*z^2 + 0.4142111*y*z*u - 4.7994029*z^3 - 0.6979054*z^2*u
  #dy <- -1.3175430*y + 0.9415805*z + 0.4901105*y^3
  #dz <- 2.606444*x + 2.487481*z

  # From generated data(not normalized) with Milank forcing (normalized), times = t_star * 10, rev(orbital_forcing) and ode generated using t_star
  #dx <- -0.9701782*x - 0.9924768*y - 0.1751359*z - 0.4886877*u
  #dy <- 1.2909528*y - 0.9966243*z - 0.5949453*y^2 - 0.9899890*y^3
  #dz <- -2.487969*x - 2.487968*z

  # From generated data (normalized) with Milank forcing (normalized), times = t_star * 10, rev(orbital_forcing) and ode generated using t_star
  #dx <- -0.9690131*x - 0.9738830*y - 0.1655041*z - 0.6649171*u
  #dy <- 1.3959684*y - 0.9655704*z + 0.1338047*y^2 - 0.5138604*y^3
  #dz <- -2.614082*x - 2.487804*z

  # Smoothed dataset data, with milankovich forcing and lambda = 0.01
  dx <- 0.03305706 - 0.17166919*x - 0.13689997*y - 0.07031300*u - 0.07580152*x^2 - 0.21686702*x*y - 0.05408039*x*u - 0.12018713*y^2 + 0.03549356*y*z - 0.04060199*y*u - 0.04235663*z*u - 0.05413009*u^2 + 0.08059330*x^3 + 0.15171564*x^2*y - 0.01803044*x^2*z + 0.04005045*x^2*u + 0.14292149*x*y^2 - 0.04684976*x*y*z - 0.03980929*x*u^2 + 0.05880614*y^3 - 0.05746835*y^2*u + 0.04778819*y*z*u - 0.01720171*y*u^2 - 0.02944868*z*u^2 - 0.02320989*u^3
  dy <- -0.05905496 + 0.02361021*x + 0.01017174*y + 0.03494736*u + 0.06048077*x^2 + 0.02646227*x*y + 0.04259619*x*u + 0.01284813*y^2 + 0.02064518*y*u + 0.03162277*z*u + 0.03007478*u^2 - 0.03498883*x^2*y - 0.06482912*x^2*z + 0.02624559*x*y^2 - 0.19042663*x*y*z - 0.03187537*x*y*u + 0.08353988*x*z^2 + 0.10855532*x*z*u + 0.01832256*x*u^2 + 0.03745223*y^3 - 0.10590102*y^2*z + 0.04895868*y*z^2 + 0.01506577*z^2*u + 0.03687272*z*u^2 + 0.02200117*u^3
  dz <- -0.01321882 + 0.07883262*x + 0.05775988*z + 0.04033800*u + 0.08873214*x^2 + 0.15122189*x*y + 0.01481510*x*u + 0.06328353*y^2 + 0.04142469*y*u - 0.03524741*x^3 - 0.03957674*x^2*y - 0.02386801*x^2*z + 0.03105348*x*y^2 - 0.03868658*x*y*z + 0.03525951*x*y*u - 0.03087339*x*z^2 - 0.04984465*x*z*u + 0.04047085*y^3 + 0.03785207*y^2*u - 0.03477240*y*z^2 - 0.06716840*y*z*u + 0.02831880*y*u^2 - 0.01474965*z^3 - 0.04435001*z*u^2 - 0.01291093*u^3
  
  list(c(dx, dy, dz))
}

get_equations <- function() {

  # Equations for dataset data
  #dx_equation <- "dx = -0.3408959*x - 0.6553861*y + 0.5914701*z + 0.9908404*x^2 + 1.2947956*x*y - 1.4918978*x*z + 1.6883207*y^2 - 1.1672142*y*z - 0.6813752*x^3 - 0.7988121*x^2*y + 0.8126251*x^2*z - 1.6405513*x*y^2 + 1.8816830*x*y*z - 0.8566249*y^3 + 0.3341154*y^2*z"
  #dy_equation <- "dy = -0.3710603 + 1.1563477*x + 1.0282703*y - 0.9894505*x^2 - 1.8122942*x*y - 0.5600973*y^2 - 1.0615566*y*z + 0.4749107*z^2 + 0.7124899*x^2*y + 0.4846215*x^2*z + 0.4238264*x*y^2 + 0.8508179*x*y*z - 0.8996630*x*z^2 + 0.6048002*y^2*z"
  #dz_equation <- "dz = -0.5859008 + 1.6880933*x - 0.5573774*y + 1.7447657*z - 1.4640817*x^2 + 0.9334177*x*y - 3.3018644*x*z + 1.8605015*y^2 - 1.5307828*y*z - 0.8136639*z^2 + 0.3576680*x^3 - 0.5424893*x^2*y + 1.4410937*x^2*z - 1.4551431*x*y^2 + 1.3245459*x*y*z + 0.9100749*x*z^2 - 0.8440160*y^3 + 0.6907300*y*z^2"

  # From generated normalized data with milankovitch forcing, age was not handled properly
  #dx_equation <- "-0.9743054*x - 0.9905567*y - 0.1695750*z - 0.6530760*u"
  #dy_equation <-  "1.3918707*y - 0.9470937*z - 0.5413422*y**3"
  #dz_equation <- "-2.626116*x - 2.489054*z"

  # From generated normalized data with milankovitvh forcing, times = 500 - t_star * 10 and oribat forcing rev
  #dx_equation <- "-0.1050817 + 0.9452510*x - 1.2341000*y + 0.8407446*z + 0.2892553*u - 0.1836959*x^2 - 0.3210831*x*z + 0.1994932*x*u + 0.2323529*y*u - 0.8115401*x^3 - 1.6594451*x^2*y - 0.5806966*x^2*z - 0.3151262*x*y^2 - 0.9783775*x*y*z + 0.1663541*x*y*u + 0.8925698*x*z^2 - 0.2343162*x*z*u + 0.4264578*y^3 + 0.3605472*y^2*z + 0.2126087*y*z^2 + 1.0944890*z^3 - 0.2747702*z^2*u"
  #dy_equation <- "1.3918752*y - 0.9470934*z - 0.5413448*y^3"
  #dz_equation <- "-2.626120*x - 2.489058*z"

  # From generated normalized data with milank forcing, times = t_star * 10 and rev(orbital_forcing)
  #dx_equation <- "-0.7521633*x - 0.9463480*y - 0.6415690*u - 0.2516033*x^3 - 0.6365085*x^2*z - 0.6159260*x*z^2 - 0.2146560*z^3"
  #dy_equation <- "1.1438281*y - 0.9578595*z + 0.3521423*x^3 + 0.3638186*x^2*y + 1.3898780*x^2*z + 1.0405844*x*y*z + 1.5228534*x*z^2 - 0.4361185*y^3 + 0.7635863*y*z^2 + 0.4940146*z^3"
  #dz_equation <- "-2.606444*x - 2.487481*z"

  # From generated normalized data with Milank forcing, times = t_star * 10, rev(orbital_forcing) and reveresed xs_gen_normalized
  #dx_equation <- "-0.2002370 - 1.2094880*x + 0.5577688*y - 0.3590866*z - 0.3800624*x^2 - 0.4604083*x*y - 0.5044614*x*z - 0.4966272*x*u - 0.5560662*y*z - 0.4612516*z*u - 2.3729734*x^3 - 0.1825837*x^2*y - 8.6882477*x^2*z + 0.5426606*x*y^2 - 1.9028461*x*y*z + 0.6833494*x*y*u - 10.3207904*x*z^2 - 0.9344692*x*z*u - 0.4061986*y^2*z - 1.0780707*y*z^2 + 0.4142111*y*z*u - 4.7994029*z^3 - 0.6979054*z^2*u"
  #dy_equation <- "-1.3175430*y + 0.9415805*z + 0.4901105*y^3"
  #dz_equation <- "2.606444*x + 2.487481*z"

  # From generated data with Milank forcing (normalized), times = t_star * 10, rev(orbital_forcing) and ode generated using t_star
  #dx_equation <- "-0.9701782*x - 0.9924768*y - 0.1751359*z - 0.4886877*u"
  #dy_equation <- "1.2909528*y - 0.9966243*z - 0.5949453*y^2 - 0.9899890*y^3"
  #dz_equation <- "-2.487969*x - 2.487968*z"

  # From generated data (normalized) with Milank forcing (normalized), times = t_star * 10, rev(orbital_forcing) and ode generated using t_star
  #dx_equation <- "-0.9690131*x - 0.9738830*y - 0.1655041*z - 0.6649171*u"
  #dy_equation <- "1.3959684*y - 0.9655704*z + 0.1338047*y^2 - 0.5138604*y^3"
  #dz_equation <- "-2.614082*x - 2.487804*z"

  dx_equation <- "in progress"
  dy_equation <- "in progress"
  dz_equation <- "in progress"

  c(dx_equation, dy_equation, dz_equation)
}

parms <- list(R_interp = R_interp)

# Solve ODEs for each initial condition
solutions <- ode(
  y = initial_conditions,
  # use t_star for generated data
  #times = t_star,
  times = model_time,
  func = system,
  parms = parms,
  method = "lsoda"
)

# Solve system using recovered ODEs
df.b <- data.frame(solutions)
colnames(df.b) <- c("time", "x", "y", "z")

#df.c <- data.frame(solutions_generated_data)
#colnames(df.c) <- c("time", "x", "y", "z")

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

plot_stacked_result <- function(common_timsecale, xs_normalized, recovered_equations, recovered_data, external_forcing, title = NULL){

  par(mfrow = c(2, 2))

  plot(
    x = common_timsecale,
    y = xs_normalized$x,
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    xlim = c(0, 500),
    main = "Ice Extent (X)"
  )

  lines(x = common_timsecale,
        recovered_data$x,
        col = "orange",
        lty = 2, lwd = 2)

  legend(
    "top",
    c("Ice Data", "Recovered"),
    col = c("black", "orange"),
    lty = c(1, 3),
    lwd = c(1, 4),
    bty = "n"
  )

  legend(
    "bottom",
    legend = recovered_equations[1],
    bty = "n",
    cex = 0.6
  )


  plot(
    x = common_timsecale,
    y = xs_normalized$y,
    col = "red",
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    xlim = c(0, 500),
    main = "CO2 Concentration (Y)"
  )

  lines(x = common_timsecale,
        recovered_data$y,
        col = "orange",
        lty = 2, lwd = 2)

  legend(
    "top",
    c("CO2 Dataset", "Recovered"),
    col = c("red", "orange"),
    lty = c(1, 3),
    lwd = c(1, 4),
    bty = "n"
  )

  legend(
    "bottom",
    legend = recovered_equations[2],
    bty = "n",
    cex = 0.6
  )

  plot(
    x = common_timsecale,
    y = xs_normalized$z,
    col = "blue",
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    xlim = c(0, 500),
    main = "Ocean Temperature (Z)"
  )
  lines(x = common_timsecale,
        recovered_data$z,
        col = "orange",
        lty = 2, lwd = 2)

  legend(
    "top",
    c("Ocean Temp Dataset", "Recovered"),
    col = c("blue", "orange"),
    lty = c(1, 3),
    lwd = c(1, 4),
    bty = "n"
  )

  legend(
    "bottom",
    legend = recovered_equations[3],
    bty = "n",
    cex = 0.6
  )

  plot(
    x = common_timsecale,
    y = external_forcing,
    col = c("orange"),
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    xlim = c(0, 500),
    main = " External forcing factor (M)"
  )

  legend(
    "top",
    c("Insolation"),
    col = c("orange"),
    lty = c(1),
    lwd = c(1),
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

plot_milank <- function(common_timsecale, generated_data, external_focring, title = NULL){

  par(mfrow = c(2, 2))

  plot(
    x = common_timsecale,
    y = generated_data$x,
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    #xlim = c(11, 400),
    xlim = rev(range(common_timsecale)),
    main = "Ice Extent (X)"
  )

  legend(
    "top",
    c("Ice Data"),
    col = c("black"),
    lty = c(1),
    lwd = c(1),
    bty = "n"
  )

  abline(h = 0, lty = 2)

  plot(
    x = common_timsecale,
    y = generated_data$y,
    col = "red",
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    #xlim = c(11, 400),
    xlim = rev(range(common_timsecale)),
    main = "CO2 concentration (Y)"
  )

  legend(
    "top",
    c("CO2 Dataset"),
    col = c("red"),
    lty = c(1),
    lwd = c(1),
    bty = "n"
  )

  abline(h = 0, lty = 2)

  plot(
    x = common_timsecale,
    y = generated_data$z,
    col = "blue",
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    #xlim = c(11, 400),
    xlim = rev(range(common_timsecale)),
    main = "Ocean temp (Z)"
  )

  legend(
    "top",
    c("Ocean Temp Dataset"),
    col = c("blue"),
    lty = c(1),
    lwd = c(1),
    bty = "n"
  )

  abline(h = 0, lty = 2)

  plot(
    x = common_timsecale,
    y = external_focring,
    col = c("orange"),
    type = "l",
    ylim = c(-2.5, 2.5),
    ylab = "(-)",
    xlab = "Time (ky)",
    #xlim = c(11, 400),
    xlim = rev(range(common_timsecale)),
    main = " External forcing factor (M)"
  )

  legend(
    "top",
    c("Insolation"),
    col = c("orange"),
    lty = c(1),
    lwd = c(1),
    bty = "n"
  )

  abline(h = 0, lty = 2)

  par(mfrow = c(1, 1))
}


#plot_scatter(ice_volume_clean$Age, ice_volume_clean$Ice_Volume, smooth_ice, common_timsecale, "ice_volume", "Raw Ice data vs smoothed curve")
#plot_scatter(co2_clean$Age, co2_clean$CO2, smooth_co2, common_timsecale, "co2", "Raw co2 data vs smoothed curve")
#plot_scatter(ocean_temp_clean$Age, ocean_temp_clean$Ocean_Temp, smooth_ocean_temp, common_timsecale, "Ocean Temp", "Raw ocean temp data vs smoothed curve")


plot_stacked_result(common_timsecale = common_timsecale, xs_normalized = xs_normalized, recovered_data = df.b, recovered_equations = get_equations(), external_forcing = isl_normalized)


#plot_milank(common_timsecale = common_timescale, generated_data = xs_normalized, external_focring = isl_normalized)

print("i finishied running")

# Todo reverse age to have first row as oldest time then rename colums to start at 0
# Double check milankovich is the right way and initial conditions