source("Solving_SINDy.R")
source("Plot.R")

recovered_data <- data.frame(solutions_3D)
colnames(recovered_data) <- c("time", "x", "y", "z")
recovered_data <- recovered_data[, c("x", "y", "z")]
recovered_data <- na.omit(recovered_data)

original_data <- xs_normalized

n <- min(nrow(original_data), nrow(recovered_data))

original_data <- original_data[1:(n), ]
recovered_data <- recovered_data[1:(n), ]

rmse_x <- sqrt(mean((original_data$x - recovered_data$x)^2))
rmse_y <- sqrt(mean((original_data$y - recovered_data$y)^2))
rmse_z <- sqrt(mean((original_data$z - recovered_data$z)^2))

print(rmse_x)
print(rmse_y)
print(rmse_z)