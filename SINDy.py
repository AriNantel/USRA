import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import pysindy as ps
from scipy.integrate import solve_ivp

def data_cleaning():
    # Import the CO2 dataset from an excel file, specify which sheet to look at and where the data starts
    co2_concentration_original = pl.read_excel(
        source = "antarctica2015co2.xls", 
        sheet_name="CO2 Composite",
        read_options={"header_row":14}
    )
    
    # Import the deep ocean dataset from a text file, define how the data is sepreated and filter out the comments at the beginning of the file
    deep_ocean_temp_original = pl.read_csv(
        source= "miller2024-sealevel.txt",
        separator="\t",
        comment_prefix="#"
    )

    # Import the gloabl ocean mean dataset from an excel file and specify where the data starts
    ocean_temp_original = pl.read_excel(
        source = "edc2021noblegastemp.xlsx", 
        read_options={"header_row":116}
        )
    
    # Import the ice volume dataset from a text file, define how the data is sepreated, filter out the comments at the beginning of the file and the non-existing values
    global_ice_volume_original = pl.read_csv(
        source="rohling2021-wh-main.txt",
        separator="\t",
        comment_prefix="#",
        null_values="NA",
        infer_schema_length=100000
    )

# Clean the CO2 dataset, we keep only the age and CO2 column. Ensuring the age is kyr and only between 0 and 400 000 years
    co2_concentration_clean = (
        co2_concentration_original
        .select(
            (pl.col("Gasage (yr BP)")/1000).alias("Age"),
            pl.col("CO2 (ppmv)").alias("CO2"))
        .filter(
           pl.col("Age").is_between(0,400))
        .sort("Age")
    )    

    # Clean the deep_ocean dataset, we keep only the age and smoothed ocean temp columns. Ensuring the age is kyr and only between 0 and 400 000 years
    deep_ocean_temp_clean = (
        deep_ocean_temp_original
        .select(
            pl.col("Age"),
            pl.col("pT_0 ").alias("Ocean Temp"))
        .filter(
            pl.col("Age").is_between(0, 400))
        .sort("Age")
    )

    # Clean the ocean temp dataset, we keep only the age and ocean temp columns. Ensuring the age is kyr and only between 0 and 400 000 years. We also filter out any column where the values are non-existing
    ocean_temp_clean = (
        ocean_temp_original
        .select(
            (pl.col("gas_age_ky")).alias("Age"),
            pl.col("MOT_KrN2").alias("Ocean Temp"))
        .filter(
            ~pl.any_horizontal(pl.all() == -999)
            & pl.col("Age").is_between(0, 400))
        .sort("Age")
    )

    # Clean the ice colume dataset, we keep only the age and global ice volume columns. Ensuring the age is kyr and only between 0 and 400 000 years.
    global_ice_volume_clean = (
        global_ice_volume_original
        .select(
            (pl.col("t(ka)") * -1).alias("Age"),
            pl.col("V_NH+SH").alias("Ice_Volume"))
        .filter(
            pl.col("Age").is_between(0, 400))
        .sort("Age")
        )

    # Return the 3 clean datasets we want to use
    return co2_concentration_clean, ocean_temp_clean, global_ice_volume_clean

def interpolate(common_timescale, df, x_var, y_var):

    # Interpolate a dataset so it is on a common timescale
    return np.interp(common_timescale["Age"], df[x_var], df[y_var])

def normalize_data(xs, col1, col2, col3):
    # normalize the data. For each column take a value substract the minimum value in the column and divide it by the biggest difference in the column
    xs_norm = (
        xs.with_columns(
            ((pl.col(col1) - pl.col(col1).min()) /
            (pl.col(col1).max() - pl.col(col1).min())
            ).alias("x_norm"),

            ((pl.col(col2) - pl.col(col2).min()) /
            (pl.col(col2).max() - pl.col(col2).min())
            ).alias("y_norm"),

            ((pl.col(col3) - pl.col(col3).min()) /
            (pl.col(col3).max() - pl.col(col3).min())
            ).alias("z_norm")
        )
        .select("x_norm", "y_norm", "z_norm")
    )

    # Return the normalized data
    return xs_norm

def generate_data(dt = 0.01, initial_conditions = [0.1, 0.8, -0.2], start=0, end=100): 
    
    # Create a timeseries
    times = np.arange(start, end+dt, dt)

    # Generate data for the initial condition given
    data = solve_ivp(system_original,[0,100],initial_conditions,t_eval=times, method="LSODA")
    xs = data.y.T

    # Return the data and the timeseries
    return xs, times

def generate_data_with_noise(noise_level, dt = 0.01, initial_conditions = [0.1, 0.8, -0.2], start=0, end=100):
    # generate clean data
    xs, times = generate_data(dt, initial_conditions, start, end)

    # Create gaussian noise based on the standard deviation of the data in xs.
    noise = noise_level * np.std(xs, axis=0) * np.random.randn(*xs.shape)
    
    # Return the noisy data and the timeseries
    return xs+noise, times

def plot_data_and_derivative(X, dt):
    # CHAT WROTE THIS
    fd = ps.FiniteDifference()
    Xdot = fd._differentiate(X, t=dt)

    fig, ax = plt.subplots(3, 2, figsize=(12, 8))

    names = ["x", "y", "z"]

    for i in range(3):
        # Data
        ax[i, 0].plot(X[:, i])
        ax[i, 0].set_title(names[i])

        # Derivative
        ax[i, 1].plot(Xdot[:, i])
        ax[i, 1].set_title(f"{names[i]}'")

    plt.tight_layout()
    plt.show()

def plot(dt=0.01,initial_conditions=[0.1, 0.8, -0.2]):
    # Get the time scale for the plot
    times = np.arange(0,30 + dt, dt)

    # Numerically solve the original equations (the Saltzman equations)
    sol_clean = solve_ivp(
        system_original,
        [0, 30],
        initial_conditions,
        t_eval=times,
        method="LSODA"
    )

    # Numerically solve the recovered ODEs from SINDy
    sol_recovered = solve_ivp(
        system_recovered,
        [0, 30],
        initial_conditions,
        t_eval=times,
        method="LSODA"
    )

    # Transpose the two solution arrays so they are the right shape
    df_clean = sol_clean.y.T
    df_recovered = sol_recovered.y.T

    plt.figure(figsize=(10, 8))

    # Recovered SINDy lines
    plt.plot(times, df_recovered[:, 0], color="orange", linestyle=":", linewidth=4, label="x, y, z from Weak SINDy")
    plt.plot(times, df_recovered[:, 1], color="orange", linestyle=":", linewidth=4)
    plt.plot(times, df_recovered[:, 2], color="orange", linestyle=":", linewidth=4)

    # Original lines
    plt.plot(times, df_clean[:, 0], color="black", linewidth=1.5, label="x = Ice")
    plt.plot(times, df_clean[:, 1], color="red", linewidth=1.5, label="y = CO2")
    plt.plot(times, df_clean[:, 2], color="blue", linewidth=1.5, label="z = T_Ocean")

    plt.axhline(0, color="black", linestyle="--", linewidth=1)

    plt.ylim(-2, 2)
    plt.xlim(0, 30)

    plt.xlabel("Time (10 ky)")
    plt.ylabel("(-)")
    plt.title("SINDy ODE extraction for simulated data with 5% noise", fontweight="bold")

    plt.legend(loc="upper center", ncol=3, frameon=False)

    plt.show()

def system_original(state):
    # define parameters from the Saltzman equations
    p = 1
    q = 2.5
    r = 1.3
    s = 0.6
    v = 0.2

    x,y,z = state

    # Define the saltzman equations
    dx = -x - y - v * z
    dy = -p * z + r * y - s * y**2 - y**3
    dz = -q * (x + z)

    return [dx, dy, dz]

def system_recovered(state):
    x,y,z = state

    # The recovered sindy equations
    dx = -0.968*x - 0.993*y - 0.172*z
    dy = 0.384*x + 1.342*y - 0.683*z - 0.555*y**2 - 0.916*y**3 + 0.012*z**3
    dz = 0.515*y - 0.423*z + 0.143*y**2 + 0.240*y**3

    return [dx, dy, dz]

def normal_sindy(xs, lam, alpha, timescale,  library = ps.PolynomialLibrary(degree=3)):
    
    # Define an optimizer to use
    opt = ps.STLSQ(threshold = lam, alpha = alpha)
    
    # Define the SINDy model
    model = ps.SINDy(feature_library=library, optimizer=opt)

    # Fit the data to the model
    model.fit(xs, t=timescale, feature_names=["x", "y", "z"])
    
    # Print the model (prints the recovered equtions)
    model.print()

def weak_sindy(xs, lam, alpha, timescale, library = ps.PolynomialLibrary(degree=3, include_bias=False)):    
    # Define a library for SINDy to use
    ode_library = ps.WeakPDELibrary(
        function_library= library,
        spatiotemporal_grid=timescale,
        include_bias=True)

    # Define an optimizer to use
    opt = ps.STLSQ(threshold=lam, alpha=alpha)
    #opt = ps.SR3(reg_weight_lam=0.1, regularizer="L1", relax_coeff_nu=0.01,max_iter=100 )

    # Define the SINDy model
    model = ps.SINDy(feature_library=ode_library, optimizer=opt)

    # Fit the data to the model
    model.fit(xs, t=timescale, feature_names=["x", "y", "z"])
    
    # Print the model (prints the recovered equtions)
    model.print()

def main():
    # Get the data from the datasets
    co2, temp, ice_vol = data_cleaning()
    
    # Define a common timescale for the datasets going from 0 to 400 000 with jumps of 5000 years
    dt = 5
    common_timescale = pl.DataFrame({"Age": range(0, 401, dt)})

    # Interpolate the data to a common timsecale
    co2_inter = interpolate(common_timescale, co2, "Age", "CO2")
    deep_temp_inter = interpolate(common_timescale, temp, "Age", "Ocean Temp")
    ice_vol_inter = interpolate(common_timescale, ice_vol, "Age", "Ice_Volume")

    # Combine the data into one dataframe
    xs = pl.DataFrame({"Age": common_timescale["Age"], "CO2": co2_inter, "temp": deep_temp_inter, "ice_vol": ice_vol_inter})

    # normalize the data
    xs_norm = normalize_data(xs, "CO2", "temp", "ice_vol")

    # Convert the dataframe to a numpy array to ensure compatibility with sindy librart
    xs_norm = xs_norm.to_numpy()

main()

data, timescale = generate_data_with_noise(0.05)
#data, timescale = generate_data()

print("normal sindy")
normal_sindy(data, lam=0.1, alpha=1, timescale=timescale)

print("weak sindy")
weak_sindy(data, lam=0.1, alpha=1, timescale=timescale)


