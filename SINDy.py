import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import pysindy as ps
from scipy.integrate import solve_ivp
from sklearn.metrics import mean_squared_error
from scipy.interpolate import UnivariateSpline
from rpy2.robjects import FloatVector
from rpy2.robjects.packages import importr
from rpy2.robjects import r

np.random.seed(42)

def data_cleaning():
    """
    This function imports the datasets used for SINDy and cleans the data so it can be used at a later time

    Args:
        None
    
    Returns:
         tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]: A tuple containing:
            - global_ice_volume_clean: Cleaned global ice volume data.
            - co2_concentration_clean: Cleaned CO2 concentration data.
            - ocean_clean: Cleaned ocean data.
    """

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

    # Import the ocean dataset starting at the year 150 ka BP, define the sheet name and the row at which the header starts
    ocean_150_original = pl.read_excel(
        source="sosdian2009.xls",
        sheet_name="607 MgCa BWT",
        read_options={"header_row":5}
    )

    # Import the ocean dataset starting at the year 10 ka BP, define the sheet name and the row at which the header starts
    ocean_10_original = pl.read_excel(
        source="sosdian2009.xls",
        sheet_name="23-24 MgCa BWT",
        read_options={"header_row":5}
    )

    # Concatenate the two ocean datasets, stacking them one on the other vertically. Starting with the earlier data
    ocean_original = pl.concat([ocean_10_original, ocean_150_original], how="vertical")

    # Cleaning the ocean data by selecting the age and ocean temperaature columns, filtering the age range from 10.37 to 400 ka bp
    ocean_clean = (
        ocean_original
        .select(
            pl.col("Age (ka)").alias("Age"),
            pl.col("BWT (°C)").alias("Ocean Temp"))
        .filter(
            pl.col("Age").is_between(10.37, 400))
        .sort("Age")
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
    return global_ice_volume_clean, co2_concentration_clean, ocean_clean

def interpolate(common_timescale, df, x_var, y_var):
    """
    This function takes in a dataframe and linearly interpolates it to a common timescale

    Args:
        common_timescale (np.ndarray): An array providing the timescale the data needs to be converted to
        df (pl.DataFrame): A datafram containing the data that needs to be interpolated to the common timescale
        x_var (string): Name of the X column for the data being interpolated
        y_var (String): Name of the Y column for the data being interpolated

    
    Returns:
         np.ndarray: An array containing the interpolated data and the common timescale
    """

    # Interpolate a dataset so it is on a common timescale
    return np.interp(common_timescale["Age"], df[x_var], df[y_var])

def smoothing_spline(common_timescale, df, x_var, y_var, s):
    """
    This function takes in a dataframe with an x and y column and creates a smoothing spline (fits a smoothed line to the data). 
    The smoothing spline is evaluated on a common timescale

    Args:
        common_timescale (np.ndarray): An array providing the timescale the data needs to be converted to
        df (pl.DataFrame): A datafram containing the data that needs to be smoothed and evaluated on a common timescale
        x_var (string): Name of the X column for the data being smoothed, containing the original age values
        y_var (String): Name of the Y column for the data being smoothed, containing the variable that requires smothing
        s (int): The smoothing factor. The larger the value the more smooth the curve will be

    
    Returns:
         np.ndarray: An array containing the smoothed data evaluated on the common timescale
    """

    # Convert the x and y columns to numpy so they can be evaluated
    x = df[x_var].to_numpy()
    y = df[y_var].to_numpy()

    # Create the smothed spline based on the provided smoothing factor
    smooth = UnivariateSpline(x, y)
    smooth.set_smoothing_factor(s)

    #return smoothed data evaluated on the common timescale
    return smooth(common_timescale).ravel()

def normalize_data(xs, col1, col2, col3):
    """
    This function takes in a dataframe with three column and performs min-max normalization so that the values vary between 0 and 1.

    Args:
        xs (pl.DataFrame): A datafram containing the data that needs to be normalized
        col1 (string): Name of the first column that will be normalized
        col2 (String): Name of the second column that will be normalized
        col3 (String): Name of the third column that will be normalized
    
    Returns:
         pl.DataFrame: A polars dataframe containing the normalized data for the three paleoclimate variables
    """

    # normalize the data. For each column take a value substract the minimum value in the column and divide it by the biggest difference in the column. Once normalized, keep only the normalized columns
    """xs_norm = (
        xs.with_columns(
            ((pl.col(col1) - pl.col(col1).min()) /
            (pl.col(col1).max() - pl.col(col1).min())
            ).alias("ice_norm"),

            ((pl.col(col2) - pl.col(col2).min()) /
            (pl.col(col2).max() - pl.col(col2).min())
            ).alias("co2_norm"),

            ((pl.col(col3) - pl.col(col3).min()) /
            (pl.col(col3).max() - pl.col(col3).min())
            ).alias("temp_norm")
        )
        .select("ice_norm", "co2_norm", "temp_norm")
    )"""

    xs_norm = (
        xs.with_columns(
            ((pl.col(col1) - pl.col(col1).min()) /
            (pl.col(col1).max() - pl.col(col1).min())
            ).alias("ice_norm"),

            ((pl.col(col2) - pl.col(col2).min()) /
            (pl.col(col2).max() - pl.col(col2).min())
            ).alias("co2_norm")
        )
        .select("ice_norm", "co2_norm")
    )

    # Return the normalized data
    return xs_norm

def generate_data(dt = 0.01, initial_conditions = [0.1, 0.8, -0.2], start=0, end=100): 
    """
    This function generates synthetic data for a defined timeseries with a given initial condition, by solving the original Barry Saltzman system of equations.

    Args:
        dt (float): The timestep between consecutive data points, if not specified the timestep will be 0.01
        initial_conditions (list[float]): The initial conditions for the 3 state variables, if not specified the initial conditions will be [0.1, 0.8, -0.2]
        start (int): The age at which the timeseries is starting, if not specified the timeseries will be 0
        end (int): The age at which the timeseries is ending, if not specified the timeseries will end at 100
    
    Returns:
         tuple [np.ndarray, np.ndarray]: A tuple containing:
            - xs: Array of shape (n_samples, 3) containing the simulated state variables.
            - times: One-dimensional array of time values corresponding to each row of xs.
    """ 
    # Create a timeseries
    times = np.arange(start, end+dt, dt)
    timespan = [start, end]
    
    # Generate data starting with the initial condition given
    xs = solve_system(system_original, timespan, times=times, initial_conditions=initial_conditions).y.T

    # Return the data and the timeseries
    return xs, times

def generate_data_with_noise(noise_level, dt = 0.01, initial_conditions = [0.1, 0.8, -0.2], start=0, end=100):
    """
    This function generates noisy synthetic data for a defined timeseries with a given initial condition, by solving the original Barry Saltzman system of equations.

    Args:
        noise_level (float): The percentage of noise added to the generated data
        dt (float): The timestep between consecutive data points, if not specified the timestep will be 0.01
        initial_conditions (list[float]): The initial conditions for the 3 state variables, if not specified the initial conditions will be [0.1, 0.8, -0.2]
        start (int): The age at which the timeseries is starting, if not specified the timeseries will be 0
        end (int): The age at which the timeseries is ending, if not specified the timeseries will end at 100
    
    Returns:
         tuple [np.ndarray, np.ndarray]: A tuple containing:
            - xs: Array of shape (n_samples, 3) containing the simulated state variables with the added noise.
            - times: One-dimensional array of time values corresponding to each row of xs.
    """
    
    # generate clean data
    xs, times = generate_data(dt, initial_conditions, start, end)

    # Create gaussian noise based on the standard deviation of the data in xs.
    noise = noise_level * np.std(xs, axis=0) * np.random.randn(*xs.shape)
    
    # Return the noisy data and the timeseries
    return xs+noise, times

def plot_graph(x, y, df_xaxis, df_yaxis, df_title):
    """
    This function plots a timeseries

    Args:
        x (np.ndarray): Values for the x-axis (typically the time array).
        y (np.ndarray): Values for the y-axis
        df_xaxis (String): Label for the x axis
        df_yaxis (String): Label for the y axis
        df_title (String): Title for the plot
    
    Returns:
        matplotlib.figure.Figure: The generated figure.
    """

    # Initialize the figure and the axes
    fig, ax = plt.subplots()

    # Plot the x and y data
    ax.plot(x, y)

    # invert the x axis the data is plotted forward in time
    ax.invert_xaxis()

    # Set the labels for the x and y axis and the title
    ax.set_xlabel(df_xaxis)
    ax.set_ylabel(df_yaxis)
    ax.set_title(df_title)

    # Return the figure
    return fig

def plot_data_and_derivative(xs, dt, labels = ["ice_vol", "co2", "ocean temp"]):
    """
    This function plots timeseries data and their derivatives

    Args:
        xs (np.ndarray): An array of data that needs to be plotted
        dt (float): The timestep between consecutive data points
        labels list[String]: labels for the columns in the xs array
    
    Returns:
        matplotlib.figure.Figure: The generated figure.
    """

    # CHAT WROTE THIS
    fd = ps.SmoothedFiniteDifference()

    # Get the derivatives for the data
    Xdot = fd._differentiate(xs, t=dt)

    # Initialize the figure to have 2 columns and 3 rows
    fig, ax = plt.subplots(3, 2, figsize=(12, 8))

    for i in range(3):
        # Data
        ax[i, 0].plot(xs[:, i])
        ax[i, 0].set_title(labels[i])

        # Derivative
        ax[i, 1].plot(Xdot[:, i])
        ax[i, 1].set_title(f"{labels[i]}'")

    plt.tight_layout()

    # Invert the axes for each graph so the data is plotted forward in time
    for row in ax:
        for a in row:
            a.invert_xaxis()

    # return the generated figure
    return fig

def plot_forward_in_time_noisy_data(noisy_data, dt=0.01,initial_conditions=[0.1, 0.8, -0.2], start = 0, end = 100):
    """
    This function plots timeseries data intergrated forward in time

    Args:
        noisy_data (np.ndarray): An array of the noisy data that needs to be plotted
        dt (float): The timestep between consecutive data points, if not specified the timestep will be 0.01
        initial_conditions (list[float]): The initial conditions for the 3 state variables, if not specified the initial conditions will be [0.1, 0.8, -0.2]
        start (int): The age at which the timeseries is starting, if not specified the timeseries will be 0
        end (int): The age at which the timeseries is ending, if not specified the timeseries will end at 100
    
    Returns:
        
    """
    
    # Get the time scale for the plot
    times = np.arange(start, end + dt, dt)
    timespan = [start, end]

    # Transpose the two solution arrays so they are the right shape
    df_clean = solve_system(system_original,timespan, times=times, initial_conditions=initial_conditions).y.T
    df_recovered = solve_system(system_recovered_gen_data, timespan, times = times, initial_conditions=initial_conditions).y.T

    plt.figure(figsize=(10, 8))

    # Noisy lines
    plt.plot(times, noisy_data[:, 0], color="purple", linestyle="--", linewidth=1, label="x = Noisy Ice")
    plt.plot(times, noisy_data[:, 1], color="purple", linestyle="--", linewidth=1, label="y = Noisy CO2")
    plt.plot(times, noisy_data[:, 2], color="purple", linestyle="--", linewidth=1, label="z = Noisy T_Ocean")

    # Recovered SINDy lines
    plt.plot(times, df_recovered[:, 0], color="orange", linestyle=":", linewidth=4, label="x, y, z from Weak SINDy")
    plt.plot(times, df_recovered[:, 1], color="orange", linestyle=":", linewidth=4)
    plt.plot(times, df_recovered[:, 2], color="orange", linestyle=":", linewidth=4)

    # Original lines
    plt.plot(times, df_clean[:, 0], color="black", linewidth=2, label="x = Ice")
    plt.plot(times, df_clean[:, 1], color="red", linewidth=2, label="y = CO2")
    plt.plot(times, df_clean[:, 2], color="blue", linewidth=2, label="z = T_Ocean")

    # plot the origin line
    plt.axhline(0, color="black", linestyle="--", linewidth=1)

    # set limits for the x and y axis
    plt.ylim(-2, 2)
    plt.xlim(0, 30)

    # set labels for the x-axis, y-axis and title
    plt.xlabel("Time (10 ky)")
    plt.ylabel("(-)")
    plt.title("SINDy ODE extraction for datasets", fontweight="bold")
    
    # Invert the axis so the data is plotted forward in time
    plt.gca().invert_xaxis()

    # Add a legend to the graph
    plt.legend(loc="upper center", ncol=3, frameon=False)

    plt.show()

def plot_forward_in_time_dataset(xs_norm, timescale = None, title = "Recovered SINDy ODE vs real dataset"):
    """
    This function plots timeseries data intergrated forward in time

    Args:
        xs_norm (np.ndarray): An array of the noisy data that needs to be plotted
        timescale (np.ndarray): The timescale associated with the data in xs_norm
        title (string): The title for the figure
    
    Returns:
        
    """
    
    # Get the time scale for the plot
    #times = np.arange(timespan[0], timespan[-1] + dt, dt)

    # Solve recovered equations, the solution array needs to be transposed so they are the right shape
    df_recovered = solve_system(system_recovered, times = [timescale[0], timescale[-1]], timescale = timescale, initial_conditions=xs_norm[0,:]).y.T

    plt.figure(figsize=(10, 8))

    # Recovered SINDy lines
    plt.plot(timescale, df_recovered[:, 0], color="orange", linestyle=":", linewidth=4, label="x, y, z from Weak SINDy")
    plt.plot(timescale, df_recovered[:, 1], color="orange", linestyle=":", linewidth=4)
    #plt.plot(timescale, df_recovered[:, 2], color="orange", linestyle=":", linewidth=4)

    # Original lines
    plt.plot(timescale, xs_norm[:, 0], color="black", linewidth=1.5, label="x = Ice Extent")
    plt.plot(timescale, xs_norm[:, 1], color="red", linewidth=1.5, label="y = CO2 Concentration")
    #plt.plot(timescale, xs_norm[:, 2], color="blue", linewidth=1.5, label="z = T_Ocean")

    plt.axhline(0, color="black", linestyle="--", linewidth=1)

    # Set the limits for the x axis
    plt.xlim(timescale[0], timescale[-1])
    #plt.xlim(0, 30)

    # Add labels for the x and y axis and a title to the figure
    plt.xlabel("Time (in ky BP)")
    plt.ylabel("Normalized values")
    plt.title(title, fontweight="bold")
    
    # Invert the xaxis
    #plt.gca().invert_xaxis()
    
    # Add a legend
    plt.legend(loc="upper center", ncol=3, frameon=False)

    plt.show()

def plot_stacked_timeseries(time, series, labels, xlabel="Time (kyr BP)", title=None):
    """
    time: 1D array of time values
    series: list of 1D arrays, one per variable
    labels: list of y-axis labels
    colors: optional list of colors
    """

    n = len(series)

    colors = ["tab:blue", "tab:red", "tab:green", "black"]

    fig, axes = plt.subplots(
        n,
        1,
        figsize=(9, 7),
        sharex=True
    )

    if n == 1:
        axes = [axes]

    for i in range(n):
        axes[i].plot(
            time,
            series[i],
            color=colors[i % len(colors)],
            linewidth=1.8
        )

        axes[i].set_ylabel(labels[i])
        axes[i].spines["top"].set_visible(False)
        axes[i].spines["right"].set_visible(False)

    axes[-1].set_xlabel(xlabel)

    if title is not None:
        fig.suptitle(title, fontweight="bold")

    plt.tight_layout()
    for ax in axes:
        ax.invert_xaxis()
    #plt.show()
    return fig

def plot_scatter_smooth(raw_x, raw_y, smooth_y, timescale, y_label, title=None):
    """
    This function plots the original data against the smoothed data

    Args:
        raw_x (np.ndarray): An array of the raw x data, typically the age of the sample
        raw_y: An array of the raw y data
        smooth_y: An array of the smoothed y data
        timescale (np.ndarray): The timescale associated with the data in xs_norm
        y_label: The label for the y axis
        title (string): The title for the figure
    
    Returns:
        matplotlib.figure.Figure: The generated figure.
    """

    # Initializing the plot
    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot the smoothed data
    ax.plot(timescale, smooth_y, linewidth=2, zorder=1, label= "smoothed data")
    
    # Plot the raw data
    ax.scatter(raw_x, raw_y, s=8, color="red", alpha=0.4, zorder=2, label= "raw data")

    # Add labels to the axes and the title
    ax.set_xlabel("Time (kyr BP)")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    
    # Create the legend
    ax.legend()
    ax.grid(alpha=0.3)

    # Invert the x-axis
    ax.invert_xaxis()

    # Return the generated figure
    return fig

def solve_system(system, times, timescale, initial_conditions=[0.1, 0.8, -0.2]):
    """
    This function solves a system of equations using the solve_ivp methode

    Args:
        system (function name): Which equations the system should use to solve the system
        times: the begining and end of the timescale
        timescale: The time series data associated with the data from the recovered ODEa
        inital_conditions: (list[float]): The initial conditions for the 3 state variables, if not specified the initial conditions will be [0.1, 0.8, -0.2]
    
    Returns:
        scipy.integrate._ivp.ivp.OdeResult: The solution returned by solve_ivp, containing the simulated state variables, time values, and additional solver information.
    """

    df_recovered = solve_ivp(system, times, initial_conditions, t_eval=timescale, method="LSODA")
    return df_recovered

def system_original(t, state):
    """
    This function defines the barry saltzman system

    Args:
        t
        state

    Returns:

    """

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

def system_recovered(t, state):
    '''
    
    '''
    #x,y,z = state
    x,y = state

    # The recovered sindy equations
    """dx = -1.013*x - 1.008*y - 0.209*z
    dy = 0.512*x + 1.284*y - 0.587*z - 0.494*y**2 - 0.797*y**3 + 0.047*z**3
    dz = 0.544*y - 0.430*z + 0.130*y**2 + 0.205*y**3"""

    """dx = -0.200*y + 0.060*z - 0.077*x*z + 0.626*y*z - 0.100*z**3
    dy = -0.294*y + 0.824*y*z
    dz = 1.159*y - 2.445*y*z"""

    """dx = (-0.057*x + 0.441*x**2*y + 0.023*x**2*z - 0.354*x*y**2 + 0.099*y**2*z)
    dy = (0.058*x - 0.081*x**3 - 0.164*x**2*y)
    dz = (0.043*z - 0.186*x**2 - 0.089*y*z + 0.184*x**3 + 0.176*x*y**2)"""
    
    """#weak sindy equations
    dx = 0.929 - 2.263*x - 5.078*y - 9.554*z + 1.858*x**2 + 7.390*x*y + 8.192*y**2 + 84.325*y*z - 3.231*z**2 - 0.469*x**3 - 2.503*x**2*y - 6.433*x*y**2 - 2.851*y**3 - 156.612*y**2*z
    dy = 0.038*x - 0.094*y + 9.365*z - 0.013*x*y - 12.391*x*z - 8.564*y*z + 0.291*z**2 - 0.034*x**3 - 0.077*x**2*y + 0.541*x*y**2 - 0.051*y**3
    dz = -0.110*z + 0.379*x*z + 0.499*z**2 - 0.217*x**2*z"""

    # weak sindy for co2_concentration_clean, ocean_clean, global_ice_volume_clean, lam = 0.03 and alpha = 1.25 with degree 3 library
    """dx = -0.028*x + 0.340*x**2 - 2.810*x*y + 1.343*y*z + 2.867*x**2*y - 1.332*x**2*z + 0.894*x*y**2 + 1.262*x*z**2 - 2.313*y**3 + 3.937*y**2*z - 2.989*y*z**2 - 0.088*z**3
    dy = -0.847*x + 1.384*y + 0.458*x**2 + 2.390*x*z + 0.225*y**2 - 4.157*y*z - 0.466*x**3 - 0.375*x**2*y + 0.052*x*y**2 - 1.709*x*z**2 - 0.300*y**3 + 2.975*y*z**2
    dz = 0.110*x*y - 0.144*x*z**2 - 0.021*y**3 - 0.108*y**2*z + 0.098*z**3"""

    # weak sindy for co2_concentration_clean, ocean_clean, global_ice_volume_clean, lam = 0.03 and alpha = 1.25 with degree 2 library
    """dx = -0.249*x + 0.108*z + 0.136*x**2 + 0.224*x*z + 0.051*y**2 - 0.150*z**2
    dy = 0.304*x - 0.027*y - 0.142*z - 0.172*x**2 - 0.258*x*z - 0.208*y**2 + 0.231*y*z + 0.137*z**2
    dz = -0.111 + 0.099*x + 0.061*y + 0.250*z - 0.109*x**2 - 0.063*y*z - 0.171*z**2"""

    # weak sindy for co2_concentration_clean, ocean_clean, global_ice_volume_clean, lam = 0.01 and alpha = 1.5 with degree 2 library
    """dx = -0.279*x - 0.112*y + 0.152*z + 0.251*x**2 + 0.334*x*y - 0.064*x*z + 0.239*y**2 - 0.266*y*z
    dy = 0.005 + 0.236*x - 0.127*y - 0.182*x**2 - 0.187*x*z - 0.122*y**2 + 0.238*y*z
    dz = -0.130 + 0.255*x - 0.021*y + 0.268*z - 0.258*x**2 - 0.075*x*y + 0.144*y**2 - 0.109*y*z - 0.168*z**2"""

    # weak sindy for co2_concentration_clean, global_ice_volume_clean, lam = 0.01 and alpha = 1.5 with degree 3 library
    #dx = 0.283 - 0.625*x - 0.863*y + 1.680*x*y + 0.672*y**2 + 0.321*x**3 - 0.103*x**2*y - 1.254*x*y**2
    #dy = -0.268*x - 0.101*y + 0.817*x**2 + 1.094*x*y - 0.606*x**3 - 1.230*x**2*y - 0.628*x*y**2 + 0.129*y**3

    # weak sindy for co2_concentration_clean, global_ice_volume_clean, lam = 0.03 and alpha = 1.25 with degree 3 library
    #dx = -0.061*x**2*y + 0.047*y**3
    #dy = 0.086*x - 0.167*x**3 - 0.110*x*y**2 - 0.042*y**3

    # weak sindy for co2_concentration_clean, global_ice_volume_clean, lam = 0.05 and alpha = 1.25 with degree 3 library
    #dx = -0.082*x + 0.075*y**2 + 0.043*x**3 + 0.520*x**2*y - 0.374*x*y**2
    #dy = 0.084*x - 0.121*x**3 - 0.211*x**2*y

    # Weak sindy with milankovich
    #dx = -0.054*x**3 + 0.343*x**2*y - 0.432*x*y**2 + 0.116*y**3
    #dy = -0.027*x**3 + 0.174*x*y**2 - 0.073*y**3

    # Weak sindy with milankovich cycle, lam=0.01, alpha=1.25, library of degree 3 (milankovich stuff was wrong)
    #dx = 0.398*x**3 + 0.519*x**2*y - 1.432*x**2*z - 1.174*x*y**2 + 1.421*x*z**2 - 0.073*y**3 + 0.730*y**2*z - 0.122*y*z**2 - 0.455*z**3
    #dy = 0.225*x**2*y + 0.013*x*y**2 - 0.186*x*z**2 - 0.323*y**2*z + 0.235*y*z**2 + 0.062*z**3
    #dz = -0.011*x**2*z + 0.841*x*y**2 - 0.937*x*y*z + 0.228*x*z**2 + 1.503*y**3 - 3.391*y**2*z + 2.266*y*z**2 - 0.427*z**3

    # Weak sindy with milankovich cycle, lam=0.01, alpha=1.25, library of degree 3 (milankovich stuff was wrong)
    #dx = -0.068*x**3 + 0.098*x**2*y - 1.712*x*y**2 + 1.670*x*y*z - 0.222*x*z**2 - 1.716*y**3 + 4.218*y**2*z - 2.798*y*z**2 + 0.412*z**3
    #dy = 0.030*x**3 - 0.002*x**2*y + 0.042*x*y*z - 0.121*x*z**2 + 0.036*y**3 - 0.610*y**2*z + 0.577*y*z**2 - 0.043*z**3
    #dz = -0.425*x**3 - 1.331*x**2*y + 1.860*x**2*z - 0.181*x*y**2 + 2.651*x*y*z - 2.378*x*z**2 + 1.498*y**3 - 2.308*y**2*z + 0.731*z**3

    # Weak Sindy with Milankovich focing as control for co2_concentration_clean, ocean_clean, global_ice_volume_clean,, lam=0.01, alpha=1.25, library of degree 3
    #dx = -0.131*y**2 - 0.167*x**3 + 0.428*x**2*y - 0.670*x*y**2 + 0.294*x*z**2 - 0.904*y**3 + 2.222*y**2*z - 0.989*y*z**2 - 0.042*z**3
    #dy = 0.111*x**2 + 2.049*x*y - 0.346*x*z - 0.897*y**2 - 0.158*y*z - 1.118*x**3 - 2.445*x**2*y + 3.010*x**2*z + 0.064*x*y*z - 2.133*x*z**2 + 0.317*y**3 + 0.307*y**2*z + 0.582*z**3
    #dz = 0.443*x**2 - 1.740*x*y + 0.878*y**2 - 0.441*x**3 + 0.811*x**2*y + 1.851*x*y**2 - 0.429*x*y*z + 0.246*x*z**2 + 1.996*y**3 - 5.860*y**2*z + 3.395*y*z**2 - 0.619*z**3

    # Just milankovich system normalized for lam = 0 and alpha = 0
    #dx = 0.01099069 - 0.07648081*x - 0.07648081*y + 0.13706249*x**2 + 0.13706249*x*y + 0.13706249*y**2 - 0.07141420*x**3 - 0.07141420*x**2*y - 0.07141420*x*y**2 - 0.07141420*y**3
    #dy = 0.01099069 - 0.07648081*x - 0.07648081*y + 0.13706249*x**2 + 0.13706249*x*y + 0.13706249*y**2 - 0.07141420*x**3 - 0.07141420*x**2*y - 0.07141420*x*y**2 - 0.07141420*y**3

    # Weak sindy, s = 10 000 for ice, s = 50 000 for co2, lam = 0 alpha = 0
    #dx =  1.072 - 2.764 * x - 4.252 * y +  1.813 *x**2 +  8.199 * x * y +  5.062 * y**2 - 0.107 * x**3 - 3.282 * x**2 * y - 5.553 * x * y**2 - 1.737 * y**3
    #dy = -0.734 +  1.742 * x + 3.133 * y - 0.765 * x**2 - 5.808 * x * y - 3.954 * y**2 -0.311 * x**3 +  2.215 * x**2 * y +  4.100 * x * y**2 + 1.458 * y**3
    
    # Weak sindy, s = 50 000 for ice, s = 100 000 for co2, lam = 0 alpha = 0
    #dx = 2.257 - 6.518*x - 6.227*y + 6.133*x**2 + 11.537*x*y + 5.803*y**2 - 1.885*x**3 - 5.122*x**2*y - 5.311*x*y**2 - 1.775*y**3
    #dy = -2.093 + 5.540*x + 6.134*y - 4.494*x**2 - 10.639*x*y - 6.027*y**2 + 1.030*x**3 + 4.281*x**2*y + 5.188*x*y**2 + 1.959*y**3

    # Weak sindy, s = 50 000 for ice, s = 100 000 for co2, lam = 0.1 alpha = 0.05
    #dx = 0.140 - 0.358*x - 0.194*y + 0.199*x**2 + 0.127*y**2 + 0.391*x**2*y
    #dy = -0.238 + 0.439*x + 0.504*y - 0.640*y**2 - 0.212*x**3 - 0.629*x**2*y + 0.327*y**3

    # Weak sindy fpr normalized CO2 and normalized co2 *0.5, s = 100 000, lam = 0.01, alpha = 1.25
    #dx =  0.024 * x +  0.012 * y - 0.045 * x**2 - 0.023 * x * y - 0.011 * y**2
    #dy =  0.015 * x - 0.024 * x**2 - 0.012 * x * y

    dx = -0.046*x**2 - 0.023*x*y - 0.012*y**2 + 0.069*x**3 + 0.034*x**2*y + 0.017*x*y**2
    dy = -0.024*x**2 - 0.012*x*y + 0.036*x**3 + 0.018*x**2*y

    return [dx, dy]
    #return [dx, dy, dz]

def system_recovered_gen_data(t, state):
    # Recovered equations for sindy using generated data with noise
    x,y,z = state

    # Weak sindy lam = 0.1, alpha = 1.2, 0.1 noise
    #dx = -0.921*x - 0.979*y - 0.135*z
    #dy = 0.085*x + 1.220*y - 0.925*z - 0.529*y**2 - 0.861*y**3 + 0.024*z**3
    #dz = 0.525*y - 0.428*z + 0.131*y**2 + 0.222*y**3

    # Weak sindy lam = 0.1, alpha = 1.25,0.1 noise
    #dx = -0.921*x - 0.979*y - 0.135*z
    #dy = 0.079*x + 1.261*y - 0.918*z - 0.553*y**2 - 0.907*y**3
    #dz = 0.525*y - 0.428*z + 0.131*y**2 + 0.222*y**3

    # Weak sindy lam = 0.1, alpha = 1.25,0.2 noise
    #dx = -0.679*x - 0.911*y + 0.060*z
    #dy = -0.105 - 0.267*x + 0.626*y - 1.071*z - 0.117*y**2 - 0.243*y**3
    #dz = 0.572*y - 0.442*z + 0.103*y**2 + 0.161*y**3

    # Weak sindy lam = 0.1, alpha = 1.5,0.2 noise
    #dx = -0.679*x - 0.911*y + 0.060*z
    #dy = 0.750*x - 0.052*z +  0.586*y**3
    #dz = 0.572*y - 0.442*z + 0.103*y**2 + 0.161*y**3

    # Weak sindy lam = 0.1, alpha = 1.25,0.3 noise
    dx = -0.439*x - 0.842*y + 0.253*z
    dy = 0.941*x + 0.024*y + 0.099*z + 0.575*y**3
    dz = 0.708*y - 0.497*z

    # Weak sindy lam = 0.1, alpha = 1.5,0.3 noise
    #dx = 0.220*x - 0.532*y + 0.772*z - 0.093*y**3
    #dy = 0.941*x + 0.024*y + 0.099*z + 0.575*y**3
    #dz = 0.708*y - 0.497*z

    # Weak sindy lam = 0.05, alpha = 1.25,0.3 noise
    #dx = 0.345*x - 0.537*y + 0.940*z - 0.067*y**3 + 0.029*y**2*z - 0.069*z**3
    #dy = -0.529 - 0.657*x + 0.313*y - 1.870*z + 0.309*y**2 + 0.039*y*z + 0.872*z**2 + 0.096*x**2*z - 0.140*y**3 + 0.622*y**2*z + 0.689*z**3
    #dz = -0.443*x + 0.634*y - 0.782*z + 0.007*y**2 - 0.035*y**3 - 0.061*z**3

    # Weak sindy lam = 0.075, alpha = 1.25,0.3 noise
    #dx = 0.220*x - 0.532*y + 0.772*z - 0.093*y**3
    #dy = -0.137 - 0.590*x + 0.401*y - 1.638*z - 0.020*y**2 + 0.334*y*z - 0.184*y**3 + 0.487*y**2*z + 0.238*z**3
    #dz = 0.632*y - 0.458*z + 0.072*y**2 + 0.090*y**3

    # Weak sindy lam = 0.075, alpha = 1,0.3 noise
    #dx = 0.220*x - 0.532*y + 0.772*z - 0.093*y**3
    #dy = -0.137 - 0.590*x + 0.401*y - 1.638*z - 0.020*y**2 + 0.334*y*z - 0.184*y**3 + 0.487*y**2*z + 0.238*z**3
    #dz = 0.632*y - 0.458*z + 0.072*y**2 + 0.090*y**3

    return [dx, dy, dz]

def normal_sindy(xs, lam, alpha, timescale,  library = ps.PolynomialLibrary(degree=3),u=None):
    """
    Recover a system of ordinary differential equations (ODEs) using the standard SINDy algorithm.

    Args:
        xs (np.ndarray): Array containing the data passed to sindy to recover ODEs
        lam (float): Sparsity threshold used by the STLSQ optimizer
        alpha (float): Regularization parameter used by the STLSQ optimizer
        timescale (np.ndarray): One-dimensional array of time values corresponding to each row of xs
        libray:  Feature library containing the candidate functions that SINDy can use to construct the recovered equations. Defaults to a third-degree polynomial library
        u (np.ndarray): External forcing (control input) corresponding to each time point. Defaults to None

    Returns:
        pysindy.SINDy: The fitted SINDy model
    """
        
    # Define an optimizer to use
    opt = ps.STLSQ(threshold = lam, alpha = alpha)
    
    # Define the SINDy model
    model = ps.SINDy(feature_library=library, optimizer=opt)

    # Fit the data to the model
    model.fit(xs, t=timescale, feature_names=["x", "y", "z"], u=u)
    
    # Print the model (prints the recovered equtions)
    model.print()

    print(model.coefficients())

    return model

def weak_sindy(xs, lam, alpha, timescale, library = ps.PolynomialLibrary(degree=3, include_bias=False), u = None):    
    """
    Recover a system of ordinary differential equations (ODEs) using the weak SINDy algorithm.

    Args:
        xs (np.ndarray): Array containing the data passed to sindy to recover ODEs
        lam (float): Sparsity threshold used by the STLSQ optimizer
        alpha (float): Regularization parameter used by the STLSQ optimizer
        timescale (np.ndarray): One-dimensional array of time values corresponding to each row of xs
        libray:  Feature library containing the candidate functions that SINDy can use to construct the recovered equations. Defaults to a third-degree polynomial library
        u (np.ndarray): External forcing (control input) corresponding to each time point. Defaults to None

    Returns:
        pysindy.SINDy: The fitted SINDy model
    """

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
    model.fit(xs, t=timescale, feature_names=["x", "y", "z"],u=u)
    
    # Print the model (prints the recovered equtions)
    model.print()

    #print(model.coefficients())

    #return model

def prior_knowledge_sindy(xs, timescale, library = ps.PolynomialLibrary(degree=3, include_bias=False)):
    """
    Recover a system of ordinary differential equations (ODEs) using the SINDy algorithm with underlying knowledge of the underlying behaviours.

    Args:
        xs (np.ndarray): Array containing the data passed to sindy to recover ODEs
        timescale (np.ndarray): One-dimensional array of time values corresponding to each row of xs
        libray:  Feature library containing the candidate functions that SINDy can use to construct the recovered equations. Defaults to a third-degree polynomial library
       
    Returns:
        pysindy.SINDy: The fitted SINDy model
    """

    ode_library = ps.WeakPDELibrary(
        function_library= library,
        spatiotemporal_grid=timescale,
        include_bias=True)
    
    n_targets = 3
    n_features = 20

    allowed = {
        0: [1, 2, 3],        # x' allows x, y, z
        1: [2, 3, 7, 16],    # y' allows y, z, y^2, y^3
        2: [1, 3]            # z' allows x, z
    }

    constraint_rows = []
    constraint_rhs = []

    for target in range(n_targets):
        for feature in range(n_features):
            if feature not in allowed[target]:
                row = np.zeros(n_targets * n_features)

                # feature-order indexing
                index = feature * n_targets + target

                row[index] = 1
                constraint_rows.append(row)
                constraint_rhs.append(0)

    constraint_lhs = np.asarray(constraint_rows)
    constraint_rhs = np.asarray(constraint_rhs)

    #opt = ps.SR3(reg_weight_lam=lam,relax_coeff_nu=alpha, max_iter=10000,regularizer="L1")
    
    opt = ps.ConstrainedSR3(
    reg_weight_lam=0.00001,
    regularizer="L2",
    #relax_coeff_nu=0.01,
    max_iter=1000,
    constraint_lhs=constraint_lhs,
    constraint_rhs=constraint_rhs,
    equality_constraints=True,
    constraint_order="feature"
    )
    
    # Define the SINDy model
    model = ps.SINDy(feature_library=ode_library, optimizer=opt)

    # Fit the data to the model
    model.fit(xs, t=timescale, feature_names=["x", "y", "z"])

    model.print()

    return model

def rmse (start, end):
    """
    This function calculates the root mean square error between the original data and the recovered ODEs

    Args:
        start (int): The age at which the timeseries is starting
        end (int): The age at which the timeseries is ending
    
    Returns:
        list[float]: a list of floats containing the RMSE for each column
    """

    # get the original data
    sol_clean = solve_system(system_original,[start,end], np.arange(start, end + 0.01,0.01)).y.T
    
    # get the recovered data
    sol_recovered = solve_system(system_recovered,[start, end], np.arange(start, end + 0.01,0.01)).y.T

    # calculate the rmse for each column
    rmse_x = np.sqrt(mean_squared_error(sol_clean[:, 0], sol_recovered[:, 0]))
    rmse_y = np.sqrt(mean_squared_error(sol_clean[:, 1], sol_recovered[:, 1]))
    rmse_z = np.sqrt(mean_squared_error(sol_clean[:, 2], sol_recovered[:, 2]))

    return [float(rmse_x), float(rmse_y), float(rmse_z)]

def milankovich_forcing(timescale):
    """
    This function produces timeseries data derived from the milankovich forcing

    Args:
        timescale (np.ndarray): The timescale for which the insolation data will be created
    
    Returns:
        np.ndarray: The insolation data
    """

    # import the palisol library from r
    palinsol = importr("palinsol")

    # conver the timescale to a value that r can understand
    timescale_r = FloatVector(timescale)
    
    # generate the insolation data at latitude N65
    r("""
    insolation_function <- function(times) {
        latitude <- 65 * pi / 180

        result <- sapply(times, function(tt) {
            Insol(
                orbit = ber78(tt),
                lon = pi / 2,
                lat = latitude
            )
        })

        return(result)
    }
    """)
    
    # evaluate the insolatiion data on the provided timescale
    isl_r = r["insolation_function"](timescale_r)

    # Convert back to NumPy and return the insolation data
    return np.array(isl_r)

def main():
    # Get the data from the datasets
    ice_vol, co2, temp = data_cleaning()
    
    # Define a common timescale for the datasets going from 0 to 400 ka
    dt = 0.01
    common_timescale = np.arange(0, 401, dt)
    
    # Interpolate the data to a common timsecale
    #co2_inter = interpolate(common_timescale, co2, "Age", "CO2")
    #deep_temp_inter = interpolate(common_timescale, temp, "Age", "Ocean Temp")
    #ice_vol_inter = interpolate(common_timescale, ice_vol, "Age", "Ice_Volume")

    ice_vol_smooth = smoothing_spline(common_timescale, ice_vol, "Age", "Ice_Volume",50000)
    co2_smooth = smoothing_spline(common_timescale, co2, "Age", "CO2",100000)
    deep_temp_smooth = smoothing_spline(common_timescale, temp, "Age", "Ocean Temp",0)


    # Combine the data into one dataframe
    #xs = pl.DataFrame({"Age": common_timescale["Age"], "CO2": co2_inter, "temp": deep_temp_inter, "ice_vol": ice_vol_inter})
    xs = pl.DataFrame({"Age": common_timescale, "ice_vol": ice_vol_smooth, "CO2": co2_smooth, "temp": deep_temp_smooth})   

    # normalize the data
    xs_norm = normalize_data(xs, "ice_vol", "CO2", "temp")

    # Get a datafram with only the ice extent and CO2 concentration data
    #xs_norm_ice_co2 = xs_norm.select("x_norm", "y_norm")

    # Convert the dataframe to a numpy array to ensure compatibility with sindy librart
    #xs_norm = xs_norm.to_numpy()
    #xs_norm_ice_co2 = xs_norm_ice_co2.to_numpy()

    # Generating synthethic data
    #data, timescale = generate_data_with_noise(0.3)
    #data, timescale = generate_data()

    #milank = milankovich_forcing(-common_timescale*1000)

    xs_test = pl.DataFrame({"co2": xs_norm.select("co2_norm"), "co2_half": xs_norm.select("co2_norm") * 0.5})
    xs_test = xs_test.to_numpy()

    #print(xs_test.shape)
    #print(xs_test[:5])
    #print(np.min(xs_test, axis=0))
    #print(np.max(xs_test, axis=0))

    #xs_test_norm = normalize_data(xs_test, "milank", "milank_half")
    #xs_test_norm = xs_test_norm.to_numpy()

    #print("normal sindy")
    #normal_sindy(xs_norm, lam=0.01, alpha=1.25, timescale=timescale)
    #normal_sindy(xs_test_norm, lam=0, alpha=0, timescale=common_timescale)


    xs_test = xs_test[::-1, :]
    common_timescale = common_timescale[::-1]
    timescale = np.arange(0, len(xs_test)) * dt

    print("weak sindy")
    #weak_sindy(xs_norm_ice_co2, lam=0.1, alpha=0.05, timescale=common_timescale)
    #weak_sindy(xs_test_norm, lam=0, alpha=0, timescale=common_timescale)
    weak_sindy(xs_test, lam=0.01, alpha=0, timescale=timescale)

    #print("constrained SINDy")
    #prior_knowledge_sindy(xs_norm, lam=0.04, alpha=1, timescale=timescale)

    #print("RMSE 0-30")
    #print(rmse(0,30))
    #print("RMSE 0-100")
    #print(rmse(0,100))
    #plot_forward_in_time()
    
    #fig1 = plot_graph(timescale,milank, "time", "milank", "plottng milankovich")

    #fig1 = plot_data_and_derivative(xs_norm,timescale)
    #fig2 = plot_stacked_timeseries(timescale,[xs_norm[:,0],xs_norm[:,1], xs_norm[:,2]],labels=["ice_vol", "co2", "ocean temp"], title="datasets plotted through time")
    #fig3 = plot_scatter_smooth(ice_vol["Age"].to_numpy().ravel(), ice_vol["Ice_Volume"].to_numpy().ravel(), ice_vol_smooth, common_timescale, "Ice volume", "Raw vs smoothed ice volume data, s = 50000")
    #fig4 = plot_scatter_smooth(co2["Age"].to_numpy().ravel(), co2["CO2"].to_numpy().ravel(), co2_smooth, common_timescale, "CO2 concentration", "Raw vs smoothed co2 concentration data, s = 100000")
    #fig5 = plot_scatter_smooth(temp["Age"].to_numpy().ravel(), temp["Ocean Temp"].to_numpy().ravel(), deep_temp_smooth, timescale, "Ocean temp", "Raw vs smoothed ocean temperature data")

    #fig = plot_forward_in_time_dataset(xs_norm_ice_co2, dt, xs_norm_ice_co2[0,:], common_timescale, "Recovered SINDy ODE vs real dataset, lam = 0.01, alpha = 1.5, library = degree 3")

    #fig = plot_forward_in_time_dataset(xs_norm, dt, xs_norm[0,:], common_timescale)

    fig = plot_forward_in_time_dataset(xs_test, common_timescale)

    #fig = plot_forward_in_time_noisy_data(data)

    plt.show()

main()