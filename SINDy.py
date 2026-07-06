import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import pysindy as ps
from scipy.integrate import solve_ivp
from sklearn.metrics import mean_squared_error
from scipy.interpolate import UnivariateSpline

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

    ocean_150_original = pl.read_excel(
        source="sosdian2009.xls",
        sheet_name="607 MgCa BWT",
        read_options={"header_row":5}
    )

    ocean_10_original = pl.read_excel(
        source="sosdian2009.xls",
        sheet_name="23-24 MgCa BWT",
        read_options={"header_row":5}
    )

    ocean_original = pl.concat([ocean_10_original, ocean_150_original], how="vertical")

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
    return co2_concentration_clean, ocean_clean, global_ice_volume_clean

def interpolate(common_timescale, df, x_var, y_var):

    # Interpolate a dataset so it is on a common timescale
    return np.interp(common_timescale["Age"], df[x_var], df[y_var])

def smoothing_spline(common_timescale, df, x_var, y_var, s):
    x = df[x_var].to_numpy()
    y = df[y_var].to_numpy()

    smooth = UnivariateSpline(x, y)
    smooth.set_smoothing_factor(s)

    #return smooth
    return smooth(common_timescale).ravel()

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
    timespan = [start, end]
    
    # Generate data for the initial condition given
    xs = solve_system(system_original, timespan, times=times, initial_conditions=initial_conditions).y.T

    # Return the data and the timeseries
    return xs, times

def generate_data_with_noise(noise_level, dt = 0.01, initial_conditions = [0.1, 0.8, -0.2], start=0, end=100):
    # generate clean data
    xs, times = generate_data(dt, initial_conditions, start, end)

    # Create gaussian noise based on the standard deviation of the data in xs.
    noise = noise_level * np.std(xs, axis=0) * np.random.randn(*xs.shape)
    
    # Return the noisy data and the timeseries
    return xs+noise, times

def plot_graph(x, y, df_xaxis, df_yaxis, df_title):
    fig, ax = plt.subplots()

    ax.plot(x, y)
    ax.invert_xaxis()


    ax.set_xlabel(df_xaxis)
    ax.set_ylabel(df_yaxis)

def plot_data_and_derivative(X, dt):
    # CHAT WROTE THIS
    fd = ps.SmoothedFiniteDifference()
    Xdot = fd._differentiate(X, t=dt)

    fig, ax = plt.subplots(3, 2, figsize=(12, 8))

    names = ["ice_vol", "co2", "ocean temp"]

    for i in range(3):
        # Data
        ax[i, 0].plot(X[:, i])
        ax[i, 0].set_title(names[i])

        # Derivative
        ax[i, 1].plot(Xdot[:, i])
        ax[i, 1].set_title(f"{names[i]}'")

    plt.tight_layout()

    for row in ax:
        for a in row:
            a.invert_xaxis()

    #plt.show()
    return fig

def plot_forward_in_time_clean_data(dt=0.01,initial_conditions=[0.1, 0.8, -0.2], start = 0, end = 30):
    # Get the time scale for the plot
    times = np.arange(start, end + dt, dt)
    timespan = [start, end]

    # Transpose the two solution arrays so they are the right shape
    df_clean = solve_system(system_original,timespan, times=times, initial_conditions=initial_conditions).y.T
    df_recovered = solve_system(system_recovered, timespan, times = times, initial_conditions=initial_conditions).y.T

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
    plt.title("SINDy ODE extraction for datasets", fontweight="bold")
    plt.gca().invert_xaxis()

    plt.legend(loc="upper center", ncol=3, frameon=False)

    plt.show()

def plot_forward_in_time_dataset(xs_norm, dt=0.01,initial_conditions=[0.1, 0.8, -0.2], timescale = None, title = "Recovered SINDy ODE vs real dataset"):
    # Get the time scale for the plot
    #times = np.arange(timespan[0], timespan[-1] + dt, dt)

    # Solve recovered equations, the solution array needs to be transposed so they are the right shape
    df_recovered = solve_system(system_recovered, [timescale[0], timescale[-1]], times = timescale, initial_conditions=xs_norm[0,:]).y.T

    plt.figure(figsize=(10, 8))

    # Recovered SINDy lines
    plt.plot(timescale, df_recovered[:, 0], color="orange", linestyle=":", linewidth=4, label="x, y, z from Weak SINDy")
    plt.plot(timescale, df_recovered[:, 1], color="orange", linestyle=":", linewidth=4)
    plt.plot(timescale, df_recovered[:, 2], color="orange", linestyle=":", linewidth=4)

    # Original lines
    plt.plot(timescale, xs_norm[:, 0], color="black", linewidth=1.5, label="x = Ice")
    plt.plot(timescale, xs_norm[:, 1], color="red", linewidth=1.5, label="y = CO2")
    plt.plot(timescale, xs_norm[:, 2], color="blue", linewidth=1.5, label="z = T_Ocean")

    plt.axhline(0, color="black", linestyle="--", linewidth=1)

    plt.xlim(timescale[0], timescale[-1])
    #plt.xlim(0, 30)

    plt.xlabel("Time (in ky BP)")
    plt.ylabel("Normalized values")
    plt.title(title, fontweight="bold")

    plt.gca().invert_xaxis()
    plt.legend(loc="upper center", ncol=3, frameon=False)

    plt.show()

def plot_stacked_timeseries(
    time,
    series,
    labels,
    colors=None,
    xlabel="Time (kyr BP)",
    title=None,
    figsize=(9, 7)
):
    """
    time: 1D array of time values
    series: list of 1D arrays, one per variable
    labels: list of y-axis labels
    colors: optional list of colors
    """

    n = len(series)

    if colors is None:
        colors = ["tab:blue", "tab:red", "tab:green", "black"]

    fig, axes = plt.subplots(
        n,
        1,
        figsize=figsize,
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
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(
        timescale,
        smooth_y,
        linewidth=2,
        zorder=1,
        label= "smoothed data"
    )
    
    ax.scatter(
        raw_x,
        raw_y,
        s=8,
        color="red",
        alpha=0.4,
        zorder=2,
        label= "raw data"
    )

    ax.set_xlabel("Time (kyr BP)")
    ax.set_ylabel(y_label)

    if title is not None:
        ax.set_title(title)

    ax.legend()
    ax.grid(alpha=0.3)
    ax.invert_xaxis()

    return fig

def solve_system(system, timespan, times, initial_conditions=[0.1, 0.8, -0.2]):
    df_recovered = solve_ivp(system, timespan, initial_conditions, t_eval=times, method="LSODA")
    return df_recovered

def system_original(t, state):
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
    x,y,z = state
    #x,y = state

    # The recovered sindy equations
    """dx = -1.013*x - 1.008*y - 0.209*z
    dy = 0.512*x + 1.284*y - 0.587*z - 0.494*y**2 - 0.797*y**3 + 0.047*z**3
    dz = 0.544*y - 0.430*z + 0.130*y**2 + 0.205*y**3"""

    """dx = -0.200*y + 0.060*z - 0.077*x*z + 0.626*y*z - 0.100*z**3
    dy = -0.294*y + 0.824*y*z
    dz = 1.159*y - 2.445*y*z"""

    """dx = (-0.057*x
          + 0.441*x**2*y
          + 0.023*x**2*z
          - 0.354*x*y**2
          + 0.099*y**2*z)

    dy = (0.058*x
          - 0.081*x**3
          - 0.164*x**2*y)

    dz = (0.043*z
          - 0.186*x**2
          - 0.089*y*z
          + 0.184*x**3
          + 0.176*x*y**2)"""
    
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
    dx = -0.279*x - 0.112*y + 0.152*z + 0.251*x**2 + 0.334*x*y - 0.064*x*z + 0.239*y**2 - 0.266*y*z

    dy = 0.005 + 0.236*x - 0.127*y - 0.182*x**2 - 0.187*x*z - 0.122*y**2 + 0.238*y*z

    dz = -0.130 + 0.255*x - 0.021*y + 0.268*z - 0.258*x**2 - 0.075*x*y + 0.144*y**2 - 0.109*y*z - 0.168*z**2

    # weak sindy for co2_concentration_clean, global_ice_volume_clean, lam = 0.01 and alpha = 1.5 with degree 3 library
    #dx = 0.283 - 0.625*x - 0.863*y + 1.680*x*y + 0.672*y**2 + 0.321*x**3 - 0.103*x**2*y - 1.254*x*y**2
    #dy = -0.268*x - 0.101*y + 0.817*x**2 + 1.094*x*y - 0.606*x**3 - 1.230*x**2*y - 0.628*x*y**2 + 0.129*y**3

    # weak sindy for co2_concentration_clean, global_ice_volume_clean, lam = 0.03 and alpha = 1.25 with degree 3 library
    #dx = -0.061*x**2*y + 0.047*y**3
    #dy = 0.086*x - 0.167*x**3 - 0.110*x*y**2 - 0.042*y**3

    # weak sindy for co2_concentration_clean, global_ice_volume_clean, lam = 0.05 and alpha = 1.25 with degree 3 library
    #dx = -0.082*x + 0.075*y**2 + 0.043*x**3 + 0.520*x**2*y - 0.374*x*y**2
    #dy = 0.084*x - 0.121*x**3 - 0.211*x**2*y

    #return [dx, dy]
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

def prior_knowledge_sindy(xs, timescale, library = ps.PolynomialLibrary(degree=3, include_bias=False)):

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

def rmse (start, end):

    sol_clean = solve_system(system_original,[start,end], np.arange(start, end + 0.01,0.01)).y.T
    sol_recovered = solve_system(system_recovered,[start, end], np.arange(start, end + 0.01,0.01)).y.T

    rmse_x = np.sqrt(mean_squared_error(sol_clean[:, 0], sol_recovered[:, 0]))
    rmse_y = np.sqrt(mean_squared_error(sol_clean[:, 1], sol_recovered[:, 1]))
    rmse_z = np.sqrt(mean_squared_error(sol_clean[:, 2], sol_recovered[:, 2]))

    return [float(rmse_x), float(rmse_y), float(rmse_z)]

def main():
    # Get the data from the datasets
    co2, temp, ice_vol = data_cleaning()
    
    # Define a common timescale for the datasets going from 0 to 400 000 with jumps of 5000 years
    dt = 1
    #common_timescale = pl.DataFrame({"Age": range(0, 401, dt)})
    common_timescale = np.arange(11, 401, 0.01)
    
    # Interpolate the data to a common timsecale
    #co2_inter = interpolate(common_timescale, co2, "Age", "CO2")
    #deep_temp_inter = interpolate(common_timescale, temp, "Age", "Ocean Temp")
    #ice_vol_inter = interpolate(common_timescale, ice_vol, "Age", "Ice_Volume")

    ice_vol_smooth = smoothing_spline(common_timescale, ice_vol, "Age", "Ice_Volume",0)
    co2_smooth = smoothing_spline(common_timescale, co2, "Age", "CO2",0)
    deep_temp_smooth = smoothing_spline(common_timescale, temp, "Age", "Ocean Temp",0)


    # Combine the data into one dataframe
    #xs = pl.DataFrame({"Age": common_timescale["Age"], "CO2": co2_inter, "temp": deep_temp_inter, "ice_vol": ice_vol_inter})
    xs = pl.DataFrame({"Age": common_timescale, "ice_vol": ice_vol_smooth, "CO2": co2_smooth, "temp": deep_temp_smooth})   

    # normalize the data
    xs_norm = normalize_data(xs, "ice_vol", "CO2", "temp")

    xs_norm_ice_co2 = xs_norm.select("x_norm", "y_norm")

    # Convert the dataframe to a numpy array to ensure compatibility with sindy librart
    xs_norm = xs_norm.to_numpy()

    xs_norm_ice_co2 = xs_norm_ice_co2.to_numpy()

    #data, timescale = generate_data_with_noise(0.1)
    #data, timescale = generate_data()

    timescale = xs["Age"].to_numpy()

    #print("normal sindy")
    #normal_sindy(xs_norm, lam=0.01, alpha=1.25, timescale=timescale)

    #print("weak sindy")
    #weak_sindy(xs_norm_ice_co2, lam=0.05, alpha=1.25, timescale=timescale)

    #print("constrained SINDy")
    #prior_knowledge_sindy(xs_norm, lam=0.04, alpha=1, timescale=timescale)

    #print("RMSE 0-30")
    #print(rmse(0,30))
    #print("RMSE 0-100")
    #print(rmse(0,100))
    #plot_forward_in_time()
    
    #fig1 = plot_graph(timescale,xs_norm[:,2], "time", "temp", "ocean temp normalized")

    #fig1 = plot_data_and_derivative(xs_norm,timescale)
    #fig2 = plot_stacked_timeseries(timescale,[xs_norm[:,0],xs_norm[:,1], xs_norm[:,2]],labels=["ice_vol", "co2", "ocean temp"], title="datasets plotted through time")
    #fig3 = plot_scatter_smooth(ice_vol["Age"].to_numpy().ravel(), ice_vol["Ice_Volume"].to_numpy().ravel(), ice_vol_smooth, timescale, "Ice volume", "Raw vs smoothed ice volume data")
    #fig4 = plot_scatter_smooth(co2["Age"].to_numpy().ravel(), co2["CO2"].to_numpy().ravel(), co2_smooth, timescale, "CO2 concentration", "Raw vs smoothed co2 concentration data")
    #fig5 = plot_scatter_smooth(temp["Age"].to_numpy().ravel(), temp["Ocean Temp"].to_numpy().ravel(), deep_temp_smooth, timescale, "Ocean temp", "Raw vs smoothed ocean temperature data")

    #fig = plot_forward_in_time_dataset(xs_norm_ice_co2, dt, xs_norm_ice_co2[0,:], common_timescale, "Recovered SINDy ODE vs real dataset, lam = 0.01, alpha = 1.5, library = degree 3")

    fig = plot_forward_in_time_dataset(xs_norm, dt, xs_norm[0,:], common_timescale, "Recovered SINDy ODE vs real dataset, lam = 0.01, alpha = 1.5, library = degree 2")

    plt.show()

main()