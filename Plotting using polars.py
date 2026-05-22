import polars as pl
import matplotlib.pyplot as plt

def two_yaxis_graph(df1, df1_xaxis, df1_yaxis, df1_label, df2, df2_xaxis, df2_yaxis, df2_label, title):
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Left axis: ocean temperature
    ax1.plot(
        df1[df1_xaxis],
        df1[df1_yaxis],
        color="red",
        label=df1_label
    )

    # For Epica Dome C 700 000 years ocean temp data set
    """# Left axis: ocean temperature
    ax1.plot(
        ocean_temp_v2_clean["Age"],
        ocean_temp_v2_clean["MOT_KrN2"],
        color="red",
        marker= "o",
        label="Ocean temperature"
    )"""

    ax1.set_xlabel(df1_xaxis)
    ax1.set_ylabel(df1_yaxis)

    # Reverse x-axis for paleoclimate convention
    ax1.invert_xaxis()

    # Right axis: CO2
    ax2 = ax1.twinx()

    ax2.plot(
        df2[df2_xaxis],
        df2[df2_yaxis],
        linestyle="--",
        label= df2_label
    )

    ax2.set_ylabel(df2_yaxis)

    # Title
    plt.title(title)

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()

    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    #plt.show()
    return fig

def plot_graph(df, df_xaxis, df_yaxis, df_title):
    fig, ax = plt.subplots()
    ax.plot(
        df[df_xaxis],
        df[df_yaxis]
    )
    ax.set_title(df_title)
    ax.set_xlabel(df_xaxis)
    ax.invert_xaxis()
    ax.set_ylabel(df_yaxis)

    #plt.show()
    return fig

def main():
    ocean_temp_original = pl.read_excel(
        source = "grimmer2024-mot-ohc-spline.xlsx", 
        read_options={"header_row":106},
        infer_schema_length=0
        )

    co2_concentration_original = pl.read_excel(
        source = "antarctica2015co2.xls", 
        sheet_name="CO2 Composite",
        read_options={"header_row":14}
        )

    ocean_temp_v2_original = pl.read_excel(
        source = "EPICA Dome C 700 000 Year Noble Gas Data and Mean Ocean Temp reconstructions.xlsx", 
        read_options={"header_row":116}
        )
    
    icecore_data_co2_composite_original = pl.read_excel(
        source = "aicc2012icecore-data.xls",
        sheet_name="CO2 composite",
        read_options={"header_row":1}
    )

    ice_accumulation_original = pl.read_excel(
        source="aicc2012official.xls",
        sheet_name="Vostok",
        read_options={"header_row":1}
    )

    air_temp_original = pl.read_csv(
    "friedrich2016temp.txt",
    separator="\t",
    comment_prefix="#"
    )

    ocean_temp_clean = ocean_temp_original.with_columns(pl.all().cast(pl.Float64)).select(pl.col("age_ka_spline").alias("Age") * 1000, pl.col("MOT_spline").alias("Variation in ocean temp"))

    #For grimmer dataset
    co2_concentration_clean = co2_concentration_original.select(pl.col("Gasage (yr BP)").alias("Gas Age"), "CO2 (ppmv)").filter(pl.col("Gas Age").is_between(350.168,798569))

    # For Epica Dome C 700 000 years dataset
    #co2_concentration_clean = co2_concentration_original.select(pl.col("Gasage (yr BP)").alias("Age"), "CO2 (ppmv)").filter(pl.col("Age").is_between(1001,699622))

    ocean_temp_v2_clean = ocean_temp_v2_original.select(pl.col("gas_age_ky").alias("Age") * 1000, "MOT_KrN2")

    icecore_data_co2_composite_clean = icecore_data_co2_composite_original.select(pl.col("AICC2012 gas age (a B1950)").alias("Gas Age"), "CO2 (ppmv)").filter(pl.col("Gas Age").is_between(3168.82, 798569))

    ice_accumulation_clean = ice_accumulation_original.select(pl.col("Ice age (a B1950)").alias("Ice Age"), pl.col("Gas age (a B1950)").alias("Gas Age"), pl.col("Accu. Rate (m ie)").alias("Accumulation Rate")).filter(pl.col("Gas Age").is_between(0,500000))

    air_temep_clean = air_temp_original.select(pl.col("age_calkaBP").alias("Age") * 1000, pl.col("tempanom-pxy").alias("Anomalies from proxy data")).filter(pl.col("Anomalies from proxy data").is_between(-900,900))

    print(ocean_temp_clean.head())
    print(co2_concentration_clean.head())
    print(ocean_temp_v2_clean.head())
    print(icecore_data_co2_composite_clean.head())
    print(ice_accumulation_clean.head())
    print(air_temep_clean.head())

    #figure1 = plot_graph(icecore_data_co2_composite_clean, "Gas Age", "CO2 (ppmv)", "Gas age vs CO2 concentration")
    #figure2 = plot_graph(ice_accumulation_clean, "Ice Age", "Accumulation Rate", "Ice Age vs Accumulation Rate for Vostok Ice core")
    figure3 = two_yaxis_graph(air_temep_clean, "Age", "Anomalies from proxy data", "Air Temp Anomalies", ocean_temp_v2_clean, "Age", "MOT_KrN2", "Water temp anomally", "Comparing Air and Ocean anomalies")

    plt.show()

main()