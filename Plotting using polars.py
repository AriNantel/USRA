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
    #ax1.invert_yaxis()

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
        sheet_name="EDC",
        read_options={"header_row":1}
    )

    air_temp_original = pl.read_csv(
    "friedrich2016temp.txt",
    separator="\t",
    comment_prefix="#"
    )

    accumulation_EPICA_C_original = pl.read_excel(
        source="Accum Rate EPICA C.xlsx",
        read_options={"header_row":127}
    )

    melt_rate_original = pl.read_excel(
        source="pcmeltco_50yr.xlsx",
        read_options={"header_row":17}
    )

    d18o_original = pl.read_excel(
        source="lisiecki2005-d18o-stack-noaa.xlsx",
        read_options={"header_row":126}
    )

    # Grimmer dataset
    ocean_temp_clean = (
        ocean_temp_original
        .with_columns(pl.all().cast(pl.Float64))
        .select(
            (pl.col("age_ka_spline")*1000).alias("Age"),
            pl.col("MOT_spline").alias("Variation in ocean temp"))
        .filter(
            ~pl.any_horizontal(pl.all() == -999) &
            pl.col("Age").is_between(0, 10500))
    )

    # Antartica 2015 co2 set
    co2_concentration_clean = (
        co2_concentration_original
        .select(
            pl.col("Gasage (yr BP)").alias("Gas Age"),
            "CO2 (ppmv)")
        .filter(
            pl.col("Gas Age").is_between(350.168,798569))
    )

    # EPICA DOM C Ocean temp dataset
    ocean_temp_v2_clean = (
        ocean_temp_v2_original
        .select(
            (pl.col("gas_age_ky")*1000).alias("Age"),
            "MOT_KrN2")
        .filter(
            ~pl.any_horizontal(pl.all() == -999),
            pl.col("Age").is_between(0, 10500))
    )
    
    # aicc2012 ice core dataset
    icecore_data_co2_composite_clean = (
        icecore_data_co2_composite_original
        .select(
            pl.col("AICC2012 gas age (a B1950)").alias("Gas Age"),
            "CO2 (ppmv)")
        .filter(
            pl.col("Gas Age").is_between(3168.82, 798569))
    )

    # aicc2012 official dataset
    ice_accumulation_clean = (
        ice_accumulation_original
        .select(
            pl.col("Ice age (a B1950)").alias("Ice Age"),
            pl.col("Gas age (a B1950)").alias("Gas Age"),
            pl.col("Accu. Rate (m ie)").alias("Accumulation Rate"))
        .filter(
            pl.col("Ice Age").is_between(0,40000) & 
           ~pl.any_horizontal(pl.all() == -1))
    )

    # friedrich 2016 temp dataset
    air_temp_clean = (
        air_temp_original
        .select(
            (pl.col("age_calkaBP")*1000).alias("Age"),
            pl.col("tempanom-pxy").alias("Anomalies from proxy data"))
        .filter(
            ~pl.any_horizontal(pl.all() == -999))
    )

    # Accum Rate EPICA C dataset
    accumulation_EPICA_C_clean = (
        accumulation_EPICA_C_original
        .select(
            pl.col("age_yrBP1950").alias("Age"),
            pl.col("A_EDC").alias("Accumulation Rate"))
    )

    # pcmeltco_50yry dataset
    melt_rate_clean = (
        melt_rate_original
        .select(
            (1950 - pl.col("Year")).alias("Age"),
            "Melt percentage")
    )

    # lisiecki 2005 dataset
    d18o_clean = (
        d18o_original
        .select(
            (pl.col("age_calkaBP")*1000).alias("Age"),
            "d18O_benthic")
        .filter(
            pl.col("Age").is_between(3168.82, 798569))
    )

    #print(ocean_temp_clean.head())
    #print(co2_concentration_clean.head())
    #print(ocean_temp_v2_clean.head())
    #print(icecore_data_co2_composite_clean.head())
    #print(ice_accumulation_clean.head())
    #print(air_temep_clean.head())
    #print(accumulation_EPICA_C_clean.head())
    #print(d18o_clean.head())

    #figure1 = plot_graph(icecore_data_co2_composite_clean, "Gas Age", "CO2 (ppmv)", "Gas age vs CO2 concentration")
    #figure2 = plot_graph(ice_accumulation_clean, "Ice Age", "Accumulation Rate", "Ice Age vs Accumulation Rate for Vostok Ice core")
    figure3 = two_yaxis_graph(melt_rate_clean, "Age", "Melt percentage", "melt percentage", ocean_temp_clean, "Age", "Variation in ocean temp", "Ocean Temp anomalies", "Comparing melt percentage to ocean temperature anomalies")

    plt.show()

main()