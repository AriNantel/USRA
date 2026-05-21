import polars as pl
import matplotlib.pyplot as plt


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

ocean_temp_clean = ocean_temp_original.with_columns(pl.all().cast(pl.Float64)).select(pl.col("age_ka_spline").alias("Age") * 1000, pl.col("MOT_spline").alias("Variation in ocean temp"))

co2_concentration_clean = co2_concentration_original.select(pl.col("Gasage (yr BP)").alias("Age"), "CO2 (ppmv)").filter(pl.col("Age").is_between(500,343500))
#co2_concentration_clean = co2_concentration_original.select(pl.col("Gasage (yr BP)").alias("Age"), "CO2 (ppmv)").filter(pl.col("Age").is_between(500,699622))

ocean_temp_v2_clean = ocean_temp_v2_original.select(pl.col("gas_age_ky").alias("Age") * 1000, "MOT_KrN2")

print(ocean_temp_clean.head())
print(co2_concentration_clean.head())
print(ocean_temp_v2_clean.head())

fig, ax1 = plt.subplots(figsize=(12, 6))

# Left axis: ocean temperature
ax1.plot(
    ocean_temp_clean["Age"],
    ocean_temp_clean["Variation in ocean temp"],
    color="red",
    label="Ocean temperature"
)

"""# Left axis: ocean temperature
ax1.plot(
    ocean_temp_v2_clean["Age"],
    ocean_temp_v2_clean["MOT_KrN2"],
    color="red",
    marker= "o",
    label="Ocean temperature"
)"""

ax1.set_xlabel("Age (yr BP)")
ax1.set_ylabel("Ocean temperature anomaly")

# Reverse x-axis for paleoclimate convention
ax1.invert_xaxis()

# Right axis: CO2
ax2 = ax1.twinx()

ax2.plot(
    co2_concentration_clean["Age"],
    co2_concentration_clean["CO2 (ppmv)"],
    linestyle="--",
    label="CO2"
)

ax2.set_ylabel("CO2 (ppmv)")

# Title
plt.title("Ocean Temperature and Atmospheric CO2")

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

plt.show()