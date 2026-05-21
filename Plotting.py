import pandas as pd
import matplotlib.pyplot as plt

# Read Excel data using pandas
data = pd.read_excel("antarctica2015co2.xls", sheet_name="CO2 Composite", header=14)

# Use bracket access for columns with spaces/parentheses
x = data['Gasage (yr BP)']
y = data['CO2 (ppmv)']

plt.figure(1)
plt.plot(x, y)
plt.grid(True)
plt.xlabel('Gas age (yr BP)')
plt.ylabel('CO2 Concentration (ppmv)')
plt.title('CO2 Composite data over the years')
plt.tight_layout()

# Read Excel data using pandas
data2 = pd.read_excel("grimmer2024-mot-ohc-spline.xlsx", header=106)

# Use bracket access for columns with spaces/parentheses
x2 = data2['age_ka_spline']
y2 = data2['MOT_spline']

plt.figure(2)
plt.plot(x2, y2)
plt.grid(True)
plt.xlabel('age (ka yr)')
plt.ylabel('Sea Water temperature (C)')
plt.title('Evolution of sea water temperature')
plt.tight_layout()

# Read Excel data using pandas
data3 = pd.read_excel("EPICA Dome C 700 000 Year Noble Gas Data and Mean Ocean Temp reconstructions.xlsx", header=116)

# Use bracket access for columns with spaces/parentheses
x3 = data3['gas_age_ky']
y3 = data3['MOT_KrN2']

plt.figure(3)
plt.plot(x3, y3)
plt.grid(True)
plt.xlabel('age (ka yr)')
plt.ylabel('Reconstructed avergae ocean temperature (C)')
plt.title('Evolution of sea water temperature')
plt.tight_layout()

plt.show()

