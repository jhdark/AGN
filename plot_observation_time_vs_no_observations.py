from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import csv

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

maintable = Table.read(
    "COSMOS_PSdetection_filter2.fits", format="fits"
)  # Reads in the table directly from the .fits file

# maintable.write("COSMOS_PSdetection_filter2.csv", format="ascii.csv", overwrite=True)

unique_obervation_times = np.unique(np.array(maintable["obsTime"]))
obervation_times = np.array(maintable["obsTime"])
number_of_occurances = []

for time in unique_obervation_times:
    n = 0
    for obs_time in obervation_times:
        if obs_time == time:
            n += 1
    number_of_occurances.append(n)

plt.figure()
plt.scatter(unique_obervation_times, number_of_occurances)
plt.ylim(0, 60000)
plt.xlabel("Observation time (days)")
plt.ylabel("Number of occurances")
ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.show()
