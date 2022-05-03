from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import csv

maintable = Table.read(
    "COSMOS_PSdetection_filter2.fits", format="fits"
)  # Reads in the table directly from the .fits file

# maintable.write("COSMOS_PSdetection_filter2.csv", format="ascii.csv", overwrite=True)

# retrive object ids with more than 5/10 datapoints
obj_IDS = np.array(maintable["objID"])
unique_ids, counts = np.unique(obj_IDS, return_counts=True)
more_than_5_datapoints = unique_ids[np.where(counts >= 5)]
more_than_10_datapoints = unique_ids[np.where(counts >= 10)]

# plot
rms = []
avg_mags = []
m_r = 1
F_r = 1
m = len(more_than_10_datapoints)

maintable = np.array(maintable)

for unique_id in more_than_10_datapoints:
    m -= 1
    row_data = []
    for row in maintable:
        if row[0] == unique_id:
            row_data.append(row[11])
    np.delete(maintable, np.where(row[0] == unique_id))
    # avg mags
    magnitudes = []
    for i in row_data:
        flux_value = i
        magnitude = -2.5 * np.log10(flux_value / F_r) + m_r
        magnitudes.append(magnitude)
    average_mag = np.mean(magnitudes)
    avg_mags.append(average_mag)

    # rms
    diffs = []
    for mag in magnitudes:
        diff = mag - average_mag
        diffs.append(diff)
    rms.append(np.mean(np.array(diffs)))
    print(m)

print(avg_mags)
print(rms)
plt.figure()
plt.scatter(np.array(avg_mags), np.array(rms))

data = zip(avg_mags, rms)
with open("plotting_data.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["avg_mags"], ["rms"])
    writer.writerows(data)

plt.show()
