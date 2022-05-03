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

maintable_structured = np.lib.recfunctions.structured_to_unstructured(
    maintable.as_array()
)

for unique_id in more_than_10_datapoints:
    m -= 1
    data = maintable_structured[np.where(maintable_structured[:, 0] == unique_id)]
    # avg mags
    magnitudes = []
    flux_values = data[:, 11]
    # this seems to only slow it down?
    # np.delete(maintable_structured, np.where(maintable_structured[:, 0] == unique_id))
    magnitudes = -2.5 * np.log10(flux_values / F_r) + m_r
    average_mag = np.mean(magnitudes)
    avg_mags.append(average_mag)
    # rms
    diffs = magnitudes - average_mag
    rms.append(np.mean(np.array(diffs)))
    print(m, "datapoints to go")

print(avg_mags)
print(rms)
plt.figure()
plt.scatter(np.array(avg_mags), np.array(rms), marker=".")

# data = zip(avg_mags, rms)
# with open("plotting_data.csv", "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(["avg_mags"], ["rms"])
#     writer.writerows(data)

plt.show()
