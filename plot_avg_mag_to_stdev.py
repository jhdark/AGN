from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import csv

maintable = Table.read(
    "COSMOS_PSdetection_filter2.fits", format="fits"
)  # Reads in the table directly from the .fits file

# maintable.write("COSMOS_PSdetection_filter2.csv", format="ascii.csv", overwrite=True)

data_5_points = np.genfromtxt("ids_5_data_points.txt", delimiter=",", names=True)
data_10_points = np.genfromtxt("ids_10_data_points.csv", names=True)

data_5 = data_5_points["objID"]
data_10 = data_5_points["objID"]

# retrive object ids with more than 5/10 datapoints
id_list = np.unique(np.array(maintable["objID"]))
obj_IDS = np.array(maintable["objID"])
more_than_5_datapoints = []
more_than_10_datapoints = []

# m = len(id_list)
# for unique_id in id_list:
#     m -= 1
#     n = 0
#     for id in obj_IDS:
#         if id == unique_id:
#             n += 1
#     if n > 4:
#         more_than_5_datapoints.append([unique_id])
#         if n > 9:
#             more_than_10_datapoints.append([unique_id])
#     obj_IDS = np.delete(obj_IDS, np.where(obj_IDS == unique_id))
#     print(m)
# print("Number of id's with more than 5 data points: ", len(more_than_5_datapoints))
# print("writing data 5 points")
# with open("ids_5_data_points.csv", "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(["objID"])
#     writer.writerows(more_than_5_datapoints)
# print("Number of id's with more than 10 data points: ", len(more_than_10_datapoints))
# print("writing data 10 points")
# with open("ids_10_data_points.csv", "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(["objID"])
#     writer.writerows(more_than_10_datapoints)

rms = []
avg_mags = []
m_r = 1
F_r = 1
m = len(data_5)
for unique_id in data_5:
    m -= 1
    row_data = []
    for row in np.array(maintable):
        if row[0] == unique_id:
            row_data.append(row[11])
            maintable.remove_row(0)
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
