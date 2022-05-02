from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import csv

maintable = Table.read(
    "COSMOS_PSdetection_filter2.fits", format="fits"
)  # Reads in the table directly from the .fits file

# maintable.write("COSMOS_PSdetection_filter2.csv", format="ascii.csv", overwrite=True)

id_list = np.unique(np.array(maintable["objID"]))
obj_IDS = np.array(maintable["objID"])
more_than_5_datapoints = []
more_than_10_datapoints = []

m = len(id_list)
for unique_id in id_list:
    m -= 1
    n = 0
    for id in obj_IDS:
        if id == unique_id:
            n += 1
    if n > 4:
        more_than_5_datapoints.append([unique_id])
        if n > 9:
            more_than_10_datapoints.append([unique_id])
    obj_IDS = np.delete(obj_IDS, np.where(obj_IDS == unique_id))
    print(m)
print("No id's with more than 5 data points: ", len(more_than_5_datapoints))
print("writing data 5 points")
with open("ids_5_data_points.txt", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["objID"])
    writer.writerows(more_than_5_datapoints)

print("No id's with more than 10 data points: ", len(more_than_10_datapoints))
print("writing data 10 points")
with open("ids_10_data_points.txt", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["objID"])
    writer.writerows(more_than_10_datapoints)
print("hello")
