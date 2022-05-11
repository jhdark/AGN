from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import csv

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

maintable = Table.read("../COSMOS_PSdetection_filter2.fits", format="fits")

maintable_structured = np.lib.recfunctions.structured_to_unstructured(
    maintable.as_array()
)


def plot_data():
    maintable = Table.read(
        "../COSMOS_PSdetection_filter2.fits", format="fits"
    )  # Reads in the table directly from the .fits file

    # maintable.write("COSMOS_PSdetection_filter2.csv", format="ascii.csv", overwrite=True)

    # maintable numpy conversion
    maintable_structured = np.lib.recfunctions.structured_to_unstructured(
        maintable.as_array()
    )

    # retrive object ids with more than a define number datapoints
    obj_IDS = np.array(maintable["objID"])
    unique_ids, counts = np.unique(obj_IDS, return_counts=True)
    min_data_points = 1
    datapoint_data = unique_ids[np.where(counts >= min_data_points)]

    mag_rms = []
    avg_mags = []
    n = len(datapoint_data)
    for unique_id in datapoint_data:
        n -= 1
        data = maintable_structured[np.where(maintable_structured[:, 0] == unique_id)]
        # avg mags
        e_1, e_2, e_3, e_4, e_5, e_6 = [], [], [], [], [], []
        e_7, e_8, e_9, e_10, e_11 = [], [], [], [], []
        for row in data:
            if row[5] >= 55189 and row[5] <= 55196:
                e_1.append(row[11])
            elif row[5] >= 55245 and row[5] <= 55252:
                e_2.append(row[11])
            elif row[5] >= 55265 and row[5] <= 55272:
                e_3.append(row[11])
            elif row[5] >= 55613 and row[5] <= 55620:
                e_4.append(row[11])
            elif row[5] >= 55968 and row[5] <= 55975:
                e_5.append(row[11])
            elif row[5] >= 55313 and row[5] <= 56320:
                e_6.append(row[11])
            elif row[5] >= 56331 and row[5] <= 56338:
                e_7.append(row[11])
            elif row[5] >= 56342 and row[5] <= 56349:
                e_8.append(row[11])
            elif row[5] >= 56666 and row[5] <= 56673:
                e_9.append(row[11])
            elif row[5] >= 56678 and row[5] <= 56685:
                e_10.append(row[11])
            elif row[5] >= 57005 and row[5] <= 57012:
                e_11.append(row[11])
        avg_fluxes = [
            np.mean(e_1),
            np.mean(e_2),
            np.mean(e_3),
            np.mean(e_4),
            np.mean(e_5),
            np.mean(e_6),
            np.mean(e_7),
            np.mean(e_8),
            np.mean(e_9),
            np.mean(e_10),
            np.mean(e_11),
        ]
        avg_fluxes = [x for x in avg_fluxes if np.isnan(x) == False]
        # print(len(avg_fluxes))
        if len(avg_fluxes) <= 3:
            continue
        magnitudes = -2.5 * np.log10(np.array(avg_fluxes) / 3631)
        average_mag = np.mean(magnitudes)
        avg_mags.append(average_mag)
        # rms
        diffs = (magnitudes - average_mag) ** 2
        rms = (np.mean(np.array(diffs))) ** 0.5
        mag_rms.append(rms)
        print(n, "datapoints to go")

    # export data
    data = zip(datapoint_data, avg_mags, mag_rms)
    header = ["objID", "avg_mags", "mag_rms"]
    with open("plotting_data.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)

    return mag_rms, avg_mags


def plot_data_from_file():
    data_file = np.genfromtxt("plotting_data.csv", delimiter=",", names=True)
    object_id = data_file["objID"]
    avg_mags = data_file["avg_mags"]
    mag_rms = data_file["mag_rms"]
    return object_id, avg_mags, mag_rms


# mag_rms, avg_mags = plot_data()
object_id, avg_mags, mag_rms = plot_data_from_file()

dataset = np.array(list(zip(object_id, avg_mags, mag_rms)))
# remove values where the rms is greater than 5
dataset = dataset[np.where(dataset[:, 2] < 5)]
dataset = dataset[np.where(dataset[:, 1] > 12)]

# sort the arrays
# inds = avg_mags.argsort()
# mag_rms = mag_rms[inds]
# avg_mags = np.sort(avg_mags)

#################################################
# binning

b_1, b_2, b_3, b_4 = [], [], [], []
b_5, b_6, b_7, b_8 = [], [], [], []
b_9, b_10, b_11, b_12 = [], [], [], []
b_13, b_14, b_15, b_16 = [], [], [], []


for row in dataset:
    if row[1] >= 14 and row[1] < 14.5:
        b_1.append(row)
    elif row[1] >= 14.5 and row[1] < 15:
        b_2.append(row)
    elif row[1] >= 15 and row[1] < 15.5:
        b_3.append(row)
    elif row[1] >= 15.5 and row[1] < 16:
        b_4.append(row)
    elif row[1] >= 16 and row[1] < 16.5:
        b_5.append(row)
    elif row[1] >= 16.5 and row[1] < 17:
        b_6.append(row)
    elif row[1] >= 17 and row[1] < 17.5:
        b_7.append(row)
    elif row[1] >= 17.5 and row[1] < 18:
        b_8.append(row)
    elif row[1] >= 18 and row[1] < 18.5:
        b_9.append(row)
    elif row[1] >= 18.5 and row[1] < 19:
        b_10.append(row)
    elif row[1] >= 19 and row[1] < 19.5:
        b_11.append(row)
    elif row[1] >= 19.5 and row[1] < 20:
        b_12.append(row)
    elif row[1] >= 20 and row[1] < 20.5:
        b_13.append(row)
    elif row[1] >= 20.5 and row[1] < 21:
        b_14.append(row)
    elif row[1] >= 21 and row[1] < 21.5:
        b_15.append(row)
    elif row[1] >= 21.5 and row[1] < 22:
        b_16.append(row)

b_1 = np.array(b_1)
b_2 = np.array(b_2)
b_3 = np.array(b_3)
b_4 = np.array(b_4)
b_5 = np.array(b_5)
b_6 = np.array(b_6)
b_7 = np.array(b_7)
b_8 = np.array(b_8)
b_9 = np.array(b_9)
b_10 = np.array(b_10)
b_11 = np.array(b_11)
b_12 = np.array(b_12)
b_13 = np.array(b_13)
b_14 = np.array(b_14)
b_15 = np.array(b_15)
b_16 = np.array(b_16)

b_1_avg = np.mean(b_1[:, 2])
b_2_avg = np.mean(b_2[:, 2])
b_3_avg = np.mean(b_3[:, 2])
b_4_avg = np.mean(b_4[:, 2])
b_5_avg = np.mean(b_5[:, 2])
b_6_avg = np.mean(b_6[:, 2])
b_7_avg = np.mean(b_7[:, 2])
b_8_avg = np.mean(b_8[:, 2])
b_9_avg = np.mean(b_9[:, 2])
b_10_avg = np.mean(b_10[:, 2])
b_11_avg = np.mean(b_11[:, 2])
b_12_avg = np.mean(b_12[:, 2])
b_13_avg = np.mean(b_13[:, 2])
b_14_avg = np.mean(b_14[:, 2])
b_15_avg = np.mean(b_15[:, 2])
b_16_avg = np.mean(b_16[:, 2])

running_avgs = [
    b_1_avg,
    b_2_avg,
    b_3_avg,
    b_4_avg,
    b_5_avg,
    b_6_avg,
    b_7_avg,
    b_8_avg,
    b_9_avg,
    b_10_avg,
    b_11_avg,
    b_12_avg,
    b_13_avg,
    b_14_avg,
    b_15_avg,
    b_16_avg,
]
x_values = np.linspace(14.25, 21.75, num=16)


running_avgs = np.array(running_avgs)
avg_value = np.mean(running_avgs)
running_avg_rms_diffs = (running_avgs - avg_value) ** 2
running_avg_rms = (np.mean(np.array(running_avg_rms_diffs))) ** 0.5

agn_candidates = []

for row in b_1:
    rms_star = (row[2] - b_1_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_2:
    rms_star = (row[2] - b_2_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_3:
    rms_star = (row[2] - b_3_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_4:
    rms_star = (row[2] - b_4_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_5:
    rms_star = (row[2] - b_5_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_6:
    rms_star = (row[2] - b_6_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_7:
    rms_star = (row[2] - b_7_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_8:
    rms_star = (row[2] - b_8_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_9:
    rms_star = (row[2] - b_9_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_10:
    rms_star = (row[2] - b_10_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_11:
    rms_star = (row[2] - b_11_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_12:
    rms_star = (row[2] - b_12_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_13:
    rms_star = (row[2] - b_13_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_14:
    rms_star = (row[2] - b_14_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_15:
    rms_star = (row[2] - b_15_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

for row in b_16:
    rms_star = (row[2] - b_16_avg) / running_avg_rms
    if rms_star >= 3:
        agn_candidates.append(row[0])

print(agn_candidates)
print(len(agn_candidates))

avg_flux_data = []
n = len(agn_candidates)
for candidate in agn_candidates:
    n -= 1
    data = maintable_structured[np.where(maintable_structured[:, 0] == candidate)]
    fluxes = data[:, 11]
    avg_flux = np.mean(fluxes)
    avg_flux_data.append(avg_flux)
    print(n, "datapoints to go")

data = zip(agn_candidates, avg_flux_data)
header = ["objID", "avg_flux"]
with open("AGN_candidates.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(data)


# rms
plt.figure()
plt.scatter(dataset[:, 1], dataset[:, 2], marker=".", color="black", s=1)
plt.plot(x_values, running_avgs, linestyle="dashed", color="red")

# compute moving average
# N = 800  # N is the window
# moving_avg = np.convolve(mag_rms, np.ones(N) / N, mode="valid")
# plt.plot(
#     np.array(avg_mags)[N - 1 :], moving_avg, label="Moving average N = {:.0f}".format(N)
# )

plt.xlim(14, 22)
plt.ylim(0, 1)
plt.xlabel("Average magnitude (mag)")
plt.ylabel("r.m.s deviation (mag)")
ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

plt.show()
