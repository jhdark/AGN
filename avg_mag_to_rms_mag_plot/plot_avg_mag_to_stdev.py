from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import csv

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)


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
    min_data_points = 10
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
        e_1_avg = np.mean(e_1)
        e_2_avg = np.mean(e_2)
        e_3_avg = np.mean(e_3)
        e_4_avg = np.mean(e_4)
        e_5_avg = np.mean(e_5)
        e_6_avg = np.mean(e_6)
        e_7_avg = np.mean(e_7)
        e_8_avg = np.mean(e_8)
        e_9_avg = np.mean(e_9)
        e_10_avg = np.mean(e_10)
        e_11_avg = np.mean(e_11)
        avg_fluxes = [
            e_1_avg,
            e_2_avg,
            e_3_avg,
            e_4_avg,
            e_5_avg,
            e_6_avg,
            e_7_avg,
            e_8_avg,
            e_9_avg,
            e_10_avg,
            e_11_avg,
        ]
        avg_fluxes = [x for x in avg_fluxes if np.isnan(x) == False]
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
    avg_mags = data_file["avg_mags"]
    mag_rms = data_file["mag_rms"]
    return mag_rms, avg_mags


# mag_rms, avg_mags = plot_data()
mag_rms, avg_mags = plot_data_from_file()

# sort the arrays
# inds = avg_mags.argsort()
# mag_rms = mag_rms[inds]
# avg_mags = np.sort(avg_mags)

# rms
plt.figure()
plt.scatter(np.array(avg_mags), np.array(mag_rms), marker=".", color="black", s=1)

# compute moving average
# N = 5  # N is the window
# moving_avg = np.convolve(mag_rms, np.ones(N) / N, mode="valid")
# plt.plot(
#     np.array(avg_mags)[N - 1 :], moving_avg, label="Moving average N = {:.0f}".format(N)
# )

# plt.xlim(10, 25)
plt.ylim(0, 0.5)
# plt.legend()
plt.xlabel("Average magnitude (mag)")
plt.ylabel("r.m.s deviation (mag)")
ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

plt.show()
