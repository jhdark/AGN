from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import csv

# plt.rc("text", usetex=True)
# plt.rc("font", family="serif", size=12)

epoch_time_1 = 55189 + 3.5
epoch_time_2 = 55245 + 3.5
epoch_time_3 = 55265 + 3.5
epoch_time_4 = 55613 + 3.5
epoch_time_5 = 55968 + 3.5
epoch_time_6 = 56313 + 3.5
epoch_time_7 = 56331 + 3.5
epoch_time_8 = 56342 + 3.5
epoch_time_9 = 56666 + 3.5
epoch_time_10 = 56678 + 3.5
epoch_time_11 = 57005 + 3.5

mags_1 = np.array(
    [
        21.31870254,
        21.32152425,
        21.36189678,
        21.1857959,
        21.02380393,
        21.27319379,
        21.37390231,
    ]
)

obs_times_1 = np.array(
    [
        epoch_time_1,
        epoch_time_3,
        epoch_time_4,
        epoch_time_5,
        epoch_time_6,
        epoch_time_9,
        epoch_time_10,
    ]
)

plt.figure()
plt.scatter(obs_times_1, mags_1)


# all on one graph
# plt.figure()
# for id in ids:
#     data = maintable_structured[np.where(maintable_structured[:, 0] == id)]
#     mags = -2.5 * np.log10(np.array(data[:, 11]) / 3631)
#     obs_times = data[:, 5]
#     plt.scatter(obs_times, mags)

# plt.xlabel("Observation time (days)")

plt.show()
