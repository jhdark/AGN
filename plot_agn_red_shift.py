import matplotlib.pyplot as plt
import numpy as np

# lum_distance = cosmo.luminosity_distance(new_red_shift_values)
# luminosity_values = new_flux_values * (4 * np.pi * (lum_distance) ** 2)
# plt.figure()
# plt.scatter(new_red_shift_values, luminosity_values, s=5, marker=".", color="black")
# plt.ylim(0, 0.8e7)
# plt.xlim(left=0)
# plt.xlabel("Redshift")
# plt.ylabel("luminosity")
# plt.show()

red_shift_data_file = np.genfromtxt("new_redshift_data.csv", delimiter=",", names=True)
object_id = red_shift_data_file["objID"]
red_shifts = red_shift_data_file["red_shift"]
fluxes = red_shift_data_file["flux"]

agn_data = np.genfromtxt("AGN_candidates.csv", delimiter=",", names=True)
AGN_candidates = agn_data["AGN_cadidates"]

data = np.array(list(zip(object_id, red_shifts, fluxes)))
plot_data = []

for candidiate in AGN_candidates:
    for row in data:
        if row[0] == candidiate:
            plot_data.append(row)

print(len(AGN_candidates))
print(len(plot_data))

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

plot_data = np.array(plot_data)
print(plot_data)
plt.figure()
plt.hist(
    plot_data[:, 1], bins=10, color="red"
)  # Still contains bad redshift. Another round of filtering to get the subset with good redshifts.
plt.xlabel("Redshift")
plt.ylabel("No. of AGN")
plt.xlim(left=0)
ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.show()
