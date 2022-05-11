from cProfile import label
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import csv

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

tab = Table.read("../COSMOS_zcat_forProject.fits", format="fits")  #  redshift table
tab.colnames

tab["Z_BEST"]  # contains your best estimate of the redshift
(index,) = np.where(
    tab["Z_BEST"] > 0
)  # but not all redshifts are good. Some are less than 0....


def find_unique_obj_ra_dec_values():

    tab1 = Table.read(
        "../COSMOS_PSdetection_filter2.fits"
    )  # example table for crossmatching. You may have a different version, where you only use the positions of the unique objects from your full table.
    tab1.colnames
    tab1_structured = np.lib.recfunctions.structured_to_unstructured(tab1.as_array())
    unique_ids = np.unique(np.array(tab1["objID"]))

    obj_IDS = np.array(tab1["objID"])
    unique_ids, counts = np.unique(obj_IDS, return_counts=True)
    min_data_points = 1
    datapoint_data = unique_ids[np.where(counts >= min_data_points)]

    ra_values = []
    dec_values = []
    flux_values = []
    n = len(datapoint_data)
    for unique_id in datapoint_data:
        n -= 1
        data = tab1_structured[np.where(tab1_structured[:, 0] == unique_id)]
        # ra and dec values
        all_ra_values = data[:, 6]
        avg_ra_values = np.mean(all_ra_values)
        ra_values.append(avg_ra_values)
        all_dec_values = data[:, 7]
        avg_dec_values = np.mean(all_dec_values)
        dec_values.append(avg_dec_values)
        # flux values
        flux_data = data[:, 11]
        avg_flux = np.mean(flux_data)
        flux_values.append(avg_flux)
        print(n, "datapoints to go")

    data = zip(datapoint_data, ra_values, dec_values, flux_values)
    header = ["objID", "ra", "dec", "flux"]
    with open("plotting_data.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)

    return ra_values, dec_values, flux_values


def plot_data_from_file():
    tab_alt = np.genfromtxt("plotting_data.csv", delimiter=",", names=True)
    ra_values = tab_alt["ra"]
    dec_values = tab_alt["dec"]
    flux_values = tab_alt["flux"]

    return ra_values, dec_values, flux_values


# ra_values, dec_values = find_unique_obj_ra_dec_values()
ra_values, dec_values, flux_values = plot_data_from_file()


# plot the coordinate of the objects on the sky
plt.figure()
plt.title("Panstarrs dataset")
# plt.plot(tab1["ra"], tab1["dec"], "k.", markersize=0.1)
plt.plot(ra_values, dec_values, "k.", markersize=0.1)
plt.xlabel(r"Ra ($\mathrm{^{\circ}}$)")
plt.ylabel(r"Dec ($\mathrm{^{\circ}}$)")
ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)


plt.figure()
plt.title("Panstarrs dataset compared to redshift catalogue")
plt.plot(
    tab["RA"], tab["DEC"], "r.", label="Redshift catalogue", markersize=0.1, alpha=0.5
)
# plt.plot(tab1["ra"], tab1["dec"], "k.", markersize=0.1)
plt.plot(ra_values, dec_values, "k.", label="Panstarrs dataset", markersize=0.1)
plt.xlabel(r"Ra ($\mathrm{^{\circ}}$)")
plt.ylabel(r"Dec ($\mathrm{^{\circ}}$)")
# plt.xlim(149, 153)
# plt.ylim(1, 4)
ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# Plot legend.
lgnd = plt.legend(loc="upper right", numpoints=1, fontsize=10)
lgnd.legendHandles[0]._legmarker.set_markersize(6)
lgnd.legendHandles[1]._legmarker.set_markersize(6)
plt.tight_layout()
plt.show()
# redshift table already has units for RA, DEC

# coords1 = SkyCoord(
#     tab["RA"], tab["DEC"]
# )  # redshift table already has units for RA, DEC
# coords2 = SkyCoord(ra_values, dec_values, unit=u.deg)

# match, dist2d, dist3d = coords1.match_to_catalog_sky(coords2)
# print(np.array(match))
# print(dist2d)
# coords1[0]
# match[0]
# coords2[624568]
# dist2d[0].arcsec

# plt.figure()
# plt.hist(
#     dist2d.arcsec, bins=100, range=(0, 3)
# )  # decided that optimal maximum matching distance is 0.5" for this table
# (index,) = np.where(dist2d.arcsec < 0.5)
# coords1[index]
# coords2[match[index]]

# header = ["Z_BEST", "dist2d"]
# rows = zip(tab["Z_BEST"], dist2d.degrees)
# with open("redshift_data.csv", "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerow(header)
#     for row in rows:
#         writer.writerow(row)


# (index,) = np.where(dist2d.arcsec < 0.5)
# (index,) = np.where(tab["Z_BEST"] < 8)
# (index,) = np.where(tab["Z_BEST"] > 0)

# plt.figure()
# plt.hist(tab["Z_BEST"][index], bins=100)  # Still contains bad redshift. Another round of filtering to get the subset with good redshifts.
# plt.show()
