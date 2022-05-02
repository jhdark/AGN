from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


tab = Table.read("COSMOS_zcat_forProject.fits", format="fits")  #  redshift table
tab.colnames

tab["Z_BEST"]  # contains your best estimate of the redshift
(index,) = np.where(
    tab["Z_BEST"] > 0
)  # but not all redshifts are good. Some are less than 0....
# plt.figure()
# plt.hist(
#     tab["Z_BEST"][index], bins=30
# )  # and some are around 10. Filter these cases out after cross-matching

tab1 = Table.read(
    "COSMOS_PSdetection_filter2.fits"
)  # example table for crossmatching. You may have a different version, where you only use the positions of the unique objects from your full table.

tab1.colnames
# plt.figure()
# plt.plot(
#     tab1["ra"], tab1["dec"], "k.", markersize=1
# )  # plot the coordinate of the objects on the sky
# plt.plot(tab1["ra"], tab1["dec"], "k.", markersize=1)
# plt.plot(tab["RA"], tab["DEC"], "r.", markersize=1)

coords1 = SkyCoord(
    tab["RA"], tab["DEC"]
)  # redshift table already has units for RA, DEC
coords2 = SkyCoord(
    tab1["ra"], tab1["dec"], unit=u.deg
)  # panstarrs table does not, so units have to be specified


match, dist2d, dist3d = coords1.match_to_catalog_sky(coords2)
len(match)
len(dist2d)
coords1[0]
match[0]
coords2[624568]
dist2d[0].arcsec

# plt.figure()
# plt.hist(
#     dist2d.arcsec, bins=100, range=(0, 3)
# )  # decided that optimal maximum matching distance is 0.5" for this table
(index,) = np.where(dist2d.arcsec < 0.5)
coords1[index]
coords2[match[index]]

# plt.figure()
# plt.hist(
#     tab["Z_BEST"][index], bins=100
# )  # Still contains bad redshift. Another round of filtering to get the subset with good redshifts.
# plt.show()
