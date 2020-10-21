# FITS file code extracted from this tutorial:
# https://learn.astropy.org/FITS-images.html

# Flux data obtained from here:
# https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/

# NFW profile formula obtained from here:
# https://arxiv.org/pdf/1911.12369.pdf

import numpy as np
from astropy.utils.data import download_file
from astropy.io import fits
import matplotlib.pyplot as plt

MAX_SENS = np.inf

image_file = "detthresh_P8R3_source_10years_PL22.fits"

hdu_list = fits.open(image_file)
image_data = hdu_list[0].data
hdu_list.close()
pi = 3.1415926535

f = open("sensitivity_10.txt", 'w')
x = 0
for line in image_data:
    f.write(','.join([str(item) for item in line]))
    if x != image_data.shape[0] - 1:
        f.write("\n")
    x+= 1
f.close()


#################################################################
# Now for masked sensitivity
#################################################################
hdul = fits.open('../flux-histograms/data/gll_psc_v16.fit')# 3FGL
fgl_3 = hdul[1].data
hdul.close()

fgl3_poses = []
fgl3_mask_radius_squared = []
for row in fgl_3:
    if 20 < row["GLON"] < 340:
        continue
    if row["GLAT"] > 20 or row["GLAT"] < -20:
        continue
    if -2 < row["GLAT"] < 2:
        continue
    TS = row["Sqrt_TS100_300"]**2 + row["Sqrt_TS300_1000"]**2 + \
        row["Sqrt_TS1000_3000"]**2 + row["Sqrt_TS3000_10000"]**2 + row["Sqrt_TS10000_100000"]**2
    if 9 < TS < 49:
        fgl3_mask_radius_squared.append(0.3**2)
    elif TS >= 49:
        fgl3_mask_radius_squared.append(1.0**2)
    else:
        continue
    lon = row["GLON"]
    if lon > 180:
        lon -= 360
    fgl3_poses.append((lon, row["GLAT"]))


def is_cut(ilat, ilon):
    for i, (lon, lat) in enumerate(fgl3_poses):
        if (lon - ilon)**2 + (lat - ilat)**2 < fgl3_mask_radius_squared[i]:
            return True
    return False


image_file = "detthresh_P8R3_source_8years_PL22.fits"

hdu_list = fits.open(image_file)
image_data = hdu_list[0].data
hdu_list.close()
deltaLat = pi / image_data.shape[0]
deltaLon = 2 * pi / image_data.shape[1]

f = open("sensitivity_8_mask.txt", 'w')
for x, line in enumerate(image_data):
    lat = x * deltaLat - pi / 2
    write_line = []
    for y, item in enumerate(line):
        lon = y * deltaLon - pi
        sensitivity = item
        if abs(lat) < 20 * np.pi / 180 and abs(lon) < 20 * np.pi / 180:
            if is_cut(lat * 180 / np.pi, lon * 180 / np.pi):
                sensitivity = str(MAX_SENS)
        write_line.append(str(sensitivity))
    f.write(','.join(write_line))
    if x != image_data.shape[0] - 1:
        f.write("\n")
f.close()



######################################################################
# Display
######################################################################

f = open("sensitivity_8_mask.txt", 'r')
data = []
for x, line in enumerate(f.readlines()):
    if x == '': continue
    lat = x * deltaLat - pi / 2
    data_line = []
    for y, item in enumerate(line.split(',')):
        lon = y * deltaLon - pi
        if abs(lon) < 20 * pi / 180 and abs(lat) < 20 * pi / 180:
            data_line.append(float(item))
    if len(data_line) > 0:
        data.append(data_line)
f.close()

c = plt.imshow(data)
plt.colorbar(c)
plt.show()