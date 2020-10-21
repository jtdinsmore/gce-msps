# FITS file code extracted from this tutorial:
# https://learn.astropy.org/FITS-images.html

# Flux data obtained from here:
# https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/

# NFW profile formula obtained from here:
# https://arxiv.org/pdf/1911.12369.pdf

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.utils.data import download_file
from astropy.io import fits
from math import pi, cos, sqrt

plt.style.use('jcap')

DIST_TO_CENTER = 8.5
KAPPA = 0.5 # dist to center / r_s
GAMMA = 1.2 # Shape facter of the NFW profile
A = 1 # Something like the number of objects per cubic kiloparsec; in fact, it's the coefficient of the NFW number density
CM_PER_KPC = 3.086e21
FLUX_TO_LUM = 1.1095246594108431e-46

DISPLAY_SIZE = 20.0
THRESHOLD = 1e34
NUM_X_LABELS = 8
NUM_Y_LABELS = 8

hdu_list_10 = fits.open("detthresh_P8R3_source_10years_PL22.fits")
image_data_10 = hdu_list_10[0].data
image_data_10 = np.flip(image_data_10, axis=0)
print(hdu_list_10[0].header)
hdu_list_10.close()

print('\n'*10)

hdu_list_8 = fits.open("detthresh_P8R3_source_8years_PL22.fits")
print(hdu_list_8[0].header)
image_data_8 = hdu_list_8[0].data
image_data_8 = np.flip(image_data_8, axis=0)
hdu_list_8.close()

assert image_data_8.shape == image_data_10.shape
deltaLat = pi / image_data_10.shape[0]
deltaLon = 2 * pi / image_data_10.shape[1]

trimmed_flux_data_10 = []
trimmed_flux_data_8 = []
for x in range(len(image_data_10)):
    fluxLine_10 = []
    fluxLine_8 = []
    for y in range(len(image_data_10[x])):
        flux_10 = image_data_10[x][y]
        flux_8 = image_data_8[x][y]
        lat = x * deltaLat - pi / 2
        lon = y * deltaLon - pi
        if abs(lat) < 2 * pi / 180:
            flux_10 = np.nan
            flux_8 = np.nan
        if abs(lat) > DISPLAY_SIZE * pi / 180 or abs(lon) > DISPLAY_SIZE * pi / 180:
            pass
        else:
            fluxLine_10.append(flux_10)
            fluxLine_8.append(flux_8)
    if len(fluxLine_10) > 0:
        trimmed_flux_data_10.append(np.asarray(fluxLine_10))
        trimmed_flux_data_8.append(np.asarray(fluxLine_8))

trimmed_flux_data_10 = np.asarray(trimmed_flux_data_10)
trimmed_flux_data_8 = np.asarray(trimmed_flux_data_8)

plt.figure(figsize=(6, 5))
#mpl.rcParams["font.size"]=12

# Display flux data:
value = trimmed_flux_data_10
cf = plt.imshow(value,
                   vmin=np.nanmin(value),
                   vmax=np.nanpercentile(value, 99))
cbar1 = plt.colorbar(cf)
cbar1.set_label("$F_\\textrm{th}^\\textrm{10 year}(b, \\ell)$ (erg/cm$^2$/s)")

xPositions = np.arange(0, value.shape[1], value.shape[1]//NUM_X_LABELS) # pixel count at label position
xLabels = np.linspace(start=DISPLAY_SIZE, stop=-DISPLAY_SIZE, num=NUM_X_LABELS+1) # labels you want to see
np.append(xPositions, value.shape[1])
plt.xticks(xPositions, [int(i) for i in xLabels])

yPositions = np.arange(0, value.shape[0], value.shape[0]//NUM_Y_LABELS) # pixel count at label position
yLabels = np.linspace(start=DISPLAY_SIZE, stop=-DISPLAY_SIZE, num=NUM_Y_LABELS+1) # labels you want to see
np.append(yPositions, value.shape[0])
plt.yticks(yPositions, [int(i) for i in yLabels])
plt.xlabel("$\\ell$ (deg)")
plt.ylabel("$b$ (deg)")
plt.tight_layout()
plt.savefig("flux-thresholds.pdf")
plt.show()


plt.figure(figsize=(6, 5))
quotient = trimmed_flux_data_10 / trimmed_flux_data_8
print(quotient)
c = plt.imshow(quotient, vmin=np.nanmin(quotient), vmax=np.nanpercentile(quotient, 99))
cbar = plt.colorbar(c)
cbar.set_label("$F_\\textrm{th}^\\textrm{10 year}(b, \\ell) / F_\\textrm{th}^\\textrm{8 year}(b, \\ell)$")

xPositions = np.arange(0, quotient.shape[1], quotient.shape[1]//NUM_X_LABELS) # pixel count at label position
xLabels = np.linspace(start=DISPLAY_SIZE, stop=-DISPLAY_SIZE, num=NUM_X_LABELS+1) # labels you want to see
np.append(xPositions, quotient.shape[1])
plt.xticks(xPositions, [int(i) for i in xLabels])

yPositions = np.arange(0, quotient.shape[0], quotient.shape[0]//NUM_Y_LABELS) # pixel count at label position
yLabels = np.linspace(start=DISPLAY_SIZE, stop=-DISPLAY_SIZE, num=NUM_Y_LABELS+1) # labels you want to see
np.append(yPositions, quotient.shape[0])
plt.yticks(yPositions, [int(i) for i in yLabels])
plt.xlabel("$\\ell$ (deg)")
plt.ylabel("$b$ (deg)")
plt.tight_layout()
plt.savefig("flux-increase.pdf")
plt.show()