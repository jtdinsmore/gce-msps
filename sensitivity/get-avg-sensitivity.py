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
from multiprocessing import Pool

plt.style.use('jcap')

DIST_TO_CENTER = 8.5
KAPPA = 0.5 # dist to center / r_s
GAMMA = 1.2 # Shape facter of the NFW profile
A = 1 # Something like the number of objects per cubic kiloparsec; in fact, it's the coefficient of the NFW number density
CM_PER_KPC = 3.086e21
FLUX_TO_LUM = 1.1095246594108431e-46
D_S = 0.001

DISPLAY_SIZE = 20.0
NUM_X_LABELS = 8
NUM_Y_LABELS = 8

R_S = 20
GAMMA = 1.2
R_CUT = None

def nfw(r):
    return ((r / R_S)**(-GAMMA) * (1 + (r/R_S))**(-3 + GAMMA))**2

image_file = "detthresh_P8R3_source_10years_PL22.fits"

hdu_list = fits.open(image_file)
#hdu_list.info()
image_data = np.flip(hdu_list[0].data, axis=1)
hdu_list.close()


trimmed_data = []
for x in range(len(image_data)):
    line = []
    for y in range(len(image_data[x])):
        flux = image_data[x][y]
        lat = x * pi / image_data.shape[0] - pi / 2
        lon = y * 2 * pi / image_data.shape[1] - pi
        if abs(lat) < 2 * pi / 180:
            flux = 0#1e-12
        if abs(lat) < DISPLAY_SIZE * pi / 180 and abs(lon) < DISPLAY_SIZE * pi / 180:
            line.append(flux)
    if len(line) > 0:
        trimmed_data.append(np.asarray(line))

trimmed_data = np.asarray(trimmed_data)


DTHETA = 40.0 / trimmed_data.shape[0] * pi / 180
DPHI = 40.0 / trimmed_data.shape[1] * pi / 180

print(trimmed_data.shape, DTHETA, DPHI)

fluxes = np.zeros_like(trimmed_data)

def count(index):
    i, j = index
    b = -20 * pi / 180 + j * DPHI
    l = -20 * pi / 180 + i * DTHETA
    if abs(b) <= 2 * pi / 180:
        return 0

    flux = 0
    for s in np.arange(D_S, 18 * DIST_TO_CENTER, D_S):
        r = sqrt(s**2 + DIST_TO_CENTER**2 - 2 * s * DIST_TO_CENTER * cos(l) * cos(b))
        if R_CUT is not None and r > R_CUT: continue
        flux += nfw(r) * cos(b)

    return flux

with Pool() as pool:
    indices = []
    for j in range(trimmed_data.shape[0]):
        for i in range(trimmed_data.shape[1]):
            indices.append((i, j))
    fluxes = pool.map(count, indices)

fluxes = np.array(fluxes).reshape(trimmed_data.shape[0], trimmed_data.shape[1])

plt.imshow(trimmed_data)
cbar = plt.colorbar()
cbar.set_label("$F_\\textrm{th}(b, \\ell)$ (erg/cm$^2$/s)")
plt.show()

plt.imshow(fluxes)
cbar = plt.colorbar()
cbar.set_label("$F_\\textrm{th}(b, \\ell)$ (erg/cm$^2$/s)")
plt.show()

average_sensitivity = np.sum(trimmed_data * fluxes) / np.sum(fluxes)
print("Average threshold (erg/cm^2 / s):", average_sensitivity)
print("Average threshold (erg/s):", average_sensitivity / FLUX_TO_LUM)
