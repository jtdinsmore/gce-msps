import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, sin, cos, pi
from astropy.utils.data import download_file
from astropy.io import fits

# See the point source map in fig 8 of Ajello 2017.

PIXEL_WIDTH = 0.46 # degrees
ROI_EDGE = 10 # degrees


NUM_ANGLE_STEPS = 1000
NUM_S_STEPS = 10000
R_C=8.5
R_S=20
INTEGRAL_LIMIT=20 * R_C
R_CUT = 2

side_length = int((2 * ROI_EDGE) // PIXEL_WIDTH)

coords_horiz = np.linspace(-10, 10, side_length)
coords_vert = np.linspace(-10, 10, side_length)

hdu_list = fits.open("gll_psc_v16.fit")
# COLUMN NAME INFO HERE: https://arxiv.org/pdf/1501.02003.pdf

threefgl = hdu_list[1].data
hdu_list.close()

ps_poses = []
fluxes = [f for f in threefgl["Flux1000"]]
fluxes.sort()

mask = np.ones((side_length, side_length))

for row in threefgl:
    if row["Flux1000"] >= fluxes[-200]:
        l = row["GLON"]
        b = row["GLAT"]
        if l > 180:
            l-= 360
        ps_poses.append((l, b))



        # Exclude for disk
for i in range(side_length):
    b = coords_vert[i]
    for j in range(side_length):
        l = coords_horiz[j]
        if sqrt(l**2 + b**2) > 10:
            mask[i][j] = 0;
            continue
empty_sum = np.sum(mask)


# Exclude for point sources
for i in range(side_length):
    b = coords_vert[i]
    for j in range(side_length):
        l = coords_horiz[j]
        for ps_l, ps_b in ps_poses:
            if sqrt((b - ps_b)**2 + (l - ps_l)**2) < 1 + PIXEL_WIDTH * sin(pi / 4):
                mask[i][j] = 0;
                break

plt.pcolor(coords_horiz, coords_vert, mask)
plt.xlabel("$l$ (deg)")
plt.ylabel("$b$ (deg)")
plt.axis('square')

print(np.sum(mask), "pixels in the ROI")
#plt.show()




def NFWSquared(r, gamma):
    return ((r/R_S)**(-gamma) * (1 + r/R_S)**(-3+gamma))**2

def my_roi_deny(l, b):
    return abs(b) < 2 * pi / 180

def ajello_roi_deny(l, b):
    if sqrt(l*l+b*b) > 10 * pi / 180:
        return True
    i = 0
    sign = coords_vert[0] * pi / 180 < b
    while (coords_vert[i+1] * pi / 180 < b) == sign:
        i += 1
    j = 0
    sign = coords_horiz[0] * pi / 180 < l
    while (coords_horiz[j+1] * pi / 180 < l) == sign:
        j += 1
    return mask[i][j] == 0

def get_flux(deny):
    integral = 0
    for b in np.linspace(-20 * pi / 180, 20 * pi / 180, NUM_ANGLE_STEPS):
        for l in np.linspace(-20 * pi / 180, 20 * pi / 180, NUM_ANGLE_STEPS):
            if deny(l, b):
                continue
            for s in np.linspace(0, INTEGRAL_LIMIT, NUM_S_STEPS):
                r = sqrt(R_C*R_C + s * s - 2 * s * R_C * cos(b) * cos(l))
                if R_CUT is not None and r > R_CUT: continue
                integral += NFWSquared(r, 1.2)
    return integral

my_roi_flux = get_flux(my_roi_deny)
ajello_flux = get_flux(ajello_roi_deny)

print("Ajello flux ratio:", ajello_flux / my_roi_flux) # (1000, 1000) gets 0.5555527292865902
plt.show()
