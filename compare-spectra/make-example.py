import numpy as np
import parse_spectrum as ps
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtrans

plt.style.use("jcap")

NAME = ps.DI_MAURO
SPECTRUM_RANGE = [0.1, 100]# GeV
NUM_PLOTS = 7
SHAPES=['o', 's', '^', '*', 'd', '+', 'x']

ERGS_PER_GEV = 0.00160218

fig, ax = plt.subplots()

f = ps.Spectrum(NAME)
f.display_data(ax)
f.display_calore(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
f.display_power_law(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
f.label_axes(ax)

print("Lbreak:", f.fit_l_break)
print("alpha above:", f.fit_alpha_above)
print("alpha below:", f.fit_alpha_below)
print("norm:", f.fit_norm)

ax.plot([], [], color='red', linestyle='--', label="Ref. [5] fit")
ax.plot([], [], color='blue', linestyle='-.', label="Power law fit")
ax.set_ylim(2e-8, 5e-7)

plt.title("")
fig.tight_layout()
ax.legend()
fig.savefig("example.pdf")
plt.show()
