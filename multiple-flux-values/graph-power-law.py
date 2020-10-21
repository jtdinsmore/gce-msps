from matplotlib import pyplot as plt
import matplotlib.colors as colors
from math import log
import numpy as np
from matplotlib.lines import Line2D

plt.style.use('jcap')


PLOT_SIZE = 50

DI_MAURO_FLUX = 1.794113925439598e-09

ALPHA = 1.94
L_MIN_RANGE = [1.0e28, 1.0e34]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
minPowerStep = (L_MIN_RANGE[1] / L_MIN_RANGE[0]) ** (1.0 / PLOT_SIZE)
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 1

paperPoint = [1e35, 1e29]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.17

str_flux_quot = ["3", "2", "1/2", "1/3", "1/4", "1/5"]
FLUX_QUOT = np.asarray([1/3, 1/2, 2, 3, 4, 5])
N_COLORS = np.flip([[l, l, l] for l in np.linspace(0, 1, len(FLUX_QUOT))])
R_COLORS = np.flip([[l, l, l] for l in np.linspace(0, 1, len(FLUX_QUOT))])

DRAW_EXTRA_CONTOURS = True
SHOW_NUMBERS = False
LINE_COLOR = "C0"
PATH_TO_FILE = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-1x/power-law/"
SHADE_SCALE=25

def shade(field, threshold, xs, ys, off=False):
    px = []
    py = []
    for x in range(0 if off else 1, SHADE_SCALE, 1):
        inx = int(float(x) / SHADE_SCALE * field.shape[1])
        for y in range(0 if off else 1, SHADE_SCALE, 1):
            iny = int(float(y) / SHADE_SCALE * field.shape[0])
            if field[iny][inx] < threshold:
                fracx = float(x) / SHADE_SCALE * field.shape[1] - inx
                fracy = float(y) / SHADE_SCALE * field.shape[0] - iny
                px.append(xs[inx] + fracx * (xs[inx+1] - xs[inx]))
                py.append(ys[iny] + fracy * (ys[iny+1] - ys[iny]))
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=1)

# ========================== Load data ===========================

totalNum = []
numSeen = []
lumSeen = []

f = open(PATH_TO_FILE + "total-num.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    totalNum.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "num-seen.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    numSeen.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "lum-seen.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item) / DI_MAURO_FLUX)
    lumSeen.append(np.asarray(enterLine))

totalNum = np.stack(totalNum, axis=0)
numSeen = np.stack(numSeen, axis=0)
lumSeen = np.stack(lumSeen, axis=0)

# ========================= Display data =========================

dimMin = len(totalNum)
dimMax = len(totalNum[0])

lMinVals = [L_MIN_RANGE[0] * minPowerStep**i for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

fig, ax = plt.subplots(figsize=(6, 4))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("$L_\\mathrm{max}$ [erg / s]")
plt.ylabel("$L_\\mathrm{min}$ [erg / s]")

if SHOW_NUMBERS:
    c1 = plt.pcolor(lMaxVals, lMinVals, totalNum,
                    norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                    vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
    cbar = plt.colorbar(c1, extend='max')
    cbar.set_label("$N_\\text{GCE}$")



n_contours = plt.contour(lMaxVals, lMinVals, numSeen, NUM_PULSARS_ABOVE_THRESHOLD * FLUX_QUOT, cmap="Greys",
    linestyles="solid", linewidths=1, vmax=np.max(NUM_PULSARS_ABOVE_THRESHOLD * FLUX_QUOT) / 10)
r_contours = plt.contour(lMaxVals, lMinVals, lumSeen, FRAC_ABOVE_THRESHOLD * FLUX_QUOT, cmap="Greys",
    linestyles="dashed", linewidths=1, vmax=np.max(FRAC_ABOVE_THRESHOLD * FLUX_QUOT) / 10)

plt.contour(lMaxVals, lMinVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=LINE_COLOR, linestyles="solid", linewidths=2)
plt.contour(lMaxVals, lMinVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2)

# Final points

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor='aquamarine', marker='^', markersize=6)


# Final points

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor='aquamarine', marker='^', markersize=6)

custom_lines = [Line2D([0], [0], color=LINE_COLOR, lw=2),
                Line2D([0], [0], color=LINE_COLOR, linestyle='dashed', lw=2, dashes=(4, 2))]\
             + [Line2D([], [], markeredgecolor='black', markerfacecolor='aquamarine', marker='^', linestyle='none', markersize=6),]
plt.legend(custom_lines, ['$N_r$', '$R_r$'] + ["Wavelet 1"], loc="lower right", prop={'size': 14})
plt.xlim(lMaxVals[0], lMaxVals[-1])
plt.ylim(lMinVals[0], lMinVals[-1])

fmt_n = {}
fmt_r = {}
for i, item in enumerate(n_contours.levels):
    fmt_n[item] = "$\\times$"+str_flux_quot[i]
for i, item in enumerate(r_contours.levels):
    fmt_r[item] = "$\\times$"+str_flux_quot[i]

n_labels = plt.clabel(n_contours, n_contours.levels, fmt=fmt_n, inline=True, fontsize=10, colors='k', manual=True)
for l in n_labels:
    pass
    #l.set_rotation(0)

r_labels = plt.clabel(r_contours, r_contours.levels, fmt=fmt_r, inline=True, fontsize=10, colors='k', manual=True)
for l in r_labels:
    pass
    #l.set_rotation(0)
plt.tight_layout()


plt.savefig("overlay-power-law.pdf")

plt.show()
