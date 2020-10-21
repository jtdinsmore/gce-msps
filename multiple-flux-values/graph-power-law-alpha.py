from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')

DI_MAURO_FLUX = 1.794113925439598e-09

PLOT_SIZE = 50

L_MIN = 1e29
ALPHA_RANGE = [1.1, 2.5]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 1

paperPoint = [1.94, 1e35]
bartels15Point = [1.5, 7e34]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.17

str_flux_quot = ["3", "2", "1/2", "1/3", "1/4", "1/5"]
FLUX_QUOT = np.asarray([1/3, 1/2, 2, 3, 4, 5])
N_COLORS = np.flip([[l, l, l] for l in np.linspace(0, 1, len(FLUX_QUOT))])
R_COLORS = np.flip([[l, l, l] for l in np.linspace(0, 1, len(FLUX_QUOT))])

PATH_TO_FILE = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-1x/power-law-alpha/"
SHADE_SCALE=25
LINE_COLOR = "C2"
SHOW_NUMBERS = False
DRAW_EXTRA_CONTOURS = True

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

totalNum = np.transpose(np.stack(totalNum, axis=0))
numSeen = np.transpose(np.stack(numSeen, axis=0))
lumSeen = np.transpose(np.stack(lumSeen, axis=0))

# ========================= Display data =========================

alphaVals = [ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / PLOT_SIZE) for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

fig, ax = plt.subplots(figsize=(6, 4))

#plt.xscale("log")
plt.yscale("log")
plt.xlabel("$\\alpha$")
plt.ylabel("$L_\\mathrm{max}$ [erg / s]")

if SHOW_NUMBERS:
    c1 = plt.pcolor(alphaVals, lMaxVals, totalNum,
                    norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                    vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
    cbar = plt.colorbar(c1, extend='max')
    cbar.set_label("$N_\\text{GCE}$")


n_contours = plt.contour(alphaVals, lMaxVals, numSeen, NUM_PULSARS_ABOVE_THRESHOLD * FLUX_QUOT, cmap="Greys",
    linestyles="solid", linewidths=1, vmax=np.max(NUM_PULSARS_ABOVE_THRESHOLD * FLUX_QUOT) / 10)
r_contours = plt.contour(alphaVals, lMaxVals, lumSeen, FRAC_ABOVE_THRESHOLD * FLUX_QUOT, cmap="Greys",
    linestyles="dashed", linewidths=1, vmax=np.max(FRAC_ABOVE_THRESHOLD * FLUX_QUOT) / 10)


plt.contour(alphaVals, lMaxVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=LINE_COLOR, linestyles="solid", linewidths=2)
plt.contour(alphaVals, lMaxVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2)
    #l.set_rotation(0)

# Final points

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor='aquamarine', marker='^', markersize=6)
plt.plot(bartels15Point[0], bartels15Point[1], markeredgecolor='black', markerfacecolor='darkgreen', marker='v', markersize=6)


custom_lines = [Line2D([0], [0], color=LINE_COLOR, linestyle='solid', lw=2),
                Line2D([0], [0], color=LINE_COLOR, linestyle='dashed', lw=2, dashes=(4, 2))] \
             + [Line2D([], [], markeredgecolor='black', markerfacecolor="aquamarine", marker='^', linestyle='none', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="darkgreen", marker='v', linestyle='none', markersize=6)]
plt.legend(custom_lines, ['$N_r$', '$R_r$'] + ["Wavelet 1", "Wavelet 2"], loc="lower right", prop={'size': 14})
plt.ylim(lMaxVals[0], lMaxVals[-1])
plt.xlim(alphaVals[0], alphaVals[-1])

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
plt.tight_layout()


plt.savefig("overlay-power-law-alpha.pdf")

plt.show()
