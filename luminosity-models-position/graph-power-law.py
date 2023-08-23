from matplotlib import pyplot as plt
import matplotlib.colors as colors
from math import log
import numpy as np
from matplotlib.lines import Line2D

plt.style.use('jcap')


PLOT_SIZE = 50

ALPHA = 1.94
L_MIN_RANGE = [1.0e28, 1.0e34]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
minPowerStep = (L_MIN_RANGE[1] / L_MIN_RANGE[0]) ** (1.0 / PLOT_SIZE)
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 1

paperPoint = [1e35, 1e29]

TOTAL_FGL_NUMBER = 265 # 109 unflagged
PERCENTAGES = [0.05, 0.1, 0.2, 0.4, 1]
FRACS_ABOVE_THRESHOLD=[0.05, 0.1, 0.2, 0.4]
TOTAL_FLUX = 1.794113925439598e-09#7.494712733226778e-10

DRAW_EXTRA_CONTOURS = False
LINE_COLOR ="C0"
PATH_TO_FILE = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-" \
+str(MULTIPLIER) + "x/power-law/"
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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.7)

# ========================== Load data ===========================

totalNum = []
numSeen = []
lumSeen = []

f = open(PATH_TO_FILE + "total-num-DR2.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    totalNum.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "num-seen-DR2.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    numSeen.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "lum-seen-DR2.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item) / TOTAL_FLUX)
    lumSeen.append(np.asarray(enterLine))

totalNum = np.stack(totalNum, axis=0)
numSeen = np.stack(numSeen, axis=0)
lumSeen = np.stack(lumSeen, axis=0)

print(numSeen)

# ========================= Display data =========================

dimMin = len(totalNum)
dimMax = len(totalNum[0])

lMinVals = [L_MIN_RANGE[0] * minPowerStep**i for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

fig, ax = plt.subplots()

plt.xscale("log")
plt.yscale("log")
plt.xlabel("$L_\\mathrm{max}$ [erg / s]")
plt.ylabel("$L_\\mathrm{min}$ [erg / s]")

c1 = plt.contourf(lMaxVals, lMinVals, totalNum,
                   norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N_\\mathrm{GCE}$")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, numSeen, [10*i for i in range(1, 20)],
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)])
n_contours = plt.contour(lMaxVals, lMinVals, numSeen, TOTAL_FGL_NUMBER * np.array(PERCENTAGES), colors=[LINE_COLOR], linewidths=[1])

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, lumSeen, [0.1*i for i in range(1, 15)],
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)])
r_contours = plt.contour(lMaxVals, lMinVals, lumSeen, FRACS_ABOVE_THRESHOLD, colors=[LINE_COLOR], linestyles='dashed', linewidths=[1])

fmt_n = {}
for i, item in enumerate(n_contours.levels):
    fmt_n[item] = str(int(100*PERCENTAGES[i])) + "\%"
fmt_r = {}
for i, item in enumerate(r_contours.levels):
    fmt_r[item] = str(int(100*PERCENTAGES[i])) + "\%"
n_labels = plt.clabel(n_contours, n_contours.levels, fmt=fmt_n, inline=True, fontsize=10, colors='k', manual=True)
r_labels = plt.clabel(r_contours, r_contours.levels, fmt=fmt_r, inline=True, fontsize=10, colors='k', manual=True)

# Observation
shade(numSeen, TOTAL_FGL_NUMBER * PERCENTAGES[0], lMaxVals, lMinVals)
shade(lumSeen, FRACS_ABOVE_THRESHOLD[0], lMaxVals, lMinVals, True)


# Final points

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor="aquamarine", marker='^', markersize=6)

custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor="aquamarine", marker='^', linestyle='None', markersize=6),]
plt.legend(custom_lines, ["$N_r$", "$R_r$", "Wavelet 1"], loc='lower right', prop={'size': 14})
plt.xlim(lMaxVals[0], lMaxVals[-1])
plt.ylim(lMinVals[0], lMinVals[-1])
plt.tight_layout()

plt.savefig("power-law-pos-x" + str(MULTIPLIER) + ".pdf")
plt.savefig("power-law-pos-x" + str(MULTIPLIER) + ".png")

plt.show()
