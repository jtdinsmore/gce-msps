from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')


PLOT_SIZE = 50

L_MIN = 1e29
ALPHA_RANGE = [1.1, 2.5]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 1

paperPoint = [1.94, 1e35]
bartels15Point = [1.5, 7e34]

TOTAL_FGL_NUMBER = 265 # 109 unflagged
PERCENTAGES = [0.05, 0.1, 0.2, 0.4, 1]
FRACS_ABOVE_THRESHOLD=[0.05, 0.1, 0.2, 0.4]
TOTAL_FLUX = 1.794113925439598e-09#7.494712733226778e-10

DRAW_EXTRA_CONTOURS = False
LINE_COLOR = "C2"
PATH_TO_FILE = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-"\
+str(MULTIPLIER)+"x/power-law-alpha/"
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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.5)

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

totalNum = np.transpose(np.stack(totalNum, axis=0))
numSeen = np.transpose(np.stack(numSeen, axis=0))
lumSeen = np.transpose(np.stack(lumSeen, axis=0))

# ========================= Display data =========================

dimMin = len(totalNum)
dimMax = len(totalNum[0])

alphaVals = [ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / PLOT_SIZE) for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

fig, ax = plt.subplots()

#plt.xscale("log")
plt.yscale("log")
plt.xlabel("$\\alpha$")
plt.ylabel("$L_\\mathrm{max}$ [erg / s]")

c1 = plt.contourf(alphaVals, lMaxVals, totalNum,
                   norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N_\\mathrm{GCE}$")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(alphaVals, lMaxVals, numSeen, [10*i for i in range(1, 20)],
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)])
n_contours = plt.contour(alphaVals, lMaxVals, numSeen, TOTAL_FGL_NUMBER * np.array(PERCENTAGES), colors=[LINE_COLOR], linewidths=[1])

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(alphaVals, lMaxVals, lumSeen, [0.1*i for i in range(1, 15)],
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)])
r_contours = plt.contour(alphaVals, lMaxVals, lumSeen, FRACS_ABOVE_THRESHOLD, colors=[LINE_COLOR], linestyles='dashed', linewidths=[1])

fmt_n = {}
for i, item in enumerate(n_contours.levels):
    fmt_n[item] = str(int(100*PERCENTAGES[i])) + "\%"
fmt_r = {}
for i, item in enumerate(r_contours.levels):
    fmt_r[item] = str(int(100*PERCENTAGES[i])) + "\%"
n_labels = plt.clabel(n_contours, n_contours.levels, fmt=fmt_n, inline=True, fontsize=10, colors='k', manual=True)
r_labels = plt.clabel(r_contours, r_contours.levels, fmt=fmt_r, inline=True, fontsize=10, colors='k', manual=True)


# Observation
shade(numSeen, TOTAL_FGL_NUMBER * PERCENTAGES[0], alphaVals, lMaxVals)
shade(lumSeen, FRACS_ABOVE_THRESHOLD[0], alphaVals, lMaxVals, True)


# Final points

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor="aquamarine", marker='^', markersize=6)
plt.plot(bartels15Point[0], bartels15Point[1], markeredgecolor='black', markerfacecolor="darkgreen", marker='v', markersize=6)

custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor="aquamarine", marker='^', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="darkgreen", marker='v', linestyle='None', markersize=6),]
plt.legend(custom_lines, ["$N_r$", "$R_r$", "Wavelet 1", "Wavelet 2"], loc='upper right', prop={'size': 14})
plt.xlim(alphaVals[0], alphaVals[-1])
plt.ylim(lMaxVals[0], lMaxVals[-1])
plt.tight_layout()

plt.savefig("power-law-alpha-pos-x" + str(MULTIPLIER) + ".pdf")
plt.savefig("power-law-alpha-pos-x" + str(MULTIPLIER) + ".png")

plt.show()
