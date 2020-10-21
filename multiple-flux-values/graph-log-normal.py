from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors
import numpy as np

from matplotlib.lines import Line2D

plt.style.use('jcap')


PLOT_SIZE = 50
L_0_RANGE=[1.0e30, 2.0e36]#[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]

MULTIPLIER = 1

lOPowerStep = (L_0_RANGE[1] / L_0_RANGE[0]) ** (1.0 / PLOT_SIZE)

DI_MAURO_FLUX = 1.794113925439598e-09

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.17

DRAW_EXTRA_CONTOURS = True
DRAW_PLOEG_POINT = True

paperPoint = [0.88e34, 0.62]
ploegPoint = [1.3023e+32, 0.69550156]
gautamPoint = [4.2970e+30, 0.93936155]
SHOW_NUMBERS = False

str_flux_quot = ["i", "3", "2", "1/2", "1/3", "1/4", "1/5"]
FLUX_QUOT = np.asarray([0, 1/3, 1/2, 2, 3, 4, 5])
N_COLORS = np.flip([[l, l, l] for l in np.linspace(0, 1, len(FLUX_QUOT))])
R_COLORS = np.flip([[l, l, l] for l in np.linspace(0, 1, len(FLUX_QUOT))])

PATH_TO_FILE = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-1x/log-normal/"
SHADE_SCALE=25
LINE_COLOR = "C1"

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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=1.0)

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

xVals = [L_0_RANGE[0] * lOPowerStep**i for i in range(PLOT_SIZE)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / PLOT_SIZE * j for j in range(PLOT_SIZE)]


fig, ax = plt.subplots()
plt.xlim(left=L_0_RANGE[0], right=L_0_RANGE[1])
plt.ylim(bottom=SIGMA_L_RANGE[0], top=SIGMA_L_RANGE[1])

plt.xscale("log")
plt.ylabel("$\sigma$")
plt.xlabel("$L_0$ [erg / s]")


#plt.contourf(xVals, yVals, numSeen, NUM_PULSARS_ABOVE_THRESHOLD * FLUX_QUOT, colors=COLORS_NUM, alpha=0.4)
#plt.contourf(xVals, yVals, lumSeen, FRAC_ABOVE_THRESHOLD * FLUX_QUOT, colors=COLORS_FRAC, alpha=0.5)

n_contours = plt.contour(xVals, yVals, numSeen, NUM_PULSARS_ABOVE_THRESHOLD * FLUX_QUOT, cmap="Greys",
    linestyles="solid", linewidths=1, vmax=np.max(NUM_PULSARS_ABOVE_THRESHOLD * FLUX_QUOT) / 10)
r_contours = plt.contour(xVals, yVals, lumSeen, FRAC_ABOVE_THRESHOLD * FLUX_QUOT, cmap="Greys",
    linestyles="dashed", linewidths=1, vmax=np.max(FRAC_ABOVE_THRESHOLD * FLUX_QUOT) / 10)

plt.contour(xVals, yVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=LINE_COLOR, linestyles="solid", linewidths=2)
plt.contour(xVals, yVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2)
#

# Observation
#shade(numSeen, NUM_PULSARS_ABOVE_THRESHOLD, xVals, yVals)
#shade(lumSeen, FRAC_ABOVE_THRESHOLD, xVals, yVals, True)


# Final points

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', markersize=6)
plt.errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[0.41e34], [0.79e34]], yerr=[[0.16], [0.15]], linewidth=1, color=LINE_COLOR)
plt.plot(ploegPoint[0], ploegPoint[1], markeredgecolor='black', markerfacecolor="fuchsia", marker='s', markersize=6)
plt.errorbar([ploegPoint[0]], [ploegPoint[1]], xerr=[0.0044191e+32], yerr=[0.00121058], linewidth=1, color="fuchsia")
plt.plot(gautamPoint[0], gautamPoint[1], markeredgecolor='black', markerfacecolor="red", marker='*', markersize=8)
plt.errorbar([gautamPoint[0]], [gautamPoint[1]], xerr=[0.19540e+30], yerr=[0.01028865], linewidth=1, color="red")


custom_lines = [Line2D([0], [0], color=LINE_COLOR, linestyle="solid", lw=2),
                Line2D([0], [0], color=LINE_COLOR, linestyle='dashed', lw=2, dashes=(4, 2))]\
             + [Line2D([], [], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='o', linestyle='none', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="fuchsia", marker='s', linestyle='none', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="red", marker='*', linestyle='none', markersize=8),]
plt.legend(custom_lines, ['$N_r$', '$R_r$'] + ["GLC", "GCE", "AIC"], loc="lower left", prop={'size': 14})
plt.xlim(xVals[0], xVals[-1])
plt.ylim(yVals[0], yVals[-1])

fmt_n = {}
fmt_r = {}
for i, item in enumerate(n_contours.levels):
    fmt_n[item] = "$\\times$"+str_flux_quot[i]
for i, item in enumerate(r_contours.levels):
    fmt_r[item] = "$\\times$"+str_flux_quot[i]

n_labels = plt.clabel(n_contours, n_contours.levels, fmt=fmt_n, inline=True, fontsize=10, colors='k', manual=True)
for l in n_labels:
    l.set_rotation(0)

r_labels = plt.clabel(r_contours, r_contours.levels, fmt=fmt_r, inline=True, fontsize=10, colors='k', manual=True)
for l in r_labels:
    l.set_rotation(0)

plt.tight_layout()

# Save
plt.savefig("overlay-log-normal.pdf")

plt.show()
