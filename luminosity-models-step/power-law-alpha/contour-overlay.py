from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')

POWER_STEP = 1.1 # 1 is the minimum

L_EXCESS = 1.794113925439598e-09 / 1.1095246594108431e-46# 6.756e36  # All units are in erg per second
L_THRESH = 1.0e34
L_MIN = 1e29
ALPHA_RANGE = [1.1, 2.5]
L_MAX_RANGE = [1.0e34, 1.0e38]

TOTAL_FGL_NUMBER = 265 # 109 unflagged
PERCENTAGES = [0.05, 0.1, 0.2, 0.4, 1]
FRACS_ABOVE_THRESHOLD=[0.05, 0.1, 0.2, 0.4]

DRAW_EXTRA_CONTOURS = False
LINE_COLOR = "C2"
SHADE_SCALE=25

dimAlpha= 50
dimMax= int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))

paperPoint = [1.94, 1e35]
bartels15Point = [1.5, 7e34]

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


def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def getNumPulsars(alpha, lMax):
    lumExp = lMax**(2-alpha) * Gamma(2-alpha, L_MIN / lMax)
    Aexp = L_EXCESS / lumExp
    Nexp = Aexp *  lMax**(1-alpha) * Gamma(1-alpha, L_MIN / lMax)

    return Nexp

def getNumPulsarsAboveThreshold(alpha, lMax):
    lumExp = lMax**(2-alpha) * Gamma(2-alpha, L_MIN / lMax)
    Aexp = L_EXCESS / lumExp
    nAbove = Aexp *  lMax**(1-alpha) * Gamma(1-alpha, L_THRESH / lMax)
    return nAbove

def getFracLumAboveThreshold(alpha, lMax):
    # return (lum above thresh) / totalLum
    fracAbove = Gamma(2-alpha, L_THRESH / lMax) / (Gamma(2-alpha, L_MIN / lMax))
    return fracAbove


numPulsars = []
for i in range(dimAlpha):
    line = []
    for j in range(dimMax):
        nExp = getNumPulsars(ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha), L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    numPulsars.append(np.asarray(line))

numAboveThreshold = []
for i in range(dimAlpha):
    line = []
    for j in range(dimMax):
        nExp = getNumPulsarsAboveThreshold(ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha), L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    numAboveThreshold.append(np.asarray(line))

fracAboveThreshold = []
for i in range(dimAlpha):
    line = []
    for j in range(dimMax):
        nExp = getFracLumAboveThreshold(ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha), L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    fracAboveThreshold.append(np.asarray(line))

numPulsars = np.transpose(np.stack(numPulsars))
numAboveThreshold = np.transpose(np.stack(numAboveThreshold))
fracAboveThreshold = np.transpose(np.stack(fracAboveThreshold))



alphaVals = [ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha) for i in range(dimAlpha)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]


print("""Bartels 15 point info:
    Coordinates: L0 = {0}, sigma = {1}
    Total num pulsars: {2}
    Num pulsars above threshold: {3}
    Fraction luminosity above threshold: {4}""".format(bartels15Point[0], bartels15Point[1], getNumPulsars(bartels15Point[0], bartels15Point[1]),
    getNumPulsarsAboveThreshold(bartels15Point[0], bartels15Point[1]), getFracLumAboveThreshold(bartels15Point[0], bartels15Point[1])))




fig, ax = plt.subplots()

plt.yscale("log")
plt.xlabel("$\\alpha$")
plt.ylabel("$L_\\mathrm{max}$ [erg / s]")

c1 = plt.contourf(alphaVals, lMaxVals, numPulsars,
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars])), cmap='Greys_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N_\\mathrm{GCE}$")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(alphaVals, lMaxVals, numAboveThreshold, [10*i for i in range(1, 20)],
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)])
n_contours = plt.contour(alphaVals, lMaxVals, numAboveThreshold, TOTAL_FGL_NUMBER * np.array(PERCENTAGES), colors=[LINE_COLOR], linewidths=[1])

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(alphaVals, lMaxVals, fracAboveThreshold, [0.1*i for i in range(1, 15)],
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)])
r_contours = plt.contour(alphaVals, lMaxVals, fracAboveThreshold, FRACS_ABOVE_THRESHOLD, colors=[LINE_COLOR], linestyles='dashed', linewidths=[1])

fmt_n = {}
for i, item in enumerate(n_contours.levels):
    fmt_n[item] = str(int(100*PERCENTAGES[i])) + "\%"
fmt_r = {}
for i, item in enumerate(r_contours.levels):
    fmt_r[item] = str(int(100*PERCENTAGES[i])) + "\%"
n_labels = plt.clabel(n_contours, n_contours.levels, fmt=fmt_n, inline=True, fontsize=10, colors='k', manual=True)
r_labels = plt.clabel(r_contours, r_contours.levels, fmt=fmt_r, inline=True, fontsize=10, colors='k', manual=True)

# Observation
shade(numAboveThreshold, TOTAL_FGL_NUMBER * PERCENTAGES[0], alphaVals, lMaxVals)
shade(fracAboveThreshold, FRACS_ABOVE_THRESHOLD[0], alphaVals, lMaxVals, True)



plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor='aquamarine', marker='^', markersize=6)
plt.plot(bartels15Point[0], bartels15Point[1], markeredgecolor='black', markerfacecolor='darkgreen', marker='v', markersize=6)

custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor='aquamarine', marker='^', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor='darkgreen', marker='v', linestyle='None', markersize=6),]
plt.legend(custom_lines, ["$N_r$", "$R_r$", "Wavelet 1", "Wavelet 2"], loc='upper right', prop={'size': 14})
plt.xlim(alphaVals[0], alphaVals[-1])
plt.ylim(lMaxVals[0], lMaxVals[-1])
plt.tight_layout()

plt.savefig("power-law-alpha-step.pdf")

plt.show()
