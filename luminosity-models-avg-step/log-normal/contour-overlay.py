from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')

DIM_TRIALS=100

L_EXCESS = 1.794113925439598e-09 / 1.1095246594108431e-46#  All units are in erg per second
L_THRESH = 3.4345358338912146e+34
L_0_RANGE=[1.0e30, 2.0e36]#[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]
powerStep =(L_0_RANGE[1] / L_0_RANGE[0])**(1/DIM_TRIALS)


TOTAL_FGL_NUMBER = 265 # 109 unflagged
PERCENTAGES = [0.05, 0.1, 0.2, 0.4, 1]
FRACS_ABOVE_THRESHOLD=[0.05, 0.1, 0.2, 0.4]

LINE_COLOR = "C1"
HALF_COLOR = "#cf5f0c"
DRAW_EXTRA_CONTOURS =False
DRAW_PLOEG_POINT = True

paperPoint = [0.88e34, 0.62]
ploegPoint = [1.3023e+32, 0.69550156]
gautamPoint = [4.2970e+30, 0.93936155]
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

def getTotalLum(L0, sigma):
    return exp(0.5 * sigma**2 * log(10)**2) * L0

def getNumPulsars(L0, sigma):
    scale = L_EXCESS / getTotalLum(L0, sigma)
    return scale * 1 # scale times number of pulsars per unit

def getNumPulsarsAboveThreshold(L0, sigma):
    scale = L_EXCESS / getTotalLum(L0, sigma)
    return scale * 0.5 * erfc((log(L_THRESH) - log(L0)) / (sqrt(2) * sigma * log(10)))

def getFracLumAboveThreshold(L0, sigma):
    # return (lum above thresh) / (total lum)
    erfArgument = sqrt(2) * (sigma**2 * log(10)**2 - log(L_THRESH) + log(L0)) / (sigma * log(100))
    amountAbove = 0.5 * (1 + erf(erfArgument))
    return amountAbove # See the mathematica notebook for a derivation


def getMinNumPulsarsInTriangle():
    EXTRA_DENSITY = 2
    minSoFar = -1
    powerStepNew =(L_0_RANGE[1] / L_0_RANGE[0])**(1/(DIM_TRIALS * EXTRA_DENSITY))
    for j in range(DIM_TRIALS * EXTRA_DENSITY):
        sigma= SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / (DIM_TRIALS * EXTRA_DENSITY) * j
        for i in range(DIM_TRIALS * EXTRA_DENSITY):
            L0=L_0_RANGE[0] * powerStepNew**i
            nAbove = getNumPulsarsAboveThreshold(L0, sigma)
            if nAbove > NUM_PULSARS_ABOVE_THRESHOLD:
                continue
            fracAbove = getFracLumAboveThreshold(L0, sigma)
            if fracAbove > FRAC_ABOVE_THRESHOLD:
                continue
            n = getNumPulsars(L0, sigma)
            if(n < minSoFar or minSoFar < 0):
                minSoFar = n

    print("Minimum number of total pulsars in allowed triangle: ", minSoFar)
    return minSoFar

def getMinPulsarsWithinOneStdevOfSigma():
    SHORT_SIGMA_RANGE = [0.62 - 0.16, 0.62 + 0.15]
    EXTRA_DENSITY = 2
    minPulsars = -1
    minL0 = 0
    minSigma = 0
    powerStepNew =(L_0_RANGE[1] / L_0_RANGE[0])**(1/(DIM_TRIALS * EXTRA_DENSITY))
    for j in range(DIM_TRIALS * EXTRA_DENSITY):
        sigma= SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / (DIM_TRIALS * EXTRA_DENSITY) * j
        if sigma < SHORT_SIGMA_RANGE[0] or SHORT_SIGMA_RANGE[1] < sigma:
            continue
        tooFewPulsarsAboveThreshold = False
        for i in range(DIM_TRIALS * EXTRA_DENSITY):
            L0=L_0_RANGE[0] * powerStepNew**i

            if sigma > 0.7 and L0 > 2e33:
                continue # Cut out top green line

            tooFewPulsarsAboveThresholdNow = getNumPulsarsAboveThreshold(L0, sigma) < NUM_PULSARS_ABOVE_THRESHOLD
            if tooFewPulsarsAboveThresholdNow ==(not tooFewPulsarsAboveThreshold) and i > 0:
                # We have just crossed over the green line
                numPulsarsNow = getNumPulsars(L0, sigma)
                if numPulsarsNow < minPulsars or minPulsars < 0:
                    minPulsars = numPulsarsNow
                    minL0 = L0
                    minSigma = sigma
            tooFewPulsarsAboveThreshold = tooFewPulsarsAboveThresholdNow

    print("Fewest possible pulsars required to hit the green line with sigma in 1 stdev of paper values: {0} at coordinates L_0={1}, sigma={2}".format(minPulsars, minL0, minSigma))
    return (minL0, minSigma)

def getPaperPointInfo(name, L0 = paperPoint[0], sigma = paperPoint[1]):
    print("""{0} point info:
    Coordinates: L0 = {1}, sigma = {2}
    Total num pulsars: {3}
    Num pulsars above threshold: {4}
    Fraction luminosity above threshold: {5}""".format(name, L0, sigma, getNumPulsars(L0, sigma),
    getNumPulsarsAboveThreshold(L0, sigma), getFracLumAboveThreshold(L0, sigma)))


#getMinNumPulsarsInTriangle()
#minPoint = getMinPulsarsWithinOneStdevOfSigma()
getPaperPointInfo("Paper (GLC)")
getPaperPointInfo("Ploeg (GCE)", ploegPoint[0], ploegPoint[1])
getPaperPointInfo("Gautam (AIC)", gautamPoint[0], gautamPoint[1])


numPulsars = []
numAboveThreshold = []
fracAboveThreshold = []
for j in range(DIM_TRIALS):
    lineNumPulsars = []
    lineNumAboveThreshold = []
    lineFracAboveThreshold = []
    for i in range(DIM_TRIALS):
        L0 = L_0_RANGE[0] * powerStep**i
        sigma =  SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j

        numNow = getNumPulsars(L0, sigma)
        numAbove = getNumPulsarsAboveThreshold(L0, sigma)
        fracAbove = getFracLumAboveThreshold(L0, sigma)

        lineNumPulsars.append(numNow)
        lineNumAboveThreshold.append(numAbove)
        lineFracAboveThreshold.append(fracAbove)

    numPulsars.append(lineNumPulsars)
    numAboveThreshold.append(np.asarray(lineNumAboveThreshold))
    fracAboveThreshold.append(np.asarray(lineFracAboveThreshold))

numAboveThreshold = np.stack(numAboveThreshold)
fracAboveThreshold = np.stack(fracAboveThreshold)


xVals = [L_0_RANGE[0] * powerStep**i for i in range(DIM_TRIALS)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j for j in range(DIM_TRIALS)]


fig, ax = plt.subplots()
plt.xlim(left=L_0_RANGE[0], right=L_0_RANGE[1])
plt.ylim(bottom=SIGMA_L_RANGE[0], top=SIGMA_L_RANGE[1])

plt.xscale("log")
plt.ylabel("$\sigma$")
plt.xlabel("$L_0$ [erg / s]")

cols = colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars]))
c1 = plt.contourf(xVals, yVals, numPulsars,
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars])), cmap='Greys_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N_\\textrm{GCE}$")

#plt.contour(xVals, yVals, numAboveThreshold, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR])
#plt.contour(xVals, yVals, fracAboveThreshold, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed')

n_contours = plt.contour(xVals, yVals, numAboveThreshold, np.array(PERCENTAGES) * TOTAL_FGL_NUMBER, colors=[LINE_COLOR], linewidths=[1])
r_contours = plt.contour(xVals, yVals, fracAboveThreshold, FRACS_ABOVE_THRESHOLD, colors=[LINE_COLOR], linestyles='dashed', linewidths=[1])

fmt_n = {}
for i, item in enumerate(n_contours.levels):
    fmt_n[item] = str(int(100*PERCENTAGES[i])) + "\%"
fmt_r = {}
for i, item in enumerate(r_contours.levels):
    fmt_r[item] = str(int(100*PERCENTAGES[i])) + "\%"

n_labels = plt.clabel(n_contours, n_contours.levels, fmt=fmt_n, inline=True, fontsize=10, colors='k', manual=True)
r_labels = plt.clabel(r_contours, r_contours.levels, fmt=fmt_r, inline=True, fontsize=10, colors='k', manual=True)
#[txt.set_backgroundcolor('white') for txt in n_labels]
#[txt.set_backgroundcolor('white') for txt in r_labels]


plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', markersize=6)
plt.errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[0.41e34], [0.79e34]], yerr=[[0.16], [0.15]], linewidth=1, color=LINE_COLOR)
plt.plot(ploegPoint[0], ploegPoint[1], markeredgecolor='black', markerfacecolor="fuchsia", marker='s', markersize=6)
plt.errorbar([ploegPoint[0]], [ploegPoint[1]], xerr=[1.3878e+30], yerr=[0.00191187], linewidth=1, color="fuchsia")
plt.plot(gautamPoint[0], gautamPoint[1], markeredgecolor='black', markerfacecolor="red", marker='*', markersize=8)
plt.errorbar([gautamPoint[0]], [gautamPoint[1]], xerr=[1.9543e+29], yerr=[0.01028865], linewidth=1, color="red")
#plt.scatter(minPoint[0], minPoint[1], c='cyan')


# Observation
shade(numAboveThreshold, PERCENTAGES[0] * TOTAL_FGL_NUMBER, xVals, yVals)
shade(fracAboveThreshold, FRACS_ABOVE_THRESHOLD[0], xVals, yVals, True)



custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="fuchsia", marker='s', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="red", marker='*', linestyle='None', markersize=8),]
plt.legend(custom_lines, ['$N_r$', '$R_r$', "GLC", "GCE", "AIC"], loc="lower left", prop={'size': 14})
plt.xlim(xVals[0], xVals[-1])
plt.ylim(yVals[0], yVals[-1])

plt.tight_layout()

# Save
plt.savefig("log-normal-step.pdf")

plt.show()
