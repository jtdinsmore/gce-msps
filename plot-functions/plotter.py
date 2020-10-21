from matplotlib import pyplot as plt
from math import pi
from scipy.special import gammainc, gamma
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * np.exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

LUM_TO_FLUX = 1.1095246594108431e-46
TOTAL_FLUX = 1.794113925439598e-09

NUM_PLOT_POINTS = 200
PLOT_LOG_MIN = 29.0
PLOT_LOG_MAX = 38.0

ALPHA = 1.94
L_MIN = 1.0e29
L_MAX = 1.0e35
L_THRESH = 10e34

L0_LOG_NORMAL = 8.8e33
SIGMA_LOG_NORMAL = 0.62

L0_PLOEG_FIT = 1.3023e+32
SIGMA_PLOEG_FIT = 0.69550156

L0_GAUTAM_FIT = 4.2970e+30
SIGMA_GAUTAM_FIT = 0.93936155

N_BELOW_NFW = -0.66
N_ABOVE_NFW = 18.2
L_BREAK_NFW = 2.5389429e+34

def powerLaw(x, lmin, lmax, alpha):
    return x**-alpha * np.exp(-x / lmax) / (Gamma(1-alpha, lmin/lmax) * lmax**(1-alpha))
def logNormal(x, L0, sigma):
    return np.log10(np.exp(1)) / (sigma * np.sqrt(2 * np.pi) * x)* np.exp(-(np.log10(x) - np.log10(L0))**2 / (2 * sigma**2))
def nptf(x, n1, n2, LBreak):
    if x < LBreak:
        return (1 - n1) * (1 - n2) / (LBreak * (n1 - n2)) * (x/LBreak)**(-n1)
    else:
        return (1 - n1) * (1 - n2) / (LBreak * (n1 - n2)) * (x/LBreak)**(-n2)
def lum_to_flux(lum):
    return lum * LUM_TO_FLUX
def flux_to_lum(flux):
    return flux / LUM_TO_FLUX

fig, ax1=plt.subplots()

ax1.set_xlabel("Luminosity $L$ [erg / s]")
ax1.set_ylabel("$dN/dL$")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(1e-45, 1e-27)

x = 10**np.linspace(PLOT_LOG_MIN, PLOT_LOG_MAX, NUM_PLOT_POINTS)

ys = [
    powerLaw(x, L_MIN, L_MAX, ALPHA),
    powerLaw(x, 1e29, 7e34, 1.5),
    logNormal(x, L0_LOG_NORMAL, SIGMA_LOG_NORMAL),
    logNormal(x, L0_PLOEG_FIT, SIGMA_PLOEG_FIT),
    logNormal(x, L0_GAUTAM_FIT, SIGMA_GAUTAM_FIT),
    [nptf(l, N_BELOW_NFW, N_ABOVE_NFW, L_BREAK_NFW) for l in x],
    [nptf(l, 0.97, 2.6, 1.7e33) for l in x],
    #[nptf(l, N_BELOW_DISK, N_ABOVE_DISK, L_BREAK_DISK) for l in x],
]
names = ["Wavelet 1", "Wavelet 2", "GLC", "GCE", "AIC", "NPTF", "Disk"]
styles = ["solid", "solid", "dashed", "dashed", "dashed", "dotted", "dotted",]
colors = ["aquamarine", "darkgreen", "C1", "fuchsia", "red", "black", "gray"]

for i in range(len(ys)):
    ax1.plot(x, np.abs(ys[i]), label = names[i], linestyle=styles[i],
        c=colors[i])
#ax1.axvline(x=L_THRESH, color='k')

ax2 = ax1.secondary_xaxis("top", functions=(lum_to_flux, flux_to_lum))
ax2.set_xlabel("Flux $F$ [erg / cm$^2$ / s]", labelpad=10)

plt.tight_layout()
ax1.legend(loc='upper right', ncol=2)
ax1.set_xlim(10**PLOT_LOG_MIN, 10**PLOT_LOG_MAX)
ax2.set_xlim(10**PLOT_LOG_MIN * LUM_TO_FLUX, 10**PLOT_LOG_MAX * LUM_TO_FLUX)
plt.savefig("lum-funcs.pdf")
plt.show()



fig, ax1 = plt.subplots()

ax1.set_xlabel("Luminosity $L$ [erg / s]")
ax1.set_ylabel("$LdN/dL$ [erg / s]")
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1e-7, 1e+3)

ys = [x * y for y in ys]

for i in range(len(ys)):
    ax1.plot(x, np.abs(ys[i]), label = names[i], linestyle=styles[i],
    c=colors[i])
#ax1.axvline(x=L_THRESH, color='k')

ax2 = ax1.secondary_xaxis('top', functions=(lum_to_flux, flux_to_lum))
ax2.set_xlabel("Flux $F$ [erg / cm$^2$ / s]", labelpad=10)

ax1.legend(loc='upper left', ncol=3)
plt.tight_layout()
ax1.set_xlim(10**PLOT_LOG_MIN, 10**PLOT_LOG_MAX)
ax2.set_xlim(10**PLOT_LOG_MIN * LUM_TO_FLUX, 10**PLOT_LOG_MAX * LUM_TO_FLUX)

plt.savefig("l-lum-funcs.pdf")
plt.show()





fig, ax1 = plt.subplots()


ax1.set_xlabel("Luminosity $L$ [erg / s]")
ax1.set_ylabel("$L^2dN/dL$ [erg$^2$ / s$^2$]")
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1e26, 1e36)

ys = [x * y for y in ys]

for i in range(len(ys)):
    ax1.plot(x, np.abs(ys[i]), label = names[i], linestyle=styles[i],
    c=colors[i])
#ax1.axvline(x=L_THRESH, color='k')
ax2 = ax1.secondary_xaxis('top', functions=(lum_to_flux, flux_to_lum))
ax2.set_xlabel("Flux $F$ [erg / cm$^2$ / s]", labelpad=10)

ax1.legend(loc='upper left')
plt.tight_layout()
ax1.set_xlim(10**PLOT_LOG_MIN, 10**PLOT_LOG_MAX)
ax2.set_xlim(10**PLOT_LOG_MIN * LUM_TO_FLUX, 10**PLOT_LOG_MAX * LUM_TO_FLUX)

plt.savefig("l2-lum-funcs.pdf")
plt.show()
