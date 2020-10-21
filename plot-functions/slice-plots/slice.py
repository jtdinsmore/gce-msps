import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gammainc, gamma
import matplotlib.colors as mc
import colorsys


NUM_PLOTS = 20
PLOT_LOG_MIN = 29.0
PLOT_LOG_MAX = 38.0
L_THRESH = 3.56e34

plt.style.use('jcap')

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * np.exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def powerLaw(x, lmin, lmax, alpha):
    return x**-alpha * np.exp(-x / lmax) / (Gamma(1-alpha, lmin/lmax) * lmax**(1-alpha))

def logNormal(x, L0, sigma):
    return np.log10(np.exp(1)) / (sigma * np.sqrt(2 * np.pi) * x)* np.exp(-(np.log10(x) - np.log10(L0))**2 / (2 * sigma**2))

def nptf(x, n1, n2, LBreak):
    if x < LBreak:
        return (1 - n1) * (1 - n2) / (LBreak * (n1 - n2)) * (x/LBreak)**(-n1)
    else:
        return (1 - n1) * (1 - n2) / (LBreak * (n1 - n2)) * (x/LBreak)**(-n2)
        

def lighten_color(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def plot_linear(fn, param_low, param_high, param_label, color):
    fig, ax = plt.subplots(figsize=(4,6))
    lums = 10**np.linspace(PLOT_LOG_MIN, PLOT_LOG_MAX, 200)
    param = 10**np.linspace(param_low, param_high, NUM_PLOTS)
    values = np.array([fn(lums, p) for p in param])
    
    high_values = values[:,lums>L_THRESH]
    plot_height = np.nanmax(high_values) * NUM_PLOTS / 5

    values *= (np.max(param) - np.min(param)) / plot_height
    values += np.min(param)
    high_values = values[:,lums>L_THRESH]

    # Scale the output so that the image has proper units for param.

    for i, p in enumerate(param):
        offset = plot_height * (i / NUM_PLOTS)
        line = values[i] + offset
        ax.plot(lums, line, lw=1, c=color, alpha=1, zorder=-i/NUM_PLOTS)
        ax.fill_between(lums, np.min(param), line, facecolor='white', alpha = 0.15, zorder=-i/NUM_PLOTS)
        ax.fill_between(lums[lums>L_THRESH], offset, line[lums>L_THRESH],
        facecolor=lighten_color(color), alpha=1, zorder=-i/NUM_PLOTS)
    
    #ax.set_yticks([])
    ax.set_ylabel(param_label)
    ax.set_xlabel("$L$ (erg/s)")
    ax.set_xscale('log')
    ax.set_ylim(np.min(param) / plot_height ** (1 / NUM_PLOTS), np.max(param) * (np.max(high_values[-1]) - np.max(param)))
    ax.set_xlim(np.min(lums), np.max(lums))
    ax.axvline(L_THRESH, color='k', linestyle='dashed', lw=1)
    fig.tight_layout()

def plot_log(ax, fn, param_low, param_high, param_label, color):
    lums = 10**np.linspace(PLOT_LOG_MIN, PLOT_LOG_MAX, 200)
    param = 10**np.linspace(param_low, param_high, NUM_PLOTS)
    values = np.array([fn(lums, p) for p in param])
    zero_value = 1 / np.nanmean(values)
    values /= np.nanmean(values)
    values = np.exp(values)
    values -= np.nanmin(values)
    
    high_values = values[:,lums>L_THRESH]
    plot_height = np.max(param) / np.min(param)

    #values /= np.nanmax(high_values) / plot_height
    # Eventual difference will be plot_height ** (i / NUM_PLOTS)
    # Want this to be equal to np.nanmax(high_values)
    # Set 
    factor = None
    for i, row in enumerate(high_values):
        space = plot_height ** ((i+1) / NUM_PLOTS) / plot_height ** (i / NUM_PLOTS)
        new_factor = (np.min(param) * (space - 1)) / np.nanmax(row)
        if factor is None:
            factor = new_factor
        else: 
            factor = min(factor, new_factor)
    values *= factor * 5
    values += np.min(param)

    # Scale the output so that the image has proper units for param.

    for i, _ in enumerate(param):
        offset = plot_height ** (i / NUM_PLOTS)
        line = values[i] * offset
        ax.plot(lums, line, lw=1, c=color, alpha=1, zorder=-i/NUM_PLOTS)
        ax.fill_between(lums, np.min(param), line, facecolor='white', alpha = 0.17, zorder=-i/NUM_PLOTS)
        ax.fill_between(lums[lums>L_THRESH], np.min(param) * offset, line[lums>L_THRESH],
            facecolor=lighten_color(color), alpha=1, zorder=-i/NUM_PLOTS)
    
    #ax.set_yticks([])
    ax.set_ylabel(param_label)
    #ax.set_xlabel("$L$ (erg/s)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(np.min(param) / plot_height ** (1 / NUM_PLOTS), plot_height * np.max(values[-1, lums>L_THRESH]))
    ax.set_xlim(np.min(lums), np.max(lums))
    ax.set_xticks([])
    ax.axvline(L_THRESH, color='k', linestyle='dashed', lw=1)


def plot_one(ax, fn, param, color, label=False, norm=None):
    lums = 10**np.linspace(PLOT_LOG_MIN, PLOT_LOG_MAX, 200)
    values = fn(lums, 10**param)
    high_values = values[lums>L_THRESH]
    if norm is None:
        norm = np.max(high_values)

    ax.plot(lums, values, color)
    ax.set_ylabel("$dN/dL$")
    ax.set_xscale('log')
    ax.set_xlim(np.min(lums), np.max(lums))
    ax.set_ylim(0, 1.2 * norm)
    if label:
        ax.set_xlabel("$L$ (erg/s)")
    else:
        ax.set_xticks([])
    ax.fill_between(lums[lums>L_THRESH], values[lums>L_THRESH], facecolor=lighten_color(color))
    ax.set_yticks([])
    ax.axvline(L_THRESH, color='k', linestyle='dashed', lw=1)
    return norm

fig = plt.figure(constrained_layout=True, figsize=(4, 8))
grids = fig.add_gridspec(3, 1, height_ratios = [1, 10, 1])
normpl = plot_one(fig.add_subplot(grids[0]), lambda x,p: powerLaw(x, 1e31, p, 1.94), 38.0, "C0")
plot_log(fig.add_subplot(grids[1]), lambda x,p: powerLaw(x, 1e31, p, 1.94), 34, 38, "$L_\\mathrm{max}$ (erg/s)", "C0")
plot_one(fig.add_subplot(grids[2]), lambda x,p: powerLaw(x, 1e31, p, 1.94), 34.0, "C0", label=True, norm=normpl)
plt.savefig("powerlaw.png")
plt.savefig("powerlaw.pdf")

fig = plt.figure(constrained_layout=True, figsize=(4, 8))
grids = fig.add_gridspec(3, 1, height_ratios = [1, 10, 1])
normln = plot_one(fig.add_subplot(grids[0]), lambda x,p: logNormal(x, p, 0.5), 36.0, "C1")
plot_log(fig.add_subplot(grids[1]), lambda x,p: logNormal(x, p, 0.5), 30, 36, "$L_0$ (erg/s)", "C1")
plot_one(fig.add_subplot(grids[2]), lambda x,p: logNormal(x, p, 0.5), 30.0, "C1", label=True, norm=normln)
plt.savefig("lognormal.png")
plt.savefig("lognormal.pdf")
plt.show()