import numpy as np
import parse_spectrum as ps
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtrans

plt.style.use("jcap")


SPECTRUM_RANGE = [0.1, 100]# GeV
SPECTRUM_RANGE_LOW = [0.1, 10]# GeV
NUM_PLOTS = 10
NUM_ROWS = 5
NUM_COLS = 2
SHAPES = ['o', 's', '*', 'd', 'd', 'o', 's', '<', '>', 'd']
FILL_STYLE = [None, None, 'none', None, 'none', 'none', 'none', 'none', 'none', 'none']
TEXT_SIZE=12

ERGS_PER_GEV = 0.00160218

fig, axes = plt.subplots(nrows=NUM_ROWS, ncols=NUM_COLS, figsize=(8, 10), sharex='col', sharey='row')

allfig, allax = plt.subplots(figsize=(6.8, 5.5))


allax.set_ylim(1e-9, 6e-7)

def sciNot(i):
    if i == None:
        return "None"
    if i == np.nan:
        return "NaN"
    if i == 0:
        return "0"
    negative = False
    if i < 0:
        negative = True
        i *= -1
    d = int(np.log10(i))
    if np.log10(i) < 0: d -= 1
    if not negative:
        return str(round(i / 10**d, 3)) + "e"+str(d)
    return "-"+str(round(i / 10**d, 3)) + "e"+str(d)

num_fluxes = []
num_flux_uncs = []
calore_fluxes = []
fit_fluxes = []
num_fluxes_low = []
num_flux_low_uncs = []
calore_fluxes_low = []
fit_fluxes_low = []
names = []

calore_flux_uncs = []
fit_flux_uncs = []
calore_flux_low_uncs = []
fit_flux_low_uncs = []

for i in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
    ax = axes[i%NUM_ROWS][i//NUM_ROWS]
    ax.set_ylim(1e-8, 2e-6)
    f = ps.Spectrum(i)
    names.append(f.get_name())
    f.display_data(allax, color="C"+str(i), shape=SHAPES[i],
        fill_color=FILL_STYLE[i])
    f.display_data(ax, title_size=TEXT_SIZE)
    f.display_calore(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    f.display_power_law(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    numerical_flux = f.get_numerical_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    calore_flux = f.get_calore_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    fit_flux = f.get_power_law_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    num_fluxes.append(numerical_flux[0] * ERGS_PER_GEV)
    calore_fluxes.append(calore_flux[0] * ERGS_PER_GEV)
    fit_fluxes.append(fit_flux[0] * ERGS_PER_GEV)
    num_flux_uncs.append(numerical_flux[1] * ERGS_PER_GEV)
    calore_flux_uncs.append(calore_flux[1] * ERGS_PER_GEV)
    fit_flux_uncs.append(fit_flux[1] * ERGS_PER_GEV)
    if i // NUM_ROWS == 0:
        ax.set_ylabel(f.get_y_label(), fontsize=TEXT_SIZE)
    if i % NUM_ROWS == NUM_ROWS - 1:
        ax.set_xlabel(f.get_x_label(), fontsize=TEXT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=TEXT_SIZE)

    print(f.get_name(), fit_flux[0] * ERGS_PER_GEV)

    #print("{}:\tNum {},\tCalore {},\tFit {}.".format(f.get_name(), numerical_flux[0] * ERGS_PER_GEV, calore_flux[0] * ERGS_PER_GEV, fit_flux[0] * ERGS_PER_GEV))

    f.label_axes(allax)

    '''ax.annotate("Num: {}\nCalore: {}\nFit: {}".format(sciNot(numerical_flux[0] * ERGS_PER_GEV),
        sciNot(calore_flux[0] * ERGS_PER_GEV),
        sciNot(fit_flux[0] * ERGS_PER_GEV)), (0.2, 0.1), xycoords='axes fraction', size=8)'''


    numerical_flux_low = f.get_numerical_flux(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], override=True)
    calore_flux_low = f.get_calore_flux(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], override=True)
    fit_flux_low = f.get_power_law_flux(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], override=True)
    num_fluxes_low.append(numerical_flux_low[0] * ERGS_PER_GEV)
    calore_fluxes_low.append(calore_flux_low[0] * ERGS_PER_GEV)
    fit_fluxes_low.append(fit_flux_low[0] * ERGS_PER_GEV)
    num_flux_low_uncs.append(numerical_flux_low[1] * ERGS_PER_GEV)
    calore_flux_low_uncs.append(calore_flux_low[1] * ERGS_PER_GEV)
    fit_flux_low_uncs.append(fit_flux_low[1] * ERGS_PER_GEV)

    #f.display_calore(ax, SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], linestyle=':')
    #f.display_power_law(ax, SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], linestyle=':')


#axes[(NUM_PLOTS-1)%NUM_ROWS][(NUM_PLOTS-1)//NUM_ROWS].set_xlabel(f.get_x_label(), fontsize=TEXT_SIZE)

fig.delaxes(axes[4][1])

while None in fit_fluxes:
    index = fit_fluxes.index(None)
    fit_fluxes[index] = 0
    #del fit_fluxes[index]
    #del calore_fluxes[index]
    #del num_fluxes[index]
    #del names[index]

x = np.arange(len(fit_fluxes))
width = 0.25
delta = 0#width/6

hist_fig, hist_ax = plt.subplots(figsize=(8, 5))
hist_ax.bar(x, num_fluxes, width, label="Numerical ({}-{} GeV)".format(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1]), color='purple', alpha=0.3)
hist_ax.bar(x + width, calore_fluxes, width, label="Ref. [5] fit", linewidth=0.5, alpha=0.5, color='r')
hist_ax.bar(x + 2 * width, fit_fluxes, width, label="Power law fit", linewidth=0.5, alpha=0.5, color='b')
hist_ax.errorbar(x + delta, num_fluxes, yerr=num_flux_uncs, linestyle="none", linewidth=1, color='purple')
hist_ax.errorbar(x + width + delta, calore_fluxes, yerr=calore_flux_uncs, linestyle="none", linewidth=1, color='r')
hist_ax.errorbar(x + 2 * width + delta, fit_fluxes, yerr=fit_flux_uncs, linestyle="none", linewidth=1, color='b')

hist_ax.bar(x, num_fluxes_low, width, label="{}-{} GeV".format(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1]), fill=False, edgecolor='k')
hist_ax.bar(x + width, calore_fluxes_low, width, linewidth=1, fill=False, edgecolor='k')
hist_ax.bar(x + 2 * width, fit_fluxes_low, width, linewidth=1, fill=False, edgecolor='k')
#hist_ax.errorbar(x + width - delta, calore_fluxes_low, yerr=calore_flux_low_uncs, linestyle="none", linewidth=1, color='k')
#hist_ax.errorbar(x + 2 * width - delta, fit_fluxes_low, yerr=fit_flux_low_uncs, linestyle="none", linewidth=1, color='k')

hist_ax.legend()
hist_ax.set_ylabel("$F_\\mathrm{GCE}$ [erg / cm$^2$ / s]")
hist_ax.set_xticks(x + width)
hist_ax.set_xticklabels(names, rotation=70, size=11)
hist_fig.tight_layout()
hist_fig.savefig("integral-hist.pdf")



calore_fit = mlines.Line2D([], [], color='red', linestyle='--', label='Ref. [5] fit')
power_law_fit = mlines.Line2D([], [], color='blue', linestyle='--', label='Power law fit')
fig.legend(handles=[calore_fit, power_law_fit], loc='lower right')


allax.legend(loc='lower center', ncol=2)
allax.set_title("")

fig.tight_layout()
allfig.tight_layout()


fig.savefig("summary.pdf")
allfig.savefig("all-spectra.pdf")
plt.show()
