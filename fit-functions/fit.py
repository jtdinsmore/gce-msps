import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameter, report_fit
import scipy.optimize as optimization
from scipy import stats
import lmfit

plt.style.use("jcap")

LUM_TO_FLUX = 1.1095246594108431e-46
NUM_POINTS = 100

def power_law(x, c, alpha, lmax):
    return np.log10(c) + (1-alpha) * x - np.log10(np.exp(1)) * 10.0**x / lmax

def log_normal(x, c, l0, sigma):# log(dN/dE)
    return np.log10(c) - np.log10(np.exp(1)) * (x - np.log10(l0))**2 / (2 * sigma**2)

def power_law_params(x, params):
    return power_law(x, params["c"], params["alpha"], params["lmax"])

def log_normal_params(x, params):# log(dN/dE)
    return log_normal(x, params["c"], params["l0"], params["sigma"])

def chi_squared_asym(params, datax, datay, weights_up, weights_down, fn):
    expected = fn(datax, params)
    weights = np.where(datay > expected, weights_down, weights_up)
    return (expected - datay) * weights

class Data:
    def __init__(self, name, filename, bar_bottom, bar_top, convert_x_to_log=False,
        convert_y_to_log=False, convert_to_lum=False, convert_to_ergs=False):

        self.name = name
        xs, ys = self.get_data(filename, convert_x_to_log,
            convert_y_to_log, convert_to_lum, convert_to_ergs)
        xdown, ydown = self.get_data(bar_bottom, convert_x_to_log,
            convert_y_to_log, convert_to_lum, convert_to_ergs)
        xup, yup = self.get_data(bar_top, convert_x_to_log,
            convert_y_to_log, convert_to_lum, convert_to_ergs)
        min_x = max(np.min(xs), np.min(xdown), np.min(xup))
        max_x = min(np.max(xs), np.max(xdown), np.max(xup))
        if min_x >= max_x:
            raise Exception("Bars and points share no common domain")
        self.xs = np.linspace(min_x, max_x, 100)
        self.ys = np.asarray([self.at(x, xs, ys) for x in self.xs])
        self.ydown = self.ys - np.asarray([self.at(x, xdown, ydown)
            for x in self.xs])
        self.yup = np.asarray([self.at(x, xup, yup) for x in self.xs]) - self.ys
        self.log_normal_results = None
        self.power_law_results = None

    def plot_data(self, ax):
        ax.plot(self.xs, self.ys, color='k')
        ax.fill_between(self.xs, self.ys-self.ydown, self.ys + self.yup, alpha=0.3, color="k", label="Data")

    def plot_fit(self, ax):
        if self.power_law_results is not None:
            ax.plot(self.xs, power_law(self.xs, self.power_law_results[0],
            self.power_law_results[2],self.power_law_results[4]), color='C0',
            label="Power law fit", linestyle='dotted')
        if self.log_normal_results is not None:
            ax.plot(self.xs, log_normal(self.xs, self.log_normal_results[0],
            self.log_normal_results[2],self.log_normal_results[4]), color='C1',
            label="Log normal fit", linestyle='dashed')

    def label_axes(self, ax):
        ax.set_xlabel("$\\log_{10}(L / $ erg s$^{-1})$")
        ax.set_ylabel("$\\log_{10}\\left(dN / dL\\right)+C$")
        #ax.set_title(self.name)

    def plot(self):
        fig, ax = plt.subplots()
        self.plot_data(ax)
        self.plot_fit(ax)
        self.label_axes(ax)
        ax.legend(prop={'size':16})
        fig.tight_layout()
        return fig

    def fit_power_law(self, c=1, alpha=1.3, lmax=1e31, display=True):
        params = lmfit.Parameters()
        params.add(Parameter('c', value=c, min=0))
        params.add(Parameter('alpha', value=alpha))
        params.add(Parameter('lmax', value=lmax, min=0))

        weights_down = 1 / self.ydown
        weights_up= 1 / self.yup

        result = lmfit.minimize(chi_squared_asym, params, args=(self.xs, self.ys, weights_up, weights_down, power_law_params))

        if display:
            report_fit(result)

        chisq = np.sum(chi_squared_asym(result.params, self.xs, self.ys, weights_up, weights_down, power_law_params)**2)

        self.power_law_results = (
            result.params["c"].value, result.params["c"].stderr,
            result.params["alpha"].value, result.params["alpha"].stderr,
            result.params["lmax"].value, result.params["lmax"].stderr)
        self.power_law_p = 1 - stats.chi2.cdf(chisq, len(self.ys)-3)
        self.power_law_chisq = chisq

        return self.power_law_results


    def fit_log_normal(self, c=1e40, l0=1e32, sigma=1, display=True):
        params = lmfit.Parameters()
        params.add(Parameter('c', value=c, min=0))
        params.add(Parameter('l0', value=l0, min=1))
        params.add(Parameter('sigma', value=sigma, min=1e-3))

        weights_down = 1 / self.ydown
        weights_up= 1 / self.yup

        result = lmfit.minimize(chi_squared_asym, params, args=(self.xs, self.ys, weights_up, weights_down, log_normal_params))

        if display:
            report_fit(result)

        chisq = np.sum(chi_squared_asym(result.params, self.xs, self.ys, weights_up, weights_down, log_normal_params)**2)

        self.log_normal_results = (
            result.params["c"].value, result.params["c"].stderr,
            result.params["l0"].value, result.params["l0"].stderr,
            result.params["sigma"].value, result.params["sigma"].stderr)
        self.log_normal_p = 1 - stats.chi2.cdf(chisq, len(self.ys)-3)
        self.log_normal_chisq = chisq

        return self.log_normal_results


    def get_data(self, filename, convert_x_to_log, convert_y_to_log,
        convert_to_lum, convert_to_ergs):

        xs = []
        ys = []
        f = open(filename, 'r')
        for line in f.readlines():
            if line == '':
                continue
            x, y = line.split(', ')
            x = float(x)
            y = float(y)
            if convert_x_to_log:
                x = np.log10(x)
            if convert_y_to_log:
                if y <= 0: continue
                y = np.log10(y)
            if convert_to_lum:
                x -= np.log10(LUM_TO_FLUX)
            if convert_to_ergs:
                raise Exception("Unimplemented")
            xs.append(x)
            ys.append(y)
        return np.asarray(xs), np.asarray(ys)

    def at(self, x, xs, ys):
        if x < np.min(xs):
            return 0
        if x > np.max(xs):
            return 0
        i = 0
        while xs[i+1] < x:
            i += 1

        return (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i]) + ys[i]

gautam = Data("Gautam", "new-gautam/points.csv", "new-gautam/down-bands.csv", "new-gautam/up-bands.csv",
    convert_to_lum=True)
gautam.fit_power_law()
gautam.fit_log_normal()
print(gautam.power_law_chisq - gautam.log_normal_chisq)
print(gautam.log_normal_chisq / len(gautam.xs))
print(gautam.log_normal_p)
fig = gautam.plot()
fig.savefig("gautam-fit.pdf")

ploeg = Data("Ploeg", "ploeg/points.csv", "ploeg/bars-bottom.csv", "ploeg/bars-top.csv",
    convert_y_to_log=True)
ploeg.fit_power_law()
ploeg.fit_log_normal()
print(ploeg.power_law_chisq - ploeg.log_normal_chisq)
print(ploeg.log_normal_chisq / len(ploeg.xs))
print(ploeg.log_normal_p)

fig = ploeg.plot()
fig.savefig("ploeg-fit.pdf")
plt.show()
