import numpy as np
from scipy.optimize import curve_fit
import lmfit, os
from scipy import stats

import warnings
warnings.filterwarnings("ignore")

DI_MAURO = 0
CALORE = 6
FERMILAB_GNFW = 1
FERMILAB_NFW = 2
ABAZAJIAN = 7
GORDON = 8
AJELLO_2017 = 3
AJELLO_BLACK = 4
AJELLO_GREEN = 5
AJELLO_2017_DI_MAURO = 9
ARTIFICIAL = 10

EPSILON=1.0e-20

MY_ROI_SIZE = 0.428821318754
BIG_ROI_SIZE = 0.477550208734
LM_FIT_SCALE = 1.0
AJELLO_2017_ROI_SIZE = 0.09569838481
LARGE_ROI_FACTOR = 1.9451918559669763 # Compare 40x40 degree region to 40x40 with mask
GORDON_ROI_FACTOR = 0.9207648341188807 # Compare 7x7 with no mask to 40x40 with mask, with gamma=1.2
ABAZAJIAN_ROI_FACTOR = 0.5594448642944202# Compare 7x7 with no mask to 40x40 with mask, with gamma=1.1
AJELLO_ROI_FACTOR = 1.3089835019210145 # Compare 15x15 with no mask to 40x40 with mask, with gamma=1.2
AJELLO_2017_ROI_FACTOR = 0.5555527292865902




# For 2 kpc cut:
'''LARGE_ROI_FACTOR = 2.702834796393673 # Compare 40x40 degree region to 40x40 with mask
GORDON_ROI_FACTOR = 1.7361316931273842 # Compare 7x7 with no mask to 40x40 with mask, with gamma=1.2
ABAZAJIAN_ROI_FACTOR = 1.1507753162003593# Compare 7x7 with no mask to 40x40 with mask, with gamma=1.1
AJELLO_ROI_FACTOR = 2.378377491975549 # Compare 15x15 with no mask to 40x40 with mask, with gamma=1.2
AJELLO_2017_ROI_FACTOR = 0.9295394545243947
'''



Y_LIM = 1e-8# lOWERMOST Y LIMIT
ERROR_BAR_SIZE = 0.2

CALORE_L_BREAK, CALORE_ALPHA_BELOW, CALORE_ALPHA_ABOVE = 2.06, 1.42, 2.63

def calore_power_law(x, scale):
    alpha = np.full_like(x, CALORE_ALPHA_ABOVE)
    alpha[(10**x) < CALORE_L_BREAK] = CALORE_ALPHA_BELOW
    return (10**x)**2 * scale * ((10**x) / CALORE_L_BREAK)**(-alpha)

def power_law(x, scale, alpha_below, alpha_above, l_break):
    alpha = np.full_like(x, alpha_above)
    alpha[(10**x) < l_break] = alpha_below
    return (10**x)**2 * scale * (10**x / l_break)**(-alpha) * LM_FIT_SCALE

def calore_power_law_params(x, params):
    return power_law(x, params["scale"].value, CALORE_ALPHA_BELOW, CALORE_ALPHA_ABOVE, CALORE_L_BREAK)

def power_law_params(x, params):
    return power_law(x, params["scale"].value, params["alpha_below"].value, params["alpha_above"].value, params["l_break"].value)

def chi_squared_sym(params, datax, datay, weights, fn):
    return (fn(datax, params) - datay) * weights

def chi_squared_asym(params, datax, datay, weights_up, weights_down, fn):
    expected = fn(datax, params)
    weights = np.where(datay > expected, weights_down, weights_up)
    return (expected - datay) * weights

_files = {
    FERMILAB_GNFW:"fermilab-gnfw",
    CALORE:"calore",
    DI_MAURO:"di-mauro",
    AJELLO_BLACK:"ajello-black",
    AJELLO_GREEN:"ajello-green",
    AJELLO_2017:"ajello-2017",
    FERMILAB_NFW:"fermilab-nfw",
    ABAZAJIAN:"abazajian",
    GORDON:"gordon",
    AJELLO_2017_DI_MAURO:"ajello-2017-di-mauro",
    ARTIFICIAL:"artificial"
}
_quadrature_unc = {
    FERMILAB_GNFW: False,
    CALORE: False,
    DI_MAURO: True,
    AJELLO_2017_DI_MAURO: False,
    AJELLO_BLACK: False,
    AJELLO_GREEN: False,
    FERMILAB_NFW: False,
    ABAZAJIAN: False,
    GORDON: True,
}

class Spectrum:
    def __init__(self, id):
        self.id = id
        self.calore_norm = None
        self.fit_norm = None
        self.fit_alpha_below = None
        self.fit_alpha_above = None
        self.fit_l_break = None

        # Store x in log space, y in linear space.
        # Both CSVs are assumed to be stored in log space
        self.pointsx, self.pointsy = self._load_csv("spectrum-data/"+_files[id]+'/points.csv')
        if os.path.exists("spectrum-data/"+_files[id]+'/up-bands.csv')\
            and os.path.exists("spectrum-data/"+_files[id]+'/down-bands.csv'):
            up_bandsx, up_bandsy = self._load_csv("spectrum-data/"+_files[id]+'/up-bands.csv')
            down_bandsx, down_bandsy = self._load_csv("spectrum-data/"+_files[id]+'/down-bands.csv')
            if _quadrature_unc[id]:
                up_barsx, up_barsy = self._load_csv("spectrum-data/"+_files[id]+'/up-bars.csv')
                down_barsx, down_barsy = self._load_csv("spectrum-data/"+_files[id]+'/down-bars.csv')
                xmin_down = min(np.min(down_barsx), np.min(down_bandsx))
                xmax_down = max(np.max(down_barsx), np.max(down_bandsx))
                xmin_up = min(np.min(up_barsx), np.min(up_bandsx))
                xmax_up = max(np.max(up_barsx), np.max(up_bandsx))
                self.up_barsx = []
                self.up_barsy = []
                self.down_barsx = []
                self.down_barsy = []

                for i, pointx in enumerate(self.pointsx):
                    if xmin_down < pointx < xmax_down:
                        down_bar = self._get_point(pointx, down_barsx, down_barsy)
                        down_band = self._get_point(pointx, down_bandsx, down_bandsy)
                        if down_band is None:
                            down_band = 0
                        self.down_barsx.append(pointx)
                        self.down_barsy.append(-np.sqrt((down_bar - self.pointsy[i])**2 + (down_band - self.pointsy[i])**2) + self.pointsy[i])
                    if xmin_up < pointx < xmax_up:
                        up_bar = self._get_point(pointx, up_barsx, up_barsy)
                        up_band = self._get_point(pointx, up_bandsx, up_bandsy)
                        if up_band is None:
                            up_band = 0
                        self.up_barsx.append(pointx)
                        self.up_barsy.append(np.sqrt((up_bar - self.pointsy[i])**2 + (up_band - self.pointsy[i])**2) + self.pointsy[i])

                self.up_barsx = np.asarray(self.up_barsx)
                self.up_barsy = np.asarray(self.up_barsy)
                self.down_barsx = np.asarray(self.down_barsx)
                self.down_barsy = np.asarray(self.down_barsy)
            else:
                self.up_barsx = up_bandsx
                self.up_barsy = up_bandsy
                self.down_barsx = down_bandsx
                self.down_barsy = down_bandsy
        else: # There are no bands
            self.up_barsx, self.up_barsy = self._load_csv("spectrum-data/"+_files[id]+'/up-bars.csv')
            self.down_barsx, self.down_barsy = self._load_csv("spectrum-data/"+_files[id]+'/down-bars.csv')

        if id in [AJELLO_BLACK, AJELLO_GREEN, AJELLO_2017_DI_MAURO]: # Convert MeV to GeV for certain plots
            self.pointsy /= 1000
            self.up_barsy /= 1000
            self.down_barsy /= 1000

        if id in [FERMILAB_NFW, FERMILAB_GNFW]: # Convert from 2 sigma to 1 sigma
            for i, barx in enumerate(self.up_barsx):
                point = self._get_point(barx, self.pointsx, self.pointsy)
                if point is None:
                    continue
                self.up_barsy[i] = (self.up_barsy[i] - point) / 2 + point
            for i, barx in enumerate(self.down_barsx):
                point = self._get_point(barx, self.pointsx, self.pointsy)
                if point is None:
                    continue
                self.down_barsy[i] = (self.down_barsy[i] - point) / 2 + point

        if id in [ABAZAJIAN]:
            self.pointsy *= 1 / ABAZAJIAN_ROI_FACTOR
            self.up_barsy *= 1 / ABAZAJIAN_ROI_FACTOR
            self.down_barsy *= 1 / ABAZAJIAN_ROI_FACTOR

        if id in [GORDON]:
            self.pointsy *= 1 / GORDON_ROI_FACTOR
            self.up_barsy *= 1 / GORDON_ROI_FACTOR
            self.down_barsy *= 1 / GORDON_ROI_FACTOR

        if id in [AJELLO_BLACK, AJELLO_GREEN]:
            self.pointsy *= 1 / AJELLO_ROI_FACTOR
            self.up_barsy *= 1 / AJELLO_ROI_FACTOR
            self.down_barsy *= 1 / AJELLO_ROI_FACTOR

        if id in [DI_MAURO]: # Convert from square ROI to ROI without band
            self.pointsy /= LARGE_ROI_FACTOR
            self.up_barsy /= LARGE_ROI_FACTOR
            self.down_barsy /= LARGE_ROI_FACTOR

        if id in [AJELLO_2017, AJELLO_2017_DI_MAURO]:
            # Di Mauro says he doesn't convert the actual flux alue.
            self.pointsy *= 1 / AJELLO_2017_ROI_FACTOR
            self.up_barsy *= 1 / AJELLO_2017_ROI_FACTOR
            self.down_barsy *= 1 / AJELLO_2017_ROI_FACTOR


    def get_name(self):
        if self.id == ARTIFICIAL:
            return "Fake"
        elif self.id == FERMILAB_GNFW:
            return "Zhong et al. 2020, $\gamma$=1.2"
        elif self.id == CALORE:
            return "Calore et al. 2015"
        elif self.id == DI_MAURO:
            return "Di Mauro 2021"
        elif self.id == AJELLO_GREEN:
            return "Ajello et al. 2016, OB"
        elif self.id == AJELLO_BLACK:
            return "Ajello et al. 2016, PSR"
        elif self.id == FERMILAB_NFW:
            return "Zhong et al. 2020, $\gamma$=1.0"
        elif self.id == ABAZAJIAN:
            return "Abazajian et al. 2014"
        elif self.id == GORDON:
            return "Gordon et al. 2013"
        elif self.id == AJELLO_2017:
            return "Ackermann et al. 2017"
        elif self.id == AJELLO_2017_DI_MAURO:
            return "Ackermann et al. 2017 from di Mauro"
        return ""

    def _load_csv(self, filepath):
        # Return x and y values of the spectrum. Both are in logarithmic units. Y is E^2 dN/dE.
        xres = []
        yres = []
        f = open(filepath, 'r')
        for line in f.readlines():
            if line == '': continue
            x, y = line.split(",")
            xres.append(float(x))
            yres.append(10**float(y))
        f.close()

        # Return the result, scaled to the ROI in question
        if self.id in [ABAZAJIAN, AJELLO_BLACK, AJELLO_GREEN, GORDON, ARTIFICIAL, AJELLO_2017]:
            return np.asarray(xres), np.asarray(yres)
        if self.id in [DI_MAURO]:
            return np.asarray(xres), np.asarray(yres) * BIG_ROI_SIZE
        if self.id in [FERMILAB_NFW, FERMILAB_GNFW, CALORE]:
            return np.asarray(xres), np.asarray(yres) * MY_ROI_SIZE
        if self.id in [AJELLO_2017_DI_MAURO]:
            # Di mauro divides by Ajello's ROI size
            return np.asarray(xres), np.asarray(yres) * AJELLO_2017_ROI_SIZE
        raise Exception("The ROI size must be specified in _load_csv")
        return None

    def _get_point(self, logx, datax, datay):
        i = 0
        if logx < datax[0] or logx > datax[-1]:
            return None
        while datax[i+1] < logx:
            i += 1
        factor = (logx - datax[i]) / (datax[i+1] - datax[i])
        return factor * (datay[i+1] - datay[i]) + datay[i]

    def get_point_log(self, logx):
        return self._get_point(logx, self.pointsx, self.pointsy)

    def get_up_bar_log(self, logx):
        e = self._get_point(logx, self.up_barsx, self.up_barsy)
        return e - self.get_point_log(logx) if e is not None else None

    def get_down_bar_log(self, logx):
        e = self._get_point(logx, self.down_barsx, self.down_barsy)
        return self.get_point_log(logx) - e if e is not None else None

    def _fit_power_law(self, l_min, l_max, override):
        if self.fit_norm != None and not override:
            return
        log_min = np.log10(l_min)
        log_max = np.log10(l_max)

        params = lmfit.Parameters()
        params.add(lmfit.Parameter("scale", value=1e-7, min=0))
        params.add(lmfit.Parameter("alpha_below", value=CALORE_ALPHA_BELOW))
        params.add(lmfit.Parameter("alpha_above", value=CALORE_ALPHA_ABOVE))
        params.add(lmfit.Parameter("l_break", value=CALORE_L_BREAK, min=0.5, max=10))

        datax = []
        datay = []
        for i in range(len(self.pointsx)):
            if log_min < self.pointsx[i] < log_max:
                if self.pointsx[i] < np.min(self.up_barsx) and self.pointsx[i] < np.min(self.down_barsx):
                    continue
                if self.pointsx[i] > np.max(self.up_barsx) and self.pointsx[i] > np.max(self.down_barsx):
                    continue
                datax.append(self.pointsx[i])
                datay.append(self.pointsy[i] * LM_FIT_SCALE)

        datax = np.asarray(datax)
        datay = np.asarray(datay)
        aligned_up_bars = np.asarray([self.get_up_bar_log(l) for l in datax])
        aligned_down_bars = np.asarray([self.get_down_bar_log(l) for l in datax])
        aligned_up_bars[aligned_up_bars == None] = EPSILON
        aligned_down_bars[aligned_down_bars == None] = EPSILON
        weights_up = np.array([1 / b for b in aligned_up_bars])
        weights_down = np.array([1 / b for b in aligned_down_bars])

        #print(power_law_params(datax, params) - datay)
        result = lmfit.minimize(chi_squared_asym, params, args=(datax, datay, weights_up, weights_down, power_law_params))
        #lmfit.report_fit(result)

        self.fit_norm = result.params["scale"].value
        self.fit_norm_unc = result.params["scale"].stderr
        self.fit_alpha_below = result.params["alpha_below"].value
        self.fit_alpha_below_unc = result.params["alpha_below"].stderr
        self.fit_alpha_above = result.params["alpha_above"].value
        self.fit_alpha_above_unc = result.params["alpha_above"].stderr
        self.fit_l_break = result.params["l_break"].value
        self.fit_l_break_unc = result.params["l_break"].stderr

        chisq = np.sum(chi_squared_asym(result.params, datax, datay, weights_up, weights_down, power_law_params)**2)
        #print(self.get_name(), "power law fit p-value\t", 1 - stats.chi2.cdf(chisq, len(datax) - 4))
        #print(self.get_name(), "power law fit redchi\t", chisq / (len(datax) - 4))

    def _fit_calore(self, l_min, l_max, override):
        if self.calore_norm != None and not override:
            return
        log_min = np.log10(l_min)
        log_max = np.log10(l_max)

        params = lmfit.Parameters()
        params.add(lmfit.Parameter("scale", value=1e-7, min=0))

        datax = []
        datay = []
        for i in range(len(self.pointsx)):
            if log_min < self.pointsx[i] < log_max:
                if self.pointsx[i] < np.min(self.up_barsx) and self.pointsx[i] < np.min(self.down_barsx):
                    continue
                if self.pointsx[i] > np.max(self.up_barsx) and self.pointsx[i] > np.max(self.down_barsx):
                    continue
                datax.append(self.pointsx[i])
                datay.append(self.pointsy[i] * LM_FIT_SCALE)

        datax = np.asarray(datax)
        datay = np.asarray(datay)

        aligned_up_bars = np.asarray([self.get_up_bar_log(l) for l in datax])
        aligned_down_bars = np.asarray([self.get_down_bar_log(l) for l in datax])
        aligned_up_bars[aligned_up_bars == None] = EPSILON
        aligned_down_bars[aligned_down_bars == None] = EPSILON
        weights_up = np.array([1 / b for b in aligned_up_bars])
        weights_down = np.array([1 / b for b in aligned_down_bars])

        result = lmfit.minimize(chi_squared_asym, params, args=(datax, datay, weights_up, weights_down, calore_power_law_params))
        #lmfit.report_fit(result)

        self.calore_norm = result.params["scale"].value
        self.calore_norm_unc = result.params["scale"].stderr

        chisq = np.sum(chi_squared_asym(result.params, datax, datay, weights_up, weights_down, calore_power_law_params)**2)
        #print(self.get_name(), "calore fit p-value\t", 1 - stats.chi2.cdf(chisq, len(datax) - 1))
        #print(self.get_name(), "calore fit redchi\t", chisq / (len(datax) - 1))

    def get_calore_flux(self, l_min, l_max, override=False):
        self._fit_calore(l_min, l_max, override)
        a = l_min
        b = l_max
        c = self.calore_norm
        L = CALORE_L_BREAK
        a1 = CALORE_ALPHA_BELOW
        a2 = CALORE_ALPHA_ABOVE
        val = c * ((L**2 - a**(2 - a1) * L**a1) / (2 - a1) - (L**2 - b**(2 - a2) * L**a2) / (2 - a2))
        return (val, val * self.calore_norm_unc / c)

    def get_power_law_flux(self, l_min, l_max, override=False):
        self._fit_power_law(l_min, l_max, override)
        if self.fit_norm < 0:
            return None

        c = self.fit_norm
        L = self.fit_l_break
        a1 = self.fit_alpha_below
        a2 = self.fit_alpha_above
        a = l_min
        b = l_max
        val = c * ((L**2 - a**(2 - a1) * L**a1) / (2 - a1) - (L**2 - b**(2 - a2) * L**a2) / (2 - a2))

        norm_unc = val * self.fit_norm_unc / c
        break_unc = c * self.fit_l_break_unc * ((2 * L - a1 * a**(2 - a1) * L**(a1-1)) / (2 - a1) - (2 * L - b**(2 - a2) * a2 * L**(a2-1)) / (2 - a2))
        below_unc = c * self.fit_alpha_below_unc * a**(-a1) / (a1 - 2)**2 * \
             (a**(a1) * L**2 - a**2 * L**a1 - (a1 - 2)* a**2 * L**a1 * np.log(a) + (a1 -2)*a**2 * L**a1 * np.log(L))
        above_unc = c * self.fit_alpha_above_unc * b**(-a2) / (a2 - 2)**2 * \
            (-b**(a2) * L**2 + b**2 * L**a2 + (a2 - 2)* b**2 * L**a2 * np.log(b) - (a2 -2)*b**2 * L**a2 * np.log(L))

        return (val, np.sqrt(norm_unc**2 + break_unc**2 + above_unc**2 + below_unc**2))

    def get_numerical_flux(self, l_min, l_max, override=False):
        log_min = np.log10(l_min)
        log_max = np.log10(l_max)
        if log_min < self.pointsx[0]:
            # Treat all points outside the data as zero flux.
            log_min = self.pointsx[0]
        if log_max > self.pointsx[-1]:
            log_max = self.pointsx[-1]
        # Trapezoidal integration
        error_squared = 0
        integral = 0
        i = 0
        while self.pointsx[i] < log_min:
            i += 1

        integral += 0.5 * (self.pointsy[i] + self.get_point_log(log_min)) * (self.pointsx[i] - log_min)

        error_squared = ((self.get_down_bar_log(self.pointsx[i]) + self.get_up_bar_log(self.pointsx[i])) / 2 * (self.pointsx[i] - log_min))**2
        while self.pointsx[i+1] < log_max:
            integral += 0.5 * (self.pointsy[i+1] + self.pointsy[i]) * (self.pointsx[i+1] - self.pointsx[i])
            error_squared += ((self.get_down_bar_log(self.pointsx[i]) + self.get_up_bar_log(self.pointsx[i])) / 2 * (self.pointsx[i+1] - self.pointsx[i]))**2

            i += 1

        integral += 0.5 * (self.get_point_log(log_max) + self.pointsy[i]) * (log_max - self.pointsx[i])
        error_squared += ((self.get_down_bar_log(log_max) + self.get_up_bar_log(log_max)) / 2 * (log_max - self.pointsx[i]))**2

        return (integral * np.log(10), np.sqrt(error_squared) * np.log(10))

    def label_axes(self, ax):
        ax.set_xlabel(self.get_x_label())
        ax.set_ylabel(self.get_y_label())

    def display_data(self, ax, color='k', shape='o', fill_color=None, title_size=None):
        fcolor = color if None else fill_color
        thisx = []
        thisy = []
        up_err = []
        down_err = []
        lolims=[]
        uplims=[]
        for i, x in enumerate(self.pointsx):
            y = self.get_point_log(x)
            u = self.get_up_bar_log(x)
            d = self.get_down_bar_log(x)
            if u is None:
                up_err.append(0)
                lolims.append(False)
            else:
                up_err.append(u)
                lolims.append(False)

            if d is None:
                down_err.append(0)
                uplims.append(False)
            elif y - d < 2e-8:
                down_err.append(ERROR_BAR_SIZE * y)
                uplims.append(True)
            else:
                down_err.append(d)
                uplims.append(False)
            thisx.append(x)
            thisy.append(y)

        thisx = np.asarray(thisx)
        thisy = np.asarray(thisy)
        down_err = np.asarray(down_err)
        up_err = np.asarray(up_err)


        ax.errorbar(10**thisx, thisy, yerr=[down_err, up_err], color=color,
            linestyle='none', elinewidth=1, lolims=lolims, uplims=uplims)

        ax.errorbar(10**thisx, thisy, yerr=[down_err, up_err], color=color,
            linestyle='none', elinewidth=1)

        ax.scatter(10**thisx, thisy, color=color,
            label=self.get_name(), marker=shape, facecolors=fcolor)
        if title_size is not None:
            ax.set_title(self.get_name(), size=title_size)
        ax.set_xlim(0.1, 100)
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log", nonpositive='clip')

    def display_calore(self, ax, l_min, l_max, show_all=False, linestyle="--"):
        self._fit_calore(l_min, l_max, False)

        if not show_all:
            calorex = np.linspace(np.log10(l_min),
                np.log10(l_max), 100)
        else:
            calorex = np.linspace(np.min(self.pointsx),
                np.max(self.pointsx), 100)
        calorey = calore_power_law(calorex, self.calore_norm)

        ax.plot(10**calorex, calorey, color='red', linestyle=linestyle)

    def display_power_law(self, ax, l_min, l_max, show_all=False, linestyle="-."):
        self._fit_power_law(l_min, l_max, False)

        if not show_all:
            fitx = np.linspace(np.log10(l_min),
                np.log10(l_max), 100)
        else:
            fitx = np.linspace(np.min(self.pointsx),
                np.max(self.pointsx), 100)
        fity = power_law(fitx, self.fit_norm, self.fit_alpha_below,
            self.fit_alpha_above, self.fit_l_break) / LM_FIT_SCALE

        ax.plot(10**fitx, fity, color='blue', linestyle=linestyle)

    def get_x_label(self):
        return "$E_\\gamma$ [GeV]"

    def get_y_label(self):
        return "$F_\\gamma$ [GeV / cm$^2$ / s]"
