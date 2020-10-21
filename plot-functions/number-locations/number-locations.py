import numpy as np
import matplotlib.pyplot as plt

FLUX_THRESHOLD = 4e-12
R_S = 20
R_C = 8.5
GAMMA = 1.2
L_BREAK =2.4127028e+34
N_ABOVE =18.2
N_BELOW =-0.66
CM_PER_KPC_SQUARED = 9.523396e+42
LAT = 0 # radians
LON = 0.2 # radians
LOG_MAX_L = 38

def gnfw_squared(r):
    return ((r / R_S)**-GAMMA * (1 + r / R_S)**(-3 + GAMMA))**2

def lum_func(l):
    if l > L_BREAK:
        return (l / L_BREAK) ** -N_ABOVE
    else:
        return (l / L_BREAK) ** -N_BELOW

def integrand(s, mult):
    r = np.sqrt(s**2 + R_C**2 - 2 * s * R_C * np.cos(LON) * np.cos(LAT))
    integral = 0
    min_lum = 4 * np.pi * s**2 * FLUX_THRESHOLD * CM_PER_KPC_SQUARED * mult
    lums = np.linspace(np.log10(min_lum), LOG_MAX_L, 10001)
    widths = [10**lums[i+1] - 10**lums[i] for i in range(len(lums) - 1)]
    for i, log_lum in enumerate(lums[:-1]):
        integral += lum_func(10**log_lum) * widths[i]
    return s**2 * gnfw_squared(r) * integral

closest_approach = R_C * np.sqrt(1 - (np.cos(LON) * np.cos(LAT))**2)
print(closest_approach)
s_cut_halfwidth = np.sqrt(max(0, 2**2 - closest_approach**2))

ss = np.linspace(2, 11, 100)
plt.plot(ss, [integrand(s, 1.4) for s in ss], color=[1, 0, 0], label="1.4 * Threshold")
plt.plot(ss, [integrand(s, 1.2) for s in ss], color=[0.8, 0, 0])
plt.plot(ss, [integrand(s, 1.0) for s in ss], color=[0.6, 0, 0])
plt.plot(ss, [integrand(s, 0.8) for s in ss], color=[0.4, 0, 0])
plt.plot(ss, [integrand(s, 0.6) for s in ss], color=[0.2, 0, 0], label="0.6 * Threshold")
plt.xlabel("s")
plt.axvline(x=R_C, linestyle='solid', color='k')
plt.axvline(x=R_C-s_cut_halfwidth, linestyle='dotted', color='k')
plt.axvline(x=R_C+s_cut_halfwidth, linestyle='dotted', color='k')
plt.ylabel("Number seen (arbitrary units)")
frame1 = plt.gca()
frame1.axes.yaxis.set_ticklabels([])
plt.legend()
plt.savefig("nptf.png")
plt.show()
