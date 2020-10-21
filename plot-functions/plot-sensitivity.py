import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

plt.style.use("jcap")

x = 10**np.linspace(-13, -10, 100)
STEP_THRESHOLD = 3e34 * 1.1093417307914119e-46
KTH = 0.45
SIGMATH = 0.28

def step(x):
    return x > STEP_THRESHOLD

def smoothed(x):
    return 0.5 * (1 + erf((np.log10(x) - (np.log10(STEP_THRESHOLD) + KTH)) / (np.sqrt(2) * SIGMATH)))


plt.xscale('log')
plt.xlabel("Flux [erg/s]")
plt.ylabel("Detection probability")
plt.plot(x, step(x), label="Step sensitivity function")
plt.plot(x, smoothed(x), label="Smoothed sensitivity function")
plt.legend()
plt.savefig("sensitivity-functions.png")
plt.show()
