from matplotlib import pyplot as plt
import numpy as np
from math import log10

plt.style.use('jcap')

L_EXCESS = 1.794113925439598e-09 / 1.1095246594108431e-46# 6.756e36  # All units are in erg per second
L_THRESH = 3.5630913584093865e+34

L_MIN=1e29
L_MAX=1e35

global outStr
outStr = ""

class LuminosityFunction:
    def __init__(self, name, minL, maxL, alpha):
        self.name = name
        self.minL = minL
        self.maxL = maxL
        self.alpha = alpha

    def integrate(self, minL, maxL=None):
        if maxL is None: maxL = self.maxL
        return (minL**(1-self.alpha) - maxL**(1-self.alpha)) / (self.minL**(1-self.alpha) - self.maxL**(1-self.alpha))

    def lintegrate(self, minL, maxL=None):
        if maxL is None: maxL = self.maxL
        return ((minL**(2-self.alpha) - maxL**(2-self.alpha)) / (self.alpha - 2)) / ((self.minL**(1-self.alpha) - self.maxL**(1-self.alpha)) / (self.alpha - 1))

    def getValue(self, l):
        scale = 1 / ((self.minL**(1-self.alpha) - self.maxL**(1-self.alpha)) / (self.alpha - 1))
        if l < self.minL or l > self.maxL:
            return 0
        else:
            return scale * l**-self.alpha

    def printEstimate(self):
        global outStr
        unscaledNumber = self.integrate(minL=self.minL, maxL=self.maxL)
        unscaledLum = self.lintegrate(minL=self.minL, maxL=self.maxL)
        unscaledNumberAbove = self.integrate(minL=max(self.minL, L_THRESH), maxL=self.maxL)
        unscaledFluxAbove = self.lintegrate(minL=max(self.minL, L_THRESH), maxL=self.maxL)

        scale = L_EXCESS / unscaledLum

        totalNumber = unscaledNumber * scale
        totalLum = unscaledLum * scale # Should be L_EXCESS
        numberAbove = unscaledNumberAbove * scale
        R = unscaledFluxAbove / unscaledLum

        addStr = """{0} luminosity function:
    Total number of pulsars:\t\t{1}
    Total luminosity:\t\t\t{2}
    Number of pulsars above threshold:\t{3}\t(Fermi: 47)
    Fraction of lum above threshold:\t{4}\t(Fermi: 0.14)
""".format(self.name, totalNumber, totalLum, numberAbove, R)
        print(addStr)
        outStr+=addStr

    def display(self):
        x=10**np.linspace(start=(29 if (L_MIN==None or L_MIN < 1e29) else log10(L_MIN)),
            stop=(35 if (L_MAX == None or L_MAX > 1e35) else log10(L_MAX)), num=100)
        y=[self.getValue(l) for l in x]
        plt.plot(x, y, label=self.name)

lumFuncs = [LuminosityFunction("Real Wavelet Fit", 1e29, 7e34, 1.5)]

for l in lumFuncs:
    l.printEstimate()
    l.display()

f = open("get-estimates-output.txt", 'w')
f.write(outStr)
f.close()

plt.axvline(x=L_THRESH, label="Threshold", color='black')
plt.xscale('log')
#plt.yscale('log')
plt.title("wavelet fit luminosity function")
plt.savefig("get-estimates.png")
plt.xlabel("Luminosity (erg/s)")
plt.ylabel("$\\frac{dN}{dL}$ (unnormalized)")
plt.legend()
plt.show()
