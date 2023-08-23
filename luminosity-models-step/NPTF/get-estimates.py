'''
The purpose of this file is to take the NPTF luminosity function
and generate predictions of how many pulsars Fermi should have seen,
and how much flux above the threshold it should have seen. We expect
the luminosity to have a short and strong peak below the threshold,
leading to very few visible pulsars and small visible luminosity.
The paper is found here https://arxiv.org/pdf/1506.05124.pdf,
with a luminosity function described in table 1.
'''

from matplotlib import pyplot as plt
import numpy as np
from math import log10

plt.style.use('jcap')

L_EXCESS = 1.794113925439598e-09 / 1.1095246594108431e-46# 6.756e36  # All units are in erg per second
L_THRESH = 1e34

L_MIN=1e29
L_MAX=1e35

CUTOFF_MINIMUM = 1e29

global outStr
outStr = ""

def pow(a, b):
    return a**b

class LuminosityFunction:
    def __init__(self, name, nBelow, nAbove, lBreak, minL, maxL):
        self.name = name
        self.nBelow = nBelow
        self.nAbove = nAbove
        self.lBreak = lBreak
        self.minL = minL
        self.maxL = maxL

    def integrate(self, minL, maxL=None):
        if minL is None: minL = CUTOFF_MINIMUM
        if maxL is not None:
            if minL > maxL: return 0

            if minL < self.lBreak:
                return (-self.nAbove * self.lBreak + pow(maxL, 1 - self.nAbove) * pow(self.lBreak, self.nAbove) -
                    minL * pow(self.lBreak/minL, self.nBelow) + self.nAbove * minL * pow(self.lBreak/minL, self.nBelow) +
                    self.nBelow * (self.lBreak - pow(maxL, 1 - self.nAbove) * pow(self.lBreak, self.nAbove)))/((-1 +
                    self.nBelow) * (-1 + self.nAbove))
            else:
                return ((pow(maxL / (maxL * minL), self.nAbove) * minL -
                   maxL * pow(minL / (maxL * minL), self.nAbove)) * pow(self.lBreak, self.nAbove))/(-1 + self.nAbove)

        nptfPremul = (self.nBelow - self.nAbove) / (self.nBelow - self.nAbove - (1 - self.nAbove) * pow(CUTOFF_MINIMUM / self.lBreak, 1 - self.nBelow))
        if (minL < self.lBreak):
            return nptfPremul * (1 - (self.lBreak / minL) ** (self.nBelow - 1) * (self.nAbove - 1) / (self.nAbove - self.nBelow))
        else:
            return nptfPremul * (self.lBreak / minL) ** (self.nAbove - 1) * (1 - self.nBelow) / (self.nAbove - self.nBelow)

    def lintegrate(self, minL, maxL=None):
        if minL is None: minL = CUTOFF_MINIMUM
        if maxL is not None:
            if minL > maxL: return 0

            if minL < self.lBreak:
                return (-self.nAbove * pow(self.lBreak, 2) - 2 * pow(minL, 2 - self.nBelow) * pow(self.lBreak, self.nBelow) +
                    self.nAbove * pow(minL, 2 - self.nBelow) * pow(self.lBreak, self.nBelow) +
                    2 * pow(maxL, 2 - self.nAbove) * pow(self.lBreak, self.nAbove) +
                    self.nBelow * (pow(self.lBreak, 2) - pow(maxL, 2 - self.nAbove) * pow(self.lBreak, self.nAbove)))/((-2 +
                    self.nBelow) * (-2 + self.nAbove))
            else:
                return ((pow(maxL, 2 - self.nAbove) - pow(minL, 2 - self.nAbove)) * pow(self.lBreak, self.nAbove))/(2 - self.nAbove)


        nptfPremul = (self.nBelow - self.nAbove) / (self.nBelow - self.nAbove - (1 - self.nAbove) * pow(CUTOFF_MINIMUM / self.lBreak, 1 - self.nBelow))
        if minL < self.lBreak:
            return nptfPremul * self.lBreak * (1 - self.nBelow) * (1-self.nAbove) * (1 / ((self.nBelow-2) * (self.nAbove-2)) + pow(self.lBreak / minL, self.nBelow - 2) / ((self.nBelow - 2) * (self.nBelow - self.nAbove)))
        else:
            return nptfPremul * self.lBreak * (1 - self.nBelow) * (1 - self.nAbove) * (pow(self.lBreak / minL, self.nAbove - 2) / ((self.nAbove - 2) * (self.nBelow - self.nAbove)))

    def getValue(self, l):
        scale = 1
        if self.nBelow > 0:
            lMin = 1e29 if (L_MIN==None or L_MIN < 1e29) else L_MIN
            scale = (lMin / self.lBreak)**(self.nBelow) # Make the highest point in the luminosity function have value one. Affects only the plot.
        if l < self.lBreak:
            return scale * (l / self.lBreak)**(-self.nBelow)
        else:
            return scale * (l / self.lBreak)**(-self.nAbove)

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

lumFuncs = [LuminosityFunction("NFW PS", -0.66, 18.2, 2.5389429e+34, 1e29, None),
            LuminosityFunction("Bartels 18", 0.97, 2.60, 10**33.24, 1e30, 1e37)]


for l in lumFuncs:
    l.printEstimate()
    l.display()

f = open("get-estimates-output.txt", 'w')
f.write(outStr)
f.close()

plt.axvline(x=L_THRESH, label="Threshold", color='black')
plt.xscale('log')
#plt.yscale('log')
plt.title("NPTF luminosity function")
plt.savefig("get-estimates.png")
plt.xlabel("Luminosity (erg/s)")
plt.ylabel("$\\frac{dN}{dL}$ (unnormalized)")
plt.legend()

plt.savefig("luminosity-func.png")
plt.show()
