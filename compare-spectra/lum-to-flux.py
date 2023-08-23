'''Convert a GCE luminosity to a flux, assuming an NFW-squared distributed
pattern of pulsars with delta function luminosity. See Jan summary for
calculation.'''

from math import pi, cos, sqrt

DELTA_ANGLE = 0.001
DELTA_S = 0.001
R_C = 8.5
R_S = 20
GAMMA = 1.2
KPC_SQUARED_PER_CM_SQUARED = 3.24078e-22**2
R_CUT = 2

def NFWSquared(r):
    return ((r/R_S)**(-GAMMA) * (1 + r/R_S)**(-3+GAMMA))**2

def calculateIntegral(function):
    integral = 0
    b = -20 * pi / 180.0
    while b < 20 * pi / 180.0:
        print(b)
        volumeElement = cos(b) * DELTA_ANGLE * DELTA_ANGLE * DELTA_S
        l = -20 * pi / 180.0
        if abs(b) < 2 * pi / 180.0:
            b += DELTA_ANGLE
            continue
        while l < 20 * pi / 180.0:
            s=0
            while s < 2*R_S:
                r = sqrt(R_C*R_C + s * s - 2 * R_C * s * cos(b) * cos(l))
                if R_CUT is not None and r > R_CUT:
                    s += DELTA_S
                    continue
                integral += function(r, s) * volumeElement
                s += DELTA_S
            l += DELTA_ANGLE
        b += DELTA_ANGLE
    return integral

def numerator(r, s):
    return NFWSquared(r)

def denominator(r, s):
    return NFWSquared(r) * s * s

print("Naive flux to luminosity ratio:", 1 / (4 * pi * R_C**2)  * KPC_SQUARED_PER_CM_SQUARED)

num = calculateIntegral(numerator)
denom = calculateIntegral(denominator)
print("Actual flux to luminosity ratio:", 1 / (4 * pi) * num / denom * KPC_SQUARED_PER_CM_SQUARED)
# 1.1095246594108431e-46 with (0.001, 0.001) Jul 28
