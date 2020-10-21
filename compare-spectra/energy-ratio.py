import numpy as np

# Jul 2, 2023: I don't know what this code does.

ALPHA_BELOW = 1.1107351768821663
ALPHA_ABOVE = 2.575199391296956
L_BREAK = 1.3491856809517384

END_LOW = 0.1
END_HIGH  = 100

def power_law(x, alpha_below, alpha_above, l_break):
    alpha = np.full_like(x, alpha_above)
    alpha[x < l_break] = alpha_below
    return (x / l_break)**(-alpha)

start_low = 1.893#float(input ("What is the low end of your origin spectrum? (GeV)"))
start_high = 11.943#float(input ("What is the high end of your origin spectrum? (GeV)"))

def lintegrate(low, high):
    return -(L_BREAK**2 - L_BREAK ** ALPHA_BELOW * low**(2-ALPHA_BELOW)) / (ALPHA_BELOW - 2)\
        + (L_BREAK**2 - L_BREAK ** ALPHA_ABOVE * high**(2-ALPHA_ABOVE)) / (ALPHA_ABOVE - 2)

result = lintegrate(start_low, start_high) / lintegrate(END_LOW, END_HIGH)
print(result / 1.11e-46 * 1.76e-10)