#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import random
from scipy import stats
from multiprocessing import Pool

SENSITIVITIES = [1, 2, 5, 10, 20]
FILE_NAMES = ["power-law", "bartels15", "log-normal", "ploeg", "gautam", "nptf", "bartels18"]
DISPLAY_NAMES = ["Wavelet 1", "Wavelet 2", "GLC", "GCE", "AIC", "NPTF", "Disk"]
SENS_TYPE = "position"
if SENS_TYPE == "smoothing":
    N_REALS = [
        [10, 31, 49, 5, 3, 25, 10],# 1
        [27, 81, 101, 21, 12, 110, 29],# 2
        [87, 221, 210, 104, 57, 399, 121],# 5
        [194, 410, 317, 307, 169, 659, 342],# 10
        [401, 697, 425, 787, 453, 831, 881],# 20
    ]
else:
    N_REALS = [
        [31, 98, 124, 20, 12, 111, 30], # 1
        [77, 210, 220, 73, 41, 460, 89], # 2
        [230, 490, 380, 340, 180, 810, 370], # 5
        [470, 830, 490, 930, 520, 907, 1032], # 10
        [955, 1300, 570, 2200, 1400, 940, 2500], # 20
    ]


# In[71]:


loadeds = {}

#The text in the file are pdfs calculated at the exact flux values written in the files. The cds are probably inaccurate.

def load_file(filename):
    global loadeds
    if filename in loadeds:
        return loadeds[filename]
    fluxes = []
    pdfs = []
    cdfs = []
    cdf = 0
    f = open(filename)
    for i, line in enumerate(f.readlines()):
        if line == "": continue
        flux, pdf, _ = line.split()
        fluxes.append(float(flux))
        pdfs.append(float(pdf))

    integral = 0
    for i in range(len(fluxes[:-1])):
        cdfs.append(integral)
        integral += pdfs[i] * (fluxes[i+1] - fluxes[i])
    cdfs.append(integral)
    pdfs = np.array(pdfs) / integral
    cdfs = 1 - np.array(cdfs) / integral
    
    loadeds[filename] = (np.asarray(fluxes), pdfs, cdfs)
    return np.asarray(fluxes), pdfs, cdfs

def draw_fluxes(ndraw, fluxes, cdfs):
    draws = np.zeros(ndraw)
    for i in range(ndraw):
        r = random.random()
        draws[i] = np.interp(r, cdfs, fluxes)
    return draws
    
    
def get_log_likelihood(sensitivity, data, hyp_func):
    hyp_fluxes, hyp_pdfs, _ = load_file("gen-cdfs/data-{}x/{}-{}.dat".format(sensitivity, FILE_NAMES[hyp_func], SENS_TYPE))
    log_like = 0
    for flux in data:
        log_like += np.log(np.interp(flux, hyp_fluxes, hyp_pdfs))
    return log_like

def get_ratio_matrix(sensitivity, ndraws, ntrials):
    mat = []
    for i in range(len(FILE_NAMES)): # True
        line = []
        
        real_fluxes, _, real_cdfs = load_file("gen-cdfs/data-{}x/{}-{}.dat".format(sensitivity,
                                               FILE_NAMES[i], SENS_TYPE))
        drawn_fluxes = draw_fluxes(ndraws[i] * ntrials, np.flip(real_fluxes), np.flip(real_cdfs))
        
        for j in range(len(FILE_NAMES)): # Hypothetical
            hyp_ll = get_log_likelihood(sensitivity, drawn_fluxes, j)
            true_ll = get_log_likelihood(sensitivity, drawn_fluxes, i)
            line.append(2 * (true_ll - hyp_ll) / ntrials)
        mat.append(np.asarray(line))
    return np.asarray(mat)


# In[69]:


N_TRIALS = 100_000

def call_matrix_47(sens):
    return get_ratio_matrix(sens, [47]*len(FILE_NAMES), N_TRIALS) 
def call_matrix_scale(sens):
    return get_ratio_matrix(sens[1], N_REALS[sens[0]], N_TRIALS) 

#with Pool() as pool:
    #mats_47 = pool.map(call_matrix_47, SENSITIVITIES)
#    mats_scale = pool.map(call_matrix_scale, enumerate(SENSITIVITIES))
mats_scale = [call_matrix_scale(e)for e in enumerate(SENSITIVITIES)]
    
#print(np.array(mats_47).shape)
print(np.array(mats_scale).shape)
    
#np.savetxt("mats-47-{}.txt".format(SENS_TYPE), np.array(mats_47).flatten())
np.savetxt("mats-scale-{}.txt".format(SENS_TYPE), np.array(mats_scale).flatten())