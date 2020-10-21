import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameter, report_fit
import scipy.optimize as optimization
from scipy import stats
import lmfit

plt.style.use("jcap")

LUM_TO_FLUX = 1.1093417307914119e-46
NUM_POINTS = 100

def get_data(filename):
    xs = []
    ys = []
    f = open(filename, 'r')
    for line in f.readlines():
        if line == '':
            continue
        x, y = line.split(', ')
        x = float(x)
        y = float(y)

        xs.append(x)
        ys.append(y)
    return np.asarray(xs), np.asarray(ys)

def at(x, xs, ys):
    if x < np.min(xs):
        return 0
    if x > np.max(xs):
        return 0
    i = 0
    while xs[i+1] < x:
        i += 1

    return (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i]) + ys[i]

disk_pointsx, disk_pointsy = get_data("disk/points.csv")
disk_upbarsx, disk_upbarsy = get_data("disk/bars-top.csv")
disk_downbarsx, disk_downbarsy = get_data("disk/bars-bottom.csv")

bb_pointsx, bb_pointsy = get_data("bb/points.csv")
bb_upbarsx, bb_upbarsy = get_data("bb/bars-top.csv")
bb_downbarsx, bb_downbarsy = get_data("bb/bars-bottom.csv")

real_xs = []
real_ys = []
real_upbars = []
real_downbars = []
for x in disk_pointsx:
    real_xs.append(x)
    disky = at(x, disk_pointsx, disk_pointsy)
    bby = at(x, bb_pointsx, bb_pointsy)
    real_ys.append(disky + bby)
    disk_upbar = abs(at(x, disk_upbarsx, disk_upbarsy) - disky)
    disk_downbar = abs(at(x, disk_upbarsx, disk_upbarsy) - disky)
    bb_upbar = abs(at(x, bb_upbarsx, bb_upbarsy) - bby)
    bb_downbar = abs(at(x, bb_upbarsx, bb_upbarsy) - bby)
    real_upbars.append(np.sqrt(disk_upbar**2 + bb_upbar**2) + disky + bby)
    real_downbars.append(-np.sqrt(disk_downbar**2 + bb_downbar**2) + disky + bby)

points_file = open("points.csv", 'w')
upbars_file = open("bars-top.csv", 'w')
downbars_file = open("bars-bottom.csv", 'w')
for (i, x) in enumerate(real_xs):
    points_file.write("{}, {}\n".format(x, real_ys[i]))
    upbars_file.write("{}, {}\n".format(x, real_upbars[i]))
    downbars_file.write("{}, {}\n".format(x, real_downbars[i]))

points_file.close()
upbars_file.close()
downbars_file.close()
