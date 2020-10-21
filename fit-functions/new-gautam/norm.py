import numpy as np

xs = []
ys = []

f = open("points.csv", 'r')
for line in f.read().split('\n'):
    if line == '':
        continue
    x, y = line.split(', ')
    xs.append(10**float(x))
    ys.append(10**float(y))
f.close()

norm = 0
for i, y in enumerate(ys[:-1]):
    width = xs[i] + xs[i+1]
    norm += y * width

print(norm)
# 1.8449245050788346e-09
print(sum(ys))
# 122491.64085316139
