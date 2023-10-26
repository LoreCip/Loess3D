import sys

import numpy as np

import h5py as h5

def reduce_points(l, s):
    points = np.arange(l)
    p1 = points[points<3*l/4]
    p2 = points[points>=3*l/4]
    points=np.concatenate((p1[::s], p2))
    l = len(points)

    return points, l


Nth = sys.argv[1]
frac = sys.argv[2]

with h5.File('/home/lorenzo/phd/DePietri/NS_HOT_EOS/EOS/compOSE/FOP(SFHoY)/FOP(SFHoY).h5', 'r') as f:

    n, m, l = int(f['pointsye'][()]), int(f['pointstemp'][()]), int(f['pointsrho'][()])

    points_n, n = reduce_points(n, 1)
    points_m, m = reduce_points(m, 1)
    points_l, l = reduce_points(l, 2)

    logtmp = f['logtemp'][()][points_m]
    lognb  = np.log10(f['nb'][()][points_l])
    ye     = f['ye'][()][points_n]
    entropy= np.log10(f['entropy'][()])
    entropy = entropy[np.ix_(points_n, points_m, points_l)]

X = np.zeros((n,m,l))
Y = np.zeros((n,m,l))
Z = np.zeros((n,m,l))
for i in range(n):
    X[i,:,:] = ye[i]
for j in range(m):
    Y[:,j,:] = logtmp[j]
for k in range(l):
    Z[:,:,k] = lognb[k]

f_ran = entropy.copy()

with open('data.init', 'w') as f:

    f.write(f"{Nth}\n")
    f.write(f"{frac}\n")
    f.write(f"{n}\t{m}\t{l}\n")
    for arr in [X, Y, Z, f_ran]:    
        arr = arr.flatten()
        for i in range(len(arr)):
            f.write(str(arr[i]) + "\t")

        f.write("\n")
