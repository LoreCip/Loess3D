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

Nth = int(sys.argv[1])
frac = float(sys.argv[2])
degree = int(sys.argv[3])

with h5.File('./EOS/FOP(SFHoY).h5', 'r') as f:

    n, m, l = int(f['pointsye'][()]), int(f['pointstemp'][()]), int(f['pointsrho'][()])

    points_n, n = reduce_points(n, 6)
    points_m, m = reduce_points(m, 6)
    points_l, l = reduce_points(l, 20)

    logtmp = f['logtemp'][()][points_m]
    lognb  = np.log10(f['nb'][()][points_l])
    ye     = f['ye'][()][points_n]
    entropy= np.log10(f['entropy'][()])
    entropy = entropy[np.ix_(points_n, points_m, points_l)]

f_ran = entropy.copy()

with h5.File('./logentropy/completeData.h5', 'w') as f:

    # Scalars
    f['n'] = n
    f['m'] = m
    f['l'] = l
    f['Nth'] = Nth
    f['degree'] = degree
    f['frac']   = frac

    # Arrays 1D
    f.create_dataset('Yq', data=ye, compression='gzip', compression_opts=9)
    f.create_dataset('logtemp', data=logtmp, compression='gzip', compression_opts=9)
    f.create_dataset('lognb', data=lognb, compression='gzip', compression_opts=9)
    f.create_dataset('LogEntropy', data=f_ran, compression='gzip')    
