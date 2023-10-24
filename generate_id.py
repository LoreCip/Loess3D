import sys

import numpy as np

import h5py as h5

Nth = sys.argv[1]

with h5.File('/home/lorenzo/phd/DePietri/NS_HOT_EOS/EOS/compOSE/FOP(SFHoY)/FOP(SFHoY).h5', 'r') as f:

    n, m, l = int(f['pointsye'][()]), int(f['pointstemp'][()]), int(f['pointsrho'][()])

    points = np.arange(m)
    p1 = points[points<3*m/4]
    p2 = points[points>=3*m/4]
    points=np.concatenate((p1[::2], p2))
    m = len(points)

    logtmp = f['logtemp'][()]
    lognb  = np.log10(f['nb'][()][points])
    ye     = f['ye'][()]
    entropy= f['entropy'][()][:, :, points]

X, Y, Z = np.meshgrid(ye, logtmp, lognb)
f_ran = entropy.copy()

with open('data.init', 'w') as f:

    f.write(f"{Nth}\n")
    f.write(f"{n}\t{m}\t{l}\n")
    for arr in [X, Y, Z, f_ran]:
        arr = arr.flatten()
        for i in range(len(arr)):
            f.write(str(arr[i]) + "\t")

        f.write("\n")
