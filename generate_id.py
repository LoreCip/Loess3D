import numpy as np

import h5py as h5

with h5.File('/home/lorenzo/phd/DePietri/NS_HOT_EOS/EOS/compOSE/FOP(SFHoY)/FOP(SFHoY).h5', 'r') as f:

    logtmp = f['logtemp'][()]
    lognb  = np.log10(f['nb'][()])
    ye     = f['ye'][()]
    entropy= f['entropy'][()]

    n, m, l = int(f['pointsye'][()]), int(f['pointsrho'][()]), int(f['pointstemp'][()])

X, Y, Z = np.meshgrid(ye, lognb, logtmp)
f_ran = entropy.copy()

with open('data.dat', 'w') as f:

    f.write(f"{n}\t{m}\t{l}\n")
    for arr in [X, Y, Z, f_ran]:
        arr = arr.flatten()
        for i in range(len(arr)):
            f.write(str(arr[i]) + "\t")

        f.write("\n")
