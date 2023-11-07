import sys

import h5py as h5

import numpy as np
from numba import njit

@njit
def firDer(a,b,dx):
    return (a-b)/(2*dx)

@njit
def compBetaV(lognb, entropy, pointsyq, pointstemp, pointsrho):
    out = np.zeros((pointsyq, pointstemp, pointsrho))
    for k in range(pointsyq):
        for i in range(pointstemp):
            out[k, i, 0] = -10**lognb[0] * (entropy[k,i,1]-entropy[k,i,0]) / (lognb[1]-lognb[0])
            for j in range(1, pointsrho-1):
                Dy = lognb[j+1] - lognb[j]
                out[k-1,i-1,j-1] = - 10**lognb[j-1] * firDer(entropy[k,i,j+1], entropy[k,i,j-1], Dy)
            out[k, i, -1] = -10**lognb[-1] * (entropy[k,i,-1]-entropy[k,i,-2]) / (lognb[-1]-lognb[-2])
            
    return out

def compKappaT(lognb, pressure, pointsyq, pointstemp, pointsrho):
    out = np.zeros((pointsyq, pointstemp, pointsrho))
    for k in range(pointsyq):
        for i in range(pointstemp):
            out[k, i, 0] = (pressure[k,i,1]-pressure[k,i,0]) / (lognb[1]-lognb[0])
            for j in range(1, pointsrho-1):
                Dx = lognb[j] - lognb[j-1]
                out[k,i,j-1] = firDer(pressure[k,i,j+1], pressure[k,i,j-1], Dx)
            out[k, i, -1] = (pressure[k,i,-1]-pressure[k,i,-2]) / (lognb[-1]-lognb[-2])

    return np.where(kappa_T != 0, 1 / kappa_T, 0)

def compCV(logtemp, entropy, pointsyq, pointstemp, pointsrho):
    out = np.zeros((pointsyq, pointstemp, pointsrho))
    for k in range(pointsyq):
        for j in range(pointsrho):
            out[k,0,j] = firDer(entropy[k,1,j], entropy[k,0,j], Dx)
            for i in range(1, pointstemp-1):    
                Dx = logtemp[i] - logtemp[i-1]
                out[k-1,i-1,j-1] = firDer(entropy[k,i+1,j], entropy[k,i-1,j], Dx)
            out[k,-1,j] = firDer(entropy[k,-1,j], entropy[k,-2,j], Dx)
    
    return out

fname = sys.argv[1]
qty = sys.argv[2]

with h5.Open(fname, 'r') as f:
    logtmp = f['logtemp'][()]
    lognb  = f['lognb'][()]
    ye     = f['ye'][()]

    LogEntropy = f['LogEntropy'][()]
    LogPressure= f['LogPressure'][()]
    

if qty == 'betaV':
    out = compBetaV(lognb, 10**LogEntropy, len(ye), len(logtmp), len(lognb))
if qty == 'kappaT':
    out = compKappaT(lognb, 10**LogPressure, len(ye), len(logtmp), len(lognb))
if qty == 'cV':
    out == compCV(logtemp, 10**LogEntropy, len(ye), len(logtmp), len(lognb))

    
with h5.Open(fname, 'w') as f:
    f.create_dataset(qty, data=out.T, compression='gzip', compression_opts=9)    
