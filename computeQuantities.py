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
                out[k,i,j] = - 10**lognb[j-1] * firDer(entropy[k,i,j+1], entropy[k,i,j-1], Dy)
            out[k, i, -1] = -10**lognb[-1] * (entropy[k,i,-1]-entropy[k,i,-2]) / (lognb[-1]-lognb[-2])
            
    return out

@njit
def compKappaT(lognb, pressure, pointsyq, pointstemp, pointsrho):
    out = np.zeros((pointsyq, pointstemp, pointsrho))
    for k in range(pointsyq):
        for i in range(pointstemp):
            out[k, i, 0] = (pressure[k,i,1]-pressure[k,i,0]) / (lognb[1]-lognb[0])
            for j in range(1, pointsrho-1):
                Dx = lognb[j] - lognb[j-1]
                out[k,i,j] = firDer(pressure[k,i,j+1], pressure[k,i,j-1], Dx)
            out[k, i, -1] = (pressure[k,i,-1]-pressure[k,i,-2]) / (lognb[-1]-lognb[-2])

    return np.where(out != 0, 1 / out, 0)

@njit
def compCV(logtemp, entropy, pointsyq, pointstemp, pointsrho):
    out = np.zeros((pointsyq, pointstemp, pointsrho))
    for k in range(pointsyq):
        for j in range(pointsrho):
            out[k,0,j] = firDer(entropy[k,1,j], entropy[k,0,j], logtemp[1]-logtemp[0])
            for i in range(1, pointstemp-1):    
                Dx = logtemp[i] - logtemp[i-1]
                out[k,i,j] = firDer(entropy[k,i+1,j], entropy[k,i-1,j], Dx)
            out[k,-1,j] = firDer(entropy[k,-1,j], entropy[k,-2,j], logtemp[-1]-logtemp[-2])
    
    return out

@njit
def compdPdE_nb(logtemp, ye, pressure, energy, pointsyq, pointstemp, pointsrho):
    out = np.zeros((pointsyq, pointstemp, pointsrho))

    for j in range(pointsrho):
        for k in range(pointsyq):
            for i in range(pointstemp):

                if (i == 0) and (k < pointsyq-1):
                    Dt = logtemp[i+1] - logtemp[i]
                    Dy = ye[k+1] - ye[k]
                    p1 = firDer(pressure[k,i+1,j], pressure[k,i,j], Dt/2)
                    p2 = firDer(energy[k,i+1,j], energy[k,i,j], Dt/2)
                    p3 = firDer(pressure[k+1,i,j], pressure[k-1,i,j], Dy)
                    p4 = firDer(energy[k+1,i,j], energy[k-1,i,j], Dy)
                
                elif (i < pointstemp-1) and (k == 0):
                    Dt = logtemp[i+1] - logtemp[i]
                    Dy = ye[k+1] - ye[k]
                    p1 = firDer(pressure[k,i+1,j], pressure[k,i-1,j], Dt)
                    p2 = firDer(energy[k,i+1,j], energy[k,i-1,j], Dt)
                    p3 = firDer(pressure[k+1,i,j], pressure[k,i,j], Dy/2)
                    p4 = firDer(energy[k+1,i,j], energy[k,i,j], Dy/2)

                elif (i == 0) and (k == 0):
                    Dt = logtemp[i+1] - logtemp[i]
                    Dy = ye[k+1] - ye[k]
                    p1 = firDer(pressure[k,i+1,j], pressure[k,i,j], Dt/2)
                    p2 = firDer(energy[k,i+1,j], energy[k,i,j], Dt/2)
                    p3 = firDer(pressure[k+1,i,j], pressure[k,i,j], Dy/2)
                    p4 = firDer(energy[k+1,i,j], energy[k,i,j], Dy/2)

                if (i == pointstemp-1) and (k < pointsyq-1):
                    Dt = logtemp[i] - logtemp[i-1]
                    Dy = ye[k+1] - ye[k]
                    p1 = firDer(pressure[k,i,j], pressure[k,i-1,j], Dt/2)
                    p2 = firDer(energy[k,i,j], energy[k,i-1,j], Dt/2)
                    p3 = firDer(pressure[k+1,i,j], pressure[k-1,i,j], Dy)
                    p4 = firDer(energy[k+1,i,j], energy[k-1,i,j], Dy)
                
                elif (i < pointstemp-1) and (k == pointsyq-1):
                    Dt = logtemp[i+1] - logtemp[i]
                    Dy = ye[k] - ye[k-1]
                    p1 = firDer(pressure[k,i+1,j], pressure[k,i-1,j], Dt)
                    p2 = firDer(energy[k,i+1,j], energy[k,i-1,j], Dt)
                    p3 = firDer(pressure[k,i,j], pressure[k-1,i,j], Dy/2)
                    p4 = firDer(energy[k,i,j], energy[k-1,i,j], Dy/2)

                elif (i == pointstemp-1) and (k == pointsyq-1):
                    Dt = logtemp[i] - logtemp[i-1]
                    Dy = ye[k] - ye[k-1]
                    p1 = firDer(pressure[k,i,j], pressure[k,i-1,j], Dt/2)
                    p2 = firDer(energy[k,i,j], energy[k,i-1,j], Dt/2)
                    p3 = firDer(pressure[k,i,j], pressure[k-1,i,j], Dy/2)
                    p4 = firDer(energy[k,i,j], energy[k-1,i,j], Dy/2)

                else:
                    Dt = logtemp[i+1] - logtemp[i]
                    Dy = ye[k+1] - ye[k]
                    p1 = firDer(pressure[k,i+1,j], pressure[k,i-1,j], Dt)
                    p2 = firDer(energy[k,i+1,j], energy[k,i-1,j], Dt)
                    p3 = firDer(pressure[k+1,i,j], pressure[k-1,i,j], Dy)
                    p4 = firDer(energy[k+1,i,j], energy[k-1,i,j], Dy)
                
                if p2 == 0:
                    p1 = 0
                    p2 = 1
                if p4 == 0:
                    p3 = 0
                    p4 = 1
                out[k,i,j] = p1/p2 + p3/p4 # Nel caso usare i log di P ed E e moltiplicare out per P/E

    return out*10**pressure/10**energy

fname = sys.argv[1]
qty = sys.argv[2]

with h5.File(fname, 'r') as f:
    logtmp = f['logtemp'][()]
    lognb  = f['lognb'][()]
    ye     = f['Yq'][()]

    LogEntropy = f['LogEntropy'][()].T
    LogPressure= f['LogPressure'][()].T
    LogEnergy  = f['LogEnergy'][()].T

if qty == 'betaV':
    out = compBetaV(lognb, 10**LogEntropy, len(ye), len(logtmp), len(lognb))
if qty == 'kappaT':
    out = compKappaT(lognb, 10**LogPressure, len(ye), len(logtmp), len(lognb))
if qty == 'cV':
    out = compCV(logtmp, 10**LogEntropy, len(ye), len(logtmp), len(lognb))
if qty == 'dPdE_nb':
    out = compdPdE_nb(logtmp, ye, LogPressure, LogEnergy, len(ye), len(logtmp), len(lognb))

out = np.where(out > 0, np.log10(out), -15)
out[np.isinf(out)] = -15
with h5.File(fname, 'a') as f:
    f.create_dataset('O_' + qty, data=out.T, compression='gzip', compression_opts=9)    
