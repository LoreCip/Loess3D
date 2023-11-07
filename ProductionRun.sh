#!bin/bash

PYTHON=$runpath/venv/bin/python3

## Add correct imports

export HDF5_USE_FILE_LOCKING='FALSE'

OC=/nfs/scratch/lorecip/.localHDF5_1.14.1/bin/h5fc
make FCOMP=$OC

# Prepare h5 with pressure and entropys

PYTHON generate_id.py $nthreads $frac $degree

# Entropy
mpiexec -n 4 $runpath/bin/run data.h5 LogEntropy
# Pressure
mpiexec -n 4 $runpath/bin/run data.h5 LogPressure

# Compute betaV
PYTHON computeQuantities.py data.h5 betaV
# Compute kappaT
PYTHON computeQuantities.py data.h5 kappaT
# Compute cV
PYTHON computeQuantities.py data.h5 cV

# betaV
mpiexec -n 4 $runpath/bin/run data.h5 betaV

# kappaT
mpiexec -n 4 $runpath/bin/run data.h5 kappaT

# cV
mpiexec -n 4 $runpath/bin/run data.h5 cV