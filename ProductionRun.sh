#!/bin/bash

# Python
PYTHON=$runpath/../venv/bin/python3

# Path of EOS file
EOSpath="/nfs/scratch/lorecip/Loess3D/EOS/FOP(SFHoY).h5"

# Add imports
export MPI_INCLUDE=/nfs/scratch/lorecip/.localMPICH4.1.1/include
export MPI_LIB=/nfs/scratch/lorecip/.localMPICH4.1.1/lib
export MPI_HOME=/nfs/scratch/lorecip/.localMPICH4.1.1

export LD_LIBRARY_PATH=/nfs/scratch/lorecip/Loess3D/mpiTesting:/nfs/scratch/lorecip/.localzlib1.2.13/lib:/nfs/scratch/lorecip/.localHDF5_1.14.1/lib:/nfs/scratch/lorecip/.localMPICH4.1.1/lib:/nfs/scratch/lorecip/.localOpenSSL3.1.0/lib64:/nfs/scratch/lorecip/.gcc-12.3.0/lib64:$LD_LIBRARY_PATH

export LD_RUN_PATH=/nfs/scratch/lorecip/.localHDF5_1.14.1/lib

export HDF5_USE_FILE_LOCKING='FALSE'

# Compile
OC=/nfs/scratch/lorecip/.localHDF5_1.14.1/bin/h5fc
make FCOMP=$OC

# Prepare h5 with pressure and entropy
$PYTHON generate_id.py $EOSpath $nthreads $frac $degree

# Entropy
mpiexec -n 4 $runpath/bin/run data.h5 O_LogEntropy

# Pressure
mpiexec -n 4 $runpath/bin/run data.h5 O_LogPressure

# Compute betaV
$PYTHON computeQuantities.py data.h5 betaV
# Compute kappaT
$PYTHON computeQuantities.py data.h5 kappaT
# Compute cV
$PYTHON computeQuantities.py data.h5 cV

# betaV
mpiexec -n 4 $runpath/bin/run data.h5 O_betaV

# kappaT
mpiexec -n 4 $runpath/bin/run data.h5 O_kappaT

# cV
mpiexec -n 4 $runpath/bin/run data.h5 O_cV