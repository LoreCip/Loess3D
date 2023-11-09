#!/bin/bash

mkdir $runpath/$name
tracker="$runpath/$name/tracker.out"

# Python
PYTHON=$runpath/venv/bin/python3

# Path of EOS file
EOSpath="/nfs/scratch/lorecip/Loess3D_Production/EOS/FOP(SFHoY).h5"

# Add imports
export MPI_INCLUDE=/nfs/scratch/lorecip/.localMPICH4.1.1/include
export MPI_LIB=/nfs/scratch/lorecip/.localMPICH4.1.1/lib
export MPI_HOME=/nfs/scratch/lorecip/.localMPICH4.1.1
export LD_LIBRARY_PATH=/nfs/scratch/lorecip/Loess3D_Production/libs:/nfs/scratch/lorecip/.localzlib1.2.13/lib:/nfs/scratch/lorecip/.localHDF5_1.14.1/lib:/nfs/scratch/lorecip/.localMPICH4.1.1/lib:/nfs/scratch/lorecip/.localOpenSSL3.1.0/lib64:/nfs/scratch/lorecip/.gcc-12.3.0/lib64:$LD_LIBRARY_PATH
export LD_RUN_PATH=/nfs/scratch/lorecip/.localHDF5_1.14.1/lib

# Set ENVS
export HDF5_USE_FILE_LOCKING='FALSE'
export OMP_WAIT_POLICY=active
export OMP_PROC_BIND=true

dd=$(date)
echo "Starting: $dd" >> $tracker

# Copy and compile
cp -r $runpath/Code $runpath/$name
cd $runpath/$name/Code
OC=/nfs/scratch/lorecip/.localHDF5_1.14.1/bin/h5fc
make FCOMP=$OC
cd $runpath/$name

echo "Preparing initial data." >> $tracker

# Prepare h5 with pressure and entropy
$PYTHON $runpath/$name/Code/PScripts/generate_id.py $EOSpath $nthreads $frac $degree
mv data.h5 $runpath/$name/

# Entropy
echo "Smoothing entropy..." >> $tracker
date >> $tracker
mpiexec -n $nexec $runpath/$name/Code/bin/run $runpath/$name/data.h5 O_LogEntropy >> $tracker

# Pressure
echo "Smoothing pressure..." >> $tracker
date >> $tracker
mpiexec -n $nexec $runpath/$name/Code/bin/run $runpath/$name/data.h5 O_LogPressure

echo "Computing derived quantities..." >> $tracker
# Compute betaV
$PYTHON $runpath/$name/Code/PScripts/computeQuantities.py $runpath/$name/data.h5 betaV
# Compute kappaT
$PYTHON $runpath/$name/Code/PScripts/computeQuantities.py $runpath/$name/data.h5 kappaT
# Compute cV
$PYTHON $runpath/$name/Code/PScripts/computeQuantities.py $runpath/$name/data.h5 cV

# betaV
echo "Smoothing betaV..." >> $tracker
date >> $tracker
mpiexec -n $nexec $runpath/$name/Code/bin/run $runpath/$name/data.h5 O_betaV

# kappaT
echo "Smoothing kappaT..." >> $tracker
date >> $tracker
mpiexec -n $nexec $runpath/$name/Code/bin/run $runpath/$name/data.h5 O_kappaT

# cV
echo "Smoothing cV..." >> $tracker
date >> $tracker
mpiexec -n $nexec $runpath/$name/Code/bin/run $runpath/$name/data.h5 O_cV