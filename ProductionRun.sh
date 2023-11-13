#!/bin/bash

mkdir $runpath/$name
tracker="$runpath/$name/tracker.out"

# Script to compute the derived quantities
PY_SRC=$runpath/$name/Code/PScripts

# Launcher
LAUNCHER=$runpath/$name/Code/bin/run
if [ "${USE_MPI}"=true ]
then
    LAUNCHER="mpiexec -n $nexec $LAUNCHER"
fi

# Path of ID file
FILE=$runpath/$name/data.h5

# Which python to use
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
if [ "${USE_MPI}"=true ]
then
    make FCOMP=$OC
else
    make FCOMP=$OC USEMPI=
fi

cd $runpath/$name

# Prepare h5 with pressure and entropy
echo "Preparing initial data." >> $tracker
$PYTHON $PY_SRC/generate_id.py $EOSpath $nthreads $frac $degree
mv data.h5 $runpath/$name/

for qty in O_LogEntropy O_LogPressure O_LogEnergy
do

    echo "Smoothing $qty..." >> $tracker
    date >> $tracker
    $LAUNCHER $FILE $qty

done

echo "Computing derived quantities..." >> $tracker
date >> $tracker
for qty in betaV kappaT cV dPdE_nb
do

    $PYTHON $PY_SRC/computeQuantities.py $FILE $qty

done

for qty in O_betaV O_kappaT O_cV O_dPdE_nb
do

    echo "Smoothing $qty..." >> $tracker
    date >> $tracker
    $LAUNCHER $FILE $qty

done