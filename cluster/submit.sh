#!/bin/sh

# -- choices for queue, pick either fotonano or hpc
#PBS -q fotonano
#PBS -A fotonano

# embedded options to qsub - start with #(PBS) (-l feature=XeonE5-2660 can be used to request 20 new nodes)
#PBS -M tomch@fotonik.dtu.dk
#PBS -m ae
#PBS -l nodes=1:ppn=10
#PBS -l mem=24gb
#PBS -l walltime=48:00:00

# -- run in the current working (submission) directory --
cd "$PBS_O_WORKDIR/../execute"

# here follow the commands you want to execute
matlab -nodisplay -nosplash -r "runRibbonLatticeMagneto('$lattype',$numcells,$Ns,[],2,[],1)"  -logfile "../output/logs/mphi_ribbon_${lattype}_${numcells}_${Ns}.tex"
