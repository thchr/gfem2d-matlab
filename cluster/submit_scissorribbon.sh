#!/bin/sh

# -- choices for queue, pick either fotonano or hpc
#PBS -q fotonano
#PBS -A fotonano

# embedded options to qsub - start with #(PBS) (-l feature=XeonE5-2660 can be used to request 20 new nodes)
#PBS -M tomch@fotonik.dtu.dk
#PBS -m ae
#PBS -l nodes=1:ppn=10
#PBS -l mem=20gb
#PBS -l walltime=84:00:00

# -- run in the current working (submission) directory --
cd "$PBS_O_WORKDIR/../execute"

# here follow the commands you want to execute
matlab -nodisplay -nosplash -r "runRibbonScissorLatticeMagneto('$lattype',$numcells,0.401267289091325,$Ns,[],[],[],[])"  -logfile "../output/logs/sc40ribbon_${lattype}_${numcells}_${Ns}.tex"
