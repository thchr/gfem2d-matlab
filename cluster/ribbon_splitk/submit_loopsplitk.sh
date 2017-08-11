#!/bin/sh

# -- choices for queue, pick either fotonano or hpc
#PBS -q fotonano
#PBS -A fotonano

# embedded options to qsub - start with #(PBS) (-l feature=XeonE5-2660 can be used to request 20 new nodes)
#PBS -M tomch@fotonik.dtu.dk
#PBS -m ae
#PBS -l nodes=1:ppn=10
#PBS -l mem=60gb
#PBS -l walltime=75:00:00

# -- run in the current working (submission) directory --
cd "$PBS_O_WORKDIR/../../execute"

# here follow the commands you want to execute
matlab -nodisplay -nosplash -r "runRibbonScissorLatticeMagneto('triangular',$numcells,0.401267289091325,$Ns,[],[],linspace($ss*$numk-$numk,$ss*$numk-1,$numk)/($numsplit*$numk-1)/2,[],[],'split$ss')"  -logfile "../output/logs/scrib_${numcells}_ss${ss}.tex"