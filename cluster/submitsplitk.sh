#!/bin/sh

# -- choices for queue, pick either fotonano or hpc
#PBS -q hpc
#PBS -A hpc

# embedded options to qsub - start with #(PBS) (-l feature=XeonE5-2660 can be used to request 20 new nodes)
#PBS -M tomch@fotonik.dtu.dk
#PBS -m ae
#PBS -l nodes=1:ppn=20
#PBS -l mem=72gb
#PBS -l walltime=75:00:00

# -- run in the current working (submission) directory --
cd "$PBS_O_WORKDIR/../execute"

# here follow the commands you want to execute
matlab -nodisplay -nosplash -r "runRibbonScissorLatticeMagneto('triangular',$numcells,0.401267289091325,$Ns,[],[],linspace($kxa,$kxb,$Nk),[],[],'$knum')"  -logfile "../output/logs/scribbon_${numcells}_k${knum}.tex"