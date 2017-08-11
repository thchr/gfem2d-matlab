#!/bin/sh

# -- choices for queue, pick either fotonano or hpc
#PBS -q fotonano
#PBS -A fotonano

# embedded options to qsub - start with #(PBS) (-l feature=XeonE5-2660 can be used to request 20 new nodes)
#PBS -M tomch@fotonik.dtu.dk
#PBS -m ae
#PBS -l nodes=1:ppn=10
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

# -- run in the current working (submission) directory --
cd "$PBS_O_WORKDIR/../../execute"

kxlist="linspace(0,pi,${numkx})";
kylist="(${ky}-1)/(${numky}-1)*2*pi/sqrt(3)";
savename="${ky}of${numky}";
# here follow the commands you want to execute
matlab -nodisplay -nosplash -r "runLatticeMagneto('$lattype',{'${kxlist}','${kylist}'},$Ns,$Ncirc,[],[],[],[],'${savename}')"  -logfile "../output/logs/2d_manproj_${Ns}_${savename}.tex"
