lattype="triangular";numeigs="70";numcells="6";Ns="69";
echo "lattice=$lattype"; \
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v lattype=$lattype,numcells=$numcells,Ns=$Ns,numeigs=$numeigs -N "${lattype}_ribbon_${numcells}_${Ns}" submit.sh; \