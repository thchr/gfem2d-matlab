lattype="triangular";numeigs="70";numcells="12";Ns="68";
echo "lattice=$lattype"; \
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v lattype=$lattype,numcells=$numcells,Ns=$Ns,numeigs=$numeigs -N "${lattype}_ribbon_${numcells}_${Ns}" submit.sh; \