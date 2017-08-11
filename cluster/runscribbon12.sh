lattype="triangular";numcells="12";Ns="58";
echo "lattice=$lattype"; \
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v lattype=$lattype,numcells=$numcells,Ns=$Ns -N "scribbon_${numcells}_${Ns}" submit_scissorribbon.sh; \