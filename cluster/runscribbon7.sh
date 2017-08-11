lattype="triangular";numcells="7";Ns="65";
echo "lattice=$lattype"; \
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v lattype=$lattype,numcells=$numcells,Ns=$Ns -N "sc60ribbon_${numcells}_${Ns}" submit_scissorribbon.sh; \
