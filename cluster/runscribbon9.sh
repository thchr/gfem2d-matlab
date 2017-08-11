lattype="triangular";numcells="9";Ns="49";
echo "lattice=$lattype"; \
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v lattype=$lattype,numcells=$numcells,Ns=$Ns -N "sc40ribbon_${numcells}_${Ns}" submit_scissorribbon.sh; \
