lattype="triangular";numcells="11";Ns="41";
echo "lattice=$lattype"; \
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v lattype=$lattype,numcells=$numcells,Ns=$Ns -N "sc40ribbon_${numcells}_${Ns}" submit_scissorribbon.sh; \
