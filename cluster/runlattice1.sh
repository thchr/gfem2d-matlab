lattype="triangular";Ns="70";Ncirc="36";
echo "lattice=$lattype"; \
echo "Ns=${Ns} | Ncirc=$Ncirc"; \
for kspec in {"projx","irrfbz"}; 
	do 	echo "   kspec=${kspec}"; \
	qsub -z -v lattype=$lattype,kspec=$kspec,Ns=$Ns,Ncirc=$Ncirc -N "2d_${lattype}_${kspec}_${Ns}" submit_lattice.sh; \
done