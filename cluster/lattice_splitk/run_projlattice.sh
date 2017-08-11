lattype="triangular";Ns="149";Ncirc="64";numky=49;numkx=49;
echo "lattice=$lattype"; \
echo "Ns=${Ns} | Ncirc=$Ncirc"; \
for ky in $(seq 1 ${numky}); 
	do 	echo "   ky=$ky"; \
	qsub -z -v lattype=$lattype,ky=$ky,numky=$numky,numkx=$numkx,Ns=$Ns,Ncirc=$Ncirc -N "2d_manproj_${Ns}_${ky}of${numky}" submit_projlattice.sh; \
done
