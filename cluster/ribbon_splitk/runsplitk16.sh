numcells="16";Ns="53";numk="5";numsplit="18";
echo "numk=${numk}"; \
echo "numsplit=${numsplit}"; \
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
for ss in {1..18}; 
	do	echo "   split, ss = ${ss}";
	qsub -z -v numcells=$numcells,Ns=$Ns,ss=$ss,numk=$numk,numsplit=$numsplit -N "scrib_16_ss$ss" submit_loopsplitk.sh; \
done;
	