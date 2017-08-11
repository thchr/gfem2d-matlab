numcells="20";Ns="70";knum="5";kxa="21/60";kxb="25/60";Nk="5";
echo "knum=${knum}"; \
echo "kxa=${kxa}"; \
echo "kxb=${kxb}"; \
echo "Nk=${Nk}"; \
qsub -z -v numcells=$numcells,Ns=$Ns,knum=$knum,kxa=$kxa,kxb=$kxb,Nk=$Nk -N "scribbon_20_k$knum" submitsmallsplitk.sh; \