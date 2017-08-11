numcells="20";Ns="70";knum="3";kxa="11/60";kxb="15/60";Nk="5";
echo "knum=${knum}"; \
echo "kxa=${kxa}"; \
echo "kxb=${kxb}"; \
echo "Nk=${Nk}"; \
qsub -z -v numcells=$numcells,Ns=$Ns,knum=$knum,kxa=$kxa,kxb=$kxb,Nk=$Nk -N "scribbon_20_k$knum" submitsmallsplitk.sh; \