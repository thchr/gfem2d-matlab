numcells="20";Ns="70";knum="4";kxa="16/60";kxb="20/60";Nk="5";
echo "knum=${knum}"; \
echo "kxa=${kxa}"; \
echo "kxb=${kxb}"; \
echo "Nk=${Nk}"; \
qsub -z -v numcells=$numcells,Ns=$Ns,knum=$knum,kxa=$kxa,kxb=$kxb,Nk=$Nk -N "scribbon_20_k$knum" submitsmallsplitk.sh; \