numcells="20";Ns="70";knum="1";kxa="0";kxb=".15+1/60";Nk="11";
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v numcells=$numcells,Ns=$Ns,knum=$knum,kxa=$kxa,kxb=$kxb,Nk=$Nk -N "scribbon_20_k$knum" submitsplitk.sh; \