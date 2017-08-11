numcells="20";Ns="70";knum="3";kxa=".3+3/60";kxb=".5";Nk="10";
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v numcells=$numcells,Ns=$Ns,knum=$knum,kxa=$kxa,kxb=$kxb,Nk=$Nk -N "scribbon_20_k$knum" submitsplitk.sh; \