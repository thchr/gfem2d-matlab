numcells="20";Ns="70";knum="2";kxa=".15+2/60";kxb=".3+2/60";Nk="10";
echo "numcells=${numcells}"; \
echo "Ns=${Ns}"; \
qsub -z -v numcells=$numcells,Ns=$Ns,knum=$knum,kxa=$kxa,kxb=$kxb,Nk=$Nk -N "scribbon_20_k$knum" submitsplitk.sh; \