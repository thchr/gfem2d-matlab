clear all; close all; clc;
load('../output/lattice/lattice_irrfbz_Ns149_Ncirc64_savephi0.mat')
cols = flatcolors; models(1).B = 0; 

cut{1} = 1:size(eneeig_eV{1},1);
cut{2} = 1:size(eneeig_eV{2},1);
cut{3} = 1:size(eneeig_eV{3},1);
for mm = 1:3
    set_figsize(mm,35,30)
plot(kplot,real(eneeig_eV{mm}(cut{mm},:)),'.k'); hold on
plot(minmax(kplot),[1,1]*omegac_eV(models(mm).B),'-','color',cols{8},'LineWidth',2); hold off
view(2); ylim([0,.12]); xlim(minmax(kplot)); set(gca,'Layer','Top'); box on
end
%export_fig('SpuriousIllustration','-pdf')