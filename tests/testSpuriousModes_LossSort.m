clear all; close all; clc;
cols = flatcolors; 
ef_eV = .2; L = 400e-9; omegac_eV = @(B) 5.403927588936477e-04*B/ef_eV;

models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',[] ,'system',[]        ,'approx','intra','keep_eV',[0,ef_eV/2.5],'gamc',[]);
models(2) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',8  ,'system','graphene','approx',[]     ,'keep_eV',[0,ef_eV/2.5],'gamc',.001); %The lower cutoff is the cyclotron frequency
models(3) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',12 ,'system','graphene','approx',[]     ,'keep_eV',[0,ef_eV/2.5],'gamc',.001); %---||---


Ns = 70; Ncirc = 32; Nk = 75; savename = 'testspurious_loss';
savepath = runLatticeMagneto([],'irrfbz',Ns,Ncirc,[],Nk,[],models,savename);

%%
%load(savepath)
load('../output/lattice/lattice_irrfbz_Ns149_Ncirc64_savephi0.mat')

set_figsize(2,25,25)
subaxis(1,3,1)
plot(kplot,eneeig_eV{1},'.k')

for mm = 2:3
    spurio_eV{mm} = eneeig_eV{mm}; 
%    spurio_eV{mm}(-imag(eneeig_eV{mm}) < models(mm).gamc*omegac_eV(models(mm).B)*.9 & eneeig_eV{mm} > omegac_eV(models(mm).B)) = NaN;
%    eneeig_eV{mm}(-imag(eneeig_eV{mm}) > models(mm).gamc*omegac_eV(models(mm).B)*.9 & eneeig_eV{mm} < omegac_eV(models(mm).B)) = NaN;

subaxis(1,3,mm)
plot(kplot,real(spurio_eV{mm})/omegac_eV(models(mm).B),'.','color',cols{8}); hold on
plot(kplot,real(eneeig_eV{mm})/omegac_eV(models(mm).B),'.k')
plot(minmax(kplot),[1,1],'-','color',cols{4}); hold off
view(2); ylim(models(mm).keep_eV/omegac_eV(models(mm).B)); set(gca,'Layer','Top'); box on
end
%export_fig('SpuriousIllustration','-pdf')