clear all; close all; clc; 
addpath(genpath('..'));
ConstantsUnits0; 
cols=flatcolors();

N = 151; 

ef = .1; kf = abs(ef)/hbar_eV/vf;
ene = linspace(0.01,2.75,N)*abs(ef);
gam = 1e-3; 
kbT = kb_eV*1; 
sigmaLRA = conducLRA(ene,ef,gam,kbT);
q = kf*1e-3;
sigmaRPA = conducRPA(ene,q,ef,gam);
sigmaRPA_T = conducRPA_Maldague(ene,q,ef,gam,kbT);


%%
set_figsize(1,30,20)
subaxis(1,2,1)
plot(ene,real(sigmaLRA),'-','color',cols{14}); hold on
plot(ene,real(sigmaRPA),'--','color',cols{21}); 
plot(ene,real(sigmaRPA_T),'--','color',cols{3}); hold off
ylim([-.03,1]*2e-3); xlim(minmax(ene))

subaxis(1,2,2)
plot(ene,imag(sigmaLRA),'-','color',cols{14}); hold on
plot(ene,imag(sigmaRPA),'--','color',cols{21}); 
plot(ene,imag(sigmaRPA_T),'--','color',cols{3}); hold off
ylim([-.15,1]*1e-3); xlim(minmax(ene))

return
%%
set_figsize(2,20,20)

q = linspace(0,3,N)*kf;
[Q,ENE] = meshgrid(q,ene);
[~,~,CHIRPA] = conducRPA(ENE,Q,ef,gam);
pcolor(Q,ENE,real(CHIRPA))
shading interp; 
colormap(bluewhitered)