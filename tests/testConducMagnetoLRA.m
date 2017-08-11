clc; close all; clear all;
%Test case equal to Figure 1 in the publication by Gusynin et al., J. Phys. 
%Condens. Matter 19 (2007) 026222.

addpath(genpath('...'));

ConstantsUnits0; cols={'r','k','b','g'};flatcolors();
vf=1e6; %Value used by Gusynin et al.
invcm2eV = 1/8065.54;
ene = linspace(0.1,1300,750)*invcm2eV;
ef = kb_eV * [50,50,510,660]; 
B = [ 0.001, 1, 1, 1]; 
kbT = kb_eV*10; 
gam = kb_eV*15;
eneB = sqrt(2*e*B*hbar)*vf/ev2jo;

for c = 1:4
    sigmaLRA(:,c) = conducLRA(ene,ef(c),gam,kbT);
    [~,sigmaxx(:,c),sigmaxy(:,c),Nmax(c)] = conducMagnetoLRA(ene,ef(c),eneB(c),gam,kbT);
    plot(ene/invcm2eV,real(sigmaxx(:,c)*hbar/(2*pi)/e^2),'-','color',cols{c},'LineWidth',1), hold on
    plot(ene/invcm2eV,real(sigmaLRA(:,c)*hbar/(2*pi)/e^2),':','color',cols{c},'LineWidth',1), hold on

end
xlim([0 1300])
ylim([0 .9])
    drawnow
