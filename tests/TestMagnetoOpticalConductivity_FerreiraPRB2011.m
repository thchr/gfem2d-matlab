clc;clear all;close all;

addpath(genpath('...'));

ConstantsUnits0;
vf = 1e6; %Value used by Ferreira 
ene = linspace(0.01,650e-3,1500);
ef_eV = 0; 
gam_eV = 6.8e-3;
T = 0;
B= 7; eneB_eV = sqrt(2*e*B*hbar)*vf/ev2jo;

sigma =  conducMagnetoLRA(ene,ef_eV,eneB_eV,gam_eV,kb_eV*T);
sigmaxx = sigma(1,1,:);
sigmaxy = sigma(1,2,:); 

plot(minmax(ene),[0 0],'-k')
hold on; 
plot(ene,real(sigmaxy(:))/(e^2/(hbar*2*pi)),'-b')
plot(ene,imag(sigmaxy(:))/(e^2/(hbar*2*pi)),'-r')
hold off

xlim(minmax(ene))
ylim([-14,7])

