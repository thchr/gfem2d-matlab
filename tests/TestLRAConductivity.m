%Test case is Figure 17 of your PhD thesis.

clear all; close all; clc

addpath(genpath('...'));

cols=flatcolors; ConstantsUnits0;
ef = 0.4; 
gam = [0,12e-3,12e-3];
T = [0,0,0.026];
ene = linspace(0,3*ef,300);

set_figsize()
for nn=1:3;
    sig{nn} = conducLRA(ene,ef,gam(nn),T(nn));
    plot(ene,real(sig{nn}/(e^2/(4*hbar))),'-','color',cols{nn}); hold on
end
ylim([0 2])
xlim(minmax(ene))

set_figsize()
for nn=1:3;
    plot(ene,imag(sig{nn}/(e^2/(4*hbar))),'-','color',cols{nn}); hold on
end
plot(minmax(ene),[0 0],':k')
ylim([-2 3])
xlim(minmax(ene))