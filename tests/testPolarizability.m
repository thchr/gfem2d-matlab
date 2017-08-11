clear all; close all; clc; 

addpath(genpath('..'))

ConstantsUnits0;
N = 250; shift = 1+.05;
%[~,p,t] = evalc('geomDisk(N,8.5/N,1)'); 
[~,p,t] = evalc('geomDisk(N,[],1)'); 
%[~,p,t] = evalc('geomEquiTriangle([],''const'',1)'); 
%nodes = {[0,0; 1,0; 1,1; 0,1],+[0,shift; 1,+shift; 1,1+shift; 0,1+shift]};
%triangle = [cosd([-30;90;210]),sind([-30;90;210])]; trishift = 1+.0125;
%nodes = {[triangle(:,1),triangle(:,2)-trishift],[triangle(:,1),-triangle(:,2)+trishift]};
%[~,p,t] = evalc('geomPolygon(nodes,.0325,1)');

mesh.p = p; mesh.t = t;
%
ef_eV = .4; 
gam_eV = 12e-3; 
L = 30e-9; 
ene_eV = linspace(0,2,2500)*ef_eV; omega = ene_eV/hbar_eV; 
sigma =  conducLRA(ene_eV,ef_eV,gam_eV,300*kb_eV);
zeta = 2i*eps0*omega*L./sigma; 
%%
%mesh.V = calcCoulombNumQuad(p,t);
[alpha,eigstruct] = calcPolarizability(mesh,[],1);

%%
cols=flatcolors;
set_figsize([],25,15);
plot(ene_eV,omega/c.*imag(alpha.x(zeta,L))/(pi*L^2),'color',cols{1}); hold on
plot(ene_eV,omega/c.*imag(alpha.y(zeta,L))/(pi*L^2),'-','color',cols{24})

zetacheck=1./abs(zeta - eigstruct.zeta(1));
%plot(ene_eV,zetacheck/max(zetacheck)*.25,':k')
 plot(omega*hbar_eV,omega/c.*imag( 2*L^3 * ( 2.8912 ./ (1.0977 - zeta) + ...
                                             0.1120 ./ (4.9140 - zeta) + ...
                                             0.0424 ./ (8.1337 - zeta) ) )/(pi*L^2),':k'); 
hold off
%%
plotVertexData(p,t,eigstruct.phi(:,23))