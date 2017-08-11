clear all; clc; close all; 

addpath(genpath(pwd));
ConstantsUnits0;
[p,t] = geomDisk(50);drawnow;
%[p,t] = geomEquiTriangle(1,.04);drawnow;
%[p,t] = geomHoleSquareArray(1,1,0);
%p = [0,0 ; 1,0 ; 0,1]; t = [1,2,3];
    
f = [1];

L=50e-9;
theta=0;
ef_eV = 0.5; 
gam_eV = 1.08e-3;
T = 4.2;
B = 5;

sigmalra = @(omega) conducLRA(omega*hbar_eV,ef_eV,gam_eV,kb_eV*T);
sigma = @(omega) conducMagnetoLRA(omega*hbar_eV,ef_eV,sqrt(2*e*B*hbar)*vf/ev2jo,gam_eV,kb_eV*T);
%sigma = @(omega) calcMagnetoClassical(omega*hbar_eV,ef_eV,B,gam_eV);

f = @(omega) sigma(omega)./sigmalra(omega);

omega = linspace(.15,.19,100)/hbar_eV;
%omega = linspace(.19,.24,300)/hbar_eV;
zetaomega = 2i*eps0*omega*L./sigmalra(omega);

Q_abs = calcAbsorption(p,t,omega,L,theta,sigmalra,f);
Q_abs_B0 = calcAbsorption(p,t,omega,L,theta,sigmalra);

%%
figure
plot(omega*hbar_eV,Q_abs,'-r'); hold on; 
plot(omega*hbar_eV,Q_abs_B0,'-b'); hold on; 
plot(omega*hbar_eV,omega/c.*imag( 2*L^3 * ( 2.8912 ./ (1.0977 - zetaomega) + ...
                                            0.1120 ./ (4.9140 - zetaomega) + ...
                                            0.0424 ./ (8.1337 - zetaomega) ) )/(pi*L^2),':k'); 
xlim([.15 .19])
%%
%%
%[DR,DL,D] = calcDifferentialPaper(p,t,f);
%V = calcCoulomb(p,t);

%[zeta,W] = calcEigen(V,D);
%disp(zeta(1:12))
