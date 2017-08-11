clear all; clc; close all;

addpath(genpath(pwd));
ConstantsUnits0; cols=flatcolors;
%hdata.fun = @(x,y) real(.005*(cos(pi/2*(x.^2+y.^2)).^.5+2).^2);
[p,t] = geomDisk();
%[p,t] = geomDisk(30,.085);drawnow;
%[p,t] = geomAnnulus();drawnow;
%[p,t] = geomEquiTriangle(1,.04);drawnow;
%[p,t] = geomHoleSquareArray(1,1,0);
%p = [0,0 ; 1,0 ; 0,1]; t = [1,2,3];

f = [1];

R{1} = [2.2222,0]; R{2} = [0,2.2222];
L=3600e-9;
theta=0;
ef_eV = 0.2;
gam_eV = (1/1e-12)*hbar_eV; 1.08e-3;
T = 300;
B = 5;

sigma = @(omega) conducLRA(omega*hbar_eV,ef_eV,gam_eV,kb_eV*T);
%sigma = @(omega) conducMagnetoLRA(omega*hbar_eV,ef_eV,sqrt(2*e*B*hbar)*vf/ev2jo,gam_eV,kb_eV*T);
%sigma = @(omega) calcMagnetoClassical(omega*hbar_eV,ef_eV,B,gam_eV);

%f = @(omega) sigma(omega)./sigmalra(omega);

omega = linspace(.025,.35,250)*ef_eV/hbar_eV; lambda = 2*pi*c./omega;

zetaomega = 2i*eps0*omega*L./sigma(omega);

[Q_abs,~,W,DL,DR] = calcLatticeAbsorption(p,t,R,omega,L,theta,sigma,1,2);
[Q_abs0,~,V,DL0,DR0] = calcAbsorption(p,t,omega,L,theta,sigma,1);
%%
figure
semilogy(lambda*1e6,Q_abs,'-','color',cols{1}); hold on;
semilogy(lambda*1e6,Q_abs0,'-','color',cols{2}); hold on;
semilogy(lambda*1e6,omega/c.*imag( 2*L^3 * ( 2.8912 ./ (1.0977 - zetaomega) + ...
    0.1120 ./ (4.9140 - zetaomega) + ...
    0.0424 ./ (8.1337 - zetaomega) + ...
    0.0224 ./ (11.3079 - zetaomega) + ...
    0.0140 ./ (14.4675 - zetaomega) + ...
    0.0096 ./ (17.6205 - zetaomega) ) )/(pi*L^2),':k');
xlim([20 120])
%%
%%
%close all
%[DR,DL,D] = calcDifferentialPaper(p,t,f);
%V = calcCoulomb(p,t);
D0 = -DL0\DR0;
%D0 = -DL0\DR0;
[zeta0,eigV0] = calcEigen(D0,V); zeta(1:15);
set_figsize([],18,25.2)
for nn = 1:48;
    subaxis(8,6,nn,'S',0.03,'M',0.01,'MT',0.025)
    plotVertexData(p,t,V*eigV0(:,nn),'colormap(bluewhitered_mod)')
    title([num2str(nn) '| ' num2str(zeta0(nn))])
end

%%
D = -DL\DR;
%D0 = -DL0\DR0;
[zeta,eigV] = calcEigen(D,W); zeta(1:15);
set_figsize([],18,25.2)
for nn = 1:48;
    subaxis(8,6,nn,'S',0.03,'M',0.01,'MT',0.025)
    plotVertexData(p,t,V*eigV(:,nn),'colormap(bluewhitered_mod)')
    title([num2str(nn) '| ' num2str(real(zeta(nn)))])
end

%%
N = 15;
guess = .01;
for ee=1:N
    ene_eig_eV(ee) = fzero(@(omega_eV) real( 2i*eps0*(omega_eV/hbar_eV)*L./sigma(omega_eV/hbar_eV)-zeta(ee)),guess);
end

ene_eig_eV
2*pi*c./(ene_eig_eV/hbar_eV)*1e6

areas = meshArea(p,t);
vertInt = [2,1,1;1,2,1;1,1,2]/12; 
pxt = p(:,1); pxt = bsxfun(@times,areas,pxt(t)).';
pyt = p(:,2); pyt = bsxfun(@times,areas,pyt(t)).'; 
for ee=1:N
    rho=eigV0(:,ee);
dipmom(ee,1) = sum(sum(pxt.*(vertInt*rho(t.')),1),2); %Summing contributions from every triangle [to the dipole moment at frequency omega(oo)]
    dipmom(ee,2) = sum(sum(pyt.*(vertInt*rho(t.')),1),2); %with exact integration for linear basis functions
end
log10(sqrt(abs(dipmom(:,1)).^2 + abs(dipmom(:,2)).^2))
%[zeta,W] = calcEigen(V,D);
%disp(zeta(1:12))
