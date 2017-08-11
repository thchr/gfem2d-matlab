clear all; clc; close all;

addpath(genpath(pwd));
ConstantsUnits0; cols=flatcolors;




R{1} = [3,0]; R{2} = [0,3];
L=1000e-9;

%[p,t] = geomDisk();
%[p,t] = geomAnnulus();drawnow;
%[p,t] = geomEquiTriangle(1,.04);drawnow;
%[p,t] = geomHoleSquareArray(4,1,[],[],struct('fun',@(x,y) .0125.*(sqrt(x.^2 + y.^2).^1.5+2)));
%[p,t] = geomHoleSquareArray(3,1);
[p,t] = geomSquare(3,3,.04);

theta=0;
ef_eV = 0.4;
gam_eV = (1/1e-12)*hbar_eV; 1.08e-3;
T = 300;
B = 5;

sigma = @(omega) conducLRA(omega*hbar_eV,ef_eV,gam_eV,kb_eV*T);
%sigma = @(omega) conducMagnetoLRA(omega*hbar_eV,ef_eV,sqrt(2*e*B*hbar)*vf/ev2jo,gam_eV,kb_eV*T);
%sigma = @(omega) calcMagnetoClassical(omega*hbar_eV,ef_eV,B,gam_eV);

%f = @(omega) sigma(omega)./sigmalra(omega);

omega = linspace(.062,.0177,200)/hbar_eV; lambda = 2*pi*c./omega; invcm = omega*hbar_eV*8065.54;

zetaomega = 2i*eps0*omega*L./sigma(omega);

[Q_abs,~,W,DL,DR] = calcLatticeAbsorption(p,t,R,omega,L,theta,sigma,1,4);
[Q_abs0,~,V,DL0,DR0] = calcAbsorption(p,t,omega,L,theta,sigma,1);
%%
figure
plot(invcm,Q_abs,'-','color',cols{1}); hold on;
plot(invcm,Q_abs0,':','color',cols{2}); hold on;
%ylim([0,.46])
%xlim([20 70])
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
    plotVertexData(p,t,abs(V*eigV0(:,nn)),'colormap(bluewhitered_mod)')
    title([num2str(nn) '| ' num2str(zeta0(nn))])
end

%%
D = -DL\DR;
%D0 = -DL0\DR0;
[zeta,eigV] = calcEigen(D,W); zeta(1:15);
set_figsize([],18,25.2)
for nn = 1:48;
    subaxis(8,6,nn,'S',0.03,'M',0.01,'MT',0.025)
    plotVertexData(p,t,real(V*eigV(:,nn)),'colormap(bluewhitered_mod)')
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
    rho=eigV(:,ee);
dipmom(ee,1) = sum(sum(pxt.*(vertInt*rho(t.')),1),2); %Summing contributions from every triangle [to the dipole moment at frequency omega(oo)]
    dipmom(ee,2) = sum(sum(pyt.*(vertInt*rho(t.')),1),2); %with exact integration for linear basis functions
end
(sqrt(abs(dipmom(:,1)).^2 + abs(dipmom(:,2)).^2))
%[zeta,W] = calcEigen(V,D);
%disp(zeta(1:12))
