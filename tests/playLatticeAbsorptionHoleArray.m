clear all; clc; close all;

addpath(genpath(pwd));
ConstantsUnits0; cols=flatcolors;

L=1250e-9;

%[p,t] = geomDisk();
%[p,t] = geomAnnulus();drawnow;
%[p,t] = geomEquiTriangle(1,.04);drawnow;
%[p,t] = geomHoleSquareArray(4,1,[],[],struct('fun',@(x,y) .0125.*(sqrt(x.^2 + y.^2).^1.5+2)));
%[p,t] = geomHoleSquareArray(3,1);
%[p,t] = geomSquare(3,3,.04);
%blochmesh = geomPeriodicInclusion('square','disk',200,'antidot');
blochmesh = geomPeriodicInclusion('square',[],100,'const');
R = blochmesh.R;

ef_eV = 0.2;
gam_eV = (1/6.37e-12)*hbar_eV; 
T = 300;

sigma = @(omega) conducLRA(omega*hbar_eV,ef_eV,gam_eV,kb_eV*T);
%sigma = @(omega) conducMagnetoLRA(omega*hbar_eV,ef_eV,sqrt(2*e*B*hbar)*vf/ev2jo,gam_eV,kb_eV*T);
%sigma = @(omega) calcMagnetoClassical(omega*hbar_eV,ef_eV,B,gam_eV);
%f = @(omega) sigma(omega)./sigmalra(omega);

omega = linspace(.062,.0177,100)/hbar_eV; lambda = 2*pi*c./omega; invcm = omega*hbar_eV*8065.54;

zetaomega = 2i*eps0*omega*L./sigma(omega);
[zeta,eigV,W,DL,DR,ene_eV,lambda_nm] = calcLatticeEigen(blochmesh,[0,0],[],L,sigma);

%trimesh(blochmesh.remesh.t,blochmesh.remesh.p(:,1),blochmesh.remesh.p(:,2),blochmesh.remesh.values(W(:,515)))

%%
vertInt = [2,1,1;1,2,1;1,1,2]/12; totalcharge = zeros(35,1);
set_figsize([],18,25.2)
for nn = 1:35;
    eigVnn = blochmesh.remesh.values(W*eigV(:,nn+1));
    plotval = imag(eigVnn);
    subaxis(7,5,nn,'S',0.03,'M',0.01,'MT',0.025)
    plotVertexData(blochmesh.remesh.p,blochmesh.remesh.t,plotval,'colormap(bluewhitered_mod)')
    caxis(max(abs(plotval(:)))*[-1,1])
    title([num2str(nn) '| ' num2str(real(zeta(nn+1)))])
    indcharge = blochmesh.remesh.values(eigV(:,nn)); 
    totalcharge(nn) = sum(vertInt*eigVnn(blochmesh.remesh.t.'),1)*blochmesh.remesh.area;
end


