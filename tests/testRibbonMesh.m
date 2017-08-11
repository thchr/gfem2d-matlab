clc; clear all; close all; 
cols=flatcolors;
%addpath(genpath('..'))

%Inclusion boundary
adivd = 2;
Ncirc = 25; 
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
inclusion = [];1/(2*adivd)*[cos(theta);sin(theta)].';

blochmesh=geomRibbonInclusion('triangular',[1,1],inclusion,15,[],'remesh');


%Momentum list
Nk = 10; %Number of k-points
kvec = bsxfun(@times,linspace(0,1,Nk)',blochmesh.Grib);

%%
kk=3;
V=calcRibbonCoulomb(blochmesh,kvec(kk,:));
%%
clc
Vd=calcRibbonCoulombDirectSum(blochmesh,kvec(kk,:),5000);
%%
clc; close all
rho = zeros(size(blochmesh.p,1),1); rho(14)=1;
phi = V*rho;
phid = Vd*rho;

pn = []; tn = []; phin = []; phind = []; nL = 0;-1:1;
for nn = 1:numel(nL);
    n = nL(nn);
    pn = [pn; bsxfun(@plus,blochmesh.remesh.p,n*blochmesh.Rrib)];
    tn = [tn; blochmesh.remesh.t + (nn-1)*size(blochmesh.remesh.p,1)];
    phin = [phin; blochmesh.remesh.values(phi) ];
    phind = [phind; blochmesh.remesh.values(phid) ];
end


%plotVertexData(pn,tn,real(phin),{'colormap(morgenstemning)','axis normal on','axis vis3d'})
trisurf(tn,pn(:,1),pn(:,2),(imag(phin-0*min(phin))),'EdgeColor',cols{1},'FaceColor','none');
hold on
trisurf(tn,pn(:,1),pn(:,2),(imag(phind-0*min(phind))),'EdgeColor',cols{2},'FaceColor','none');
hold off

figure
trisurf(tn,pn(:,1),pn(:,2),(abs((phin-phind)./phin)),'EdgeColor',cols{1},'FaceColor','none');


