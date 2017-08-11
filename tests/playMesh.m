clear all; clc; close all; 

addpath(genpath(pwd));

[p,t] = geomDisk(100);%,.0225);
%[p,t] = geomEquiTriangle(1,.02);
%[p,t] = geomHoleSquareArray(1,1,0);
%p = [0,0 ; 1,0 ; 0,1]; t = [1,2,3];
Nvert = size(p,1); Ntri = size(t,1); area = meshArea(p,t);
f = [1,.5;.5,1];
%%
[DR,DL,D] = calcDifferentialPaper(p,t,f);
V = calcCoulomb(p,t);

%%

[zeta,W] = calcEigen(V,D);

figure
plot(real(zeta),'.k'); drawnow;
zeta(1:10)




figure
rc = meshCentroid(p,t); 
for k = 1:numel(zeta)
    trisurf(t,p(:,1),p(:,2),real(W(:,k)),'EdgeColor','none'); shading interp; colormap(bluewhitered); title(num2str(zeta(k))); view(2); axis equal off; drawnow; 
    pause(1);
end
