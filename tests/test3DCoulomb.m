clear all; close all; clc;
ConstantsUnits0;
% MESH
mesh = geomEllipse([1,1],[],.075,1);

% PLOT EVALUATION REGION
Nxy = 65; Nz = 35;
x = linspace(-1.5,1.5,Nxy);
y = x;
z = linspace(-.5,.5,Nz);
[X,Y,Z] = meshgrid(x,y,z); r = [X(:), Y(:), Z(:)];
[Xz0,Yz0] = meshgrid(x,y);   rz0 = [Xz0(:),Yz0(:)];

% SYSTEM MATRICES
V = calcCoulomb(mesh.p,mesh.t);
Vmap = calcCoulombGeneralPoints(mesh.p,mesh.t,r);
Vmapz0 = calcCoulombGeneralPoints(mesh.p,mesh.t,rz0);
[DR,DL] = calcDifferential(mesh.p,mesh.t);
%%
%MODELS
    models.ef_eV = 0.2; models.gam_eV = 1e-3; models.kT_eV = kb_eV*300;
    models.ene = 48.5e-3; models.L = 200e-9;
    models.sigma = conducLRA(models.ene,models.ef_eV,models.gam_eV,models.kT_eV);
    models.approx = 'full'; models.type = 'isotropic'; models.keep_eV = [0,.2];
%eigene_eV = calcEigenEnergiesAny(mesh,models);
%% EXTERNAL POTENTIAL

r0 = [0,0,10]; pol = [1,0,0];
phiext = calcDipoleSource([mesh.p,zeros(size(mesh.p,1),1)],r0,pol);

%% RESPONSE
%Compute response due to each dipole in each model, and find induced field
for mm = 1:numel(models)
    for ee = 1:numel(models(mm).ene)
        rho =  ( (DL+1i*models(mm).sigma(ee)/(4*pi*eps0*models(mm).L*models(mm).ene(ee)/hbar_eV)*DR*V) ) \ ... %We write it this (lengthy) way, to avoid potential memory issues of the mesh is very large.
            ( -1i*models(mm).sigma(ee)/(models(mm).ene(ee)/hbar_eV*models(mm).L^2)*DR*phiext );
        phiind = Vmap*rho;
        phiindz0 = Vmapz0*rho; 
    end
end
phiind = reshape(phiind,size(X));
phiindz0 = reshape(phiindz0,size(Xz0));


%% SURFACE PLOT
set_figsize(2,20,20)
for zz = 1:Nz;
    fill3(cos(linspace(0,2*pi,50)),sin(linspace(0,2*pi,50)),zeros(1,50),[1,1,1]*.75,'facealpha',.5); hold on; 
    surf(X(:,:,zz),Y(:,:,zz),real(phiind(:,:,zz))); 
    surf(Xz0,Yz0,real(phiindz0),'facealpha',.25,'edgecolor','none'); hold off; 
    zlim(minmax(real(phiind)));
    caxis(minmax(real(phiind)));
    colormap(bluewhitered_mod)
    drawnow; title(z(zz)); pause(.5);
end

%% VOLUMETRIC PLOT
set_figsize(3,20,20)
v = real(phiind); 
isovals = -mean(abs(v(:)))*[6,4,2,1.5]; isovals = [isovals,-fliplr(isovals)];
caxis(minmax(isovals)); cmap = bluewhitered_mod(numel(isovals));
fill3(cos(linspace(0,2*pi,50)),sin(linspace(0,2*pi,50)),zeros(1,50),[1,1,1]*.75,'facealpha',.5); hold on; 
for ii = 1:numel(isovals)
[faces,verts] = isosurface(X,Y,Z,v,isovals(ii));
patch('vertices',verts,'faces',faces,'facecolor',cmap(ii,:),'facealpha',.25,'edgecolor','none'); hold on
end
axis equal

