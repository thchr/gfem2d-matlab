clear all; clc; close all
ConstantsUnits0;

%% mesh
ntheta = 50;
[mesh.p,mesh.t,node,cnct] = geomDisk(ntheta,.05,1);

%% response operators
[DR, DL] = calcDifferential(mesh.p,mesh.t);
V = calcCoulomb(mesh.p,mesh.t);

%% external potential
r0 = [.5, 0, .1]; %position of dipole
polar = { [0,0,1] }; % orientation of dipole
phiext = zeros(size(mesh.p,1),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiext(:,pp) = calcDipoleSource([mesh.p,zeros(size(mesh.p,1),1)],r0,polar{pp});
end
hold on; plot3(r0(1), r0(2), r0(3), 'o'); axis equal ; hold off

%% finite differences for induced field
delta = 0.001;
rG(:,1) = r0(1) + [0;0];
rG(:,2) = r0(2) + [0;0];
rG(:,3) = r0(3) + [-delta; delta]/2;

VG = calcCoulombGeneralPoints(mesh.p,mesh.t,rG);

%% material
ef_eV = .4; L = 40e-9; gam_eV = 1.6e-3; %12e-3;
ene_eV_L = linspace(0.01, ef_eV, 200);

[zetaeig,rhoeig] = calcEigen(V,DL,DR);
eneeig_eV = sqrt( (e^2*ef_eV*ev2jo)/(2*pi*eps0)*(real(zetaeig)/L) ) / ev2jo;

phiext = phiext./(4*pi*eps0*L^2); % get units right for external dipole potential
Gind = zeros(size(ene_eV_L));
for ee = 1:length(ene_eV_L)
    ene_eV = ene_eV_L(ee);
    omega = ene_eV/hbar_eV;
    sigma0 = conducIntra(ene_eV,ef_eV,gam_eV);
    
    %% solving
    rhoind = ( (DL+1i*sigma0/(4*pi*eps0*L*omega)*DR*V) ) \ ...
        ( -1i*sigma0/(omega*L^2)*DR*phiext );
    %phiind = L/(4*pi*eps0) * V*rhoind;
    %plotVertexData(mesh.p,mesh.t,real(phiind))
    
    %% green function
    phiindG = L/(4*pi*eps0) * VG * rhoind;
    Eind = -diff(phiindG)/(L*delta);
    Gind(ee) = Eind/(omega^2*mu0);
    
    if mod(ee,10) == 0
        fprintf('%g/%g\n', ee, length(ene_eV_L))
    end
end


%%
imG0 = (ene_eV_L/hbar_eV)/(6*pi*c);
imG = imG0 + imag(Gind);
Fp = imG./imG0;

figure
semilogy(ene_eV_L,Fp)
hold on
for j = 1:10
    semilogy(eneeig_eV(j)*[1,1], minmax(Fp),'-k')
end
hold off


%%







%%
Nsq = 100;
xsq = linspace(-2,2,Nsq); ysq=xsq; [Xsq,Ysq] = meshgrid(xsq,ysq);
rsq = [Xsq(:), Ysq(:), ones(size(Xsq(:)))*.15];

Vsq = calcCoulombGeneralPoints(mesh.p,mesh.t,rsq);

%% square grid
set_figsize(3,20,20)
for j = 1:12
    subaxis(4,3,j)
    
    phiindsq = L * reshape(Vsq*rhoeig(:,j), Nsq, Nsq);
    
    
    contourf(Xsq, Ysq, real(phiindsq), 10); hold on;
    plot(node(:,1), node(:,2), '--k'); hold off;
    axis equal tight
    colormap(bluewhitered_mod)
    shading interp
    
    title(['#' num2str(j) ' (' num2str(eneeig_eV(j), '%.2g') ' eV)'])
end
