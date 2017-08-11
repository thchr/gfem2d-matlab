clear all; close all;  
ConstantsUnits0; 

%% SET UP A NUMBER OF MODELS
clc;
ef_eV = .2; L = 50e-9;

%models = struct('type','magneto',  'ef_eV',ef_eV,'gam_eV',6e-3,'L',L,'B',8,'system','graphene'); N = 388*2+2;
models = struct('type','isotropic', 'ef_eV',ef_eV,'gam_eV',6e-3,'L',L,      'approx','intra'); N = 1; 
%[mesh.p,mesh.t] = geomDisk(50,struct('fun',@(x,y) ones(size(x))*.085),0);
[mesh.p,mesh.t] = geomDisk(50,struct('fun',@(x,y) ones(size(x))*.045),0);
%%
%ene_eV=calcEigenEnergiesAny(mesh,models);
%%
timeenv = @(T) sin(1*T).*exp( -(T-4).^2/1.5);
spatialenv = calcDipoleSource(mesh.p,[0,0,.25],[0,1,0].');
phiext = @(T) bsxfun(@times,spatialenv,timeenv(T).');
dTphiext = @(T) -mesh.p(:,1)*cos(T);
Tmax = 60;

set_figsize(1,20,20);
plot(linspace(0,Tmax,400),timeenv(linspace(0,Tmax,400))); drawnow
[T,rho,phiind,V] = timeEvolution(mesh,models,Tmax,phiext);

 

%%
try close(3); end
set_figsize(3,20,20);
plot(T,max(abs(phiind),[],1),'.-'); hold on
if exist('ene_eV','var')
semilogy(T,abs(real(exp(-1i*T/ef_eV*(ene_eV(N))-1i*2.25)))*1.95/pi,'-r'); 
end
hold off; drawnow

%%
%{
try close(2); end
set_figsize(2,20,20);
for tt = 1:numel(T)
    
trisurf(mesh.t,mesh.p(:,1),mesh.p(:,2),phiind(:,tt))
zlim(minmax(phiind)*2); 
caxis(minmax(phiind))
colormap(bluewhitered_mod)
title(num2str(T(tt)))
drawnow; pause(3/numel(T))
end
%}
%%
x = linspace(-1.5,1.5,100);
[x,y] = meshgrid(x,x); r = [x(:),y(:)];
Vsq = calcCoulombGeneralPoints(mesh.p,mesh.t,r);
phiindsq = Vsq*rho; 
phiindsq = reshape(phiindsq,[size(x),numel(T)]);
%%
visstyle = 'reg'; %abs or reg(ular)
cols=flatcolors;
try close(4); end
set_figsize(4,20,20);
for tt = 1:numel(T)
    if strcmpi(visstyle,'reg')
        pcolor(x,y,phiindsq(:,:,tt)); hold on; 
        caxis(minmax(phiind)); 
        colormap(blueblackred)
    elseif strcmpi(visstyle,'abs')
        pcolor(x,y,abs(phiindsq(:,:,tt))); hold on; 
        caxis(minmax(abs(phiind(:))));
        colormap(morgenstemning)
    elseif strcmpi(visstyle,'logabs')
        pcolor(x,y,log(abs(phiindsq(:,:,tt)))); hold on; 
        caxis(log([1e-3,1]*max(abs(phiind(:)))));
        colormap(morgenstemning)
    end
plot(cos(linspace(0,2*pi,50)),sin(linspace(0,2*pi,50)),':w'); hold off

title(num2str(T(tt)))
drawnow; pause(1/numel(T))
end
