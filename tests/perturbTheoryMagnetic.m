clc; clear all; close all; 
ConstantsUnits0; cols = flatcolors(); 
%% Mesh
mesh = geomEllipse([1,1],[],.1,1);
L = 200e-9; 
%% System matrices
DRa = calcDifferential(mesh.p,mesh.t,[0,1;-1,0]);
[DRi,DL] = calcDifferential(mesh.p,mesh.t); 
V = calcCoulomb(mesh.p,mesh.t); 
Da = -full(DL\DRa); 
Di = -full(DL\DRi); 

%% Graphene
ef = .2;
gam = 1e-3;
ene = linspace(0.001,.1,100);
sigma = conducIntra(ene,ef,gam); 
B = 4;
[~,sigmaxx,sigmaxy] = conducMagnetoClassical(ene,ef,B,gam);
f(1,1,:) = sigmaxx./sigma;
f(1,2,:) = sigmaxy./sigma;
f(2,1,:) = -sigmaxy./sigma;
f(2,2,:) = sigmaxx./sigma;
omegac = e*B*vf^2/(ef*ev2jo)*hbar_eV;

%% Absorption
polar = [1,0;0,1]; 
Q_abs = calcAbsorption(mesh.p,mesh.t,ene/hbar_eV,L,polar,sigma,f);

%% "Bare" eigensystem
[zeta0,rhoeig0] = calcEigen(V,DL,DRi,2);
ene0 = hbar_eV/hbar*sqrt(e^2*(ef*ev2jo)/(2*pi*eps0)*zeta0/L);

%% Perturbation terms
braket_overlap = integrateMeshFunction(mesh.p,mesh.t,rhoeig0.*(V*rhoeig0))*L/(4*pi*eps0);
braket_ani = integrateMeshFunction(mesh.p,mesh.t, (V*rhoeig0).*(Da*V*rhoeig0))*L/(4*pi*eps0);
braket_ani2 = integrateMeshFunction(mesh.p,mesh.t, (Di\rhoeig0).*(Da/Di*rhoeig0)).*zeta0.^2*pi*L/eps0;
brakets = braket_ani./braket_overlap;

prefactor = brakets./zeta0/(4*pi);
prefactor = [1i;-1i]/2; 
delta = i*prefactor.*(1-1i*gam./ene0/2).*omegac;
%delta(1:2) = [-1,+1]*omegac/4; 
%% Plotting
set_figsize(2,20,15)
plot([1,1]*omegac,[0,2],':','color',cols{8}); hold on
hold on
plot(ene,Q_abs(:,1),'-.','color',cols{1},'LineWidth',2) %x
plot(ene,Q_abs(:,2),':','color',cols{2},'LineWidth',2)  %y
plot(ene0,ones(size(ene0)),'.')
plot(ene0+delta,1.25*ones(size(ene0)),'ok','MarkerFaceColor',[1,1,1]*.7)