function [Q_abs,absorb,W,DL,DR] = calcLatticeAbsorption(blochmesh,omega,L,sigma,f,nmax)
%CALL: [Q_abs,absorb,W,DL,DR] =
%              calcLatticeAbsorption(blochmesh,omega,L,sigma,f,nmax)
%Calculates the absorption in a periodic system illuminated by a in-plane 
%sine wave wave incident at normal direction, i.e. with k = [0,0] - there
%are serious problems with the physical context of this implementation: USE
%ONLY AS A ROUGH GUIDE TO RESONANCE POSITIONS AND RELEVANCE.

ConstantsUnits0;

%Unwrap and
p = blochmesh.p;
remesh = blochmesh.remesh;
R = blochmesh.R;
areas = blochmesh.area;

%Misc preparation
area_unitcell = abs(R{1}(1)*R{2}(2)-R{2}(1)*R{1}(2)); %Area of the unit cell |cross(R{1},R{2})| (the correct normalization for a lattice)
if ~exist('nmax','var'); %If unspecified, let calcLatticeEwaldCoulomb subroutine choose its own default value
    nmax = [];
end
%Calculate the lattice-summed Coloumb potential
W = calcLatticeCoulombOptimized(blochmesh,[0,0],nmax);
phiext = - sin( p(:,1)/R{1}(1) * 2*pi);  %This is totally AD HOC because it is impossible to define a plane wave in a periodic lattice via a potential 
                                         %!!!! AS A RESULT THE FOLLOWING DOES NOT CORRESPOND TO PLANE WAVE SCATTERING !!!! [%(p(:,1)*cos(theta) + p(:,2)*sin(theta) );]

if ~exist('f','var') || (isnumeric(f) && f == 1) %Isotropic scenario
    [DR,DL] = calcLatticeDifferential(blochmesh,[0;0]); %These are not lattice-versions, because at k=0, there is no influence on the D matrices
    anisotropic = 0; 
else
    [~,DL] = calcLatticeDifferential(blochmesh,[0;0]);
    anisotropic = 1;
end

vertInt = [2,1,1;1,2,1;1,1,2]/12; %Integration matrix for vertex-specified functions

fprintf('Calculates the absorption cross-section\n')
fprintf('   Energies in range %.1f eV to %.1f eV requested, over %g points\n',min(omega)*hbar_eV,max(omega)*hbar_eV,numel(omega))

absorb=zeros(numel(omega),1); %Preallocation
timeAbs = tic;
%pxt = remesh.p(:,1); pxt = bsxfun(@times,areas,pxt(remesh.t)).';
%pyt = remesh.p(:,2); pyt = bsxfun(@times,areas,pyt(remesh.t)).'; 
remesh.phiextt = remesh.values(phiext); remesh.phiextt = bsxfun(@times,areas,remesh.phiextt(remesh.t)).';
for oo = 1:numel(omega)
    if anisotropic == 1;
        [~,DR] = evalc('calcDifferential(p,t,f(omega(oo)))'); %Run without output to command line
    end
    LHS =  (DL+1i*sigma(omega(oo))/(4*pi*eps0*L*omega(oo))*DR*W);
    RHS = -1i*sigma(omega(oo))/(omega(oo)*L)*DR*phiext;
    rho = LHS\RHS;
    remesh.rho = remesh.values(rho);
    
    absorb(oo) = sum(sum(remesh.phiextt.*(vertInt*remesh.rho(remesh.t.')),1),2); %Summing contributions from every triangle [to the dipole moment at frequency omega(oo)]
                                                                                 %with exact integration for linear basis functions    
    if mod(oo,5) == 0; 
        fprintf('   %g/%g frequencies calculated (%.1f min/%.1f min)\n',oo,numel(omega),toc(timeAbs)/60,toc(timeAbs)/60*numel(omega)/oo); 
    end
end
absorb = -0.5*(omega.').*imag(absorb)*L^3; %Prefactors;
incpower = eps0*c/2; %Arbitrary NON-CORRECT scaling, reminiscent of plane-wave 
sigma_abs = absorb/incpower;
Q_abs = sigma_abs/(area_unitcell*L^2); %Absorption efficiency