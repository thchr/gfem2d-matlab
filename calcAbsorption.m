function [Q_abs,sigma_abs] = calcAbsorption(p,t,omega,L,polar,sigma,f)
%CALL: [Q_abs,sigma_abs] = calcAbsorption(p,t,omega,L,polar,sigma,f)
%DESCRIPTION: Calculates the absorption due to an incident plane wave.
%INPUT: p,t | mesh parameters
%       omega | frequencies (rad/s)
%       L | defining length scale associated with p,t (m)
%       polar | polarization vectors [x,y], where x,y can be columns to
%               request several polarizations.
%       sigma | conductivity "base line" (SI)
%       f | occupation function relative to conduc. base line (dim-less)
%OUTPUT: Q_abs | absorption efficiency (dim-less)
%        sigma_abs | absorption cross-section (SI)

ConstantsUnits0;
areas = meshArea(p,t);

numomega = numel(omega); %Number of frequencies
numpolar = size(polar,1); %Number of polarizations

V = calcCoulomb(p,t);
phiext = - ( bsxfun(@times,p(:,1),polar(:,1).') + ... %All the plane-wave realizations 
             bsxfun(@times,p(:,2),polar(:,2).') );    %requested by 'polar'

if ~exist('f','var') || (isscalar(f) && f == 1) %Isotropic scenario
    [DR,DL] = calcDifferential(p,t);
    anisotropy = 0;
else
    [~,DL] = calcDifferential(p,t);
    anisotropy = 1;
    fprintf('Anisotropic conductivy input via 2x2xN tensor form of f\n')
end

vertInt = [2,1,1;1,2,1;1,1,2]/12; %Integration matrix for vertex-specified functions

fprintf('Calculates the plane wave absorption cross-section (mesh: %g nodes and %g elements)\n',size(p,1),size(t,1))
fprintf('   Energies in range %.1f eV to %.1f eV requested, over %g points\n',min(omega)*hbar_eV,max(omega)*hbar_eV,numel(omega))
fprintf('   Incident field polarization: ')
for pp = 1:numpolar
    fprintf(['#%g, [%s,%s] * E0\n' repmat(' ',1,32)],pp,num2str(polar(pp,1),2),num2str(polar(pp,2),2)) %The awkward implementation is due to fprint not working well with complex/imaginary numbers
end
fprintf('\n')

dipmom=zeros(numomega,2,numpolar); %Preallocation
timeAbs = tic;
pxt = p(:,1); pxt = bsxfun(@times,areas,pxt(t)).';
pyt = p(:,2); pyt = bsxfun(@times,areas,pyt(t)).';
for oo = 1:numomega
    if anisotropy == 1
        [~,DR] = evalc('calcDifferential(p,t,f(:,:,oo))'); %Run without output to command line
    end
    
    rho =  ( (DL+1i*sigma(oo)/(4*pi*eps0*L*omega(oo))*DR*V) ) \ ... %We write it this (lengthy) way, to avoid potential memory issues of the mesh is very large.
           ( -1i*sigma(oo)/(omega(oo)*L)*DR*phiext ); %A factor L is absorbed in phi; hence, the numerator is not L^2 here (as it is otherwise) but L because a multiplicatition has already been done
    
    for pp = 1:numpolar %Sum dipole contributions along x and y for each polarization
        rhopp = rho(:,pp); 
        dipmom(oo,1,pp) = sum(sum(pxt.*(vertInt*rhopp(t.')),1),2); %Summing contributions from every triangle [to the dipole moment at frequency omega(oo)]
        dipmom(oo,2,pp) = sum(sum(pyt.*(vertInt*rhopp(t.')),1),2); %with exact integration for linear basis functions
    end
    
    if mod(oo,25) == 0 %Print progress
        fprintf('   %g/%g frequencies calculated (%.1f min/%.1f min)\n',oo,numomega,toc(timeAbs)/60,toc(timeAbs)/60*numel(omega)/oo);
    end
end
dipmom = dipmom*L^3; %Prefactors;

alpha = zeros(numomega,numpolar); %Preallocate polarization object
for pp = 1:numpolar
    alpha(:,pp) = -(conj(dipmom(:,1,pp))*polar(pp,1) + conj(dipmom(:,2,pp))*polar(pp,2))/eps0;
end
sigma_abs = bsxfun(@times,omega(:)/c,imag(alpha));
Q_abs = sigma_abs/(sum(areas)*L^2);