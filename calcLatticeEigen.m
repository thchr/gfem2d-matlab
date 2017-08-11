function [zeta,eigV,W,DL,DR,ene_eV,lambda_nm] = calcLatticeEigen(blochmesh,k,f,L,sigma,epsB,numeigs,nmax)
%CALL: [Q_abs,absorb,W,DL,DR] =
%        calcLatticeEigen(blochmesh,k,f,L,sigma,epsB,numeigs,nmax)
%Calculates the eigenmodes of a periodic system given specified by
%mesh 'blochmesh', wave vector 'k', occupation 'f', length scale 'L' and
%baseline conductivity 'sigma'. 
%NOTE: This function does not seem to be quite updated relative to some
%other utilities, e.g. it is probably better to simply use calcEigen() and
%calcEigenEnergiesAny().

%If unspecified; calcLatticeEwaldCoulomb subroutine chooses its own default
if ~exist('nmax','var');
    nmax = [];
end

%If unspecified; calcEigen subroutine chooses its own default for number of eigenmodes
if ~exist('numeigs','var');
    numeigs = [];
end

%Write purpose to prompt
fprintf('Calculating the eigenvalues of the system at k = [%g,%g]\n',k(1),k(2))

%Calculate the lattice-summed Coloumb potential
W = calcLatticeCoulombOptimized(blochmesh,k,nmax);

%Calculates the differential matrices
if ~exist('f','var') || isempty(f) || (isnumeric(f) && f == 1) %Isotropic scenario
    [DR,DL] =  calcLatticeDifferentialIsotropic(blochmesh,k); 
    anisotropic = 0;
else
    error('Not yet implemented\n')
    [DRiso,DL] = calcLatticeDifferentialIsotropic(blochmesh,k);
    DRani =  calcLatticeDifferentialAnisotropic(blochmesh,k);
    anisotropic = 1;
end

%Inversion of the generalized matrix system (find eigenvalues/vectors)
fprintf('Matrices assembled; solving eigenproblem by inversion\n')
[zeta,eigV] = calcEigen(W,DL,DR,numeigs,[],true);


%% Calculate eigenenergies (and wavelengths) if sigma is provided
if exist('sigma','var') && ~isempty('sigma')
    ConstantsUnits0;
    guess_eV = 0.01; %%sqrt(e^2*0.3*ev2jo*zeta/(2*pi*eps0*L))/ev2jo; %Do a random guess from intraband approx at ef_eV = 0.3. May be fragile.
    ene_eV = zeros(1,numel(zeta));
    if exist('epsilon','var') %If a background dielectric constant is given
        if ~isa(epsB,'function_handle')
            epsB = @(omega_eV) epsB*omega_eV; %Turn into a "dummy-function", which always returns a constant value
        end
    else 
        epsB = @(omega_eV) ones(size(omega_eV));
    end
            
    for ee=1:numel(zeta)
        ene_eV(ee) = fzero(@(omega_eV) real( 2i*eps0*epsB(omega_eV).*(omega_eV/hbar_eV)*L./sigma(omega_eV/hbar_eV)-zeta(ee)),guess_eV);
    end
    
    lambda_nm = 2*pi*c./(ene_eV/hbar_eV)*1e9;
end