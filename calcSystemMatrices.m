function [V,DL,DRi,DRa] = calcSystemMatrices(mesh,bloch)
%CALL:      [V,DL,DRi,DRa] = calcSystemMatrices(mesh,bloch)
%INPUT: 'mesh'  | structure with fields p (vertices) and t (connections)
%       'bloch' | structure with fields
%             type | '0D' for a finite calculation (default)
%                    '1D' for a ribbon calculation
%                    '2D' for a lattice calculation
%             k | k-vector [1x2 array] specifying the bloch momentum
%                 (necessary for '1D' and '2D')
%             n_genexpint | If bloch.type is '1D' this value must be passed
%                           to specify number of generalized exponential
%                           integral terms in the Taylor expansion of the
%                           long-range Coulomb term in the Ewald scheme
%             supercell | If bloch.type is '1D', we can alternatively
%                         calculate the Coulomb interaction by a supercell
%                         scheme; in that case, parse this value as a non-
%                         empty array with 2 elements indicating
%                         multiplication onto R{1} and R{2} (extension).
%OUTPUT: V | Coulomb operator
%        DL | Left-hand side continuity eq. operator
%        DRi | Right-hand side continuity eq. operator (isotropic part)
%        DRa | Right-hand side continuity eq. operator (anisotropic part)
%DESCRIPTION: Calculates the system matrices of a mesh (in 'mesh'). 
%Periodicity and dimensionality indicated in 'bloch' structure.

%Default values
if ~exist('bloch','var') || isempty(bloch)
    bloch.type = '0D';
end

fprintf('Calculating system matrices\n');
%Calculate Coulomb matrix
if ~isfield(mesh,'V')
    ticV = tic;
    switch bloch.type
        case '0D' %Finite system
            [~,V] = evalc('calcCoulomb(mesh.p,mesh.t)');
        case '1D' %Extended 1D lattice; generalized ribbons
            if ~isfield(bloch,'supercell') || isempty(bloch.supercell) || bloch.supercell == 0 %We use the 1D Ewald summation scheme via
                [~,V] = evalc('calcRibbonCoulombViaSeries(mesh,bloch.k,bloch.n_genexpint)');   %series summation in generalized exponentials
            else %Employ a supercell approximation
                [~,A_uc] = calcReciprocal(mesh.R); E = sqrt(pi/A_uc); %Get a useful splitting parameter before we change the lattice vectors
                mesh.R{1} = bloch.supercell(1)*mesh.R{1}; mesh.R{2} = bloch.supercell(2)*mesh.R{2};
                [~,V,nmax] = evalc('calcLatticeCoulombOptimized(meshstruct,k,[],E)');
                fprintf('      [Supercell approx calculated with Ewald parameters]\n      [%g by %g terms in R-sum | %g by %g terms in G-sum]',nmax(1,1),nmax(1,2),nmax(2,1),nmax(2,2))
            end
        case '2D' %Extended 2D lattice
            [~,V] = evalc('calcLatticeCoulombOptimized(mesh,bloch.k)');
    end
    fprintf('   Coulomb matrix V calculated in %.1f min\n',toc(ticV)/60);
else
    V = mesh.V; clear mesh.V;
    fprintf('   Coulomb matrix V supplied in mesh structure\n')
end

%Calculate differential matrices DL and DR(i,a)
if ~isfield(mesh,'DRi') && ~isfield(mesh,'DL')
    ticDi = tic; %Isotropic parts of D matrices
    switch bloch.type
        case '0D' %Finite system
            [~,DRi,DL] = evalc('calcDifferential(mesh.p,mesh.t)');
        case {'1D','2D'} %1D or 2D periodic system
            [~,DRi,DL] = evalc('calcLatticeDifferentialIsotropic(mesh,bloch.k)');
    end
    fprintf('   Isotropic differential matrices DL and DRi calculated in %.1f min\n',toc(ticDi)/60);
else
    DRi = mesh.DRi; DL = mesh.DL;
    fprintf('   Isotropic differential matrices DL and DRi supplied in mesh structure\n')
end

if nargout == 4 %Anisotropic parts of DR matrix
    if ~isfield(mesh,'DRa')
        ticDa = tic;
        switch bloch.type
            case '0D' %Finite system
                [~,DRa] = evalc('calcDifferential(mesh.p,mesh.t,[0,1;-1,0])');
            case {'1D','2D'} %1D or 2D periodic system
                [~,DRa] = evalc('calcLatticeDifferentialAnisotropic(mesh,bloch.k)');
        end
        fprintf('   Anisotropic part of differential matrix DRa calculated in %.1f min\n',toc(ticDa)/60);
    else
        DRa = mesh.DRa;
        fprintf('   Anisotropic part of differential matrix DRa supplied in mesh structure\n');
    end
end