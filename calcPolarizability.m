function [alpha,eigstruct] = calcPolarizability(mesh,numeigs,fieldsout,sympairs)
%CALL:          [alpha,eigstruct] = calcPolarizability(mesh,n,numeigs)
%DESCRIPTION: Calculates the polarizability as a function of zeta and 
%characteric length parameter L; returned as a structure of function 
%handles alpha.x(zeta,L) and alpha.y(zeta,L) [SI units]. The .x and .y 
%fields give the polarizability for perturbations with polarization along x
% and y-directions, respectively. The zeta-parameter is defined by:
%       zeta = 2i*eps0*epsB*omega*L./sigma
%The method employs an eigenexpansion approach in 'numeigs' eigenvectors.
%The absorption cross-section can subsequently be obtained from 'alpha' as:
%       sigma_abs = imag(alpha.x(zeta,L))*omega/c;
%INPUT: mesh | mesh structure containing a triangular mesh (e.g. created by
%              mesh2d), with fields 'p' (vertices) and 't' (triangulation)
%       numeigs | number of eigenvectors/values used in eigendecomposition 
%                 of polarizability (default, 48)
%       fieldsout | boolean for whether to include rho and phi fields (for 
%                   visualization with e.g. plotVertexData) (1=yes, 0=no)
%       sympairs | cell-array of two-element arrays to _force_ a symmetry
%                  between certain eigenstates (numbering must be known a 
%                  priori) while ensuring orthogonality between states.
%                  e.g. for the disk, to ensure symmetry between the n = 1
%                  and 2 dipole modes, set sympairs = {[1,2],[12,13]}
%OUTPUT: alpha | structure of two function handles alpha.x(zeta,L) (and .y)
%        eigstruct | contains fields 'zeta' and 'oscstrength' (and possibly
%                    'rho' and 'phi' if 'fieldsout' = 1) for further
%                    analysis
%REVISED: 10 July, 2017                                 Thomas Christensen


%Default values   
if ~exist('numeigs','var') || isempty(numeigs)
    numeigs=48;
end

fprintf('Calculating polarizability function for specified mesh (%g nodes and %g elements)\n',size(mesh.p,1),size(mesh.t,1))

%If V, DR, and DL are not included in the mesh structure, we calculate them.
if ~isfield(mesh,'V')
    ticV = tic; 
    [~,mesh.V] = evalc('calcCoulomb(mesh.p,mesh.t)'); %Run without textual output
    fprintf('   Coulomb matrix V calculated in %.1f min\n',toc(ticV)/60);
end
if ~isfield(mesh,'DR') || ~isfield(mesh,'DL')
    ticD = tic; 
    if isfield(mesh,'f')
        [~,mesh.DR,mesh.DL] = evalc('calcDifferential(mesh.p,mesh.t,mesh.f)'); %Run without textual output
    else
        [~,mesh.DR,mesh.DL] = evalc('calcDifferential(mesh.p,mesh.t)'); %Run without textual output
    end
    fprintf('   Differential matrices DL and DR calculated in %.1f min\n',toc(ticD)/60);
end
if ~isfield(mesh,'areavec') %Similarly, we will need the areavec for integration (calc if not given)
    [~,mesh.areavec] = meshArea(mesh.p,mesh.t);
end

%Calculate the eigenvalues and eigenmodes
[zetaeig,rhoeig] = calcEigen(mesh.V,mesh.DL,mesh.DR,numeigs);

%Force symmetry if requested; this is to overcome issues with missing
%orthogonality in very symmetric geometries; so far hand-coded to 2-dim
%degeneracies only 
if exist('sympairs','var') && ~isempty(sympairs)
    [zetaeig,rhoeig] = forceSymmetryInPairs(sympairs,mesh,zetaeig,rhoeig);
end

%Normalize eigenvectors such that <rho_n|V|rho_n> = 2*pi*zeta_n (consistent
%with choice in your thesis, in our present notation)
overlap = integrateMeshFunction(mesh.p,mesh.t,rhoeig.*(mesh.V*rhoeig),mesh.areavec);
rhoeig = bsxfun(@times,rhoeig,sqrt(2*pi*zetaeig./overlap).');
    
%Perform analytical integration over each triangular element to obtain
%normalized oscillator strengths along x and y (stored in 1st and 2nd column, respectively)
oscstrength(:,1) = abs(integrateMeshFunction(mesh.p,mesh.t,bsxfun(@times,mesh.p(:,1),rhoeig),mesh.areavec)).^2;
oscstrength(:,2) = abs(integrateMeshFunction(mesh.p,mesh.t,bsxfun(@times,mesh.p(:,2),rhoeig),mesh.areavec)).^2;

%Create function for polarizability evaluable by inserting fixed zeta and L
alpha.x = @(zeta,L) (2*L^3)*reshape( sum( bsxfun(@times, oscstrength(:,1), 1./bsxfun(@minus,zetaeig,zeta(:).') ),1 ), size(zeta) );
alpha.y = @(zeta,L) (2*L^3)*reshape( sum( bsxfun(@times, oscstrength(:,2), 1./bsxfun(@minus,zetaeig,zeta(:).') ),1 ), size(zeta) );

%Create eigstruct structure if requested as output
if nargout >= 2
    eigstruct.zeta = zetaeig;  
    eigstruct.oscstrength = oscstrength; 
    if exist('fieldsout','var') && fieldsout == 1
        eigstruct.rho = rhoeig;
        eigstruct.phi = mesh.V*rhoeig; 
    end
end

function [zetaeig,rhoeig] = forceSymmetryInPairs(sympairs,mesh,zetaeig,rhoeig)
%A Gram-Schidt step for pairs of eigenvectors that we want to force to be
%symmetric and simultaneously ensure orthogonality
%Note that since the discretization breaks strict orthogonality, the
%Gram-Schmidt'ed eigenpairs are still not strictly orthogonal after the
%procedure; one of the overlaps <a|b> or <b|a> (but small).
%This could be expanded to do stabilized Gram-Schmidt on more than two
%eigenvectors if degeneracies beyond dim=2 are involved
%(https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Algorithm)

    fprintf('   Symmetry enforced for the following mode indices:\n\t\t')
    for pair = sympairs; fprintf('[%g,%g] & ',pair{:}(1),pair{:}(2)); end; fprintf('\b\b\n\n');
    
    for paircell = sympairs
        pair = paircell{:};
        
        %Average the "symmetric-forced" eigenvalues to force equality
        zetaeig(pair) = sum(zetaeig(pair))/numel(pair); 
        overlap_ab = integrateMeshFunction(mesh.p,mesh.t,... %The first element of this array is <a|a> and the next is 
            bsxfun(@times,rhoeig(:,pair(1)),mesh.V*rhoeig(:,pair)), mesh.areavec);  %<a|b>, in the index-counting of sympairs{pp}
        
        %Gram-Schmidt step |b'> = |b> - <a|b>/<a|a>|a> (orthogonalize without normalizing step)
        rhoeig(:,pair(2)) = rhoeig(:,pair(2)) - ...
            overlap_ab(2)/overlap_ab(1)*rhoeig(:,pair(1));
        
    end
