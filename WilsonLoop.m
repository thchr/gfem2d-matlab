function [berry, wcc] = WilsonLoop(loop, mesh, ind, plotloops)
%CALL:     [berry, wcc, berryeig] = WilsonLoop(loop, mesh, ind, plotloops)
%DESCRIPTION: Calculates the (possibly non-Abelian) Berry phase along a
%loop (specified by 'loop'), for a structure specified by 'mesh'. Only
%bands that are indicated in the index list 'ind' are included. Multiplets
%can be included as array-elements in the 'ind' cell array.
%Present implementation only allows time-invariant (isotropic sigma) case.
%INPUT:
%   loop | struct, describing the loop, with the following fields
%          .kstart | 1x2 array of [kx,ky] start value of BZ loop
%          .G  | 1x2 reciprocal lattice vector giving the direction of the
%                loop (going from kstart to kstart+G)
%          .N  | number of steps taken along loop
%   mesh | periodic 2D mesh describing the unit cell (fields, p, t, remesh)
%   ind  | indices as cell array; each element of cell array interpreted as
%          a multiplet
%   plotloops | if not zero, will plot the loop in question in assoc. fig #
%OUTPUT:
%   berry | Berry phase along loop (vector if multiplet)
%   wcc   | Wannier "charge" center assoc. with loop (vector if multiplet)
%REFERENCE: See Phys. Rev. B 89, 115102 (2014) by Taherinejad, Garrity, &
%Vanderbilt.


% preprocessing
bloch = struct('type','2D','k',[NaN,NaN]);
allind = unique([ind{:}]);
if ~isfield(loop,'k')
    for xy = 1:2
        loop.k(:,xy) = linspace(loop.kstart(xy),loop.kstart(xy)+loop.G(xy),loop.N).';
    end
end
loop.Gnorm = sqrt(loop.G(1)^2+loop.G(2)^2);

% plot loops if plotloops is true
if exist('plotloops','var') && plotloops ~= 0
    figure(plotloops); cols = flatcolors;
        
    hold on; plot(loop.k(:,1),loop.k(:,2),'ok','MarkerFaceColor',cols{14}); hold off;
    axis tight equal; box on; drawnow;
end

% check if areavec calculated
if ~isfield(mesh.remesh,'areavec')
    [~,mesh.remesh.areavec] = meshArea(mesh.remesh.p,mesh.remesh.t);
end
%if k == [0,0], we retain zero-solution (in contrast to spurious zero-solutions)
loop.kzero = all(loop.k == 0,2);

% a remesh function (all integrals MUST be done with the full mesh, not the
% folded mesh!). This in-line function is slightly more general than
% mesh.remesh.values, allowing matrix args instead of only vector args.
remeshf = @(f) [f; f(mesh.remesh.partner(:,1),:)];

% preallocation
Lambda = cell(numel(ind),1); berry = Lambda; wcc = Lambda;
for ii = 1:numel(ind)
    Lambda{ii} = eye(numel(ind{ii}));
end

% Step through calculation by multiplying every k-point with its previous
% point
fprintf('   Running through k-points in current loop (%g steps)\n      Step #: ',loop.N)
tloop=tic; 
for kk = 1:loop.N
    tstep=tic;  fprintf('%g',kk)
    
    bloch.k = loop.k(kk,:);
    if kk ~= loop.N
        % construct the necessary matrices for the calculation
        [~,V,DL,DRi] = evalc('calcSystemMatrices(mesh,bloch)');
        
        
        % solve eigenvalue problem
        [~,~,rho] = evalc('calcEigen(V,DL,DRi,max(allind),[],loop.kzero(kk))'); % retain zero-solution if k = [0,0]
        rho = rho(:, allind);
        
        % get dual eigenvectors (potentials)
        phi = V*rho;
        
        % normalize the direct and dual eigenvectors st. <phi_n|rho_m> = delta_nm
        overlap = integrateMeshFunction(mesh.remesh.p,mesh.remesh.t,...
            remeshf(conj(rho).*phi),mesh.remesh.areavec);
        rho = bsxfun(@times,rho,sqrt(1./overlap).');
        phi = bsxfun(@times,phi,sqrt(1./overlap).');
        
    else %if kk == loop.N
        % gauge fixing required between first and last k-points in loop
        phi = bsxfun(@times, phi_first, exp(-1i*mesh.p*loop.G(:)));
    end
    
    % Wilson loop multiplication part (note that the loop does not "flip
    % around", i.e. there's no overlap term between first and last k-point)
    if kk ~= 1
        for ii=1:numel(ind)
            Lambda{ii} = Lambda{ii}*calcU(rho_prev, phi, ii, mesh, ind, allind);
        end
    else % gauge-fixing between first and last points in loop
        phi_first = phi;
    end
    
    rho_prev = rho;  % for next iteration step
    if toc(tstep) < 60; fprintf(' (%.1f s), ',toc(tstep));
    else;               fprintf(' (%.1f min), ',toc(tstep)/60); end
end
fprintf('\n      Time spent in loop: %.1f min\n',toc(tloop)/60)

% calculate (non-Abelian) berry phases and Wannier centers from Wilson loops
for ii = 1:numel(ind)
    lambda = eig(Lambda{ii});
    berry{ii} = sort(-imag(log(lambda)),'ascend');
    if nargout > 1
        wcc{ii} = (loop.Gnorm/(2*pi))*berry{ii};
    end
end

end

% SUBFUNCTIONS
function U = calcU(rho_prev, phi, ii, mesh, ind, allind)

M = zeros(numel(ind{ii}),numel(ind{ii})); % preallocate
for jj1 = 1:numel(ind{ii})
    local_jj1 = allind == ind{ii}(jj1);     % logical indexing
    for jj2 = 1:numel(ind{ii})
        local_jj2 = allind == ind{ii}(jj2); % logical indexing
        
        % compute <phi_{m,k}|rho_{n,k+delta}>
        M(jj1,jj2) = integrateMeshFunction(mesh.remesh.p,mesh.remesh.t,...
            mesh.remesh.values( conj(rho_prev(:, local_jj1)).*phi(:,local_jj2) ), ...
            mesh.remesh.areavec);
    end
end

% do the svd trick prescribed by Vanderbilt
[Vsvd,~,Wsvd] = svd(M);
U = Vsvd*Wsvd';
end
