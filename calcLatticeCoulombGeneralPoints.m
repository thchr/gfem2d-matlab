function [V,nmax] = calcLatticeCoulombGeneralPoints(meshstruct,k,r,nmax,E)
%Takes a Mx2 vector r of points and calculates the Coulomb matrix TO these
%points from the mesh specified by meshstruct.

fprintf('Calculates the (lattice-summed) Coulomb matrix via summation in 16 subtriangle elements V\n')

%Unwrapping the necessary elements of the blochmesh structure
t = meshstruct.t;
areas = meshstruct.area;
R = meshstruct.R;
remesh = meshstruct.remesh;

%Number of vertices and elements (triangles)
Nvert = size(meshstruct.p,1);
Ntri = size(t,1);

%The fundamental reciprocal lattice vectors and unit cell area
[G,A_uc] = calcReciprocal(R); 

%The optimal integral splitting parameter in the Ewald scheme (as
%suggested by Gallinet et al.)
if nargin < 5 || isempty(E)
    E = sqrt(pi/A_uc);
end

%If unspecified, choose default value that ensure convergence with certainty
if ~exist('nmax','var') || isempty(nmax);
    erfccut = 5.5; %Cutoff; such that erfc(erfccut) has a value less than 1e-14
    nmax(1,1) = max(ceil(erfccut/(E*norm(R{1},2))+1),2);
    nmax(1,2) = max(ceil(erfccut/(E*norm(R{2},2))+1),2);
    nmax(2,1) = max(ceil(2*E*erfccut/norm(G{1},2)) + (norm(k) == 0)*1,2);
    nmax(2,2) = max(ceil(2*E*erfccut/norm(G{2},2)) + (norm(k) == 0)*1,2);
end
if numel(nmax) == 1 %Equal number of each direct/reciprocal lat. vec.
    loopR1 = -nmax:nmax; loopR2 = loopR1; loopG1 = loopR1; loopG2 = loopR1;      
elseif numel(nmax) == 2; %Different number of direct/reciprocal lat. vec.
    loopR1 = -nmax(1):nmax(1);   loopR2 = loopR1; 
    loopG1 = -nmax(2):nmax(2);   loopG2 = loopG1;
elseif all(size(nmax) == [2,2]) %Distinction in numbers of each direct/reciprocal lat. vec. along each direction
    loopR1 = -nmax(1,1):nmax(1,1);   loopR2 = -nmax(1,2):nmax(1,2);
    loopG1 = -nmax(2,1):nmax(2,1);   loopG2 = -nmax(2,2):nmax(2,2);  
end
numRvec = numel(loopR1)*numel(loopR2);  numGvec = numel(loopG1)*numel(loopG2); 
fprintf('   Ewald summation parameters | %g by %g terms in R-sum\n',numel(loopR1),numel(loopR2))
fprintf('                              | %g by %g terms in G-sum\n',numel(loopG1),numel(loopG2))

%Direct and reciprocal lattice vectors  needed in the Ewald scheme are
%computed once and for all here
Rnm(:,1) = reshape(bsxfun(@plus,loopR1'*R{1}(1),loopR2*R{2}(1)),numRvec,1);
Rnm(:,2) = reshape(bsxfun(@plus,loopR1'*R{1}(2),loopR2*R{2}(2)),numRvec,1);
Gnm(:,1) = reshape(bsxfun(@plus,loopG1'*G{1}(1),loopG2*G{2}(1)),numGvec,1);
Gnm(:,2) = reshape(bsxfun(@plus,loopG1'*G{1}(2),loopG2*G{2}(2)),numGvec,1);

%Reocurring exponential elements needed in the Ewald scheme
expiGnmr = exp(-1i*( bsxfun(@times, r(:,1), Gnm(:,1).') + bsxfun(@times, r(:,2), Gnm(:,2).') ) );
expiRnmk = exp(1i*(k(1)*Rnm(:,1) + k(2)*Rnm(:,2)));
expikr =   exp(-1i*( k(1)*r(:,1)+k(2)*r(:,2) ) );

%Reoccuring erfc elements needed in Ewald summation scheme
abskGnm = sqrt( (k(1) - Gnm(:,1)).^2 + (k(2) - Gnm(:,2)).^2 );
erfcabskGnm = erfc(abskGnm/(2*E))./abskGnm;

%Gaussian quadrature weights
GaussQuadMat = [ 10, 7, 8, 7, 4, 5, 4, 5, 4,  1, 2, 1, 2, 1, 2,  1 ; ...
                  1, 4, 2, 1, 7, 5, 4, 2, 1, 10, 8, 7, 5, 4, 2,  1 ; ...
                  1, 1, 2, 4, 1, 2, 4, 5, 7,  1, 2, 4, 5, 7, 8, 10 ];

%Preallocation
V = zeros(size(r,1),Nvert); 

%Construct the total Coulomb matrix by summing over contributions from all
%triangle elements
timeV = tic;
for j = 1:Ntri
    if mod(j,250) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeV)/60,toc(timeV)/60/(j/Ntri)); end
    s = areas(j);
    l = remesh.t(j,1); m = remesh.t(j,2); n = remesh.t(j,3);    %This must be the remeshed positions, in order to
    rl = remesh.p(l,:); rm = remesh.p(m,:); rn = remesh.p(n,:); %get the subtriangle positions correct (they are NOT on the edge)
    
    cpnt = calcQuadraturePoints(rl,rm,rn); %Center points for 16 subtriangles in the j'th triangle
    
    %Evaluates the effective (lattice-summed) Coulomb potential for these
    %vectors for the specified lattice R at momentum k (via funV1 & funV2)
    %Each element is finally the summation of contributions from all sub-
    %triangles, added into the total Coulomb matrix at the vertex sites
    V(:,t(j,:)) = V(:,t(j,:)) + ( GaussQuadMat * ... 
                   ( funV1(cpnt,k,Gnm,expiGnmr,erfcabskGnm,A_uc) + ...
                     funV2(r,cpnt,k,Rnm,expiRnmk,expikr,E) ) ).'*s;
end
V = V/192; %Normalize by (straggling numerical factors)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%     SUBFUNTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V1 = funV1(rp,k,Gnm,expiGnmr,erfcabskGnm,A_uc)

%The rp-dependent part of the nm-sum in matrix form (vector operation)
expiGnmrp = exp(1i*( bsxfun(@times, rp(:,1), Gnm(:,1).') + bsxfun(@times, rp(:,2), Gnm(:,2).' ) ) );
erfcexpGnmrp = bsxfun(@times, expiGnmrp, erfcabskGnm.');

%Calculate V1 by a vectorized approach (faster than looping)
excludekzero = Gnm(:,1) == k(1) & Gnm(:,2) == k(2); 
if any(excludekzero) %Skip any nm components where Gnm = k , which gives the divergent part - importantly, is dr-independent and so may be excluded
    erfcexpGnmrp(:,excludekzero) = []; %Just remove the "offending" column
    expiGnmr(:,excludekzero) = [];
end
V1 = erfcexpGnmrp*expiGnmr.'; %Efficient summation approach by matrix multiplication

V1 = V1*(2*pi/A_uc); %No reshaping necessary


function V2 = funV2(r,rp,k,Rnm,expiRnmk,expikr,E)

%Number of points in r and rp (assumed in column form)
numr = size(r,1); numrp = size(rp,1);
%Calculate V2 by vector operations only (the following two calculations, i.e.
%absdrRnm and V, take >95% of all evaluation-time for calcLatticeCoulombOptimized!)
absdrRnm = pdist2([reshape(bsxfun(@minus, r(:,1).', rp(:,1)),numrp*numr,1),... %Calculates all distances between (r, rp), and Rnm 
                   reshape(bsxfun(@minus, r(:,2).', rp(:,2)),numrp*numr,1)],...%This implementation is faster than a sqrt-implementation
                   Rnm);                                                       %because it utilizes a .mex function (so it is C I guess)
V2 = (erfc(absdrRnm*E)./absdrRnm) * expiRnmk; %The matrix multiplication takes care of a summation-type operation
prefactor = bsxfun(@times, expikr.', exp(1i*(k(1)*rp(:,1) + k(2)*rp(:,2))));
V2 =  prefactor(:).*V2;

%Reshape to the size we would get from bsxfun'ing r.' with rp
V2 = reshape(V2, [numrp, numr]);