function W = calcRibbonCoulombViaSeries(meshstruct,k,n_genexpint,r)
%CALL:      W = calcRibbonCoulombViaSeries(meshstruct,k,n_genexpint)
%Lattice calculation for structure that are only 1D periodic, which
%utilizes a series-expansion of the remaining integral in the 1D Ewald
%technique. n_genexpint indicates the number of terms included in the
%truncated series [default=10, which the method is, in practice, optimized
%roughly for (accuracy of series representation is then at worst ~0.5^10)]
%If a 4th argument is passed, 'r' the value can be obtained at general
%points in the 2D plane (i.e. not necessarily the mesh vertices)
%
%AUTHOR: Thomas Christensen (28 April, 2016)
%                        (edits to allow general points on 25 August, 2016)

if ~exist('n_genexpint','var') || isempty(n_genexpint) %Typically, n_genexpint = 10 will do well.
    n_genexpint = 10;
end

fprintf('Calculates the (ribbon-summed) Coulomb matrix via summation in 16 subtriangle elements V\n')

%Unwrapping the necessary elements of the blochmesh structure
if nargin <= 3
    p = meshstruct.p;
else
    if iscell(r) == 1 %Assume input as {X,Y}
        p = [r{1}(:),r{2}(:)];
    else %Interpret as Nx2 column vector
        p = r; 
    end
end
t = meshstruct.t;
areas = meshstruct.area;
R = meshstruct.Rrib; %Lattice vector
G = meshstruct.Grib; %Reciprocal lattice vector
remesh = meshstruct.remesh;


%Number of vertices and elements (triangles)
Nvert = size(meshstruct.p,1);
Ntri = size(t,1);

%Number of points to evaluate Coulomb operator at 
Npoints = size(p,1); 

%Unit cell length
L_uc = norm(R,2);

%Find maximum xperp and, based on this, choose a Ewald splitting parameter, 
%such that the exponential integral sum converges well
twiddlefac = .5; %"By hand" adjustment of overall max magnitude of E*drperp
pardir = R./norm(R,2); 
rperpdiffmax = diff(minmax(pardir*p.')); 
E = twiddlefac/L_uc/rperpdiffmax;
fprintf('   Splitting parameter E: %.3f\n',E)

%Then we choose an nmax, for the Ewald sum, to ensure that all erf-like
%terms have converged to their asymptotic values
vcut_sr = 6;  
nmax_sr = ceil(vcut_sr/(norm(R,2)*E));  
vcut_lr = 6;  
nmax_lr = ceil(vcut_lr^2*4*E^2/norm(k+G,2).^2)+1;
fprintf('   Ewald nmax | short-range sum: %.3f\n',nmax_sr)
fprintf('              | long-range sum: %.3f\n',nmax_lr)

%Direct and reciprocal lattice vectors  needed in the Ewald scheme are
%computed once and for all here
Rn = bsxfun(@times,(-nmax_sr:nmax_sr).',R); %[(2*nmax+1) x 2] vector
Gn = bsxfun(@times,(-nmax_lr:nmax_lr).',G); %[(2*nmax+1) x 2] vector

%Reocurring exponential elements needed in the Ewald scheme
expiRnk = exp(1i* ( k(1)*Rn(:,1) + k(2)*Rn(:,2) ) );
expikr  = exp(-1i*( k(1)*p(:,1)  + k(2)*p(:,2)  ) );

%The generalized exponential integral values and associated prefactors 
EnkG = zeros(size(Gn,1),n_genexpint);
for gg=1:size(Gn,1);
    kG2E2 = ((k+Gn(gg,:))*pardir.')^2/(4*E^2);
    for nn=1:n_genexpint
        EnkG(gg,nn) = genexpint(nn,kG2E2); %gen. exp. integrals
    end
end

%Gaussian quadrature weights
GaussQuadMat = [ 10, 7, 8, 7, 4, 5, 4, 5, 4,  1, 2, 1, 2, 1, 2,  1 ; ...
                  1, 4, 2, 1, 7, 5, 4, 2, 1, 10, 8, 7, 5, 4, 2,  1 ; ...
                  1, 1, 2, 4, 1, 2, 4, 5, 7,  1, 2, 4, 5, 7, 8, 10 ];

%Preallocation
W = zeros(Npoints,Nvert); 

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
    %vectors for the specified lattice R at momentum k (via funW1 & funW2)
    %Each element is finally the summation of contributions from all sub-
    %triangles, added into the total Coulomb matrix at the vertex sites
    W(:,t(j,:)) = W(:,t(j,:)) + ( GaussQuadMat * ... 
                   ( funW1(p,cpnt,k,Gn,L_uc,E,EnkG) + ...
                     funW2(p,cpnt,k,Rn,expiRnk,expikr,E) ) ).'*s;
end
W = W/192; %Normalize by (straggling numerical factors)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%     SUBFUNTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W1 = funW1(r,rp,k,Gn,L_uc,E,EnkG)



%Number of points in r and rp (assumed in column form)
numr = size(r,1); numrp = size(rp,1);
pardir = Gn(end,:)./norm(Gn(end,:),2);  %This is the direction along the lattice vector
dr = [reshape(bsxfun(@minus,r(:,1).',rp(:,1)),numrp*numr,1), ...
      reshape(bsxfun(@minus,r(:,2).',rp(:,2)),numrp*numr,1)]; %[Nx2 column]

%Perpendicular parts of dr-vector (perp to R, G, and k)
drperp = dr - bsxfun(@times,pardir.',pardir*dr.').'; %Components [Nx2 column]
drperp2 = drperp(:,1).^2+drperp(:,2).^2;        %Squared magnitude [Nx1 column]

pref2 = E^2*drperp2; 
fGn = zeros(numel(drperp2),size(Gn,1)); %Preallocate
pref2n = ones(size(E^2*drperp2));
for nn = 1:size(EnkG,2);
    for gg = 1:size(Gn,1); %Integrate
        if ~all(Gn(gg,:) == [0,0]) || ~all(k == [0,0]) || nn > 1 %All nonsingular terms where either Gn ~= 0 or k ~= 0
            fGn(:,gg) = fGn(:,gg) + pref2n*EnkG(gg,nn);
        end
    end
    pref2n = pref2n.*(-pref2/nn);
end

expiGndr = exp(-1i*( bsxfun(@times, dr(:,1), Gn(:,1).') + bsxfun(@times, dr(:,2), Gn(:,2).') ) ); %NOTE: There is a chance that there is a sign error here - ultimately, it does not matter though, because only the G sum is even (I believe)

W1 = sum(fGn.*expiGndr,2); %Efficient summation approach by matrix multiplication

%Reshape to the size we would get from bsxfun'ing r.' with rp
W1 = reshape(W1, [numrp, numr])/L_uc; 


function W2 = funW2(r,rp,k,Rn,expiRnk,expikr,E)

%Number of points in r and rp (assumed in column form)
numr = size(r,1); numrp = size(rp,1);
%Calculate W2 by vector operations only (the following two calculations, i.e.
%absdrRnm and W2, take >95% of all evaluation-time for calcLatticeCoulombOptimized!)
absdrRn = pdist2([reshape(bsxfun(@minus, r(:,1).', rp(:,1)),numrp*numr,1),...  %Calculates all distances between (r, rp), and Rnm 
                  reshape(bsxfun(@minus, r(:,2).', rp(:,2)),numrp*numr,1)],... %This implementation is faster than a sqrt-implementation
                  Rn);                                                         %because it utilizes a .mex function (so it is C I guess)
W2 = (erfc(absdrRn*E)./absdrRn) * expiRnk; %The matrix multiplication takes care of a summation-type operation
prefactor = bsxfun(@times, expikr.', exp(1i*(k(1)*rp(:,1) + k(2)*rp(:,2))));
W2 =  prefactor(:).*W2;

%Reshape to the size we would get from bsxfun'ing r.' with rp
W2 = reshape(W2, [numrp, numr]);