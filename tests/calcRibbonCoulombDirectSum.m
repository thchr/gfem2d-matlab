function W = calcRibbonCoulombDirectSum(meshstruct,k,nmax)
%Lattice calculation for structure that are only 1D periodic

fprintf('Calculates the (lattice-summed) Coulomb matrix via summation in 16 subtriangle elements V\n')

%Unwrapping the necessary elements of the blochmesh structure
p = meshstruct.p;
t = meshstruct.t;
areas = meshstruct.area;
R = meshstruct.Rrib; %Lattice vector
remesh = meshstruct.remesh;

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

%If unspecified, choose default value
if ~exist('nmax','var') || isempty(nmax);
    nmax = 2 + (norm(k) == 0)*2;
end

%Gaussian quadrature weights
GaussQuadMat = [ 10, 7, 8, 7, 4, 5, 4, 5, 4,  1, 2, 1, 2, 1, 2,  1 ; ...
                  1, 4, 2, 1, 7, 5, 4, 2, 1, 10, 8, 7, 5, 4, 2,  1 ; ...
                  1, 1, 2, 4, 1, 2, 4, 5, 7,  1, 2, 4, 5, 7, 8, 10 ];

%Preallocation
W = zeros(Nvert,Nvert); 

%Construct the total Coulomb matrix by summing over contributions from all
%triangle elements
timeV = tic;
for j = 1:Ntri
    if mod(j,1) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeV)/60,toc(timeV)/60/(j/Ntri)); end
    s = areas(j);
    l = remesh.t(j,1); m = remesh.t(j,2); n = remesh.t(j,3);    %This must be the remeshed positions, in order to
    rl = remesh.p(l,:); rm = remesh.p(m,:); rn = remesh.p(n,:); %get the subtriangle positions correct (they are NOT on the edge)
    
    cpnt = calcQuadraturePoints(rl,rm,rn); %Center points for 16 subtriangles in the j'th triangle
    
    %Evaluates the effective (lattice-summed) Coulomb potential for these
    %vectors for the specified lattice R at momentum k (via funW1 & funW2)
    %Each element is finally the summation of contributions from all sub-
    %triangles, added into the total Coulomb matrix at the vertex sites
    W(:,t(j,:)) = W(:,t(j,:)) + ( GaussQuadMat * sumterms(p,cpnt,R,k,nmax ) ).'*s;
end
W = W/192; %Normalize by (straggling numerical factors)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%     SUBFUNTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Wd = sumterms(p,cpnt,R,k,nmax)

expikRn = exp(1i*(k*R(:))*(-nmax:nmax));
Wd = zeros(size(p,1),size(cpnt,1));
for nn=-nmax:nmax;
    Wd = Wd + expikRn(nn+nmax+1)./pdist2(p,bsxfun(@plus,cpnt,nn*R));
end
Wd=Wd.'.*exp(-1i* ( k(1)*bsxfun(@minus,p(:,1).',cpnt(:,1)) + k(2)*bsxfun(@minus,p(:,2).',cpnt(:,2)) ) );