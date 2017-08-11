function mesh = geomEllipse(r,ntheta,hdata,plotmesh)
%Create the triangulation data 'p' (points) and 't' (connectedness) for an
%annulus of outer radius 1 and inner radius 'rIn' (<1) with 'nThetaIn' and 
%'nThetaOut' boundary points (in the small and large circles, respectively)
%and (possibly specified) mesh characteristics specified by 'hdata'.

% Default values
if ~exist('r','var') || isempty(r); r = [1,.75]; end
if ~exist('ntheta','var') || isempty(ntheta); ntheta = 80; end
if ~exist('plotmesh','var') || isempty(plotmesh); plotmesh = 0; end

% Ellipse boundary
theta  = linspace(-pi,pi,ntheta+1); theta = theta(1:end-1).';
node   = [r(1)*cos(theta) r(2)*sin(theta)];

% Connectedness of node points
cnct = [(1:ntheta).',[2:ntheta,1].'];

% Mesh density via hdata
if ~exist('hdata','var') || isempty(hdata) || (ischar(hdata) && strcmpi(hdata,'const'))
    hdata.fun = @(x,y) .15*ones(size(x));
elseif isnumeric(hdata)
    hdataval = hdata; clear hdata; 
    hdata.fun = @(x,y) hdataval*ones(size(x));
end

% Make mesh
[p,t,stats] = mesh2d(node,cnct,hdata,struct('plot',false)); disp(stats);
[mesh.p,mesh.t] = smoothmesh(p,t); 

% If requested, plot the mesh
if plotmesh == 1
    plotWireMesh(p,t)
end