function [p,t] = geomDisk(ntheta,hdata,plotmesh)
%CALL:          [p,t] = geomDisk(ntheta,hdata,plotmesh)
%DESCRIPTION: Create the triangulation data 'p' (points) and 't' 
%(connectedness) for a disk with 'ntheta' boundary points and (possibly
%specified) mesh characteristics specified by 'hdata'. Set 'plotmesh' = 1
%to depict mesh.

% Geometry boundary
if nargin == 0 || isempty(ntheta); ntheta = 402; end
theta  = linspace(-pi,pi,ntheta+1); theta = theta(1:end-1).';
node   = [cos(theta) sin(theta)];

%Mesh density via hdata
if ~exist('hdata') || isempty(hdata)
    hdata.fun = @(x,y) real(.035*(cos(pi/2*(x.^2+y.^2)).^.5+.45).^(.75));
elseif ischar(hdata) && strcmpi(hdata,'const')
    clear hdata
    hdata.fun = @(x,y) .15*ones(size(x));
elseif isnumeric(hdata)
    hdataval = hdata; clear hdata
    hdata.fun = @(x,y) hdataval*ones(size(x));
end

% Make mesh
[p,t,stats] = mesh2d(node,[],hdata,struct('plot',false)); disp(stats)
[p,t] = smoothmesh(p,t);

%Figure
if exist('plotmesh','var') && ~isempty(plotmesh) && plotmesh == 1
    plotWireMesh(p,t)
end