function [mesh,nodes] = geomPeanut(adivb,width,depth,steepness,ntheta,hdata,dispstats,plotmesh)
%CALL:   [mesh,nodes] = geomPeanut(adivb,width,depth,steepness,ntheta,hdata,dispstats,plotmesh)
%DESCRIPTION: Create the triangulation data 'p' (points) and 't' 
%(connectedness) for a peanut-like geometry.
%INPUT: 
%   adivb | x-axis radius divided by y-axis radius of initial ellipse
%   width, depth, steepness | parameters in a squeeze function which is
%       applied to the initial ellipse along the y-direction:
%           squeeze = exp(-abs(x).^abs(steepness)/width)
%       depth = 0 is no squeezing, depth = 1 is (singular) complete squeeze
%   ntheta | number of  boundary points (currently just distributed in equal
%            angles, rather than equal segments, which would be preferable)
%   hdata | This may be specified by the user, to indicate mesh size
%           function.
%Set 'dispstats' = 1 to print mesh attributes and 'plotmesh' = 1 to plot
%mesh.
%OUTPUT:
%   mesh | structure with nodal elements 'p' and connectivity 't'
%   nodes | boundary points of structure in [x y] array


%% DEFAULTS
% Geometry boundary
if ~exist('adivb','var') || isempty(adivb); adivb = 2; end
if ~exist('width','var') || isempty(width); width = .5; end
if ~exist('depth','var') || isempty(depth); depth = .5; end
if ~exist('steepness','var') || isempty(steepness); steepness = 2; end
if ~exist('ntheta','var') || isempty(ntheta); ntheta = 50; end
if ~exist('dispstats','var') || isempty(dispstats); dispstats = 0; end
if ~exist('plotmesh','var') || isempty(plotmesh); plotmesh = 0; end

disp([adivb,width,depth,steepness])

%% BOUNDARY 
theta  = linspace(0,2*pi,ntheta+1); theta = theta(1:end-1).';
nodes   = [cos(theta) sin(theta)/adivb];

squeeze = @(x) exp(-abs(x).^abs(steepness)/width);
nodes(:,2) = nodes(:,2).*(1-depth*squeeze(nodes(:,1))); 

%Mesh density via hdata
if ~exist('hdata') || isempty(hdata)
    dist = pdist2(nodes,nodes); dist(dist==0) = inf; meandist = mean(min(dist));
    hdata.fun = @(x,y) meandist*ones(size(x));
elseif ischar(hdata) && strcmpi(hdata,'const')
    hdata.fun = @(x,y) .1*ones(size(x));
elseif isnumeric(hdata)
    hdata.fun = @(x,y) hdata*ones(size(x));
end

%% MESHING
% Make mesh
[p,t,stats] = mesh2d(nodes,[],hdata,struct('plot',false,'output',false)); 
[mesh.p,mesh.t] = smoothmesh(p,t);
if dispstats == 1
    disp(stats)
end

%Figure
if exist('plotmesh','var') && ~isempty(plotmesh) && plotmesh == 1
    plotWireMesh(p,t)
end