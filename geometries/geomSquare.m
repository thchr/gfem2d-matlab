function [p,t] = geomSquare(a,nside,hdata,plotmesh)
%CALL:           [p,t] = geomHoleSquare(a,nside,hdata,plotmesh)
%Create the triangulation data p (points) and t (connectedness) for a
%square of side 'a' with 'nside' points along the periphery, and
%mesh-density specified by 'hdata'. Set 'plotmesh' to 1 to disply mesh.

%Default input
if ~exist('a','var') || isempty(a); a = 4; end
if ~exist('nside','var') || isempty(nside); nside = 20; end
if nside < 3; error('Number of points per side must be at least 3 in this implementation'); end

%Outer boundary (square)
step = 2/nside; nsquare = 4*(nside-1);
node(:,1) = [linspace(-1,1,nside).';
                linspace(1,1,nside-1).';
                linspace(1-step,-1,nside-1).';
                linspace(-1,-1,nside-2).'];
node(:,2) = [linspace(-1,-1,nside).';
                linspace(-1+step,1,nside-1).';
                linspace(1,1,nside-1).';
                linspace(1-step,-1+step,nside-2).'];
node = node*(a/2); 


%Connectedness of node points
cnct = [(1:nsquare).', [(2:nsquare).';1]];


%Mesh density via hdata
if ~exist('hdata','var') || isempty(hdata);
    hdata.fun = @(x,y) .1*ones(size(x));
elseif isnumeric(hdata)
    constval = hdata; hdata = [];
    hdata.fun = @(x,y) constval*ones(size(x));
end

% Make mesh
[p,t,stats] = mesh2d(node,cnct,hdata,struct('plot',false)); disp(stats)
[p,t] = smoothmesh(p,t);

%Figure
if exist('plotmesh','var') && ~isempty(plotmesh) && plotmesh == 1
    plotWireMesh(p,t)
end