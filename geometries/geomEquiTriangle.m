function [p,t,node,cnct] = geomEquiTriangle(side,nside,hdata,plotmesh)
%Create the triangulation data p (points) and t (connectedness) for an
%equilateral triangle of with specified sidelength 'side' and with mesh
%characteristics specified by hdata and nside (>=2) points on every side.

if ~exist('side','var') || isempty(side);   side = 1;   end
if ~exist('nside','var') || isempty(nside); nside = 20; end

% Geometry boundary
corners = [cos(pi/2),sin(pi/2); cos(-pi*1/6), sin(-pi*1/6); cos(pi*7/6), sin(pi*7/6)]*side/sqrt(3);
node(:,1) = [linspace(corners(1,1),corners(2,1),nside).'; 
             linspace(corners(2,1),corners(3,1),nside).';
             linspace(corners(3,1),corners(1,1),nside).'];
node(:,2) = [linspace(corners(1,2),corners(2,2),nside).'; 
             linspace(corners(2,2),corners(3,2),nside).';
             linspace(corners(3,2),corners(1,2),nside).'];
%Remove duplicates
node([nside+1,2*nside+1,3*nside],:) = []; 
%Connectedness
cnct = [(1:size(node,1)).',[2:size(node,1),1].'];

%Mesh density via hdata
if ~exist('hdata','var') || isempty(hdata)
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