function [p,t] = geomEquiTriangle(side,hdata,plotmesh)
%Create the triangulation data p (points) and t (connectedness) for an
%equilateral triangle of with specified sidelength 'side' and with mesh
%characteristics specified by hdata.

if ~exist('side','var') || isempty(side)
    side = 1;
end

% Geometry boundary
node = [cos(pi/2),sin(pi/2); cos(pi*7/6), sin(pi*7/6); cos(-pi*1/6), sin(-pi*1/6)]*side/sqrt(3);
%node = linspace(node(1,:),node(2,:),10)
%Mesh density via hdata
if ~exist('hdata','var') || isempty(hdata)
    hdata.fun = @(x,y) .05./sqrt(x.^2 + y.^2)+.01;
elseif ischar(hdata) && strcmpi(hdata,'const')
    hdata.fun = @(x,y) side*.025*ones(size(x));
elseif isnumeric(hdata)
    hdata.fun = @(x,y) hdata*ones(size(x));
end

% Make mesh
[p,t,stats] = mesh2d(node,[],hdata,struct('plot',false)); disp(stats)
[p,t] = smoothmesh(p,t);


%Figure
if exist('plotmesh','var') && ~isempty(plotmesh) && plotmesh == 1
    plotWireMesh(p,t)
end