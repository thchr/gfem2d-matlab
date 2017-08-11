function [p,t] = geomPolygon(periphery,hdata,plotmesh)
%CALL:           [p,t] = geomHoleSquareArray(a,r,nside,ntheta,hdata)
%Create the triangulation data p (points) and t (connectedness) for an
%arbitrary closed polygon as specified by 'periphery' (which may be several 
%disjoint elements; then a cell array), with mesh characteristics
%specified by 'hdata' (default; constant). Set 'plotmesh' to 1 to display
%mesh.

nodes = []; cnct = [];
if ~isempty(periphery)
    if iscell(periphery)                %For several disjoint inclusions
        for ii=1:numel(periphery)
            numpincl = size(periphery{ii},1);
            cnct = [cnct; +size(nodes,1)+[1:numpincl;2:numpincl,1].'];
            nodes = [nodes;periphery{ii}];
        end
    else                                %For a single inclusion
        numpincl = size(periphery,1);
        cnct = [cnct; +size(nodes,1)+[1:numpincl;2:numpincl,1].'];
        nodes = [nodes; periphery];
    end
end

%Mesh density via hdata
if ~exist('hdata','var') || isempty(hdata);
    hdata.fun = @(x,y) .1*ones(size(x));
elseif isnumeric(hdata)
    constval = hdata; hdata = [];
    hdata.fun = @(x,y) constval*ones(size(x));
end

% Make mesh
[p,t,stats] = mesh2d(nodes,cnct,hdata,struct('plot',false)); disp(stats)
[p,t] = smoothmesh(p,t);

%Figure
if exist('plotmesh','var') && ~isempty(plotmesh) && plotmesh == 1
    plotWireMesh(p,t)
end