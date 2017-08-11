function [p,t] = geomHoleSquareArray(a,r,nside,ntheta,hdata)
%CALL:           [p,t] = geomHoleSquareArray(a,r,nside,ntheta,hdata)
%Create the triangulation data p (points) and t (connectedness) for a hole
%array on a square lattice (the Wigner Seitz cell), for a square lattice
%with lattice constant L and for a hole with radius R (which must be <L/2)
%with (minimum) nTheta boundary points on the disk and mesh characteristics
%specified by hdata.

%Standard input

if ~exist('a','var') || isempty(a); a = 4; end
if ~exist('r','var') || isempty(r); r = 1; end
if ~exist('nside','var') || isempty(nside); nside = 20; end
if ~exist('ntheta','var') || isempty(ntheta); ntheta = 60; end
if nside < 3; error('Number of points per side must be at least 3 in this implementation'); end

%Outer boundary (square)
step = 2/nside; nsquare = 4*(nside-1);
nodeout(:,1) = [linspace(-1,1,nside).';
                linspace(1,1,nside-1).';
                linspace(1-step,-1,nside-1).';
                linspace(-1,-1,nside-2).'];
nodeout(:,2) = [linspace(-1,-1,nside).';
                linspace(-1+step,1,nside-1).';
                linspace(1,1,nside-1).';
                linspace(1-step,-1+step,nside-2).'];
nodeout = nodeout*(a/2); 

%Inner boundary (hole)
theta  = linspace(-pi,pi,ntheta+1); theta = theta(1:end-1).';
nodein   = r*[cos(theta) sin(theta)];

%Total node list
node = [nodeout; nodein];

%Connectedness of node points
cnct = [(1:nsquare).', [(2:nsquare).';1]];
cnct = [cnct; ((nsquare+1):(nsquare+ntheta)).',...
                           [((nsquare+2):(nsquare+ntheta)).';nsquare+1]  ];

%Mesh density via hdata
if ~exist('hdata','var')
    hdata.fun = @(x,y) .0125.*(sqrt(x.^2 + y.^2).^4+2);
elseif ischar(hdata) && strcmpi(hdata,'const')
    hdata.fun = @(x,y) .005*ones(size(x));
elseif isnumeric(hdata)
    constval = hdata; hdata = [];
    hdata.fun = @(x,y) constval*ones(size(x));
end

% Make mesh
[p,t,stats] = mesh2d(node,cnct,hdata,struct('output',false)); disp(stats);
[p,t] = smoothmesh(p,t);

%Display mesh
figure; set(gcf,'color','w');
trimesh(t,p(:,1),p(:,2),'Color',[.1,.1,.1])
title('Mesh'); axis equal tight;
drawnow;