function mesh = geomEllipseAnnulus(rout,rin,nthetaout,nthetain,hdata)
%Create the triangulation data 'p' (points) and 't' (connectedness) for an
%annulus of outer radius 1 and inner radius 'rIn' (<1) with 'nThetaIn' and 
%'nThetaOut' boundary points (in the small and large circles, respectively)
%and (possibly specified) mesh characteristics specified by 'hdata'.

%Default values
if ~exist('rout','var'); rout = [1,.75]; end
if ~exist('rin','var'); rin = [.5,.375]; end
if ~exist('nthetaout','var'); nthetaout = 80; end
if ~exist('nthetain','var'); nthetain = max(ceil(nthetaout*sqrt(rin(1)*rin(2)/(rout(1)*rout(2)))),25); end

%Outer geometry boundary
thetaout  = linspace(-pi,pi,nthetaout+1); thetaout = thetaout(1:end-1).';
nodeout   = [rout(1)*cos(thetaout) rout(2)*sin(thetaout)];

%Inner geometry boundary
thetain  = linspace(-pi,pi,nthetain+1); thetain = thetain(1:end-1).';
nodein   = [rin(1)*cos(thetain) rin(2)*sin(thetain)];

%Total node list
node = [nodeout; nodein];

%Connectedness of node points
cnctout = [(1:nthetaout).',[2:nthetaout,1].'];
cnctin = nthetaout + [(1:nthetain).',[2:nthetain,1].'];
cnct = [cnctout;cnctin];


%Mesh density via hdata
if ~exist('hdata','var') || (ischar(hdata) && strcmpi(hdata,'const'))
    hdata.fun = @(x,y) .15*ones(size(x));
elseif isnumeric(hdata)
    hdata.fun = @(x,y) hdata*ones(size(x));
end

% Make mesh
[p,t,stats] = mesh2d(node,cnct,hdata,struct('plot',false)); disp(stats);
[mesh.p,mesh.t] = smoothmesh(p,t); 