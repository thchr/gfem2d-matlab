function mesh = geomAnnulus(rin,nthetaout,nthetain,hdata)
%Create the triangulation data 'p' (points) and 't' (connectedness) for an
%annulus of outer radius 1 and inner radius 'rIn' (<1) with 'nThetaIn' and 
%'nThetaOut' boundary points (in the small and large circles, respectively)
%and (possibly specified) mesh characteristics specified by 'hdata'.

%Default values
if ~exist('rin','var'); rin = .5; end
if ~exist('nthetain','var'); nthetaout = 80; end
if ~exist('nthetain','var'); nthetain = max(ceil(nthetaout*rin),25); end

%Outer geometry boundary
thetaout  = linspace(-pi,pi,nthetaout+1); thetaout = thetaout(1:end-1).';
nodeout   = [cos(thetaout) sin(thetaout)];

%Inner geometry boundary
thetain  = linspace(-pi,pi,nthetain+1); thetain = thetain(1:end-1).';
nodein   = rin*[cos(thetain) sin(thetain)];

%Total node list
node = [nodeout; nodein];

%Connectedness of node points
cnctout = [(1:nthetaout).',[2:nthetaout,1].'];
cnctin = nthetaout + [(1:nthetain).',[2:nthetain,1].'];
cnct = [cnctout;cnctin];


%Mesh density via hdata
if ~exist('hdata','var')
    hdata.fun = @(x,y) real(0.01*((sin(pi*(sqrt(x.^2+y.^2)-rin)./(1-rin))).^1.75+2).^2);%.005*(cos(pi/2*(x.^2+y.^2))+2).^2;%.025./(x.^2 + y.^2+.5);
elseif ischar(hdata) && strcmpi(hdata,'const')
    hdata.fun = @(x,y) .15*ones(size(x));
elseif isnumeric(hdata)
    hdata.fun = @(x,y) hdata*ones(size(x));
end

% Make mesh
[p,t] = mesh2d(node,cnct,hdata);
[mesh.p,mesh.t] = smoothmesh(p,t); 